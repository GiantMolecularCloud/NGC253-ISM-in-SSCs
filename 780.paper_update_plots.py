#############################
# GAS IN SSCS: update paper #
#############################

# redo the paper plots


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
detected_species = fnunpickle(os.path.join(XCLASSdir, 'detected_species.pickle'))

all_data = fnunpickle(join(refitdir, 'all_data.pickle'))
ratios   = fnunpickle(join(refitdir, 'ratios.pickle'))


###################################################################################################
# figure 1: all spectra
###################################################################################################

paper_models = {SSC['num']: {str(num): {} for num in np.arange(0,101)} for SSC in SSCs}

# replace the species that changed names suddenly for no reason
for ssc,ds in detected_species.items():
    for i,s in enumerate(ds):
        if s=='SO;v=0;#1':
            detected_species[ssc][i] = s.replace('#1','#2')

# create model molfit file
###################################################################################################

def create_model_molfit(num, SSC):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    for specie,data in all_data[SSC['num']].items():
        if ('#1' in specie) or ('#2' in specie):
            spx = specie
        else:
            spx = specie+';'

        molfit += f"""
{spx}   {ncomp}"""

        for c in np.arange(data['components']):
            T = data['temperature']['all'][c][num]
            v = data['velocity']['all'][c][num]
            w = data['linewidth']['all'][c][num]
            N = data['column density']['all'][c][num]

            molfit += f"""
n  0.01  10  1   n  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}   n  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}   n  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""

    # save to disk
    savepath = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    paper_models[SSC['num']][str(num)]['molfit'] = savepath


for num in tqdm(np.arange(0,101)):
    for SSC in SSCs:
        create_model_molfit(num, SSC)


# create model script
###################################################################################################

def model_spectrum_XCLASS(num, SSC):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = paper_models[SSC['num']][str(num)]['molfit']
    savepath   = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'XCLASS_model.py'))
    resultpath = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'results'))
    specpath   = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'spectrum.pickle'))
    energypath = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'energy.pickle'))
    taupath    = escape_fname(os.path.join(paperdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'opacity.pickle'))
    paper_models[SSC['num']][str(num)]['scriptfile'] = savepath
    paper_models[SSC['num']][str(num)]['spectrum']   = specpath
    paper_models[SSC['num']][str(num)]['energy']     = energypath
    paper_models[SSC['num']][str(num)]['tau']        = taupath

    with open(savepath, 'w') as f:
        f.write(f"""
import pickle

FreqMin          = 342500
FreqMax          = 358500
FreqStep         = 2.5
TelescopeSize    = 0.15
Inter_Flag       = True
t_back_flag      = True
tBack            = 0.0
tslope           = 0.0
nH_flag          = False
N_H              = 1e+24
beta_dust        = 0.1
kappa_1300       = 0.01
MolfitsFileName  = '{modelpath}'
iso_flag         = False
IsoTableFileName = ' '
RestFreq         = 0.0
vLSR             = 0.0

modeldata, log, TransEnergies, IntOptical, JobDir = myXCLASS()

os.system('mkdir -p {resultpath}')
os.system('mv '+JobDir+'/* {resultpath}/')

pickle.dump(modeldata, open('{specpath}', 'w'))
pickle.dump(TransEnergies, open('{energypath}', 'w'))
pickle.dump(IntOptical, open('{taupath}', 'w'))
""")


for num in tqdm(np.arange(0,101)):
    for SSC in SSCs:
        model_spectrum_XCLASS(num, SSC)


# execute in CASA
###################################################################################################

def CASA_command(SSC):
    XCL_files = [paper_models[SSC['num']][str(num)]['scriptfile'] for num in np.arange(0,101)]
    Xcommand = f"""
import time,datetime

start = time.time()
for f in {XCL_files}:
\texecfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
"""
    return Xcommand

commands = [CASA_command(SSC) for SSC in SSCs]
pool = Pool(len(SSCs))
pool.map(run_in_casa, commands)
pool.close()
pool.join()

fnpickle(paper_models, os.path.join(paperdir, 'paper_models.pickle'))


# plot figure 1
###################################################################################################

def plot_spectrum_overview():

    def get_spectrum(SSC, band):
        spectrum = np.genfromtxt(escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_0',band+'.dat')), dtype=None)
        return spectrum[:,0]/1000, spectrum[:,1]

    def get_models(nums, SSC, band):
        models = [ paper_models[SSC['num']][str(num)]['spectrum'] for num in nums ]
        spectra = []
        for model in models:
            with open(model, 'rb') as f:
                m = pickle.load(f, encoding="latin1")
                f = (m[:,0]*u.MHz).to(u.GHz).value
                i = (m[:,1]*u.K).value
                spectra.append([f,i])
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,0], 50, axis=0)
        p16,median,p84 = np.percentile(spectra[:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def set_up_figure():
        fig,axes = plt.subplots(nrows=len(SSCs), ncols=2, squeeze=True, sharex='col', sharey='row', figsize=A4_inches)
        fig.subplots_adjust(hspace=0, wspace=0, top=0.80, bottom=0.04, left=0.04, right=0.93)
        return fig, axes

    def plot_spectrum(ax, frequency, intensity):
        ax.plot(frequency, intensity, lw=0.5, ls='-', color='k')
        ax.fill_between(frequency, intensity, [0. for f in frequency], color='k', alpha=0.5, zorder=2)

    def plot_median_model(ax, frequency, p16, median, p84):
        ax.plot(frequency, median, lw=0.5, ls='-', color='r')
        # ax.fill_between(frequency, p16, p84, color='r', alpha=0.5, zorder=2)

    def add_SSC_overlay(ax):
        ax.text(0.5, 0.9, 'SSC '+str(SSC['no']), color='k', fontsize=10, transform=ax.transAxes, ha='center', va='top', bbox=props)

    def format_panel(ax):
        if band=='LSB':
            ax.set_xlim([342.35, 346.20])
        elif band=='USB':
            ax.set_xlim([354.25, 358.10])
        ax.set_ylim(-15,120)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, top=True, bottom=True, left=True, right=True)
        ax.set_axisbelow(True)
        ax.grid(axis='y', ls=':', c='grey')

    def get_detected_lines():
        # get detected species
        all_species = []
        for SSC in SSCs:
            for specie in detected_species[str(SSC['no'])]:
                if not specie in all_species:
                    all_species.append(specie)
        # get all lines of the detected species
        all_lines = []
        for specie in all_species:
            slines = [l for l in lines if l['XCLASS']==specie]
            for sl in slines:
                all_lines.append(sl)
        return sorted(all_lines, key=lambda k: k['restfreq'])

    def label_lines_line(ax, lines):
        for line in lines:
            restfreq = line['restfreq'].to(u.GHz).value
            ax.axvline(x=restfreq, ymin=0, ymax=1, color='dimgrey', ls='--', lw=0.5, zorder=1)

    def label_lines_text(axes, lines):
        # get lines in LSB
        lines_LSB = []
        for line in lines:
            if line['restfreq']<350*u.GHz:
                lines_LSB.append(line)
        # get lines in USB
        lines_USB = []
        for line in lines:
            if line['restfreq']>350*u.GHz:
                lines_USB.append(line)
        # label lines
        for ib, ll in enumerate([lines_LSB, lines_USB]):
            for idx,line in enumerate(ll):
                ax = axes[0][ib]
                restfreq = line['restfreq'].to(u.GHz).value
                xloc = ax.get_xlim()[0] +((idx+0.5)/len(ll))*(ax.get_xlim()[1]-ax.get_xlim()[0])
                ax.plot([restfreq,xloc], [120,150], color='dimgrey', ls='--', lw=0.5, zorder=1, clip_on=False)
                ax.text(xloc, 155, line_tex(line), color='dimgrey', fontsize=9, rotation=90, ha='center', va='bottom')

    def format_figure(axes):
        axes[-1][0].tick_params(labelbottom=True, labelleft=True, labelright=False)
        axes[-1][1].tick_params(labelbottom=True, labelleft=False, labelright=True)
        axes[-1][1].yaxis.set_label_position('right')
        for i_band,band in enumerate(['LSB','USB']):
            axes[-1][i_band].xaxis.set_visible(True)
            axes[-1][i_band].set_xlabel(r'rest frequency $\nu_\mathrm{rest}$ [GHz]', fontsize=10)
            axes[-1][i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=10)

    def save_figure(fig):
        savepath = escape_fname(os.path.join(plotdir, '14.paper', 'all_spectra.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    # figure making happens here:
    fig, axes = set_up_figure()
    dlines = get_detected_lines()
    for i_SSC,SSC in tqdm(enumerate(SSCs)):
        for i_band,band in enumerate(['LSB','USB']):
            ax = axes[i_SSC,i_band]

            frequency, intensity = get_spectrum(SSC, band)
            mfrequency, p16, median, p84 = get_models(np.arange(0,101), SSC, band)

            plot_spectrum(ax, frequency, intensity)
            plot_median_model(ax, mfrequency, p16, median, p84)

            add_SSC_overlay(ax)
            label_lines_line(ax, dlines)
            format_panel(ax)

    label_lines_text(axes, dlines)
    format_figure(axes)
    save_figure(fig)


plot_spectrum_overview()


###################################################################################################
# figure 2: sample high-res spectra
###################################################################################################

def plot_combined_variation(nums, SSC, band):
    """
    run the functions that build up the figure
    """

    def get_spectra(nums, SSC, band):
        spectra = [ np.genfromtxt(escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat')), dtype=None) for num in nums ]
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,:,0], 50, axis=0)/1000
        p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def get_models(nums, SSC, band):
        models = [ paper_models[SSC['num']][str(num)]['spectrum'] for num in nums ]
        spectra = []
        for model in models:
            with open(model, 'rb') as f:
                m = pickle.load(f, encoding="latin1")
                f = (m[:,0]*u.MHz).to(u.GHz).value
                i = (m[:,1]*u.K).value
                spectra.append([f,i])
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,0], 50, axis=0)
        p16,median,p84 = np.percentile(spectra[:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def get_models2(nums, SSC, band):
        models = [ paper_models[SSC['num']][str(num)]['spectrum'] for num in nums ]
        spectra = []
        for model in models:
            with open(model, 'rb') as f:
                m = pickle.load(f, encoding="latin1")
                f = (m[:,0]*u.MHz).to(u.GHz).value
                i = (m[:,1]*u.K).value
                spectra.append([f,i])
        spectra = np.array(spectra)
        spectra[:,0] = spectra[:,0]
        return spectra

    def set_up_figure(SSC, band):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_spectra(ax, frequency, d16, dmed, d84):
        ax.plot(frequency, dmed, lw=1, ls='-', color='k', zorder=3)
        ax.fill_between(frequency, d16, d84, color='k', alpha=0.4, zorder=2)

    def plot_fitted_spectra(ax, frequency, m16, mmed, m84):
        ax.plot(frequency, mmed, lw=1, ls='-', color='r', zorder=5)
        ax.fill_between(frequency, m16, m84, color='r', alpha=0.4, zorder=4)

    def plot_fitted_spectra2(ax, spectra):
        for spectrum in spectra:
            ax.plot(spectrum[0], spectrum[1], lw=1, ls='-', color='r', alpha=2./len(nums), zorder=5)

    def get_detected_lines(band=None):
        # get detected species
        all_species = []
        for SSC in SSCs:
            for specie in detected_species[str(SSC['no'])]:
                if not specie in all_species:
                    all_species.append(specie)
        # get all lines of the detected species
        all_lines = []
        for specie in all_species:
            slines = [l for l in lines if l['XCLASS']==specie]
            for sl in slines:
                all_lines.append(sl)
        # keep only lines of given band
        if not band==None:
            bandlines = []
            for line in all_lines:
                if band=='LSB':
                    if line['restfreq']<350*u.GHz:
                        bandlines.append(line)
                elif band=='USB':
                    if line['restfreq']>350*u.GHz:
                        bandlines.append(line)
            return sorted(bandlines, key=lambda k: k['restfreq'])
        else:
            return sorted(all_lines, key=lambda k: k['restfreq'])

    def label_lines(ax, spectrum, band):
        detected_lines = get_detected_lines(band=band)
        for idx,line in enumerate(detected_lines):
            restfreq = line['restfreq'].to(u.GHz).value
            if (restfreq>frequency[0] and restfreq<frequency[-1]):
                if band=='LSB':
                    xlim = [342.4, 346.2]
                elif band=='USB':
                    xlim = [354.3, 358.1]
                xloc = xlim[0] +((idx+0.5)/len(detected_lines))*(xlim[1]-xlim[0])
                ax.axvline(x=restfreq, ymin=0, ymax=1, color='dimgrey', ls='--', lw=0.5, zorder=1)
                ax.plot([restfreq,xloc], [1.05*np.nanmax(spectrum), 1.05*1.05*np.nanmax(spectrum)], color='dimgrey', ls='--', lw=0.5, zorder=1, clip_on=False)
                ax.text(xloc, 1.06*1.05*np.nanmax(spectrum), line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')

    def format_figure(ax, frequency, spectrum, band):
        if band=='LSB':
            ax.set_xlim([342.4, 346.2])
        elif band=='USB':
            ax.set_xlim([354.3, 358.1])
        ax.set_ylim(-0.05*np.nanmax(spectrum), 1.05*np.nanmax(spectrum))
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_axisbelow(True)
        ax.grid(axis='y', ls=':', c='grey')
        ax.set_xlabel(r'rest frequency $\nu_\mathrm{rest}$ [GHz]', fontsize=12)
        ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
        fig.set_tight_layout(True)

    def save_figure(fig, band):
        savepath = escape_fname(os.path.join(plotdir, '14.paper', 'spectra', 'SSC_'+str(SSC['no'])+'.'+band+'.spectrum.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    frequency, d16,dmed,d84  = get_spectra(nums, SSC, band)
    mfrequency, m16,mmed,m84 = get_models(nums, SSC, band)
    # spectra = get_models2(nums, SSC, band)
    fig,ax = set_up_figure(SSC, band)
    plot_spectra(ax, frequency, d16,dmed,d84)
    plot_fitted_spectra(ax, mfrequency, m16,mmed,m84)
    # plot_fitted_spectra2(ax, spectra)
    label_lines(ax, dmed, band)
    format_figure(ax, frequency, dmed, band)
    save_figure(fig, band)


for SSC in tqdm(SSCs):
    for band in ['LSB','USB']:
        plot_combined_variation(np.arange(0,101), SSC, band)


###################################################################################################
# figure 3: XDR/PDR plots
###################################################################################################

M15_table3 = fnunpickle('M15_table3.pickle')
M15_table4 = fnunpickle('M15_table4.pickle')

def plot_XDR_PDR(type):
    """
    Line ratio plot to decide on XDR vs PDR and low density vs high density as in Baan+08
    """

    fig,axes = plt.subplots(nrows=2, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    fig.subplots_adjust(hspace=0, wspace=0) #, top=0.80, bottom=0.04, left=0.04, right=0.93)
    axes[0][1].tick_params(labelbottom=True)

    # get data
    sscs   = [SSC['no'] for SSC in SSCs]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    HCO_HCN, HNC_HCN, HNC_HCO = [],[],[]
    HCO_HCN_err, HNC_HCN_err, HNC_HCO_err = [],[],[]
    for SSC in SSCs:
        try:
            hco_hcn_med = ratios['HCO+/HCN'][type][SSC['num']]['median']
            hco_hcn_p16 = ratios['HCO+/HCN'][type][SSC['num']]['16th']
            hco_hcn_p84 = ratios['HCO+/HCN'][type][SSC['num']]['84th']
            hco_hcn_low = hco_hcn_med-hco_hcn_p16
            hco_hcn_hig = hco_hcn_p84-hco_hcn_med
            HCO_HCN.append(     np.log10(hco_hcn_med) )
            HCO_HCN_err.append( [0.434*hco_hcn_low/hco_hcn_med,0.434*hco_hcn_hig/hco_hcn_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )
        try:
            hnc_hcn_med = ratios['HNC/HCN'][type][SSC['num']]['median']
            hnc_hcn_p16 = ratios['HNC/HCN'][type][SSC['num']]['16th']
            hnc_hcn_p84 = ratios['HNC/HCN'][type][SSC['num']]['84th']
            hnc_hcn_low = hnc_hcn_med-hnc_hcn_p16
            hnc_hcn_hig = hnc_hcn_p84-hnc_hcn_med
            HNC_HCN.append(     np.log10(hnc_hcn_med) )
            HNC_HCN_err.append( [0.434*hnc_hcn_low/hco_hcn_med,0.434*hnc_hcn_hig/hco_hcn_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )
        try:
            hnc_hco_med = ratios['H15NC/HCO+'][type][SSC['num']]['median']*ratios['14N/15N'][type][SSC['num']]['median']
            hnc_hco_p16 = ratios['H15NC/HCO+'][type][SSC['num']]['16th']*ratios['14N/15N'][type][SSC['num']]['median']
            hnc_hco_p84 = ratios['H15NC/HCO+'][type][SSC['num']]['84th']*ratios['14N/15N'][type][SSC['num']]['median']
            hnc_hco_low = hnc_hco_med-hnc_hco_p16
            hnc_hco_hig = hnc_hco_p84=hnc_hco_med
            HNC_HCO.append(     np.log10(hnc_hco_med) )
            HNC_HCO_err.append( [0.434*hnc_hco_low/hnc_hco_med,0.434*hnc_hco_hig/hnc_hco_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )

    # comparison from Baan+08
    B_hcn = [318.2, 14]
    B_hnc = [234.0, 7]
    B_hco = [276.1, 14]
    B_hco_hcn = [B_hco[0]/B_hcn[0], B_hco[0]/B_hcn[0]*np.sqrt((B_hco[1]/B_hco[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hcn = [B_hnc[0]/B_hcn[0], B_hnc[0]/B_hcn[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hco = [B_hnc[0]/B_hco[0], B_hnc[0]/B_hco[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hco[1]/B_hco[0])**2)]
    B_HCO_HCN = [np.log10(B_hco_hcn[0]), 0.434*B_hco_hcn[1]/B_hco_hcn[0]]
    B_HNC_HCN = [np.log10(B_hnc_hcn[0]), 0.434*B_hnc_hcn[1]/B_hnc_hcn[0]]
    B_HNC_HCO = [np.log10(B_hnc_hco[0]), 0.434*B_hnc_hco[1]/B_hnc_hco[0]]

    # comparison from Meier+15
    def get_Meier15(type):
        HCO_HCN_all = []
        HNC_HCN_all = []
        HNC_HCO_all = []
        for i in np.arange(10):
            if type == 'integrated intensity':
                try:
                    hcn = [M15_table3['HCN(1-0)'][i]['value'], M15_table3['HCN(1-0)'][i]['error']]
                    hco = [M15_table3['HCO+(1-0)'][i]['value'], M15_table3['HCO+(1-0)'][i]['error']]
                    hnc = [M15_table3['HN13C(1-0)'][i]['value'] *M15_table3['HCN(1-0)'][i]['value']/M15_table3['H13CN(1-0)'][i]['value'], M15_table3['HN13C(1-0)'][i]['error'] *M15_table3['HCN(1-0)'][i]['value']/M15_table3['H13CN(1-0)'][i]['value']]
                    hco_hcn = [hco[0]/hcn[0], hco[0]/hcn[0]*np.sqrt((hco[1]/hco[0])**2+(hcn[1]/hcn[0])**2)]
                    hnc_hcn = [hnc[0]/hcn[0], hnc[0]/hcn[0]*np.sqrt((hnc[1]/hnc[0])**2+(hcn[1]/hcn[0])**2)]
                    hnc_hco = [hnc[0]/hco[0], hnc[0]/hco[0]*np.sqrt((hnc[1]/hnc[0])**2+(hco[1]/hco[0])**2)]
                    HCO_HCN = [np.log10(hco_hcn[0]), 0.434*hco_hcn[1]/hco_hcn[0]]
                    HNC_HCN = [np.log10(hnc_hcn[0]), 0.434*hnc_hcn[1]/hnc_hcn[0]]
                    HNC_HCO = [np.log10(hnc_hco[0]), 0.434*hnc_hco[1]/hnc_hco[0]]
                except:
                    HCO_HCN = [np.nan, np.nan]
                    HNC_HCN = [np.nan, np.nan]
                    HNC_HCO = [np.nan, np.nan]
            elif type == 'column density':
                try:
                    hcn = [M15_table4['HCN'][i]['value min'], M15_table4['HCN'][i]['value max']]
                    hco = [M15_table4['HCO+'][i]['value min'], M15_table4['HCO+'][i]['value max']]
                    hn13c = [M15_table4['HN13C'][i]['value min'], M15_table4['HN13C'][i]['value max']]
                    h13cn = [M15_table4['H13CN'][i]['value min'], M15_table4['H13CN'][i]['value max']]

                    hcn = [(hcn[0]+hcn[1])/2., (hcn[1]-hcn[0])/2.]
                    hco = [(hco[0]+hco[1])/2., (hco[1]-hco[0])/2.]
                    hn13c = [(hn13c[0]+hn13c[1])/2., (hn13c[1]-hn13c[0])/2.]
                    h13cn = [(h13cn[0]+h13cn[1])/2., (h13cn[1]-h13cn[0])/2.]
                    hnc   = [hn13c[0]*hcn[0]/h13cn[0], hn13c[1]*hcn[0]/h13cn[0]]
                    hco_hcn = [hco[0]/hcn[0], hco[0]/hcn[0]*np.sqrt((hco[1]/hco[0])**2+(hcn[1]/hcn[0])**2)]
                    hnc_hcn = [hnc[0]/hcn[0], hnc[0]/hcn[0]*np.sqrt((hnc[1]/hnc[0])**2+(hcn[1]/hcn[0])**2)]
                    hnc_hco = [hnc[0]/hco[0], hnc[0]/hco[0]*np.sqrt((hnc[1]/hnc[0])**2+(hco[1]/hco[0])**2)]
                    HCO_HCN = [np.log10(hco_hcn[0]), 0.434*hco_hcn[1]/hco_hcn[0]]
                    HNC_HCN = [np.log10(hnc_hcn[0]), 0.434*hnc_hcn[1]/hnc_hcn[0]]
                    HNC_HCO = [np.log10(hnc_hco[0]), 0.434*hnc_hco[1]/hnc_hco[0]]
                except:
                    HCO_HCN = [np.nan, np.nan]
                    HNC_HCN = [np.nan, np.nan]
                    HNC_HCO = [np.nan, np.nan]
            HCO_HCN_all.append(HCO_HCN)
            HNC_HCN_all.append(HNC_HCN)
            HNC_HCO_all.append(HNC_HCO)
        return np.array(HCO_HCN_all), np.array(HNC_HCN_all), np.array(HNC_HCO_all)

    def format_panel(ax):
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_axisbelow(True)
        ax.grid(axis='both', which='both')

    def label_regions(ax):
        ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
        ax.text(0.05, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='left', va='bottom', weight='bold', fontsize=16)

    def plot_panel(ax, divide_line, divide_shade, x, y, x_err, y_err, Bx, By, Bx_err, By_err, Mx, My, Mx_err, My_err, xlim, ylim, xlabel, ylabel):
        ax.plot(divide_line[0],divide_line[1], ls='-', lw=1, c='grey', zorder=2)
        ax.fill_between(divide_shade[0],divide_shade[1],divide_shade[2], color='lightgrey', alpha=0.5, zorder=1)
        label_regions(ax)
        for a,b,a_err,b_err,c,s in zip(x, y, x_err, y_err, colors, SSCs):
            if np.isfinite(a) and np.isfinite(b):
                ax.errorbar(a,b, xerr=[[a_err[0]],[a_err[1]]], yerr=[[b_err[0]],[b_err[1]]], marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, label='SSC '+str(s['no']), zorder=5)
        ax.errorbar(Bx, By, xerr=Bx_err, yerr=Bx_err, marker='o', ms=3, lw=0, color='lime', elinewidth=1, ecolor='lime', label=r'central 440\,pc (Baan+08)', zorder=6)
        ax.errorbar(Mx, My, xerr=Mx_err, yerr=Mx_err, marker='o', ms=3, lw=0, color='darkgrey', elinewidth=1, ecolor='darkgrey', label='selected regions (40\,pc)\n(Meier+15)', zorder=4)
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])
        format_panel(ax)
        if type=='column density':
            itype = 'N'
        elif type=='integrated intensity':
            itype = 'I'
        ax.set_xlabel(xlabel.replace('X',itype), fontsize=12)
        ax.set_ylabel(ylabel.replace('X',itype), fontsize=12)

    M_HCO_HCN, M_HNC_HCN, M_HNC_HCO = get_Meier15(type)

    # panel 1: HCO+/HCN over HNC/HCO+ (top left)
    plot_panel(ax = axes[0][0],
               divide_line= [[-10,10],[10,-10]],
               divide_shade=[[-10,10],[10,-10],[10,10]],
               x=HNC_HCO,         y=HCO_HCN,         x_err=HNC_HCO_err,     y_err=HCO_HCN_err,
               Bx=B_HNC_HCO[0],   By=B_HCO_HCN[0],   Bx_err=B_HNC_HCO[1],   By_err=B_HCO_HCN[1],
               Mx=M_HNC_HCO[:,0], My=M_HCO_HCN[:,0], Mx_err=M_HNC_HCO[:,1], My_err=M_HCO_HCN[:,1],
               xlim = (-0.90,0.60),
               ylim = (-0.85,0.65),
               xlabel = '', #r'log N(HNC$^{**}$) / N(HCO$^+$)',
               ylabel = r'log X(HCO$^+$) / X(HCN)'
              )

    # panel 2: HNC/HCN over HCO/HCN (top right)
    plot_panel(ax = axes[0][1],
               divide_line=[[0,0],[-10,10]],
               divide_shade=[[0,10],[-10,-10],[10,10]],
               x=HNC_HCN,         y=HCO_HCN,         x_err=HNC_HCN_err,     y_err=HCO_HCN_err,
               Bx=B_HNC_HCN[0],   By=B_HCO_HCN[0],   Bx_err=B_HNC_HCN[1],   By_err=B_HCO_HCN[1],
               Mx=M_HNC_HCN[:,0], My=M_HCO_HCN[:,0], Mx_err=M_HNC_HCN[:,1], My_err=M_HCO_HCN[:,1],
               xlim = (-0.95,0.55),
               ylim = (-0.85,0.65),
               xlabel = r'log X(HNC$^{**}$) / X(HCN)',
               ylabel = '' #r'log N(HCO$^+$) / N(HCN)'
              )

    # panel 3: HNC/HCO over HNC/HCN (bottom left)
    plot_panel(ax = axes[1][0],
               divide_line=[[-10,10],[0,0]],
               divide_shade=[[-10,10],[0,0],[10,10]],
               x=HNC_HCO,         y=HNC_HCN,         x_err=HNC_HCO_err,     y_err=HNC_HCN_err,
               Bx=B_HNC_HCO[0],   By=B_HNC_HCN[0],   Bx_err=B_HNC_HCO[1],   By_err=B_HNC_HCN[1],
               Mx=M_HNC_HCO[:,0], My=M_HNC_HCN[:,0], Mx_err=M_HNC_HCO[:,1], My_err=M_HNC_HCN[:,1],
               xlim = (-0.90,0.60),
               ylim = (-1.05,0.45),
               xlabel = r'log X(HNC$^{**}$) / X(HCO$^+$)',
               ylabel = r'log X(HNC$^{**}$) / X(HCN)'
              )

    # panel 4: legend
    ax = axes[1][1]
    ax.set_axis_off()
    handles, labels = axes[0][0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc=3, bbox_to_anchor=(0.52,0.025,0.14,0.3), ncol=1, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    savepath = escape_fname(os.path.join(plotdir, '14.paper', 'XDR-PDR_'+type.replace(' ','_')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


plot_XDR_PDR('column density')
plot_XDR_PDR('integrated intensity')


###################################################################################################
# figure 4: temperature comparison
###################################################################################################

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
colors = ['orangered','darkorange']
temperature_species = ['H2CS;v=0;#1','SO2;v=0']

for SSC in SSCs:
    for a,(specie,color) in enumerate(zip(temperature_species,colors)):
        if specie in all_data[SSC['num']].keys():
            temp = all_data[SSC['num']][specie]['temperature']['all'][0]

            if np.nanmedian(temp)<250:
                boxplot = ax.boxplot(temp, patch_artist=True, positions=[SSC['no']-0.2+a*0.4], showfliers=False, zorder=2)
                for box in boxplot['boxes']:
                    box.set(color=color, linewidth=1)
                    box.set(facecolor=mpl.colors.to_rgba(color, alpha=0.5)) #+a*0.5))
                for median in boxplot['medians']:
                    median.set(color=color, linewidth=3)
                for whisker in boxplot['whiskers']:
                    whisker.set(color=color, linewidth=3)
                for cap in boxplot['caps']:
                    cap.set(color=color, linewidth=3)
                for flier in boxplot['fliers']:
                    flier.set(marker='.', size=8, color=color, alpha=0.5) #+a*0.5)
            else:
                lowlim = all_data[SSC['num']][specie]['temperature']['16th'][0]
                ax.scatter([SSC['no']-0.2+a*0.4], [lowlim], s=40, c=color, marker='^')

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
tracers = [Line2D([0], [0],
                  marker='s', color=color, markerfacecolor=color, markersize=10,
                  label=specie_tex(specie))
           for a,(specie,color) in enumerate(zip(temperature_species,colors))]
ax.legend(handles=tracers, loc=1, ncol=2) #loc=3, bbox_to_anchor=(0.,0.97,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

# plot stripes for clarity
for idx,SSC in enumerate(SSCs):
    if idx%2==0:
        ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='lightgrey', lw=0, alpha=0.5, zorder=0)
    else:
        ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='darkgrey', lw=0, alpha=0.5, zorder=0)

ax.set_xlim(0.5, len(SSCs)+0.5)
ax.set_ylim(0,300)
ax.set_xticks(np.arange(1,len(SSCs)+1))
ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
ax.yaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(50))
ax.set_axisbelow(True)
ax.grid(axis='y', which='both')
ax.set_ylabel(r'T$_\mathrm{rot}$ [K]', fontsize=12)
ax.set_xlabel('SSC')
fig.tight_layout()

savepath = escape_fname(os.path.join(plotdir, '14.paper', 'compare_temperatures.pdf'))
os.system('mkdir -p '+os.path.dirname(savepath))
fig.savefig(savepath, dpi=300, bbox_inches='tight')



###################################################################################################
# figure 5: age comparison
###################################################################################################

# SSC ages according to Rico-Villaz+19
###################################################################################################

SSC_ages = {
'1':  {'Lp*/L*': 1.14},
'2':  {'Lp*/L*': 2.22},
'3':  {'Lp*/L*': 3.51},
'4':  {'Lp*/L*': 0.11},
'5':  {'Lp*/L*': 0.02},
'6':  {'Lp*/L*': np.nan},
'7':  {'Lp*/L*': np.nan},
'8':  {'Lp*/L*': 0.56},
'9':  {'Lp*/L*': 0.01},
'10': {'Lp*/L*': 0.03},
'11': {'Lp*/L*': 0.01},
'12': {'Lp*/L*': 0.01},
'13': {'Lp*/L*': 1.52},
'14': {'Lp*/L*': 0.32}
}
for ssc,x in SSC_ages.items():
    if x['Lp*/L*'] > 0.05:
        x['age'] = 1e5 /(1+x['Lp*/L*']) *u.yr
    else:
        x['age'] = 1.0001e5*u.yr

# plot
###################################################################################################

def compare_SSC_age_tracers(type):

    def get_data(type):
        dat = []
        for SSC in SSCs:
            age = SSC_ages[SSC['num']]['age']/1.e5
            dat.append([SSC,
                        age.value,
                        ratios['HCN/H13CN'][type][SSC['num']]['median'],
                        [ratios['HCN/H13CN'][type][SSC['num']]['84th']-ratios['HCN/H13CN'][type][SSC['num']]['median'], ratios['HCN/H13CN'][type][SSC['num']]['median']-ratios['HCN/H13CN'][type][SSC['num']]['16th']],
                        ratios['HCN/HC15N'][type][SSC['num']]['median'],
                        [ratios['HCN/HC15N'][type][SSC['num']]['84th']-ratios['HCN/HC15N'][type][SSC['num']]['median'], ratios['HCN/HC15N'][type][SSC['num']]['median']-ratios['HCN/HC15N'][type][SSC['num']]['16th']],
                        ratios['HCN/HC3N'][type][SSC['num']]['median'],
                        [ratios['HCN/HC3N'][type][SSC['num']]['84th']-ratios['HCN/HC3N'][type][SSC['num']]['median'], ratios['HCN/HC3N'][type][SSC['num']]['median']-ratios['HCN/HC3N'][type][SSC['num']]['16th']]
                       ])
        return dat

    def plot_age_ratios(type):
        for SSC, age, hcn_h13cn, hcn_h13cn_err, hcn_hc15n, hcn_hc15n_err, hcn_hc3n, hcn_hc3n_err in get_data(type):
            if age==1.0001:
                marker = '>'
            else:
                marker = 'o'
            if not np.isnan(hcn_h13cn):
                ax1.errorbar(age-0.01, hcn_h13cn, yerr=[hcn_h13cn_err],
                             marker=marker, ms=3, lw=0, color='r',
                             elinewidth=1, ecolor='r'
                             )
            if not np.isnan(hcn_hc15n):
                ax1.errorbar(age+0.00, hcn_hc15n, yerr=[hcn_hc15n_err],
                             marker=marker, ms=3, lw=0, color='g',
                             elinewidth=1, ecolor='g'
                             )
            if not np.isnan(hcn_hc3n):
                ax2.errorbar(age+0.01, hcn_hc3n, yerr=[hcn_hc3n_err],
                             marker=marker, ms=3, lw=0, color='b',
                             elinewidth=1, ecolor='b'
                             )

    def label_SSCs(type):
        for SSC, age, hcn_h13cn, hcn_h13cn_err, hcn_hc15n, hcn_hc15n_err, hcn_hc3n, hcn_hc3n_err in get_data(type):
            ax1.text(age, 22.5, SSC['num'], color='k', fontsize=12, ha='center', va='center')

    def fit_correlations(type):
        from scipy.optimize import curve_fit

        def linear(x,a,c):
            return a*x+c

        def filter_nan_old(age,val,err):
            fage,fval,ferr = [],[],[]
            for a,v,e in zip(age,val,err):
                if not a==1.0001:
                    if not np.isnan(v):
                        fage.append(a)
                        fval.append(v)
                        ferr.append(e)
            return fage,fval,ferr

        ages          = [x[1] for x in get_data(type)]
        ratios_13     = [x[2] for x in get_data(type)]
        ratios_13_err = [x[3] for x in get_data(type)]
        ratios_15     = [x[4] for x in get_data(type)]
        ratios_15_err = [x[5] for x in get_data(type)]
        ratios_c3     = [x[6] for x in get_data(type)]
        ratios_c3_err = [x[7] for x in get_data(type)]

        age13,r13,r13err = filter_nan_old(ages,ratios_13,ratios_13_err)
        age15,r15,r15err = filter_nan_old(ages,ratios_15,ratios_15_err)
        agec3,rc3,rc3err = filter_nan_old(ages,ratios_c3,ratios_c3_err)

        # x = np.array([0.075,0.475])
        x = np.array([min(ages),max(ages)])

        # HCN/H13CN
        coeff, covar = curve_fit(linear, age13, r13, sigma=np.mean(r13err,axis=1))
        ax1.plot(x, linear(x, *coeff), 'r', lw=1, ls='--', alpha=0.5)

        # HCN/HC15N
        coeff, covar = curve_fit(linear, age15, r15, sigma=np.mean(r15err,axis=1))
        ax1.plot(x, linear(x, *coeff), 'g', lw=1, ls='--', alpha=0.5)

        # HCN/HC3N
        coeff, covar = curve_fit(linear, agec3, rc3, sigma=np.mean(rc3err,axis=1))
        ax2.plot(x, linear(x, *coeff), 'b', lw=1, ls='--', alpha=0.5)


    fig,ax1 = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    ax2 = ax1.twinx()

    plot_age_ratios(type)
    label_SSCs(type)
    fit_correlations(type)

    ax1.set_xlim(0.,1.2)
    ax1.set_ylim(0.,25.)
    # ax2.set_ylim(0.,10.)
    ax.set_axisbelow(True)
    ax.grid(axis='y', which='both')
    ax1.set_xlabel(r'SSC age [$10^5$\,yr] (Rico-Villaz+19)', fontsize=12)
    ax1.set_ylabel(r'HCN/H$^{13}$CN and HCN/HC$^{15}$N', fontsize=12)
    ax2.set_ylabel(r'HCN/HC$_3$N', fontsize=12)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '14.paper', 'SSC_age_comparison.'+type+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_SSC_age_tracers('integrated intensity')
compare_SSC_age_tracers('column density')


###################################################################################################
#
###################################################################################################
