#######################
# GAS IN SSCS: XCLASS #
#######################

# Merge free fitting, temperature line fitting and ro-vib line fitting.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
detected_species = fnunpickle(os.path.join(XCLASSdir, 'detected_species.pickle'))

data_XCLASS = fnunpickle(os.path.join(resultsdir, 'data_XCLASS.pickle'))
data_Gauss  = fnunpickle(os.path.join(resultsdir, 'data_Gauss.pickle'))

model_files = {SSC['num']: {str(num): {} for num in np.arange(0,101)} for SSC in SSCs}


###################################################################################################
# create model molfit file
###################################################################################################

def create_model_molfit(num, SSC):

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset"""

    for specie,data in data_XCLASS[SSC['num']].items():
        if '#1' in specie:
            spx = specie
        else:
            spx = specie+';'

        ncomp = len(data['velocity']['median'])
        molfit += f"""
{spx}   {ncomp}"""

        for c in np.arange(ncomp):
            T = data['temperature']['all'][c][num]
            v = data['velocity']['all'][c][num]
            w = data['linewidth']['all'][c][num]
            N = data['column density']['all'][c][num]

            molfit += f"""
n  0.01  10  1   n  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}  {'{:.2E}'.format(T)}   n  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}  {'{:.2E}'.format(N)}   n  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}  {'{:.1f}'.format(w)}   n  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}  {'{:.1f}'.format(v)}   c"""

    # save to disk
    savepath = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)

    # keep savepath
    model_files[SSC['num']][str(num)]['molfit'] = savepath


for num in tqdm(np.arange(0,101)):
    for SSC in SSCs:
        create_model_molfit(num, SSC)


###################################################################################################
# create model script
###################################################################################################

def model_spectrum_XCLASS(num, SSC):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = model_files[SSC['num']][str(num)]['molfit']
    savepath   = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'XCLASS_model.py'))
    resultpath = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'results'))
    specpath   = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'spectrum.pickle'))
    energypath = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'energy.pickle'))
    taupath    = escape_fname(os.path.join(resultsdir,'models','SSC_'+str(SSC['no']),'run_'+str(num),'opacity.pickle'))
    model_files[SSC['num']][str(num)]['scriptfile'] = savepath
    model_files[SSC['num']][str(num)]['spectrum']   = specpath
    model_files[SSC['num']][str(num)]['energy']     = energypath
    model_files[SSC['num']][str(num)]['tau']        = taupath

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


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSC):
    XCL_files = [model_files[SSC['num']][str(num)]['scriptfile'] for num in np.arange(0,101)]
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

fnpickle(model_files, os.path.join(resultsdir, 'model_files.pickle'))


###################################################################################################
# plot model/data variation
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
        models = [ model_files[SSC['num']][str(num)]['spectrum'] for num in nums ]
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
        models = [ model_files[SSC['num']][str(num)]['spectrum'] for num in nums ]
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
        ax.set_xlabel(r'$\nu_\mathrm{rest}$ [GHz]', fontsize=12)
        ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
        fig.set_tight_layout(True)

    def save_figure(fig, band):
        savepath = escape_fname(os.path.join(plotdir, '10.results', 'spectra', 'SSC_'+str(SSC['no'])+'.'+band+'.spectrum.pdf'))
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
# plot spectrum overview for paper
###################################################################################################

def plot_spectrum_overview():

    def get_spectrum(SSC, band):
        spectrum = np.genfromtxt(escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_0',band+'.dat')), dtype=None)
        return spectrum[:,0]/1000, spectrum[:,1]

    def get_models(nums, SSC, band):
        models = [ model_files[SSC['num']][str(num)]['spectrum'] for num in nums ]
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
            axes[-1][i_band].set_xlabel(r'$\nu_\mathrm{rest}$ [GHz]', fontsize=10)
            axes[-1][i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=10)

    def save_figure(fig):
        savepath = escape_fname(os.path.join(plotdir, '10.results', 'all_spectra.pdf'))
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
#
###################################################################################################
