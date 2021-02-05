#######################
# GAS IN SSCS: XCLASS #
#######################

# address referee report: derive intensities from XCLASS fits

# refit a couple of CO lines because the results were inconsistent


###################################################################################################
# create model spectra
###################################################################################################

def create_model_molfit(num, SSCnum):

    # fit molfit files & model molfit files
    fit   = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'fit','results','molecules__LM__call_1.out.molfit'))
    model = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','molecules.molfit'))

    # read molfit file
    with open(fit, 'r') as file :
        filedata = file.read()

    # replace the fit switch
    filedata = filedata.replace(' y ', ' n ')

    # write the model molfit file
    os.system('mkdir -p '+os.path.dirname(model))
    with open(model, 'w') as file:
        file.write(filedata)

    refitfiles[SSCnum][str(num)]['model'] = {'molfit': model}


def model_spectrum_XCLASS(num, SSCnum):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = refitfiles[SSCnum][str(num)]['model']['molfit']
    savepath   = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','XCLASS_model.py'))
    resultpath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','results'))
    specpath   = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','spectrum.pickle'))
    energypath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','energy.pickle'))
    taupath    = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','opacity.pickle'))

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
    refitfiles[SSCnum][str(num)]['model']['file']     = savepath
    refitfiles[SSCnum][str(num)]['model']['spectrum'] = specpath
    refitfiles[SSCnum][str(num)]['model']['energy']   = energypath
    refitfiles[SSCnum][str(num)]['model']['tau']      = taupath


for num in tqdm(nums):
    for SSCnum in refitfiles.keys():
        create_model_molfit(num, SSCnum)
        model_spectrum_XCLASS(num, SSCnum)


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSCnum):
    XCL_files = [refitfiles[SSCnum][str(num)]['model']['file'] for num in nums]
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

commands = [CASA_command(SSCnum) for SSCnum in refitfiles.keys()]
pool = Pool(len(list(refitfiles.keys())))
pool.map(run_in_casa, commands)
pool.close()
pool.join()

fnpickle(refitfiles, escape_fname(join(refitdir, refitspecies, 'refitfiles.pickle')))


###################################################################################################
# plot model/data variation
###################################################################################################

def plot_combined_variation(nums, SSCnum, band):
    """
    run the functions that build up the figure
    """

    def get_spectra(nums, SSCnum,band):
        spectra = [ np.genfromtxt(refitfiles[SSCnum][str(num)]['dat'][band], dtype=None) for num in nums ]
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,:,0], 50, axis=0)/1000
        p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def get_models(nums, SSCnum,band):
        models = [ refitfiles[SSCnum][str(num)]['model']['spectrum'] for num in nums ]
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

    def get_correct_models(nums, SSCnum, band):
        spectra = []
        for num in nums:
            mfiles = glob.glob( refitfiles[SSCnum][str(num)]['model']['spectrum'].replace('spectrum.pickle','results/')+'intensity_*.dat' )
            mdats  = []
            for mf in mfiles:
                mdats.append( np.genfromtxt(mf, dtype=None, skip_header=4) )
            mdats = np.array(mdats)
            spectra.append( np.sum(mdats, axis=0) )
        spectra = np.array(spectra)
        frequency      = mdats[0][:,0]/1000
        p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def set_up_figure(SSCnum,band):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+SSCnum+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_spectra(ax, frequency, d16, dmed, d84):
        ax.plot(frequency, dmed, lw=1, ls='-', color='k', zorder=3)
        ax.fill_between(frequency, d16, d84, color='k', alpha=0.4, zorder=2)

    def plot_fitted_spectra(ax, frequency, m16, mmed, m84):
        ax.plot(frequency, mmed, lw=1, ls='-', color='r', zorder=5)
        ax.fill_between(frequency, m16, m84, color='r', alpha=0.4, zorder=4)

    def get_detected_lines(band=None):
        # get detected species
        all_species = []
        for SSCnum in refitfiles.keys():
            for specie in detected_species[SSCnum]:
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
        savepath = escape_fname(join(plotdir, '11.refit', refitspecies, 'spectra', 'SSC_'+SSCnum+'.'+band+'.spectrum.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    frequency, d16,dmed,d84  = get_spectra(nums, SSCnum, band)
    mfrequency, m16,mmed,m84 = get_models(nums, SSCnum, band)
    # mfrequency, m16,mmed,m84 = get_correct_models(nums, SSCnum, band)
    fig,ax = set_up_figure(SSCnum,band)
    plot_spectra(ax, frequency, d16,dmed,d84)
    plot_fitted_spectra(ax, mfrequency, m16,mmed,m84)
    label_lines(ax, dmed, band)
    format_figure(ax, frequency, dmed, band)
    save_figure(fig, band)


for SSCnum in tqdm(refitfiles.keys()):
    for band in ['LSB','USB']:
        plot_combined_variation(nums, SSCnum, band)


###################################################################################################
#
###################################################################################################


###################################################################################################
# check if components are fitted correctly
###################################################################################################

for SSCnum in tqdm(refitfiles.keys()):

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+SSCnum, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    # get spectra
    spectra = [ np.genfromtxt(refitfiles[SSCnum][str(num)]['dat'][refitband], dtype=None) for num in nums ]
    spectra = np.array(spectra)
    spec_frequency = np.percentile(spectra[:,:,0], 50, axis=0)/1000
    spec_p16,spec_median,spec_p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)

    ax.plot(spec_frequency, spec_median, lw=1, ls='-', color='k', zorder=3, label='data')
    ax.fill_between(spec_frequency, spec_p16, spec_p84, color='k', alpha=0.4, zorder=2)

    # get model spectra
    import matplotlib.colors as colors
    norm = colors.Normalize(vmin=0, vmax=len(nums)-1)
    for num in nums:
        color = cm.get_cmap('rainbow_r')(norm(num))
        with open(refitfiles[SSCnum][str(num)]['model']['spectrum'], 'rb') as f:
            m = pickle.load(f, encoding="latin1")
            mod_freq = (m[:,0]*u.MHz).to(u.GHz)
            mod_intensity = (m[:,1]*u.K)
            ax.plot(mod_freq, mod_intensity, lw=1, ls='--', color=color, zorder=3, label=str(num))

    fig.legend(loc='top right')
    # ax.set_xlim([345.4, 346.1])
    ax.set_ylim(-0.05*np.nanmax(spec_median), 1.05*np.nanmax(spec_median))
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
    savepath = escape_fname(join(plotdir, '11.refit', refitspecies, 'components', 'SSC_'+SSCnum+'.components.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
