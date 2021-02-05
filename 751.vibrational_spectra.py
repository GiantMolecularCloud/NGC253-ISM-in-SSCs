#######################
# GAS IN SSCS: XCLASS #
#######################

# Fitting of ro-vibational lines
# Keep other lines fixed and fit only for ro-vibrational lines.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))

nums      = fnunpickle(os.path.join(vibdir, 'nums.pickle'))
vibfiles = fnunpickle(os.path.join(vibdir, 'vibfiles.pickle'))
vibrational_species = fnunpickle(os.path.join(vibdir, 'vibrational_species.pickle'))


###################################################################################################
# create model spectra
###################################################################################################

def create_model_molfit(num, SSC):

    # fit molfit files & model molfit files
    fit   = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'fit','results','molecules__LM__call_1.out.molfit'))
    model = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','molecules.molfit'))

    # read molfit file
    with open(fit, 'r') as file :
        filedata = file.read()

    # replace the fit switch
    filedata = filedata.replace(' y ', ' n ')

    # write the model molfit file
    os.system('mkdir -p '+os.path.dirname(model))
    with open(model, 'w') as file:
        file.write(filedata)

    vibfiles[SSC['num']][str(num)]['model'] = {'molfit': model}


def model_spectrum_XCLASS(num, SSC):
    """
    A script to run the create_model_spectrum function within CASA.
    """

    modelpath  = vibfiles[SSC['num']][str(num)]['model']['molfit']
    savepath   = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','XCLASS_model.py'))
    resultpath = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','results'))
    specpath   = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','spectrum.pickle'))
    energypath = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','energy.pickle'))
    taupath    = escape_fname(os.path.join(vibdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','opacity.pickle'))

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
    vibfiles[SSC['num']][str(num)]['model']['file']     = savepath
    vibfiles[SSC['num']][str(num)]['model']['spectrum'] = specpath
    vibfiles[SSC['num']][str(num)]['model']['energy']   = energypath
    vibfiles[SSC['num']][str(num)]['model']['tau']      = taupath


for num in tqdm(nums):
    for SSC in SSCs:
        create_model_molfit(num, SSC)
        model_spectrum_XCLASS(num, SSC)


###################################################################################################
# execute in CASA
###################################################################################################

def CASA_command(SSC):
    XCL_files = [vibfiles[SSC['num']][str(num)]['model']['file'] for num in nums]
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

fnpickle(vibfiles, os.path.join(vibdir, 'vibfiles.pickle'))


###################################################################################################
# plot model/data variation
###################################################################################################

def plot_combined_variation(nums, SSC, band):
    """
    run the functions that build up the figure
    """

    def get_spectra(nums, SSC, band):
        spectra = [ np.genfromtxt(vibfiles[SSC['num']][str(num)]['dat'][band], dtype=None) for num in nums ]
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,:,0], 50, axis=0)/1000
        p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    def get_models(nums, SSC, band):
        models = [ vibfiles[SSC['num']][str(num)]['model']['spectrum'] for num in nums ]
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
        models = [ vibfiles[SSC['num']][str(num)]['model']['spectrum'] for num in nums ]
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
        detected_species = fnunpickle(os.path.join(XCLASSdir, 'detected_species.pickle'))
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
        savepath = escape_fname(os.path.join(plotdir, '06.vibrations', 'spectra', 'SSC_'+str(SSC['no'])+'.'+band+'.spectrum.pdf'))
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
        plot_combined_variation(nums, SSC, band)


###################################################################################################
#
###################################################################################################
