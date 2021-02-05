#######################
# GAS IN SSCS: XCLASS #
#######################

# Prepare the data in python3 in a way that keeps the XCLASS code simple and manageable.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

Xfiles      = fnunpickle(os.path.join(Xfinaldir, 'Xfiles.pickle'))
temperature = fnunpickle(os.path.join(Xfinaldir, 'temperature.pickle'))
nums        = fnunpickle(os.path.join(Xfinaldir, 'nums.pickle'))

execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))


###################################################################################################
# overview figure
###################################################################################################

def spectrum_overview():

    def get_spectrum(SSC, band):
        spectrum  = spectra[str(SSC['no'])][band]
        frequency = spectrum['frequency'].to(u.GHz)
        intensity = spectrum['spectrum'].to(u.K)
        # shift spectrum to rest frequency
        velshift  = SSC['velshift']
        frequency = [(-vsys-velshift).to(u.GHz, equivalencies=u.doppler_optical(f)).value for f in frequency]*u.GHz
        # remove NaNs
        frequency, intensity = crossmatch(frequency.to(u.GHz).value, intensity.to(u.K).value)
        return frequency, intensity

    def get_models(nums, SSC, band):
        with open(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_0','combined_model.spectrum.pickle')), 'rb') as f:
            m = pickle.load(f, encoding="latin1")
        frequency = (m[:,0]*u.MHz).to(u.GHz).value

        models = []
        # files = glob.glob(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_*','combined_model.spectrum.pickle')))
        files = [ Xfiles[SSC['num']][str(num)]['model']['spectrum'] for num in nums ]
        for f in files:
            with open(f, 'rb') as f:
                m = pickle.load(f, encoding="latin1")
                model = (m[:,1]*u.K).value
                models.append(model)
        m16,mmed,m84 = np.percentile(np.array(models), (16,50,84), axis=0)
        return frequency,m16,mmed,m84

    def set_up_figure():
        fig,axes = plt.subplots(nrows=len(SSCs), ncols=2, squeeze=True, sharex='col', sharey='row', figsize=A4_inches)
        fig.subplots_adjust(hspace=0, wspace=0, top=0.80, bottom=0.04, left=0.04, right=0.93)
        return fig, axes

    def plot_spectrum(ax, frequency, intensity):
        ax.plot(frequency, intensity, lw=1., ls='-', color='k')
        ax.fill_between(frequency, intensity, [0. for f in frequency], color='k', alpha=0.5, zorder=2)

    def plot_median_model(ax, frequency, m16, mmed, m84):
        ax.plot(frequency, mmed, lw=1., ls='-', color='r')
        # ax.fill_between(frequency, m16, m84, color='r', alpha=0.5, zorder=2)

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
                ax.text(xloc, 155, line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')

    def format_figure(axes):
        axes[-1][0].tick_params(labelbottom=True, labelleft=True, labelright=False)
        axes[-1][1].tick_params(labelbottom=True, labelleft=False, labelright=True)
        axes[-1][1].yaxis.set_label_position('right')
        for i_band,band in enumerate(['LSB','USB']):
            axes[-1][i_band].xaxis.set_visible(True)
            axes[-1][i_band].set_xlabel(r'$\nu_\mathrm{rest}$ [GHz]', fontsize=10)
            axes[-1][i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=10)

    def save_figure(fig):
        savepath = escape_fname(os.path.join(plotdir, '04.XCLASS_final', 'all_spectra.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    # figure making happens here:
    fig, axes = set_up_figure()
    dlines = get_detected_lines()
    for i_SSC,SSC in tqdm(enumerate(SSCs)):
        for i_band,band in enumerate(['LSB','USB']):
            ax = axes[i_SSC,i_band]

            frequency, intensity = get_spectrum(SSC, band)
            mfrequency, m16, mmed, m84 = get_models(nums, SSC, band)

            plot_spectrum(ax, frequency, intensity)
            plot_median_model(ax, mfrequency, m16, mmed, m84)

            add_SSC_overlay(ax)
            label_lines_line(ax, dlines)
            format_panel(ax)

    label_lines_text(axes, dlines)
    format_figure(axes)
    save_figure(fig)


spectrum_overview()


###################################################################################################
#
###################################################################################################
