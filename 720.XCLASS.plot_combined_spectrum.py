###################################################################################################
# plot fit results
###################################################################################################

def plot_combined_spectrum(SSC, band):
    """
    run the functions that build up the figure
    """

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

    def get_model(SSC, band):
        with open(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'combined_model.spectrum.pickle')), 'rb') as f:
            m = pickle.load(f, encoding="latin1")
        frequency = (m[:,0]*u.MHz).to(u.GHz)
        model     = m[:,1]*u.K
        return frequency.value,model.value

    def set_up_figure(SSC, band):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_spectrum(ax, frequency, spectrum):
        ax.plot(frequency, spectrum, lw=1, ls='-', color='k', zorder=3)
        ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5, zorder=2)

    def plot_fitted_spectrum(ax, frequency, model):
        ax.plot(frequency, model, lw=1, ls='-', color='r', zorder=5)
        # ax.fill_between(frequency, model, [0. for f in frequency], color='r', alpha=0.5, zorder=4)

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
        savepath = escape_fname(os.path.join(plotdir, '03.XCLASS_fit', 'combined_spectra', 'SSC_'+str(SSC['no'])+'.'+band+'.combined_spectrum.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    frequency, spectrum  = get_spectrum(SSC, band)
    mfrequency, model    = get_model(SSC, band)
    fig,ax = set_up_figure(SSC, band)
    plot_spectrum(ax, frequency, spectrum)
    plot_fitted_spectrum(ax, mfrequency, model)
    label_lines(ax, spectrum, band)
    format_figure(ax, frequency, spectrum, band)
    save_figure(fig, band)


def plot_combined_variation(nums, SSC, band, rms):
    """
    run the functions that build up the figure
    """

    def get_spectra(nums, SSC, band, rms):
        spectrum  = spectra[str(SSC['no'])][band]
        frequency = spectrum['frequency'].to(u.GHz)
        intensity = spectrum['spectrum'].to(u.K)
        # shift spectrum to rest frequency
        velshift  = SSC['velshift']
        frequency = [(-vsys-velshift).to(u.GHz, equivalencies=u.doppler_optical(f)).value for f in frequency]*u.GHz
        # remove NaNs
        frequency, intensity = crossmatch(frequency.to(u.GHz).value, intensity.to(u.K).value)
        # add noise
        intensities = []
        for num in nums:
            if not num==0:
                randstate = np.random.RandomState(num)
                noise =  np.random.normal(loc=0., scale=rms.to(u.K).value, size=len(frequency))
                int_noise = intensity+noise
                intensities.append(int_noise)
            else:
                intensities.append(intensity)
        # get percentiles
        d16,dmed,d84 = np.percentile(np.array(intensities), (16,50,84), axis=0)
        return frequency,d16,dmed,d84

    def get_models(nums, SSC, band):
        with open(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_0','combined_model.spectrum.pickle')), 'rb') as f:
            m = pickle.load(f, encoding="latin1")
        frequency = (m[:,0]*u.MHz).to(u.GHz).value

        models = []
        for num in nums:
            with open(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),'model_spectrum','run_'+str(num),'combined_model.spectrum.pickle')), 'rb') as f:
                m = pickle.load(f, encoding="latin1")
            model     = (m[:,1]*u.K).value
            models.append(model)
        m16,mmed,m84 = np.percentile(np.array(models), (16,50,84), axis=0)
        return frequency,m16,mmed,m84

    def set_up_figure(SSC, band):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_spectra(ax, frequency, d16, dmed, d84):
        ax.plot(frequency, dmed, lw=1, ls='-', color='k', zorder=3)
        ax.fill_between(frequency, d16, d84, color='k', alpha=0.5, zorder=2)

    def plot_fitted_spectra(ax, frequency, m16, mmed, m84):
        ax.plot(frequency, mmed, lw=1, ls='-', color='r', zorder=5)
        ax.fill_between(frequency, m16, m84, color='r', alpha=0.5, zorder=4)

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
        savepath = escape_fname(os.path.join(plotdir, '03.XCLASS_fit', 'combined_spectra', 'SSC_'+str(SSC['no'])+'.'+band+'.combined_spectra.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    frequency, d16,dmed,d84  = get_spectra(nums, SSC, band, rms)
    mfrequency, m16,mmed,m84 = get_models(nums, SSC, band)
    fig,ax = set_up_figure(SSC, band)
    plot_spectra(ax, frequency, d16,dmed,d84)
    plot_fitted_spectra(ax, mfrequency, m16,mmed,m84)
    label_lines(ax, dmed, band)
    format_figure(ax, frequency, dmed, band)
    save_figure(fig, band)


###################################################################################################
#
###################################################################################################
