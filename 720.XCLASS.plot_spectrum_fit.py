###################################################################################################
# plot fit results
###################################################################################################

def plot_spectrum_fit(SSC, specie, band):
    """
    run the functions that build up the figure
    """

    def get_spectrum(SSC, band):
        spec = np.genfromtxt(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_0',band+'.dat')))
        return spec[:,0]/1000, spec[:,1]

    def get_bestfit(SSC, specie, band):
        spec = np.genfromtxt(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_0','results',band+'.LM__call_1.out.dat')))
        return spec[:,0]/1000, spec[:,1]

    def get_fit_range(SSC, specie, band):
        # load fit files
        fit_spec_files = glob.glob(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_*','results',band+'.LM__call_1.out.dat')))
        fit_spec_files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
        fit_specs = []
        for f in fit_spec_files:
            spec = np.genfromtxt(f)
            fit_specs.append(spec[:,1])
        # reorder fits in array
        fit_specs = np.transpose(fit_specs)
        # get percentiles per channel
        lower,median,upper = np.percentile(fit_specs, (16,50,84), axis=1)
        return lower, median, upper

    def set_up_figure(SSC, specie, band):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+specie.replace('#','\#')+' '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_spectrum(ax, frequency, spectrum):
        ax.plot(frequency, spectrum, ls='-', color='k', zorder=3)
        ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5, zorder=2)

    def plot_fitted_spectrum(ax, frequency, lower, best, upper):
        ax.plot(frequency, best, ls='-', color='r', zorder=5)
        ax.fill_between(frequency, lower, upper, color='r', alpha=0.5, zorder=4)

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

    def velocity_indicator(ax, frequency, spectrum):
        x = np.nanmedian(frequency)
        y = 0.95*np.nanmax(spectrum)
        ax.plot([x, x+0.05836398],[y,y], c='k', ls='-', lw='3', zorder=11)
        ax.text(x, y, r'-                   $\sim 50$\,km\,s$^{-1}$', color='k', ha='left', va='center', weight='bold', bbox=props, zorder=10)

    def format_figure(ax, frequency, spectrum):
        if band=='LSB':
            xlim = [342.4, 346.2]
        elif band=='USB':
            xlim = [354.3, 358.1]
        ax.set_xlim(xlim)
        ax.set_ylim(-0.05*np.nanmax(spectrum), 1.05*np.nanmax(spectrum))
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='grey')
        ax.set_xlabel(r'$\nu_\mathrm{rest}$ [GHz]', fontsize=12)
        ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
        fig.set_tight_layout(True)

    def save_figure(fig, specie, band):
        savepath = escape_fname(os.path.join(plotdir, '03.XCLASS_fit', 'SSC_'+str(SSC['no']), specie+'.'+band+'.spectrum.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    frequency, spectrum  = get_spectrum(SSC, band)
    fit_freq, fit_spec   = get_bestfit(SSC, specie, band)
    lower, median, upper = get_fit_range(SSC, specie, band)
    fig,ax = set_up_figure(SSC, specie, band)
    plot_spectrum(ax, frequency, spectrum)
    plot_fitted_spectrum(ax, fit_freq, lower, fit_spec, upper)
    label_lines(ax, spectrum, band)
    velocity_indicator(ax, frequency, spectrum)
    format_figure(ax, frequency, spectrum)
    save_figure(fig, specie, band)


###################################################################################################
#
###################################################################################################
