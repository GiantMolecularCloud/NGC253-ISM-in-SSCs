#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Plot fitted spectra.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
bandfit = fnunpickle(os.path.join(mandir, 'updated_band_fit_parameters.pickle'))


###################################################################################################
# fitting function definitions
###################################################################################################

def filter_lines(lines):
    filtered_lines = []
    for line in lines:
        if line['confidence']!='0' and line['confidence']!='-':
            filtered_lines.append(line)
    return filtered_lines

def lines_in_band(lines, band):
    blines = []
    for line in lines:
        if band=='LSB' and line['restfreq']<350*u.GHz:
            blines.append(line)
        elif band=='USB' and line['restfreq']>350*u.GHz:
            blines.append(line)
    return blines

def get_spectrum(SSC, band):
    frequency = spectra[str(SSC['no'])][band]['frequency'].value
    spectrum  = spectra[str(SSC['no'])][band]['spectrum'].value
    return crossmatch(frequency, spectrum)

def guess_from_dict(fit_dict, SSC, band, line_id):
    amp   = fit_dict[str(SSC['no'])][band][line_id]['amp'][0]
    sigma = fit_dict[str(SSC['no'])][band][line_id]['sigma'][0]
    mu    = fit_dict[str(SSC['no'])][band][line_id]['mu'][0]
    guess = np.array([amp, sigma, mu])
    return guess

def multi_gauss(x, *params):
    """
    Sum of arbitrary many gaussians. Parameter order is: x, n*(amplitude, mu, sigma)
    """
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp   = params[i]
        sigma = params[i+1]
        mu    = params[i+2]
        y += amp*np.exp(-((x-mu)/sigma)**2)
    return y

def plot_spectrum_fit(SSC, band, blines, frequency, spectrum, fit_params, number_lines=False):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    # plot spectrum
    ax.plot(frequency, spectrum, ls='-', color='k', zorder=3)
    ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5, zorder=2)

    # plot fitted spectrum
    fit = multi_gauss(frequency, *fit_params)
    ax.plot(frequency, fit, ls='-', color='r', zorder=4)

    # label expected lines
    for idx,line in enumerate(sorted(blines, key=lambda k: k['restfreq'])):
        obsfreq = line['obsfreq'].to(u.GHz).value
        if (obsfreq>frequency[0] and obsfreq<frequency[-1]):
            obsfreq = (vsys+SSC['velshift']).to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value
            xlim = [frequency[0], frequency[-1]]
            xloc = xlim[0] +((idx+0.5)/len(blines))*(xlim[1]-xlim[0])
            if number_lines==True:
                numl = '['+str(blines.index(line))+']'
            else:
                numl = ''
            ax.axvline(x=obsfreq, ymin=0, ymax=1, color='dimgrey', ls='--', lw=0.5, zorder=1)
            ax.plot([obsfreq,xloc], [1.05*np.nanmax(spectrum), 1.05*1.05*np.nanmax(spectrum)], color='dimgrey', ls='--', lw=0.5, zorder=1, clip_on=False)
            ax.text(xloc, 1.06*1.05*np.nanmax(spectrum), numl+line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')


    ax.set_xlim([frequency[0], frequency[-1]])
    ax.set_ylim(-0.05*np.nanmax(spectrum), 1.05*np.nanmax(spectrum))
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_axisbelow(True)
    ax.grid()
    ax.set_xlabel(r'$\nu$ [GHz]', fontsize=12)
    ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(plotdir, 'SSC_'+str(SSC['no'])+'.'+band+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# re-plot fitted band spectra
###################################################################################################

plt.ioff()
for band in ['LSB','USB']:
    for SSC in tqdm(SSCs):
        blines = lines_in_band(lines, band)
        frequency, spectrum = get_spectrum(SSC, band)
        fit_params = flatten( [guess_from_dict(bandfit, SSC, band, line) for line in bandfit[str(SSC['no'])][band].keys()] )
        plot_spectrum_fit(SSC, band, blines, frequency, spectrum, fit_params, number_lines=True)
plt.clf()
plt.ion()


###################################################################################################
# single figure - all fitted spectra
###################################################################################################

fig,axes = plt.subplots(nrows=len(SSCs), ncols=2, squeeze=True, sharex='col', sharey='row', figsize=A4_inches)
fig.subplots_adjust(hspace=0, wspace=0, top=0.80, bottom=0.04, left=0.04, right=0.93)

fitted_lines = [get_line(x) for x in list(set([l for ll in [list(bandfit[str(SSC['no'])][b].keys()) for SSC in SSCs for b in ['LSB','USB']] for l in ll]))]

for i_SSC,SSC in tqdm(enumerate(SSCs)):
    for i_band,band in enumerate(['LSB','USB']):
        ax = axes[i_SSC,i_band]

        # spectrum
        frequency, spectrum = get_spectrum(SSC, band)
        ax.plot(frequency, spectrum, lw=1., ls='-', color='k')
        ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5, zorder=2)

        # plot fitted spectrum
        fit_params = flatten( [guess_from_dict(bandfit, SSC, band, line) for line in bandfit[str(SSC['no'])][band].keys()] )
        fit = multi_gauss(frequency, *fit_params)
        ax.plot(frequency, fit, lw=1., ls='-', color='r', zorder=4)

        # SSC overlay
        ax.text(0.5, 0.9, 'SSC '+str(SSC['no']), color='k', fontsize=10, transform=ax.transAxes, ha='center', va='top', bbox=props)

        # label expected lines
        for line in fitted_lines:
            obsfreq = line['obsfreq'].to(u.GHz).value
            obsfreq = (vsys+SSC['velshift']).to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value
            if (obsfreq<350 and band=='LSB') or (obsfreq>350 and band=='USB'):
                ax.axvline(x=obsfreq, ymin=0, ymax=1, color='dimgrey', lw=0.5, ls='--', zorder=1)


        # formatting
        if band=='LSB':
            ax.set_xlim([342.138, 345.835])
        elif band=='USB':
            ax.set_xlim([354.010, 357.695])
        ax.set_ylim(-15,120)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, top=True, bottom=True, left=True, right=True)
        ax.set_axisbelow(True)
        ax.grid(color='grey', alpha=0.5)

# label expected lines
for line in fitted_lines:
    obsfreq = vsys.to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value
    if obsfreq<350:
        ax = axes[0][0]
    elif obsfreq>350:
        ax = axes[0][1]
    ax.text(obsfreq, 125, line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')

axes[-1][0].tick_params(labelbottom=True, labelleft=True, labelright=False)
axes[-1][1].tick_params(labelbottom=True, labelleft=False, labelright=True)
axes[-1][1].yaxis.set_label_position('right')
for i_band,band in enumerate(['LSB','USB']):
    axes[-1][i_band].xaxis.set_visible(True)
    axes[-1][i_band].set_xlabel(r'$\nu$ [GHz]', fontsize=10)
    axes[-1][i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=10)

fig.savefig(os.path.join(plotdir, 'fitted_spectra.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
