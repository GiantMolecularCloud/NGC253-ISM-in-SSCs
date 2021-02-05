#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Fit a complete sideband together to get a more robust estimate of the line position and also
# account for blending.

###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))


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

def amp_estimate(line):
    amp_dict = {'5': 0.5, '4': 0.25, '3': 0.1, '2': 0.05, '1': 0.02, '0': 0.01, '-': 0.01}
    return amp_dict[line['confidence']]

def get_guess(lines, spectrum):
    """
    Estimate initial guesses from the spectrum and information about the spectral lines.
    """
    smax  = np.nanmax(spectrum)
    guess = []
    for line in lines:
        amp   = smax * amp_estimate(line)                                                                    # K
        sigma = 0.05 *np.sqrt(amp_estimate(line))                                                            # GHz
        mu    = (SSC['velshift']+vsys).to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value    # GHz
        guess.append([amp, sigma, mu])
    return flatten(guess)

def get_bound(lines, SSC):
    """
    Estimate the fit bounds from the spectrum and information about the spectral lines.
    """
    lower_bound = []
    upper_bound = []
    for line in lines:
        lb_amp   = 0.
        ub_amp   = 200.
        lb_sigma = 0.
        ub_sigma = 0.05
        lb_mu    = (SSC['velshift']+vsys+5*u.km/u.s).to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value
        ub_mu    = (SSC['velshift']+vsys-5*u.km/u.s).to(u.GHz, equivalencies=u.doppler_optical(line['restfreq'])).value
        lower_bound.append([lb_amp, lb_sigma, lb_mu])
        upper_bound.append([ub_amp, ub_sigma, ub_mu])
    return (flatten(lower_bound), flatten(upper_bound))

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

def fit_spectrum(frequency, spectrum, guess, bound):
    from scipy.optimize import curve_fit
    return curve_fit(multi_gauss, frequency, spectrum, p0=guess, bounds=bound)

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

def ask_repeat(SSC, band):
    return str(input('Refine fit initial guess for SSC '+str(SSC['no'])+': '+band+'? ["no" to stop]    '))

def ask_num_lines(SSC, band):
    return str(input('SSC '+str(SSC['no'])+': '+band+': which lines to plot? ("no" to stop, enter to continue)    '))

def get_new_guess(amp, sigma, mu, const):
    amp_new   = input('{:<6}'.format('amp=')  +'{:6.3f}'.format(amp)  +'{:>17}'.format('[new]  amp = '))
    sigma_new = input('{:<6}'.format('sigma=')+'{:6.3f}'.format(sigma)+'{:>17}'.format('[new]  sigma = '))
    mu_new    = input('{:<6}'.format('mu=')   +'{:6.3f}'.format(mu)   +'{:>17}'.format('[new]  mu = '))
    const_new = input('{:<6}'.format('const=')+'{:6.3f}'.format(const)+'{:>17}'.format('[new]  const = '))
    if amp_new!='':
        amp = float(amp_new)
    if sigma_new!='':
        sigma = float(sigma_new)
    if mu_new!='':
        mu = float(mu_new)
    if const_new!='':
        const = float(const_new)
    return (amp, sigma, mu, const)

def ask_line_indices(lines):
    lfitted = 'Lines fitted: '
    for line in lines:
        lfitted += str(lines.index(line))+', '
    lfitted = lfitted[:-2]
    print(lfitted)
    new_lines = input('Lines to fit: ')
    new_lines = new_lines.split(',')
    new_lines = [lines[int(l.strip())] for l in new_lines]
    return new_lines

def update_params_dict(fit_dict, SSC, band, lines, fit_params, covar):
    fit_dict[str(SSC['no'])][band] = {}
    for i,line in enumerate(lines):
        amp   = fit_params[i*3]
        sigma = fit_params[i*3+1]
        mu    = fit_params[i*3+2]
        amp_err   = covar[i*3][i*3]
        sigma_err = covar[i*3+1][i*3+1]
        mu_err    = covar[i*3+2][i*3+2]
        fit_dict[str(SSC['no'])][band][line['ID']] = {'amp': [amp, amp_err], 'sigma': [sigma, sigma_err], 'mu': [mu, mu_err]}


###################################################################################################
# interactively fit the spectra
###################################################################################################

# # filter lines to include only lines with reasonable confidence
# flines = filter_lines(lines)
#
# # line fits
# manual_bandfit = {str(SSC['no']): {} for SSC in SSCs}
#
# # plot spectrum+fit and check
# for band in ['LSB','USB']:
#     for SSC in tqdm(SSCs):
#         blines = lines_in_band(flines, band)
#         frequency, spectrum = get_spectrum(SSC, band)
#         guess = get_guess(blines, spectrum)
#         bound = get_bound(blines, SSC)
#         try:
#             fit_params, covar = fit_spectrum(frequency, spectrum, guess, bound)
#         except RuntimeError:
#             print("Fit for SSC "+str(SSC['no'])+" in "+band+" failed. Plotting guesses.")
#             fit_params = guess
#         update_params_dict(manual_bandfit, SSC, band, fit_params)
#         plot_spectrum_fit(SSC, band, frequency, spectrum, fit_params)
#
# fnpickle(manual_bandfit, os.path.join(mandir, 'manual_band_fit_parameters.pickle'))


manual_bandfit = {str(SSC['no']): {} for SSC in SSCs}

# plot spectrum+fit and check
for band in ['LSB','USB']:
    for SSC in SSCs:
        blines = lines_in_band(lines, band)
        frequency, spectrum = get_spectrum(SSC, band)
        guess = get_guess(blines, spectrum)
        bound = get_bound(blines, SSC)
        plot_spectrum_fit(SSC, band, blines, frequency, spectrum, guess, number_lines=True)

        redo = ask_num_lines(SSC, band)

        while not (redo=='n' or redo=='no'):
            new_lines = ask_line_indices(blines)
            guess = get_guess(new_lines, spectrum)
            bound = get_bound(new_lines, SSC)
            try:
                fit_params, covar = fit_spectrum(frequency, spectrum, guess, bound)
                update_params_dict(manual_bandfit, SSC, band, new_lines, fit_params, covar)
            except RuntimeError:
                print("Fit for SSC "+str(SSC['no'])+" in "+band+" failed. Plotting guesses.")
                fit_params = guess
            plot_spectrum_fit(SSC, band, new_lines, frequency, spectrum, fit_params, number_lines=True)

            redo = ask_num_lines(SSC, band)

fnpickle(manual_bandfit, os.path.join(mandir, 'manual_band_fit_parameters.pickle'))


###################################################################################################
# refit with updated code
###################################################################################################

line_idx_to_fit = {'LSB': {'1': [7, 9, 3, 8, 2, 4, 0],
                           '2': [11,7,9,3,8,14,15,16,17,2,4,0],
                           '3': [11,7,9,3,13,8,14,16,17,18,1,2,4,0],
                           '4': [7,9,3,8,2,0,4],
                           '5': [7,3,8,2,4,0],
                           '6': [7,2,4,0,21],
                           '7': [7,2,0],
                           '8': [7,9,3,8,1,2,4,0],
                           '9': [7,8,2,4,0],
                           '10': [7,8,14,2,4,0],
                           '11': [11,7,9,8,2,4,0],
                           '12': [11,7,3,8,2,4,0],
                           '13': [11,7,9,3,8,14,1,19,2,4,0],
                           '14': [11,7,9,5,10,3,8,16,17,6,1,19,2,4,20,0,21]},
                   'USB': {'1': [1,6,21,22,7,18,8,2,4,0,24,25,26,27,28,29,30],
                           '2': [1,6,21,23,7,18,8,11,4,0,16,17,24,25,26,27,28,29,30],
                           '3': [1,6,21,18,8,9,19,11,4,0,16,17,24,25,26,27,28,29,30],
                           '4': [1,6,21,7,8,9,19,11,4,0,16,17,24,25,26,27,28,29,30],
                           '5': [1,6,21,7,18,8,9,19,11,4,5,0,24,25,26,27,28,29,30],
                           '6': [1,0],
                           '7': [1,0,29,30],
                           '8': [1,6,21,7,18,8,9,19,11,4,0,16,17,24,25,26,27,28,29,30],
                           '9': [1,6,0],
                           '10': [1,6,21,0,30],
                           '11': [1,6,11,4,0,24,25,26,27,28],
                           '12': [1,6,0],
                           '13': [1,6,21,7,18,8,9,19,10,11,3,4,5,0,16,17,18,24,25,26,27,28,29,30],
                           '14': [1,6,21,7,18,9,19,11,12,13,14,2,3,4,5,0,16,17,32,24,25,26,27,28,29]}}

auto_bandfit = {str(SSC['no']): {} for SSC in SSCs}

# plot spectrum+fit and check
plt.ioff()
for band in ['LSB','USB']:
    for SSC in tqdm(SSCs):
        blines = lines_in_band(lines, band)
        new_lines = [blines[l] for l in line_idx_to_fit[band][str(SSC['no'])]]
        frequency, spectrum = get_spectrum(SSC, band)
        guess = get_guess(new_lines, spectrum)
        bound = get_bound(new_lines, SSC)
        fit_params, covar = fit_spectrum(frequency, spectrum, guess, bound)
        update_params_dict(auto_bandfit, SSC, band, new_lines, fit_params, covar)
        plot_spectrum_fit(SSC, band, blines, frequency, spectrum, fit_params, number_lines=True)

plt.clf()
plt.ion()
fnpickle(auto_bandfit, os.path.join(mandir, 'auto_band_fit_parameters.pickle'))


###################################################################################################
#
###################################################################################################
