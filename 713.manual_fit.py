#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Manually fit (by eye) lines that are not well fit by the band fitting procedure.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))


###################################################################################################
# fitting function definitions
###################################################################################################

def get_spectrum(SSC, line):
    obsfreq = line['obsfreq'].to(u.GHz).value
    band = 'LSB' if obsfreq<350. else 'USB'
    frequency = spectra[str(SSC['no'])][band]['frequency'].value
    spectrum  = spectra[str(SSC['no'])][band]['spectrum'].value
    fmask = np.logical_and(frequency>obsfreq-0.2, frequency<obsfreq+0.2)
    return frequency[fmask], spectrum[fmask]

def lines_in_band(lines, band):
    blines = []
    for line in lines:
        if band=='LSB' and line['restfreq']<350*u.GHz:
            blines.append(line)
        elif band=='USB' and line['restfreq']>350*u.GHz:
            blines.append(line)
    return blines

def gauss(x, amp=1., sigma=1., mu=0., const=0.):
    # return (amp/(sigma*np.sqrt(2*np.pi)) *np.exp(-0.5*((x-mu)/sigma)**2)) +const
    return (amp *np.exp(-1*((x-mu)/sigma)**2)) +const

def fit_spectrum(frequency, spectrum, guess):
    from scipy.optimize import curve_fit
    return curve_fit(gauss, frequency, spectrum, p0=guess, bounds=((0.,0.,-np.inf,-1e-5),(200.,0.25,np.inf,1e-5)))

def plot_spectrum_fit(SSC, band, line, frequency, spectrum, fit_params, number_lines=False):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+line_tex(line), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    # plot spectrum
    ax.plot(frequency, spectrum, ls='-', color='k')
    ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

    # plot fitted spectrum
    fit = gauss(frequency, *fit_params)
    ax.plot(frequency, fit, ls='-', color='r')

    # label expected lines
    blines = lines_in_band(lines, band)
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
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_axisbelow(True)
    ax.grid()
    ax.set_xlabel(r'$\nu$ [GHz]', fontsize=12)
    ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
    fig.tight_layout()

def ask_repeat(SSC, line):
    return str(input('Refine fit parameters for SSC '+str(SSC['no'])+': '+line['molecule']+' '+line['transition']+' '+line['vibration']+'? ["no" to stop]    '))

def get_new_guess(amp, sigma, mu, const):
    amp_new   = input('{:<6}'.format('amp=')  +'{:6.3f}'.format(amp)  +'{:>17}'.format('[new]  amp = '))
    sigma_new = input('{:<6}'.format('sigma=')+'{:6.3f}'.format(sigma)+'{:>17}'.format('[new]  sigma = '))
    mu_new    = input('{:<6}'.format('mu=')   +'{:6.3f}'.format(mu)   +'{:>17}'.format('[new]  mu = '))
    if amp_new!='':
        amp = float(amp_new)
    if sigma_new!='':
        sigma = float(sigma_new)
    if mu_new!='':
        mu = float(mu_new)
    return (amp, sigma, mu, 0.)

def update_params_dict(fit_dict, SSC, band, line, guess):
    amp   = guess[0]
    sigma = guess[1]
    mu    = guess[2]
    fit_dict[str(SSC['no'])][band][line['ID']] = {'amp': [amp, 0.], 'sigma': [sigma, 0.], 'mu': [mu, 0.]}

def guess_from_fit(fit_dict, SSC, band, line):
    amp   = fit_dict[str(SSC['no'])][band][line['ID']]['amp'][0]
    sigma = fit_dict[str(SSC['no'])][band][line['ID']]['sigma'][0]
    mu    = fit_dict[str(SSC['no'])][band][line['ID']]['mu'][0]
    guess = np.array([amp, sigma, mu, 0])
    return guess


###################################################################################################
# manually re-fit problematic lines
###################################################################################################

# lines that need to be re-fitted
refit_indeces = {'LSB': {'1': [4,0],
                         '2': [4,0],
                         '3': [0],
                         '4': [4,0],
                         '5': [4,0],
                         '6': [4,0],
                         '7': [7],
                         '8': [],
                         '9': [4,0],
                         '10': [],
                         '11': [4,0],
                         '12': [4,0],
                         '13': [],
                         '14': [0,21]},
                 'USB': {'1': [1,6,0],
                         '2': [1,6,11,4,0],
                         '3': [1,11,0,26],
                         '4': [1,11,4,0],
                         '5': [],
                         '6': [],
                         '7': [1,0],
                         '8': [1,6,16],
                         '9': [],
                         '10': [6],
                         '11': [],
                         '12': [1,6,0],
                         '13': [1,0],
                         '14': [6]}}

# get estimate from automated fitting
auto_bandfit = fnunpickle(os.path.join(mandir, 'band_fit_parameters.pickle'))

# update band fits
bandfit = copy.deepcopy(auto_bandfit)

# plot spectrum+fit and check
for band in ['LSB','USB']:
    for SSC in SSCs:
        blines = lines_in_band(lines, band)
        refit_index = refit_indeces[band][str(SSC['no'])]
        if not refit_index==[]:
            re_lines = [blines[l] for l in refit_index]
        else:
            continue
        for re_line in re_lines:
            frequency, spectrum = get_spectrum(SSC, re_line)
            guess = guess_from_fit(bandfit, SSC, band, re_line)
            plot_spectrum_fit(SSC, band, re_line, frequency, spectrum, guess, number_lines=True)

            redo = ask_repeat(SSC, re_line)
            while not (redo=='n' or redo=='no'):
                guess = get_new_guess(*guess)
                update_params_dict(bandfit, SSC, band, re_line, guess)
                plot_spectrum_fit(SSC, band, re_line, frequency, spectrum, guess, number_lines=True)
                redo = ask_repeat(SSC, re_line)

fnpickle(bandfit, os.path.join(mandir, 'band_fit_parameters.pickle'))


###################################################################################################
# re-fit again
###################################################################################################

# Fit the same peaks for the double peaks in some SSCs. This was not consistent before.

refit_indeces = {'LSB': {'3': [7]},
                 'USB': {'1': [1,0],
                         '2': [1,0]}}

for band in refit_indeces.keys():
    for s in refit_indeces[band].keys():
        SSC = SSCs[int(s)-1]
        refit_index = refit_indeces[band][s]
        blines = lines_in_band(lines, band)
        if not refit_index==[]:
            re_lines = [blines[l] for l in refit_index]
        else:
            continue
        for re_line in re_lines:
            frequency, spectrum = get_spectrum(SSC, re_line)
            guess = guess_from_fit(bandfit, SSC, band, re_line)
            plot_spectrum_fit(SSC, band, re_line, frequency, spectrum, guess, number_lines=True)

            redo = ask_repeat(SSC, re_line)
            while not (redo=='n' or redo=='no'):
                guess = get_new_guess(*guess)
                update_params_dict(bandfit, SSC, band, re_line, guess)
                plot_spectrum_fit(SSC, band, re_line, frequency, spectrum, guess, number_lines=True)
                redo = ask_repeat(SSC, re_line)

fnpickle(bandfit, os.path.join(mandir, 'band_fit_parameters.pickle'))


###################################################################################################
# remove actually undetected lines
###################################################################################################

remove_fits = [['1',  'LSB', 'HC3N 38-37'],
               ['5',  'LSB', 'HC3N 38-37'],
               ['6',  'LSB', 'HC3N 38-37'],
               ['9',  'LSB', 'HC3N 38-37'],
               ['11', 'LSB', 'HC3N 38-37'],
               ['12', 'LSB', 'HC3N 38-37']]
for (s,b,l) in remove_fits:
    del bandfit[s][b][l]

fnpickle(bandfit, os.path.join(mandir, 'band_fit_parameters.pickle'))


###################################################################################################
#
###################################################################################################
