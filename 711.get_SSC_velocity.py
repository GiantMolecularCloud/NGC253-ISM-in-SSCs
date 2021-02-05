#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# In order to fit the spectra, the expected position of the line needs to be known. This can be
# done most reliably with the CS (7-6) line as there are few nearby lines, has a nice shape and
# is quite bright.


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

def initial_guess(line, spectrum):
    amp   = np.nanmax(spectrum)                     # K
    sigma = 0.05                                    # GHz
    mu    = line['obsfreq'].to(u.GHz).value         # GHz
    const = 0.                                      # K, should be zero with the improved contsub
    return (amp, sigma, mu, const)

def gauss(x, amp=1., sigma=1., mu=0., const=0.):
    # return (amp/(sigma*np.sqrt(2*np.pi)) *np.exp(-0.5*((x-mu)/sigma)**2)) +const
    return (amp *np.exp(-1*((x-mu)/sigma)**2)) +const

def fit_spectrum(frequency, spectrum, guess):
    from scipy.optimize import curve_fit
    return curve_fit(gauss, frequency, spectrum, p0=guess, bounds=((0.,0.,-np.inf,-1e-5),(200.,0.25,np.inf,1e-5)))

def plot_spectrum_fit(SSC, line, frequency, spectrum, fit_params):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+line_tex(line), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16)

    # plot spectrum
    ax.plot(frequency, spectrum, ls='-', color='k')
    ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

    # plot fitted spectrum
    fit = gauss(frequency, *fit_params)
    ax.plot(frequency, fit, ls='-', color='r')

    # label expected lines
    for line in lines:
        obsfreq = line['obsfreq'].to(u.GHz).value
        if (obsfreq>frequency[0] and obsfreq<frequency[-1]):
            ax.axvline(x=line['obsfreq'].to(u.GHz).value, ymin=0, ymax=1, color='dimgrey', ls='--')
            ax.text(line['obsfreq'].to(u.GHz).value, 1.02*ax.get_ylim()[1], line_tex(line), color='dimgrey', rotation=90, ha='center', va='bottom')

    ax.set_xlim([frequency[0], frequency[-1]])
    ax.set_ylim(-10, 1.05*np.nanmax(spectrum))
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
    return str(input('Refine fit initial guess for SSC '+str(SSC['no'])+': '+line['molecule']+' '+line['transition']+' '+line['vibration']+'? ["no" to stop]    '))

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

def update_params_dict(SSC, line, amp, sigma, mu, const):
    CSfit[str(SSC['no'])][line['molecule']+line['transition']+line['vibration']] = {'amp': amp, 'sigma': sigma, 'mu': mu, 'const': const}


###################################################################################################
# determine velocity shift of SSCs relative to systemic velocity by manually fitting CS
###################################################################################################

# line fits
CS    = next((line for line in lines if line['molecule']=='CS'), None)
CSfit = {str(SSC['no']): {} for SSC in SSCs}

# plot spectrum+fit and check
for SSC in SSCs:

    frequency, spectrum = get_spectrum(SSC, CS)
    guess = initial_guess(CS, spectrum)
    fit_params, covar = fit_spectrum(frequency, spectrum, guess)
    update_params_dict(SSC, CS, *fit_params)
    plot_spectrum_fit(SSC, CS, frequency, spectrum, fit_params)
    redo = ask_repeat(SSC, CS)

    while not (redo=='n' or redo=='no'):
        guess = get_new_guess(*fit_params)
        fit_params, covar = fit_spectrum(frequency, spectrum, guess)
        update_params_dict(SSC, CS, *fit_params)
        plot_spectrum_fit(SSC, CS, frequency, spectrum, fit_params)
        redo = ask_repeat(SSC, CS)

fnpickle(CSfit, os.path.join(mandir, 'CS_fit_parameters.pickle'))


###################################################################################################
# calculate offsets
###################################################################################################

shifts = []
for SSC in SSCs:
    # diff_freq = CS['restfreq']-CSfit[str(SSC['no'])]['CS7-6']['mu']*u.GHz
    obsfreq   = CSfit[str(SSC['no'])]['CS7-6']['mu']*u.GHz
    diff_velo = obsfreq.to(u.km/u.s, equivalencies=u.doppler_optical(CS['restfreq']))
    shifts.append(diff_velo.value)

shifts = [s-250 for s in shifts]
SSCs.add_column(Column(name='velshift', data=shifts, unit='km/s', description='Velocity shift relative to v_sys=250 km/s'))

SSCs.write(os.path.join(subprojectdir,'SSCs.fits'), format='fits', overwrite=True)


###################################################################################################
#
###################################################################################################
