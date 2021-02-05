#######################
# GAS IN SSCS: XCLASS #
#######################

# address referee report: derive intensities from XCLASS fits


###################################################################################################
# fit model spectra & get intensities
###################################################################################################

# create structure to save data
###################################################################################################

intensity_data = {}
for SSC in SSCs:
    try:
        intensity_data[SSC['num']]
    except:
        intensity_data[SSC['num']] = {}

    for specie in data_XCLASS_refit[SSC['num']].keys():
        slines = specie_lines(specie)
        for sline in slines:
            try:
                intensity_data[SSC['num']][sline['ID']]
            except:
                intensity_data[SSC['num']][sline['ID']] = {}

            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
                try:
                    intensity_data[SSC['num']][sline['ID']][component]
                except:
                    intensity_data[SSC['num']][sline['ID']][component] = {}

                for num in nums:
                    try:
                        intensity_data[SSC['num']][sline['ID']][component][num]
                    except:
                        intensity_data[SSC['num']][sline['ID']][component][num] = {}


# fit spectrum with Gaussian
###################################################################################################

def specie_lines_in_band(specie, band):
    blines = []
    for line in lines:
        if line['XCLASS']==specie:
            if band=='LSB' and line['restfreq']<350*u.GHz:
                blines.append(line)
            elif band=='USB' and line['restfreq']>350*u.GHz:
                blines.append(line)
    return blines

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

def get_guess(SSC,specie,component,line):
    """
    Return a guess of line intensity, position and width based on the XCLASS fits.
    """
    linewidth = data_XCLASS_refit[SSC['num']][specie]['linewidth']['median'][component]*u.km/u.s
    velocity  = data_XCLASS_refit[SSC['num']][specie]['velocity']['median'][component]*u.km/u.s
    lw_freq   = linewidth.to(u.GHz, equivalencies=u.doppler_optical(line['restfreq']))
    v_freq    = velocity.to(u.GHz, equivalencies=u.doppler_optical(line['restfreq']))

    amp   = 1.0
    sigma = (line['restfreq']-lw_freq).value
    mu    = v_freq.value

    return [amp,sigma,mu]

def guess_bounds(SSC, specie, component, bline):
    guess      = []
    low_bounds = []
    up_bounds  = []
    for bline in blines:
        amp,sigma,mu = get_guess(SSC,specie,component,bline)
        guess.append([amp,sigma,mu])
        low_bounds.append( (0.01, sigma/2., mu-0.05))       # amp (K), sigma (GHz), mu (GHz)
        up_bounds.append(  (200., sigma*2., mu+0.05))       # amp (K), sigma (GHz), mu (GHz)
    guess      = flatten(guess)
    low_bounds = tuple(flatten(low_bounds))
    up_bounds  = tuple(flatten(up_bounds))
    return guess, low_bounds, up_bounds

def fit_spectrum(frequency, intensity, guess, low_bounds, up_bounds):
    from scipy.optimize import curve_fit
    return curve_fit(multi_gauss, frequency, intensity, p0=guess, bounds=(low_bounds,up_bounds))


failed = []
for num in nums:
    for SSC in SSCs:
        for specie in data_XCLASS_refit[SSC['num']].keys():
            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):
                for band in ['LSB','USB']:

                    # get lines in band
                    blines = specie_lines_in_band(specie, band)

                    if not blines==[]:

                        guess, low_bounds, up_bounds = guess_bounds(SSC, specie, component, bline)
                        frequency = intensity_spectra[SSCnum][specie][component][num][band][0].value
                        intensity = intensity_spectra[SSCnum][specie][component][num][band][1].value

                        print_overwrite("now fitting: "+'{:3d}'.format(num)+'/100, SSC '+'{:>2}'.format(SSCnum)+' '+'{:>12}'.format(specie)+'     '+str(component))
                        try:
                            fit_params, covar = fit_spectrum(frequency, intensity, guess, low_bounds, up_bounds)
                            for i,bline in enumerate(blines):
                                amp   = fit_params[i*3]*u.K
                                sigma = fit_params[i*3+1]*u.GHz
                                mu    = fit_params[i*3+2]*u.GHz
                                intensity_data_refit[SSCnum][bline['ID']][component][num] = {'amplitude': amp, 'sigma': sigma, 'mu':mu}
                        except:
                            #  break degeneracy between close lines
                            try:
                                for i,bline in enumerate(blines):
                                    guess[i*3] = np.random.random_sample()
                                fit_params, covar = fit_spectrum(frequency, intensity, guess, low_bounds, up_bounds)
                                for i,bline in enumerate(blines):
                                    amp   = fit_params[i*3]*u.K
                                    sigma = fit_params[i*3+1]*u.GHz
                                    mu    = fit_params[i*3+2]*u.GHz
                                    intensity_data_refit[SSCnum][bline['ID']][component][num] = {'amplitude': amp, 'sigma': sigma, 'mu':mu}
                            except:
                                failed.append({'SSC': SSCnum, 'specie': specie, 'component': component, 'band': band, 'num': num})
print(failed)


# get intensity
###################################################################################################

def lineID_lines(lineID):
    """
    Select all lines of a given species.
    """
    for line in lines:
        if line['ID']==lineID:
            return line


for num in tqdm(nums):
    for SSC in SSCs:
        for ID,data in intensity_data[SSC['num']].items():
            line   = lineID_lines(ID)
            specie = line['XCLASS']
            for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):

                if not data[component][num]=={}:
                    amp   = data[component][num]['amplitude']
                    sigma = data[component][num]['sigma']
                    mu    = data[component][num]['mu']

                    # convert to optical velocity
                    position  = mu.to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq']))
                    linewidth = np.abs( (mu+sigma).to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq'])) -position )
                    FWHM      = 2.355*linewidth

                    # integrated intensity
                    intensity = amp*linewidth*np.sqrt(np.pi)

                    intensity_data[SSC['num']][ID][component][num]['position']             = position
                    intensity_data[SSC['num']][ID][component][num]['linewidth']            = linewidth
                    intensity_data[SSC['num']][ID][component][num]['FWHM']                 = FWHM
                    intensity_data[SSC['num']][ID][component][num]['integrated intensity'] = intensity

                else:
                    intensity_data[SSC['num']][ID][component][num]['position']             = np.nan*u.km/u.s
                    intensity_data[SSC['num']][ID][component][num]['linewidth']            = np.nan*u.km/u.s
                    intensity_data[SSC['num']][ID][component][num]['FWHM']                 = np.nan*u.km/u.s
                    intensity_data[SSC['num']][ID][component][num]['integrated intensity'] = np.nan*u.K*u.km/u.s


# get statistics
###################################################################################################

for SSC in tqdm(SSCs):
    for ID,data in intensity_data[SSC['num']].items():
        line   = lineID_lines(ID)
        specie = line['XCLASS']
        for component in np.arange(data_XCLASS_refit[SSC['num']][specie]['components']):

            try:
                intensity_data[SSC['num']][ID][component]['median']['position']
            except:
                intensity_data[SSC['num']][ID][component]['all']    = {}
                intensity_data[SSC['num']][ID][component]['median'] = {}
                intensity_data[SSC['num']][ID][component]['16th']   = {}
                intensity_data[SSC['num']][ID][component]['84th']   = {}

            for x in ['position','linewidth','FWHM','integrated intensity']:
                xs   = [intensity_data[SSC['num']][ID][component][num][x] for num in nums]
                unit = xs[0].unit
                xs   = np.array([i.value for i in xs])

                # ignore bad fit outliers
                p16, median, p84 = np.nanpercentile(xs, (16,50,84))
                if x=='position':
                    mask_bad = [np.abs(xs)-np.abs(median)<50.]
                elif x=='linewidth' or x=='FWHM':
                    mask_bad = np.logical_and(np.abs(xs)>0.1*np.abs(median), np.abs(xs)<10.*np.abs(median))
                elif x=='integrated intensity':
                    mask_bad = np.logical_and(np.abs(xs)>0.1*np.abs(median), np.abs(xs)<10.*np.abs(median))
                else:
                    mask_bad = [True for x in xs]
                xs_masked = xs[mask_bad]
                try:
                    p16, median, p84 = np.nanpercentile(xs_masked, (16,50,84))
                except:
                    print('failed statistics: ', 'SSC', SSC['num'], ID, 'component', component)

                intensity_data[SSC['num']][ID][component]['all'][x]    = xs*unit
                intensity_data[SSC['num']][ID][component]['median'][x] = median*unit
                intensity_data[SSC['num']][ID][component]['16th'][x]   = p16*unit
                intensity_data[SSC['num']][ID][component]['84th'][x]   = p84*unit

fnpickle(intensity_data, join(intensitydir, 'intensity_data.pickle'))


###################################################################################################
# merge into large data dictionary
###################################################################################################

def get_specie_info(ID):
    line   = lineID_lines(ID)
    specie     = line['XCLASS']
    molecule   = line['molecule']
    transition = line['transition']
    vibration  = line['vibration']
    if vibration=='':
        vibration = 'v=0'
    return specie, molecule, transition, vibration


# build data structure for all_data
###################################################################################################

all_data = copy.deepcopy(data_XCLASS_refit)
for SSC in tqdm(SSCs):
    for ID in intensity_data[SSC['num']].keys():
        specie, molecule, transition, vibration = get_specie_info(ID)
        for x_I,x in zip(['position','linewidth','FWHM','integrated intensity'],['position (Gauss)','linewidth (Gauss)','FWHM (Gauss)','integrated intensity']):
            try:
                all_data[SSC['num']][specie][x]
            except:
                all_data[SSC['num']][specie][x] = {}

            for component in np.arange(all_data[SSC['num']][specie]['components']):
                try:
                    all_data[SSC['num']][specie][x][transition]
                except:
                    all_data[SSC['num']][specie][x][transition] = {}
                try:
                    all_data[SSC['num']][specie][x][transition][vibration]
                except:
                    all_data[SSC['num']][specie][x][transition][vibration] = {}
                try:
                    all_data[SSC['num']][specie][x][transition][vibration][component]
                except:
                    all_data[SSC['num']][specie][x][transition][vibration][component] = {}

                for y in ['median','16th','84th','all']:
                    try:
                        all_data[SSC['num']][specie][x][transition][vibration][component][y]
                    except:
                        all_data[SSC['num']][specie][x][transition][vibration][component][y] = {}


# fill data structure with fitted values
###################################################################################################

for SSC in tqdm(SSCs):
    for ID in intensity_data[SSC['num']].keys():
        specie, molecule, transition, vibration = get_specie_info(ID)
        for x_I,x in zip(['position','linewidth','FWHM','integrated intensity'],['position (Gauss)','linewidth (Gauss)','FWHM (Gauss)','integrated intensity']):
            for component in np.arange(all_data[SSC['num']][specie]['components']):
                for y in ['median','16th','84th','all']:
                    all_data[SSC['num']][specie][x][transition][vibration][component][y] = copy.deepcopy(intensity_data[SSC['num']][ID][component][y][x_I])

fnpickle(all_data, join(paperdir, 'all_data.pickle'))


###################################################################################################
#
###################################################################################################
