###################################################################################################
# export data to be read by XCLASS
###################################################################################################

# It is much easier to deal with a collection of datasets when the data is shifted to the rest
# frequency. That way all sources can be treated together and line misidentification cannot occur.

def export_spectrum(num, SSC, specie, rms):

    bands = bands_from_specie(specie)

    dat_files = []
    for band in bands:

        # get data
        spectrum  = spectra[str(SSC['no'])][band]
        frequency = spectrum['frequency'].to(u.MHz)
        intensity = spectrum['spectrum'].to(u.K)

        # shift spectrum to rest frequency
        velshift  = SSC['velshift']
        frequency = [(-vsys-velshift).to(u.GHz, equivalencies=u.doppler_optical(f)).value for f in frequency]*u.GHz

        # remove NaNs
        frequency, intensity = crossmatch(frequency.to(u.MHz).value, intensity.to(u.K).value)

        # add noise
        if not num==0:
            randstate = np.random.RandomState(num)
            noise =  np.random.normal(loc=0., scale=rms.to(u.K).value, size=len(frequency))
            intensity += noise

        # export to ASCII file
        savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),band+'.dat'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        np.savetxt(savepath, np.transpose([frequency, intensity]), fmt=('%10.4f', '%10.4f'))
        dat_files.append(savepath)
    return dat_files


###################################################################################################
#
###################################################################################################
