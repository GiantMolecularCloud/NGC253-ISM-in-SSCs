###################################################################################################
# figure 2: spectra
###################################################################################################

def spectrum_dat(SSC,band):

    def get_spectra(nums, SSC, band):
        spectra = [ np.genfromtxt(escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat')), dtype=None) for num in nums ]
        spectra = np.array(spectra)

        frequency      = np.percentile(spectra[:,:,0], 50, axis=0)/1000
        p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
        return frequency,p16,median,p84

    frequency, p16,med,p84  = get_spectra(np.arange(101), SSC, band)

    spec = open('spectrum.SSC_'+SSC['num']+'.CbF.dat','a+')
    spec.write('SSC '+SSC['num']+'\n')
    spec.write('    freq        int\n')
    spec.write('    [GHz]       [K]\n')
    spec.write('===================\n')

    for f,m in zip(frequency,med):
        spec.write('{:9.5f}'.format(f)+'   '+'{:7.3f}'.format(m))
        spec.write('\n')

    spec.close()


for SSC in tqdm(SSCs):
    for band in ['LSB','USB']:
        spectrum_dat(SSC,band)
