#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# plot/check each fit
###################################################################################################

def get_spectrum(SSC, band):
    spectra = [ np.genfromtxt(escape_fname(os.path.join(Xfinaldir,'SSC_'+str(SSC['no']),'run_'+str(num),band+'.dat')), dtype=None) for num in np.arange(101) ]
    spectra = np.array(spectra)
    frequency      = np.percentile(spectra[:,:,0], 50, axis=0)/1000
    p16,median,p84 = np.percentile(spectra[:,:,1], (16,50,84), axis=0)
    return frequency,p16,median,p84

def get_params(SSC, ID, component):
    amps, sigmas, mus = [],[],[]
    for num in np.arange(101):
        try:    # sometimes the fit failed, so the dict is empty
            amps.append(   intensity_data[SSC['num']][ID][component][num]['amplitude'].value )
            sigmas.append( intensity_data[SSC['num']][ID][component][num]['sigma'].value     )
            mus.append(    intensity_data[SSC['num']][ID][component][num]['mu'].value        )
        except:
            pass
    amp_p16, amp_median, amp_p84       = np.nanpercentile(amps, (16,50,84))
    sigma_p16, sigma_median, sigma_p84 = np.nanpercentile(sigmas, (16,50,84))
    mu_p16, mu_median, mu_p84          = np.nanpercentile(mus, (16,50,84))
    return [amp_p16,sigma_p16,mu_p16],[amp_median,sigma_median,mu_median],[amp_p84,sigma_p84,mu_p84]

def line_from_ID(ID):
    for line in lines:
        if line['ID']==ID:
            return line

def lines_in_band(band):
    blines = []
    for line in lines:
        if band=='LSB' and line['restfreq']<350*u.GHz:
            blines.append(line)
        elif band=='USB' and line['restfreq']>350*u.GHz:
            blines.append(line)
    return sorted(blines, key=lambda k: k['restfreq'])

def band_from_ID(ID):
    line = line_from_ID(ID)
    if line['restfreq']<350*u.GHz:
        return 'LSB'
    elif line['restfreq']>350*u.GHz:
        return 'USB'
    else:
        raise ValueError("Cannot determine band for ID: "+str(ID))

def specie_from_line(line):
    return line['XCLASS']

def specie_from_ID(ID):
    return specie_from_line( line_from_ID(ID) )


def plot_spectrum_fit(SSC, ID):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+ID.replace('_','$\_$'), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    band = band_from_ID(ID)
    frequency,_,spectrum,_ = get_spectrum(SSC, band)

    # plot spectrum
    ax.plot(frequency, spectrum, ls='-', color='k', alpha=0.5, zorder=3)
    # ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.25, zorder=2)

    # plot XCLASS model
    specie = specie_from_ID(ID)
    for component in np.arange(all_data[SSC['num']][specie]['components']):
        modfreq = model_spectra[SSC['num']][specie][component][band]['frequency']
        p16     = model_spectra[SSC['num']][specie][component][band]['p16']
        median  = model_spectra[SSC['num']][specie][component][band]['median']
        p84     = model_spectra[SSC['num']][specie][component][band]['p84']
        ax.plot(modfreq, median, ls='--', color='k', zorder=3)
        ax.fill_between(modfreq, median, [0. for f in modfreq], color='grey', alpha=0.25, zorder=2)

    # plot fitted spectrum
    norm = colors.Normalize(vmin=0, vmax=len(intensity_data[SSC['num']][ID].keys())-1)
    for component in np.arange(all_data[SSC['num']][specie]['components']):
        p16, median, p84 = get_params(SSC, ID, component)
        color            = cm.get_cmap('rainbow_r')(norm(component))
        ax.plot(frequency, multi_gauss(frequency, *median), ls='-', alpha=0.5, color=color, zorder=4)
        # ax.fill_between(frequency, multi_gauss(frequency, *p16), multi_gauss(frequency, *p84), color=color, alpha=0.25, zorder=5)

    blines = lines_in_band(band)
    for idx,line in enumerate(blines):
        restfreq = line['restfreq'].to(u.GHz).value
        if (restfreq>frequency[0] and restfreq<frequency[-1]):
            if band=='LSB':
                xlim = [342.4, 346.2]
            elif band=='USB':
                xlim = [354.3, 358.1]
            xloc = xlim[0] +((idx+0.5)/len(blines))*(xlim[1]-xlim[0])
            ax.axvline(x=restfreq, ymin=0, ymax=1, color='dimgrey', ls='--', lw=0.5, zorder=1)
            ax.plot([restfreq,xloc], [1.05*np.nanmax(spectrum), 1.05*1.05*np.nanmax(spectrum)], color='dimgrey', ls='--', lw=0.5, zorder=1, clip_on=False)
            ax.text(xloc, 1.06*1.05*np.nanmax(spectrum), line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')

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
    mkdir(escape_fname(join(plotdir, '12.intensity', 'fit_check')))
    fig.savefig(escape_fname(join(plotdir, '12.intensity', 'fit_check', 'SSC_'+str(SSC['no'])+'.'+ID+'.'+band+'.pdf')), dpi=300, bbox_inches='tight')


for SSC in tqdm(SSCs):
    # for ID,data in intensity_data[SSC['num']].items():
    ID='H13CN 4-3'
    try:
        plot_spectrum_fit(SSC, ID)
    except:
        print('failed plot: SSC', SSC['num'], ID)


###################################################################################################
# plot/check combined fit
###################################################################################################

def plot_spectrum_fit_combined(SSC, band, number_lines=False):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+band, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    frequency,_,spectrum,_ = get_spectrum(SSC, band)

    # plot spectrum
    ax.plot(frequency, spectrum, ls='-', color='k', zorder=3)
    ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5, zorder=2)

    # get fit spectrum
    p16, median, p84 = [],[],[]
    for ID,data in intensity_data[SSC['num']].items():
        for component in data.keys():
            fit_p16, fit_median, fit_p84 = get_params(SSC, ID, component)
            p16.append(fit_p16)
            median.append(fit_median)
            p84.append(fit_p84)
    p16    = flatten(p16)
    median = flatten(median)
    p84    = flatten(p84)

    # plot fitted spectrum
    ax.plot(frequency, multi_gauss(frequency, *median), ls='-', color='r', zorder=4)
    ax.fill_between(frequency, multi_gauss(frequency, *p16), multi_gauss(frequency, *p84), color='r', alpha=0.5, zorder=5)

    blines = lines_in_band(band)
    for idx,line in enumerate(blines):
        restfreq = line['restfreq'].to(u.GHz).value
        if (restfreq>frequency[0] and restfreq<frequency[-1]):
            if band=='LSB':
                xlim = [342.4, 346.2]
            elif band=='USB':
                xlim = [354.3, 358.1]
            xloc = xlim[0] +((idx+0.5)/len(blines))*(xlim[1]-xlim[0])
            ax.axvline(x=restfreq, ymin=0, ymax=1, color='dimgrey', ls='--', lw=0.5, zorder=1)
            ax.plot([restfreq,xloc], [1.05*np.nanmax(spectrum), 1.05*1.05*np.nanmax(spectrum)], color='dimgrey', ls='--', lw=0.5, zorder=1, clip_on=False)
            ax.text(xloc, 1.06*1.05*np.nanmax(spectrum), line_tex(line), color='dimgrey', fontsize=10, rotation=90, ha='center', va='bottom')


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
    fig.savefig(join(plotdir, '12.intensity', 'SSC_'+str(SSC['no'])+'.'+band+'.pdf'), dpi=300, bbox_inches='tight')


# plot
###################################################################################################

plt.ioff()
for band in ['LSB','USB']:
    for SSC in tqdm(SSCs):
        # plot_spectrum_fit_combined(SSC, band, number_lines=True)

            try:
                plot_spectrum_fit_combined(SSC, band, number_lines=True)
            except:
                print('failed plot: SSC', SSC['num'], band)


plt.clf()


###################################################################################################
#
###################################################################################################
