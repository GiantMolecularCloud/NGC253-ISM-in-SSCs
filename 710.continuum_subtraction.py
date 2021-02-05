#################################
# GAS IN SSCS: DATA PREPARATION #
#################################

# Data must be prepared. Continuum subtraction is good in LSB but must be better in USB

###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))


###################################################################################################
# update SSCs table with further info
###################################################################################################

LSB_x, LSB_y = np.round(LSB.wcs.all_world2pix(SSCs['RA'].value, SSCs['DEC'].value,1,1))[0:2]
USB_x, USB_y = np.round(LSB.wcs.all_world2pix(SSCs['RA'].value, SSCs['DEC'].value,1,1))[0:2]

SSCs.add_column(Column(name='LSB x', data=LSB_x, dtype='int'))
SSCs.add_column(Column(name='LSB y', data=LSB_y, dtype='int'))
SSCs.add_column(Column(name='USB x', data=USB_x, dtype='int'))
SSCs.add_column(Column(name='USB y', data=USB_y, dtype='int'))

SSCs.write(os.path.join(subprojectdir,'SSCs.fits'), format='fits', overwrite=True)


###################################################################################################
# plot spectra at SSC positions before continuum subtraction
###################################################################################################

# all SSCs
fig,axes = plt.subplots(nrows=len(SSCs), ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(10,20))
fig.subplots_adjust(hspace=0,wspace=0)

for i_SSC,SSC in tqdm(enumerate(SSCs)):
    axes[i_SSC][1].tick_params(labelleft=False)
    for i_band,band in enumerate([LSB,USB]):
        ax = axes[i_SSC,i_band]

        # SSC overlay
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no']), color='k', transform=ax.transAxes, ha='left', va='top')

        # spectrum
        spectrum  = band[:,SSC['LSB y'],SSC['LSB x']]
        frequency = spectrum.spectral_axis.to(u.GHz).value
        spectrum  = spectrum.value
        ax.plot(frequency, spectrum, ls='-', color='k')
        ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

        axes[i_SSC][i_band].set_xlim([frequency[0]-0.1, frequency[-1]+0.1])
        axes[i_SSC][i_band].set_ylim(-15,120)
        axes[i_SSC][i_band].xaxis.set_major_locator(MultipleLocator(0.5))
        axes[i_SSC][i_band].xaxis.set_minor_locator(MultipleLocator(0.1))
        axes[i_SSC][i_band].yaxis.set_major_locator(MultipleLocator(25))
        axes[i_SSC][i_band].yaxis.set_minor_locator(MultipleLocator(5))
        axes[i_SSC][i_band].tick_params(axis='both', which='major', labelsize=10)
        axes[i_SSC][i_band].tick_params(labelbottom=False, left=True, right=True, top=True, bottom=True)
        axes[i_SSC][i_band].set_axisbelow(True)
        axes[i_SSC][i_band].grid()

axes[-1][0].xaxis.set_visible(True)
axes[-1][1].xaxis.set_visible(True)
axes[-1][0].tick_params(labelbottom=True)
axes[-1][1].tick_params(labelbottom=True)
axes[-1][0].set_xlabel(r'$\nu$ [GHz]', fontsize=12)
axes[-1][0].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
axes[-1][0].tick_params(labelbottom=True)
axes[-1][1].tick_params(labelbottom=True)

fig.savefig(os.path.join(plotdir, 'contsub_before.pdf'), dpi=300, bbox_inches='tight')


# single SSCs
for i_SSC,SSC in tqdm(enumerate(SSCs)):
    fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    fig.subplots_adjust(hspace=0,wspace=0)
    for i_band,band in enumerate([LSB,USB]):

        axes[i_band].text(0.05, 0.9, 'SSC '+str(SSC['no']), color='k', transform=axes[i_band].transAxes, ha='left', va='top')

        spectrum  = band[:,SSC['LSB y'],SSC['LSB x']]
        frequency = spectrum.spectral_axis.to(u.GHz).value
        spectrum  = spectrum.value
        axes[i_band].plot(frequency, spectrum, ls='-', color='k')
        axes[i_band].fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

        axes[i_band].set_xlim([frequency[0]-0.1, frequency[-1]+0.1])
        axes[i_band].set_ylim(-15,120)
        axes[i_band].xaxis.set_major_locator(MultipleLocator(0.5))
        axes[i_band].xaxis.set_minor_locator(MultipleLocator(0.1))
        axes[i_band].yaxis.set_major_locator(MultipleLocator(25))
        axes[i_band].yaxis.set_minor_locator(MultipleLocator(5))
        axes[i_band].tick_params(axis='both', which='major', labelsize=10)
        axes[i_band].set_axisbelow(True)
        axes[i_band].grid()
        axes[i_band].set_xlabel(r'$\nu$ [GHz]', fontsize=12)
        axes[i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
    axes[0].tick_params(labelbottom=True, left=True, right=False, top=False, bottom=True)
    axes[1].tick_params(labelbottom=True, left=False, right=True, top=False, bottom=True)
    axes[1].yaxis.set_label_text('')
    fig.savefig(os.path.join(plotdir, 'contsub', 'SSC_'+str(SSC['no'])+'.before.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fit for continuum
###################################################################################################

# from astropy.convolution import Gaussian1DKernel
# LSB_sm = LSB.spectral_smooth(Gaussian1DKernel(3))
# USB_sm = USB.spectral_smooth(Gaussian1DKernel(3))

def linear(x,a,b):
    return a*x+b

def plot_continuum_estimate(SSC,band,a,b):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(20,15))
    fig.subplots_adjust(hspace=0,wspace=0)
    ax.text(0.05, 0.9, 'SSC '+str(SSC['no']), color='k', transform=ax.transAxes, ha='left', va='top')

    spectrum  = band[:,SSC['LSB y'],SSC['LSB x']]
    frequency = spectrum.spectral_axis.to(u.GHz).value
    spectrum  = spectrum.value

    # continuum
    continuum = linear(frequency-355.,a,b)

    ax.plot(frequency, spectrum, ls='-', color='k')
    ax.plot(frequency, continuum, ls='-', color='r')
    ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

    ax.text(0.95, 0.9, '$a='+'{:.3f}'.format(a)+'$\n$b='+'{:.3f}'.format(b)+'$', color='r', transform=ax.transAxes, ha='right', va='top')

    ax.set_axisbelow(True)
    ax.grid()


# continuum fits
continua = {'LSB': {},
            'USB': {}}

# plot continuum and check
for SSC in SSCs:
    a = 0.
    b = 0.5
    plot_continuum_estimate(SSC, USB, a, b)
    redo = str(input('Refine estimate for SSC '+str(SSC['no'])+'? ["no" to stop]    '))

    while not (redo=='n' or redo=='no'):
        a_new = input('a='+'{:6.3f}'.format(a)+'\t[new]  a = ')
        b_new = input('b='+'{:6.3f}'.format(b)+'\t[new]  b = ')
        if a_new!='':
            a = float(a_new)
        if b_new!='':
            b = float(b_new)
        continua['USB'][str(SSC['no'])] = {'a': a, 'b': b}

        plot_continuum_estimate(SSC, USB, a, b)
        redo = str(input('Refine estimate for SSC '+str(SSC['no'])+'? ["no" to stop]    '))

fnpickle(continua, os.path.join(mandir, 'continuum_fit_parameters.pickle'))


###################################################################################################
# apply continuum subtraction and write to disk
###################################################################################################

spectra = {str(SSC['no']): {'LSB': {}, 'USB': {}} for SSC in SSCs}

# LSB
for SSC in tqdm(SSCs):
    spectrum  = LSB[:,SSC['LSB y'],SSC['LSB x']]
    frequency = spectrum.spectral_axis.to(u.GHz)
    spectra[str(SSC['no'])]['LSB']['frequency']              = frequency
    spectra[str(SSC['no'])]['LSB']['spectrum (pre-contsub)'] = spectrum
    spectra[str(SSC['no'])]['LSB']['spectrum']               = spectrum

# USB
for SSC in tqdm(SSCs):
    spectrum = USB[:,SSC['USB y'],SSC['USB x']]
    frequency = spectrum.spectral_axis.to(u.GHz)
    spectra[str(SSC['no'])]['USB']['frequency']               = frequency
    spectra[str(SSC['no'])]['USB']['spectrum (pre-contsub)']  = spectrum

    a = continua['USB'][str(SSC['no'])]['a']/u.GHz*u.K
    b = continua['USB'][str(SSC['no'])]['b']*u.K
    continuum = linear(frequency-355.*u.GHz,a,b)
    spectrum_contsub = spectrum-continuum

    spectra[str(SSC['no'])]['USB']['spectrum']                = spectrum_contsub


fnpickle(spectra, os.path.join(mandir, 'spectra.pickle'))


###################################################################################################
# re-check after contsub
###################################################################################################

# all SSCs
fig,axes = plt.subplots(nrows=len(SSCs), ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(10,20))
fig.subplots_adjust(hspace=0,wspace=0)

for i_SSC,SSC in tqdm(enumerate(SSCs)):
    axes[i_SSC][1].tick_params(labelleft=False)
    for i_band,band in enumerate(['LSB','USB']):
        ax = axes[i_SSC,i_band]

        # SSC overlay
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no']), color='k', transform=ax.transAxes, ha='left', va='top')

        # spectrum
        frequency = spectra[str(SSC['no'])][band]['frequency'].value
        spectrum  = spectra[str(SSC['no'])][band]['spectrum'].value
        ax.plot(frequency, spectrum, ls='-', color='k')
        ax.fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

        axes[i_SSC][i_band].set_xlim([frequency[0]-0.1, frequency[-1]+0.1])
        axes[i_SSC][i_band].set_ylim(-15,120)
        axes[i_SSC][i_band].xaxis.set_major_locator(MultipleLocator(0.5))
        axes[i_SSC][i_band].xaxis.set_minor_locator(MultipleLocator(0.1))
        axes[i_SSC][i_band].yaxis.set_major_locator(MultipleLocator(25))
        axes[i_SSC][i_band].yaxis.set_minor_locator(MultipleLocator(5))
        axes[i_SSC][i_band].tick_params(axis='both', which='major', labelsize=10)
        axes[i_SSC][i_band].tick_params(labelbottom=False, left=True, right=True, top=True, bottom=True)
        axes[i_SSC][i_band].set_axisbelow(True)
        axes[i_SSC][i_band].grid()

axes[-1][0].xaxis.set_visible(True)
axes[-1][1].xaxis.set_visible(True)
axes[-1][0].tick_params(labelbottom=True)
axes[-1][1].tick_params(labelbottom=True)
axes[-1][0].set_xlabel(r'$\nu$ [GHz]', fontsize=12)
axes[-1][0].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
axes[-1][0].tick_params(labelbottom=True)
axes[-1][1].tick_params(labelbottom=True)

fig.savefig(os.path.join(plotdir, 'contsub_after.pdf'), dpi=300, bbox_inches='tight')


# single SSCs
for i_SSC,SSC in tqdm(enumerate(SSCs)):
    fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    fig.subplots_adjust(hspace=0,wspace=0)
    for i_band,band in enumerate(['LSB','USB']):

        axes[i_band].text(0.05, 0.9, 'SSC '+str(SSC['no']), color='k', transform=axes[i_band].transAxes, ha='left', va='top')

        frequency = spectra[str(SSC['no'])][band]['frequency'].value
        spectrum  = spectra[str(SSC['no'])][band]['spectrum'].value
        axes[i_band].plot(frequency, spectrum, ls='-', color='k')
        axes[i_band].fill_between(frequency, spectrum, [0. for f in frequency], color='grey', alpha=0.5)

        axes[i_band].set_xlim([frequency[0]-0.1, frequency[-1]+0.1])
        axes[i_band].set_ylim(-15,120)
        axes[i_band].xaxis.set_major_locator(MultipleLocator(0.5))
        axes[i_band].xaxis.set_minor_locator(MultipleLocator(0.1))
        axes[i_band].yaxis.set_major_locator(MultipleLocator(25))
        axes[i_band].yaxis.set_minor_locator(MultipleLocator(5))
        axes[i_band].tick_params(axis='both', which='major', labelsize=10)
        axes[i_band].set_axisbelow(True)
        axes[i_band].grid()
        axes[i_band].set_xlabel(r'$\nu$ [GHz]', fontsize=12)
        axes[i_band].set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
    axes[0].tick_params(labelbottom=True, left=True, right=False, top=False, bottom=True)
    axes[1].tick_params(labelbottom=True, left=False, right=True, top=False, bottom=True)
    axes[1].yaxis.set_label_text('')
    fig.savefig(os.path.join(plotdir, 'contsub', 'SSC_'+str(SSC['no'])+'.after.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
