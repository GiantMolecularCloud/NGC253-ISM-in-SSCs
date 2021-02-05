#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Get physical quantities from fits, such as integrated intensity, mass, linewidth, line shift, ...


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
bandfit = fnunpickle(os.path.join(mandir, 'updated_band_fit_parameters.pickle'))


###################################################################################################
# helper functions
###################################################################################################

def get_line(l):
    for line in lines:
        if line['ID']==l:
            return line


###################################################################################################
# calculate line properties
###################################################################################################

line_data_I = {str(SSC['no']): {} for SSC in SSCs}
for s in bandfit.keys():
    SSC = SSCs[int(s)-1]
    for b in bandfit[s].keys():
        for l in bandfit[s][b].keys():
            try:
                line = get_line(l)
                fit_params = bandfit[s][b][l]
                amp       = fit_params['amp'][0]   *u.Kelvin
                amp_err   = fit_params['amp'][1]   *u.Kelvin
                sigma     = fit_params['sigma'][0] *u.GHz
                sigma_err = fit_params['sigma'][1] *u.GHz
                mu        = fit_params['mu'][0]    *u.GHz
                mu_err    = fit_params['mu'][1]    *u.GHz

                # convert to optical velocity
                lp       = mu.to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq']))
                lp_err   = np.abs( (mu+mu_err).to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq']))-lp )
                lw       = np.abs( (mu+sigma).to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq'])) -lp )
                lw_err   = np.abs( (mu+sigma_err).to(u.km/u.s, equivalencies=u.doppler_optical(line['restfreq'])) -lp )
                FWHM     = 2.355*lw
                FWHM_err = 2.355*lw_err

                # estimate error for by eye fits without error
                if mu_err==0.:
                    lp_err=2.5*u.km/u.s
                if lp_err>100*u.km/u.s:
                    lp_err = lw
                if sigma_err==0.:
                    lw_err   = 2.5*u.km/u.s
                    FWHM_err = 2.355*lw_err

                # line shift
                shift     = lp-vsys
                shift_err = lp_err

                # integrated intensity
                int_int     = amp*lw*np.sqrt(np.pi)
                int_int_err = int_int *np.sqrt( (amp_err/amp)**2 + (lw_err/lw)**2 )
                if amp_err==0.:
                    int_int_err = int_int *np.sqrt( (0.46*u.K/amp)**2 + (lw_err/lw)**2 )

                # store data
                line_data_I[str(SSC['no'])][line['ID']] = {'line'                      : line,
                                                           'amplitude'                 : {'bestfit': amp, 'error': amp_err},
                                                           'sigma'                     : {'bestfit': sigma, 'error': sigma_err},
                                                           'mu'                        : {'bestfit': mu, 'error': mu_err},
                                                           'line position'             : {'bestfit': lp, 'error': lp_err},
                                                           'line width'                : {'bestfit': lw, 'error': lw_err},
                                                           'FWHM'                      : {'bestfit': FWHM, 'error': FWHM_err},
                                                           'integrated intensity'      : {'bestfit': int_int, 'error': int_int_err},
                                                           'line shift'                : {'bestfit': shift, 'error': shift_err}
                                                          }
            except:
                print(s,b,l)

# save SSC table
fnpickle(line_data_I, os.path.join(mandir,'line_intensity_data.pickle'))


###################################################################################################
# plot fit parameters
###################################################################################################

def plot_parameters_I(SSC):
    """
    Make an errorbar plot of all detected species for a given SSC.
    """

    for Q,ylab,ymin,ymax in [['integrated intensity',r'$\int I \mathrm{d}v$ [km\,s$^{-1}$]',1e-1,1e5],['line position',r'$v$ [km\,s$^{-1}$]',150,350],['line width',r'$\sigma$ [km\,s$^{-1}$]',0,100]]:

        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+Q, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

        colors = [plt.cm.inferno(i/(len(lines)+1)) for i,_ in enumerate(lines)]
        for idx,(line,c) in enumerate(zip(lines, colors)):
            if line['ID'] in list(line_data_I[str(SSC['no'])].keys()):
                q     = line_data_I[str(SSC['no'])][line['ID']][Q]['bestfit'].value
                q_err = line_data_I[str(SSC['no'])][line['ID']][Q]['error'].value
                ax.errorbar(idx, q, yerr=q_err, marker='o', ms=6, color=c, elinewidth=2, ecolor=c)

        ax.set_xlim(-1, len(line_data_I[str(SSC['no'])]))
        ax.set_ylim(ymin,ymax)
        ax.set_axisbelow(True)
        ax.grid(axis='y')
        ax.set_xticks(np.arange(len(lines)))
        ax.set_xticklabels([line_tex(l) for l in lines])
        ax.tick_params(axis='x', rotation=90)
        if Q=='integrated intensity':
            ax.set_yscale('log')
        ax.set_ylabel(ylab, fontsize=12)
        fig.tight_layout()

        savepath = os.path.join(plotdir, '04.fit_results', 'parameters_I', 'SSC_'+str(SSC['no'])+'.'+Q.replace(' ','_')+'.pdf')
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


###################################################################################################
# plot fit parameters
###################################################################################################

for SSC in tqdm(SSCs):
    plot_parameters_I(SSC)


###################################################################################################
# plot quantity comparison
###################################################################################################

def plot_SSC_quantity_I(lineID, quantity, label):
    """
    Make a box-and-wiskers plot of an XLCASS quantity in a given SSC.
    """

    line = get_line(lineID)

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, lineID_tex(line['ID']), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        if line['ID'] in line_data_I[str(SSC['no'])]:
            fit = line_data_I[str(SSC['no'])][line['ID']][quantity]['bestfit'].value
            err = line_data_I[str(SSC['no'])][line['ID']][quantity]['error'].value
            ax.errorbar(SSC['no'], fit, yerr=err, marker='o', ms=6, color=color, elinewidth=2, ecolor=color)

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(label, fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'SSC_'+quantity.replace(' ','_')+'_I', quantity.replace(' ','_')+'.'+line['ID'].replace(' ','_')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


def compare_linewidths_I(lineIDs):
    """
    Compare linewidths of the given lines and Leroy+18 values.
    """

    lines = [get_line(lineID) for lineID in lineIDs]

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))

    colors = [plt.cm.inferno(i/(len(lines)+0.2)) for i in np.arange(len(lines)+1)]

    # my measurements
    for idx,line in enumerate(lines):
        SSC_no = [SSC['no'] for SSC in SSCs if line['ID'] in line_data_I[str(SSC['no'])]]
        fit    = [line_data_I[str(SSC['no'])][line['ID']]['line width']['bestfit'].value for SSC in SSCs if line['ID'] in line_data_I[str(SSC['no'])]]
        err    = [line_data_I[str(SSC['no'])][line['ID']]['line width']['error'].value for SSC in SSCs if line['ID'] in line_data_I[str(SSC['no'])]]
        ax.errorbar(SSC_no, fit, yerr=err, ls='', marker='o', ms=6, color=colors[idx], elinewidth=2, ecolor=colors[idx], label=lineID_tex(line['ID']))

    # Leroy measurements
    SSC_no    = [SSC['no'] for SSC in SSCs]
    SSC_sigma = [SSC['sigma_v'].value for SSC in SSCs]
    SSC_err   = [SSC['err_sigma_v'].value for SSC in SSCs]
    ax.errorbar(SSC_no, SSC_sigma, yerr=SSC_err, ls='', marker='_', ms=6, color=colors[-1], elinewidth=2, ecolor=colors[-1], label='Leroy+18')

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'$\sigma$ [km\,s$^{-1}$]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.legend(loc=2, bbox_to_anchor=(0.02, 0.975), bbox_transform=ax.transAxes)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'linewidth_comparison_I.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


###################################################################################################
# plot quantity comparison
###################################################################################################

compare_linewidths_I(['CS 7-6','H13CN 4-3','CO 3-2','HCN 4-3','HCO+ 4-3'])

plot_SSC_quantity_I('CS 7-6', 'line width', r'$\sigma$ [km\,s$^{-1}$]')
plot_SSC_quantity_I('H13CN 4-3', 'line width', r'$\sigma$ [km\,s$^{-1}$]')

# HC3N
plot_SSC_quantity_I('HC3N 38-37', 'integrated intensity', '$\int I \mathrm{d}v$')
plot_SSC_quantity_I('HC3N 39-38', 'integrated intensity', '$\int I \mathrm{d}v$')

# sulfur chemistry
plot_SSC_quantity_I('SO 8(8)-7(7) 3Sum_v=0', 'integrated intensity', '$\int I \mathrm{d}v$')
plot_SSC_quantity_I('SO2 11(4,8)-11(3,9)',   'integrated intensity', '$\int I \mathrm{d}v$')
plot_SSC_quantity_I('CS 7-6',                'integrated intensity', '$\int I \mathrm{d}v$')
plot_SSC_quantity_I('H2CS 10(3,8)-9(3,7)',   'integrated intensity', '$\int I \mathrm{d}v$')


###################################################################################################
#
###################################################################################################
