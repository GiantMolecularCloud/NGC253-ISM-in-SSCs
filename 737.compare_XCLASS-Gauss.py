#######################
# GAS IN SSCS: XCLASS #
#######################

# Compare the simple Gaussian fits to full-scale modelling in XLCASS


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

data_N   = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))
data_I   = fnunpickle(os.path.join(mandir,'line_intensity_data.pickle'))
ratios_N = fnunpickle(os.path.join(Xfinaldir, 'ratios_column_density.pickle'))
ratios_I = fnunpickle(os.path.join(mandir, 'ratios_intensity.pickle'))


###################################################################################################
# plot opacity vs. opacity tracer ratio
###################################################################################################

def plot_opacity_check():
    """
    Plot XCLASS opacity vs. HCN/H13CN ratio as proxy for opacity.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(8,8))
    # ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        nearest = np.argmin(np.abs(data_N[SSC['num']]['HCN;v=0']['velocity']['median']))
        opt     = data_N[SSC['num']]['HCN;v=0']['integrated opacity']['median'][nearest]
        opt16   = data_N[SSC['num']]['HCN;v=0']['integrated opacity']['16th'][nearest]
        opt84   = data_N[SSC['num']]['HCN;v=0']['integrated opacity']['84th'][nearest]
        opt_l   = opt-opt16
        opt_u   = opt84-opt
        ratio   = ratios_I['HCN/HCNthin'][SSC['num']]['bestfit']
        ratio_e = ratios_I['HCN/HCNthin'][SSC['num']]['error']
        ax.errorbar(opt,ratio,
                    xerr   = [[opt_l],[opt_u]],
                    yerr   = [[ratio_e],[ratio_e]],
                    marker = 'o', ms=5, color=color, elinewidth=1, ecolor=color, ls='',
                    label  = SSC['num'],
                    zorder = 4
                   )

    ax.set_xlim(0.5,20)
    ax.set_ylim(1,100)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_ylabel(r'N (HCN) / N (H$^{13}$CN)', fontsize=12)
    ax.set_xlabel(r'$\tau$ (HCN)', fontsize=12)
    fig.legend(loc=3, bbox_to_anchor=(0.,1.,1.,0.08), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '08.compare_XCLASS-Gauss', 'opacity-HCN_ratio.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


plot_opacity_check()


###################################################################################################
#
###################################################################################################
