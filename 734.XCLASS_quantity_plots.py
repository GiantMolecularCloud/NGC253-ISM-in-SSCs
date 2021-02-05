#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

Xfiles      = fnunpickle(os.path.join(Xfinaldir, 'Xfiles.pickle'))
temperature = fnunpickle(os.path.join(Xfinaldir, 'temperature.pickle'))
nums        = fnunpickle(os.path.join(Xfinaldir, 'nums.pickle'))
final_data = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))

execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))


###################################################################################################
# compare linewidths
###################################################################################################

def compare_linewidths(species):
    """
    Compare linewidths of the given species and Leroy+18 values.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))

    colors = [plt.cm.inferno(i/(len(species)+0.2)) for i in np.arange(len(species)+1)]

    # my measurements
    for idx,specie in enumerate(species):
        spx_shift = -0.3+idx*0.6/len(species)

        SSC_no = [SSC['no'] for SSC in SSCs if specie in final_data[str(SSC['no'])]]

        med,low,upp = [],[],[]
        for SSC in SSCs:
            nearest_comp = np.argmin(np.abs(final_data[str(SSC['no'])][specie]['velocity']['median']))
            median = final_data[str(SSC['no'])][specie]['linewidth']['median'][nearest_comp]
            p16    = final_data[str(SSC['no'])][specie]['linewidth']['16th'][nearest_comp]
            p84    = final_data[str(SSC['no'])][specie]['linewidth']['84th'][nearest_comp]
            med.append(median)
            low.append(median-p16)
            upp.append(p84-median)
        ax.errorbar([s+spx_shift for s in SSC_no], med, yerr=[low,upp], ls='', marker='o', ms=6, color=colors[idx], elinewidth=2, ecolor=colors[idx], label=specie_tex(specie))
        # ax.plot([s+spx_shift for s in SSC_no], med, ls='', marker='o', ms=6, color=colors[idx], label=specie_tex(specie))

    # Leroy measurements
    SSC_no    = [SSC['no'] for SSC in SSCs]
    SSC_sigma = [SSC['sigma_v'].value for SSC in SSCs]
    SSC_err   = [SSC['err_sigma_v'].value for SSC in SSCs]
    ax.errorbar([s+0.3 for s in SSC_no], SSC_sigma, yerr=SSC_err, ls='', marker='_', ms=6, color=colors[-1], elinewidth=2, ecolor=colors[-1], label='Leroy+18')

    # plot stripes for clarity
    for idx,SSC in enumerate(SSCs):
        if idx%2==0:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(0.5, len(SSCs)+0.5)
    ax.set_ylim(0,100)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'$\sigma$ [km\,s$^{-1}$]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.legend(loc=2, bbox_to_anchor=(0.02, 0.975), bbox_transform=ax.transAxes)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.XCLASS_final', 'linewidth_comparison.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_linewidths(['CS;v=0','HC-13-N;v=0','CO;v=0','HCN;v=0','HCO+;v=0'])


###################################################################################################
#
###################################################################################################
