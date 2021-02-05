###################################################################################################
# check convergence
###################################################################################################

def plot_chi2_evolution(SSC, specie):
    """
    run the functions that build up the figure
    """

    def set_up_figure(SSC, specie):
        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+specie_tex(specie), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,ax

    def plot_chi2(ax, iterations, chi2_orig, lower, upper):
        ax.plot(iterations, chi2_orig, ls='-', color='k', zorder=3)
        ax.fill_between(iterations, lower, upper, color='grey', alpha=0.5, zorder=2)

    def format_figure(ax, iterations, chi2_orig):
        ax.set_yscale('log')
        ax.set_xlim([iterations[0]-1, iterations[-1]+1])
        # ax.set_ylim(-0.05*np.nanmax(spectrum), 1.05*np.nanmax(spectrum))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='grey')
        ax.set_xlabel(r'iterations', fontsize=12)
        ax.set_ylabel(r'$\chi^2$', fontsize=12)
        fig.set_tight_layout(True)

    def save_figure(fig, specie):
        savepath = escape_fname(os.path.join(plotdir, '03.XCLASS_fit', 'SSC_'+str(SSC['no']), specie+'.chi2.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    percentiles = parse_log(SSC, specie)
    chi2_orig  = percentiles[:,0]
    lower      = percentiles[:,1]
    median     = percentiles[:,2]
    upper      = percentiles[:,3]
    iterations = np.arange(1,len(chi2_orig)+1)
    fig,ax = set_up_figure(SSC, specie)
    plot_chi2(ax, iterations, chi2_orig, lower, upper)
    format_figure(ax, iterations, chi2_orig)
    save_figure(fig, specie)


###################################################################################################
#
###################################################################################################
