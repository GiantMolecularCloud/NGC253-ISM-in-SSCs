###################################################################################################
# plot histogram of fit parameters
###################################################################################################

def plot_fit_hist(SSC, specie):
    """
    run the functions that build up the figure
    """

    def set_up_figure(SSC, specie):
        fig,axes = plt.subplots(nrows=4, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,10))
        fig.subplots_adjust(hspace=0.5)
        axes[0].text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+specie_tex(specie), color='k', transform=axes[0].transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
        return fig,axes

    def plot_hist(ax, data, label):
        from matplotlib import cm
        colors = [cm.inferno(x) for x in np.linspace(0, 1, len(data))]
        for c,component in enumerate(data):
            if np.nanmax(component)>1e5:
                bins = np.logspace(np.log10(np.nanmin(component)),np.log10(np.nanmax(component)), 20)
                ax.set_xscale('log')
            else:
                bins = 20
            ax.hist(component, bins=bins, histtype='step',       color=colors[c],            label=str(c+1))
            ax.hist(component, bins=bins, histtype='stepfilled', color=colors[c], alpha=0.5, label=None)
            ax.set_xlabel(label, fontsize=12)

    def format_figure(axes):
        for ax in axes:
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_axisbelow(True)
            ax.grid(ls=':', c='grey')
            ax.set_ylabel(r'frequency', fontsize=12)
        fig.set_tight_layout(True)
        axes[0].legend(loc=3, bbox_to_anchor=(0.,1.1,1.,0.1), ncol=4, mode="expand", borderaxespad=0., fontsize=12)

    def save_figure(fig, specie):
        savepath = escape_fname(os.path.join(plotdir, '03.XCLASS_fit', 'SSC_'+str(SSC['no']), specie+'.hist.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    moldata, percentiles = parse_molfit(SSC, specie, return_all=True)
    fig,[ax_T, ax_N, ax_lw, ax_v] = set_up_figure(SSC, specie)
    plot_hist(ax_T, moldata['temperature'], r'T [K]')
    plot_hist(ax_N, moldata['column density'], r'$\Sigma$ [cm$^{-2}$]')
    plot_hist(ax_lw, moldata['linewidth'], r'$\sigma$ [km\,s$^{-1}$]')
    plot_hist(ax_v, moldata['velocity'], r'v$_\mathrm{rel}$ [km\,s${-1}$]')
    format_figure([ax_T, ax_N, ax_lw, ax_v])
    save_figure(fig, specie)


###################################################################################################
#
###################################################################################################
