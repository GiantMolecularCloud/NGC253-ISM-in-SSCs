###################################################################################################
# plot ratios
###################################################################################################

def plot_ratios_N(rname, rdict):
    """
    Make a box-and-wiskers plot of a line ratio in all SSCs.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    sscs   = [SSC['no'] for SSC in SSCs]
    r_all  = ratios_N[rname]['all']
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for s,r,c in zip(sscs,r_all,colors):
        try:
            boxplot = ax.boxplot(r, patch_artist=True, positions=[s], showfliers=False)
            for box in boxplot['boxes']:
                box.set(color=c, linewidth=1)
                box.set(facecolor=mpl.colors.to_rgba(c, alpha=0.5))
            for median in boxplot['medians']:
                median.set(color=c, linewidth=3)
            for whisker in boxplot['whiskers']:
                whisker.set(color=c, linewidth=3)
            for cap in boxplot['caps']:
                cap.set(color=c, linewidth=3)
            for flier in boxplot['fliers']:
                flier.set(marker='.', size=8, color=c, alpha=0.5)
        except:
            print(rname+": no boxplot in SSC "+str(s))

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_yscale('log')
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'N ('+specie_tex(rdict['a'])+r') / N ('+specie_tex(rdict['b'])+r')', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'ratios_N', 'ratio_'+rname.replace('/','-')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')
