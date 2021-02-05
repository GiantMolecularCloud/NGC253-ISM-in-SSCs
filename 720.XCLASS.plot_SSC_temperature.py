###################################################################################################
# plot SSC temperatures
###################################################################################################

def plot_SSC_temperatures(T_specie):
    """
    Make a box-and-wiskers plot of the XLCASS temperature in a given SSC.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, specie_tex(T_specie), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        if T_specie in detected_species[str(SSC['no'])]:
            try:
                temperatures = line_data_N[str(SSC['no'])][T_specie]['temperature']['all']
                boxplot = ax.boxplot(temperatures, patch_artist=True, positions=[SSC['no']], showfliers=False)
                for box in boxplot['boxes']:
                    box.set(color=color, linewidth=1)
                    box.set(facecolor=mpl.colors.to_rgba(color, alpha=0.5))
                for median in boxplot['medians']:
                    median.set(color=color, linewidth=3)
                for whisker in boxplot['whiskers']:
                    whisker.set(color=color, linewidth=3)
                for cap in boxplot['caps']:
                    cap.set(color=color, linewidth=3)
                for flier in boxplot['fliers']:
                    flier.set(marker='.', size=8, color=color, alpha=0.5)
            except:
                print('SSC '+str(SSC['no'])+', '+T_specie+': no boxplot in SSC')

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'T$_\mathrm{rot}$ [K]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'SSC_temperatures', 'temperature.'+T_specie+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')
