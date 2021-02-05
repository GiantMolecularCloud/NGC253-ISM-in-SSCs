###################################################################################################
# plot SSC values
###################################################################################################

def plot_SSC_quantity_N(specie, quantity):
    """
    Make a box-and-wiskers plot of an XLCASS quantity in a given SSC.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, specie_tex(specie), color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        if specie in detected_species[str(SSC['no'])]:
            try:
                nearest_comp = np.argmin(np.abs(line_data_N[str(SSC['no'])][specie]['velocity']['median']))
                values = line_data_N[str(SSC['no'])][specie][quantity]['all'][nearest_comp]
                if quantity=='velocity':
                    values = np.array(values) +SSC['velshift'].value +250
                boxplot = ax.boxplot(values, patch_artist=True, positions=[SSC['no']], showfliers=False)
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
                print('SSC '+str(SSC['no'])+', '+specie+' '+quantity+': no boxplot in SSC')

    ax.set_xlim(0, len(SSCs)+1)
    if quantity=='column density':
        ax.set_yscale('log')
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    if quantity=='temperature':
        label = r'T$_\mathrm{rot}$ [K]'
    elif quantity=='column density':
        label = r'N [cm$^{-2}$]'
    elif quantity=='linewidth':
        label = r'$\sigma$ [km\,s$^{-2}$]'
    elif quantity=='velocity':
        label = r'v [km\,s$^{-2}$]'
        # ax.avhline(250)
    ax.set_ylabel(label, fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'SSC_'+quantity.replace(' ','_')+'_N', quantity.replace(' ','_')+'.'+specie+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')




def compare_linewidths_N(species):
    """
    Compare linewidths of the given species and Leroy+18 values.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))

    colors = [plt.cm.inferno(i/(len(species)+0.2)) for i in np.arange(len(species)+1)]

    # my measurements
    for idx,specie in enumerate(species):
        SSC_no = [SSC['no'] for SSC in SSCs if specie in line_data_N[str(SSC['no'])]]

        med,low,upp = [],[],[]
        for SSC in SSCs:
            nearest_comp = np.argmin(np.abs(line_data_N[str(SSC['no'])][specie]['velocity']['median']))
            med.append(line_data_N[str(SSC['no'])][specie]['linewidth']['median'][nearest_comp])
            low.append(line_data_N[str(SSC['no'])][specie]['linewidth']['16th'][nearest_comp])
            upp.append(line_data_N[str(SSC['no'])][specie]['linewidth']['84th'][nearest_comp])
        # ax.errorbar(SSC_no, med, yerr=[low,upp], ls='', marker='o', ms=6, color=colors[idx], elinewidth=2, ecolor=colors[idx], label=specie_tex(specie))
        ax.plot(SSC_no, med, ls='', marker='o', ms=6, color=colors[idx], label=specie_tex(specie))

    # Leroy measurements
    SSC_no    = [SSC['no'] for SSC in SSCs]
    SSC_sigma = [SSC['sigma_v'].value for SSC in SSCs]
    SSC_err   = [SSC['err_sigma_v'].value for SSC in SSCs]
    ax.errorbar(SSC_no, SSC_sigma, yerr=SSC_err, ls='', marker='_', ms=6, color=colors[-1], elinewidth=2, ecolor=colors[-1], label='Leroy+18')

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_ylim(0,100)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'$\sigma$ [km\,s$^{-1}$]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.legend(loc=2, bbox_to_anchor=(0.02, 0.975), bbox_transform=ax.transAxes)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'linewidth_comparison_N.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')
