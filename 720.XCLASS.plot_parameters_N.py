###################################################################################################
# plot fit parameters
###################################################################################################

def plot_parameters_N(SSC):
    """
    Make a box-and-wiskers plot of all detected species for a given SSC.
    """

    for Q,ylab,ymin,ymax in [['temperature','T [K]',0,250],['column density','N [cm$^{-2}$]',1e12,1e20],['linewidth','$\sigma$ [km\,s$^{-1}$]',0,90],['velocity','v [km\,s$^{-1}$]',-50,50]]:

        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
        ax.text(0.05, 0.9, 'SSC '+str(SSC['no'])+': '+Q, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

        colors = [plt.cm.inferno(i/(len(fitable_species)+1)) for i,_ in enumerate(fitable_species)]
        for idx,(specie,c) in enumerate(zip(fitable_species, colors)):
            if specie in line_data_N[str(SSC['no'])].keys():
                
                # max_comp     = np.argmax(line_data_N[str(SSC['no'])][specie]['column density']['median'])
                nearest_comp = np.argmin(np.abs(line_data_N[str(SSC['no'])][specie]['velocity']['median']))
                q_med  = line_data_N[str(SSC['no'])][specie][Q]['median'][nearest_comp]
                q_16th = line_data_N[str(SSC['no'])][specie][Q]['16th'][nearest_comp]
                q_84th = line_data_N[str(SSC['no'])][specie][Q]['84th'][nearest_comp]
                q_all  = line_data_N[str(SSC['no'])][specie][Q]['all'][nearest_comp]
                boxplot = ax.boxplot(q_all, patch_artist=True, positions=[idx], showfliers=False)
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
                    flier.set(marker='o', color=c, alpha=0.5)


        ax.set_xlim(-1, len(fitable_species))
        ax.set_ylim(ymin,ymax)
        ax.set_axisbelow(True)
        ax.grid(axis='y')
        ax.set_xticks(np.arange(len(fitable_species)))
        ax.set_xticklabels([specie_tex(f) for f in fitable_species])
        ax.tick_params(axis='x', rotation=90)
        if Q=='column density':
            ax.set_yscale('log')
        ax.set_ylabel(ylab, fontsize=12)
        fig.tight_layout()

        savepath = os.path.join(plotdir, '04.fit_results', 'parameters_N', 'SSC_'+str(SSC['no'])+'.'+Q+'.pdf')
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')
