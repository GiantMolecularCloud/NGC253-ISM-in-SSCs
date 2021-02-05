#######################
# GAS IN SSCS: XCLASS #
#######################

# Temperature fitting
# Keep other lines fixed and fit only for temperature sensitive lines.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))

nums      = fnunpickle(os.path.join(tempdir, 'nums.pickle'))
tempfiles = fnunpickle(os.path.join(tempdir, 'tempfiles.pickle'))
temperature_species = fnunpickle(os.path.join(tempdir, 'temperature_species.pickle'))
temperature_data    = fnunpickle(os.path.join(tempdir, 'temperature_data.pickle'))
data_N   = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))


###################################################################################################
# compare fixed-free temperature
###################################################################################################

def compare_quantity(specie, quantity):
    """
    Compare a quantity (column density, opacity, linewidth) between fixed and free temperature runs.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(8,8))

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        if specie in temperature_data[SSC['num']].keys():
            nearest = np.argmin(np.abs(temperature_data[SSC['num']][specie]['velocity']['median']))
            freeq     = temperature_data[SSC['num']][specie][quantity]['median'][nearest]
            freeq16   = temperature_data[SSC['num']][specie][quantity]['16th'][nearest]
            freeq84   = temperature_data[SSC['num']][specie][quantity]['84th'][nearest]
            freeq_l   = freeq-freeq16
            freeq_u   = freeq84-freeq
            fixedq     = data_N[SSC['num']][specie][quantity]['median'][nearest]
            fixedq16   = data_N[SSC['num']][specie][quantity]['16th'][nearest]
            fixedq84   = data_N[SSC['num']][specie][quantity]['84th'][nearest]
            fixedq_l   = fixedq-fixedq16
            fixedq_u   = fixedq84-fixedq
            ax.errorbar(fixedq,freeq,
                        xerr   = [[fixedq_l],[fixedq_u]],
                        yerr   = [[freeq_l],[freeq_u]],
                        marker = 'o', ms=5, color=color, elinewidth=1, ecolor=color, ls='',
                        label  = SSC['num'],
                        zorder = 4
                       )

    if quantity in ['column density', 'integrated opacity', 'peak opacity']:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_xlim())
    ax.set_aspect('equal', 'box')
    ax.plot(ax.get_xlim(),ax.get_xlim(), ls='--', color='lightgrey', zorder=1)
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_xlabel(quantity+' ($\mathrm{T}=130$\,K, fixed)', fontsize=12)
    ax.set_ylabel(quantity+' (T fitted)', fontsize=12)
    fig.legend(loc=3, bbox_to_anchor=(0.,0.97,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '05.temperatures', 'comparison', specie+'.compare_'+quantity.replace(' ','_')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


for specie in tqdm(temperature_species):
    for quantity in ['integrated opacity','column density','linewidth']:
        compare_quantity(specie, quantity)


###################################################################################################
# compare temperature estimates
###################################################################################################

def compare_temperatures():

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    # colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    colors = ['orangered','darkorange']

    for SSC in SSCs:
        for a,(specie,color) in enumerate(zip(temperature_species,colors)):
            if specie in temperature_data[SSC['num']].keys():
                nearest = np.argmin(np.abs(temperature_data[SSC['num']][specie]['velocity']['median']))
                temp = temperature_data[SSC['num']][specie]['temperature']['all'][nearest]

                boxplot = ax.boxplot(temp, patch_artist=True, positions=[SSC['no']-0.2+a*0.4], showfliers=False, zorder=2)
                for box in boxplot['boxes']:
                    box.set(color=color, linewidth=1)
                    box.set(facecolor=mpl.colors.to_rgba(color, alpha=0.5+a*0.5))
                for median in boxplot['medians']:
                    median.set(color=color, linewidth=3)
                for whisker in boxplot['whiskers']:
                    whisker.set(color=color, linewidth=3)
                for cap in boxplot['caps']:
                    cap.set(color=color, linewidth=3)
                for flier in boxplot['fliers']:
                    flier.set(marker='.', size=8, color=color, alpha=0.5+a*0.5)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    tracers = [Line2D([0], [0],
                      marker='s', color=color, markerfacecolor=color, markersize=10,
                      label=specie_tex(specie))
               for a,(specie,color) in enumerate(zip(temperature_species,colors))]
    ax.legend(handles=tracers, loc=1, ncol=2) #loc=3, bbox_to_anchor=(0.,0.97,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    # plot upper limits for fit range
    ulim = {'H2CS;v=0': 1000, 'SO2;v=0': 500}
    ax.plot([0.5, len(SSCs)+0.5], [500,500],   ls='--', color='darkgrey', zorder=1)
    ax.plot([0.5, len(SSCs)+0.5], [1000,1000], ls='--', color='darkgrey', zorder=1)
    ax.text(0.6, 500, 'limit SO$_2$', color='darkgrey', ha='left', va='bottom', weight='bold', fontsize=12)
    ax.text(0.6, 1000, 'limit H$_2$CS', color='darkgrey', ha='left', va='bottom', weight='bold', fontsize=12)

    # plot stripes for clarity
    for idx,SSC in enumerate(SSCs):
        if idx%2==0:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(0.5, len(SSCs)+0.5)
    ax.set_ylim(0,1050)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.set_axisbelow(True)
    ax.grid(axis='y', which='both')
    ax.set_ylabel(r'T$_\mathrm{rot}$ [K]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '05.temperatures', 'comparison', 'compare_temperatures.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_temperatures()


###################################################################################################
# temperature table
###################################################################################################

def temperature_table():
    header = r'\colhead{SSC} & '
    for specie in temperature_species:
        header += r'\colhead{'+specie_tex(specie)+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # temperatures
    for SSC in SSCs:
        row = str(SSC['no'])+' & '
        for specie in temperature_species:
            try:
                nearest = np.argmin(np.abs(temperature_data[SSC['num']][specie]['velocity']['median']))
                median = temperature_data[SSC['num']][specie]['temperature']['median'][nearest]
            except:
                median = np.nan
            if np.isnan(median):
                row += ' & '
            else:
                p16    = temperature_data[SSC['num']][specie]['temperature']['16th'][nearest]
                p84    = temperature_data[SSC['num']][specie]['temperature']['84th'][nearest]
                upper  = p84-median
                lower  = median-p16
                row += '{:.0f}'.format(median)+r'$^{+'+'{:.0f}'.format(upper)+'}_{-'+'{:.0f}'.format(lower)+'}$ & '
        row = row[:-3]
        row += r'\\'
        print(row)


temperature_table()


###################################################################################################
#
###################################################################################################
