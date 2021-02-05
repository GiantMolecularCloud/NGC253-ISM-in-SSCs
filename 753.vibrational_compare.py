#######################
# GAS IN SSCS: XCLASS #
#######################

# Fitting of ro-vibational lines
# Keep other lines fixed and fit only for ro-vibrational lines.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))

nums      = fnunpickle(os.path.join(vibdir, 'nums.pickle'))
vibfiles = fnunpickle(os.path.join(vibdir, 'vibfiles.pickle'))
vibrational_species = fnunpickle(os.path.join(vibdir, 'vibrational_species.pickle'))
vibrational_data = fnunpickle(os.path.join(vibdir, 'vibrational_data.pickle'))
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
        if specie in vibrational_data[SSC['num']].keys():
            nearest = np.argmin(np.abs(vibrational_data[SSC['num']][specie]['velocity']['median']))
            freeq     = vibrational_data[SSC['num']][specie][quantity]['median'][nearest]
            freeq16   = vibrational_data[SSC['num']][specie][quantity]['16th'][nearest]
            freeq84   = vibrational_data[SSC['num']][specie][quantity]['84th'][nearest]
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
    ax.set_ylabel(quantity+' ($\mathrm{T}=300$\,K, fixed)', fontsize=12)
    fig.legend(loc=3, bbox_to_anchor=(0.,0.97,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '06.vibrations', 'comparison', specie+'.compare_'+quantity.replace(' ','_')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


for specie in tqdm(vibrational_species):
    for quantity in ['integrated opacity','column density','linewidth']:
        compare_quantity(specie, quantity)


###################################################################################################
#
###################################################################################################
