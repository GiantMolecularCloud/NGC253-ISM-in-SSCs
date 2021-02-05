############################
# GAS IN SSCS: final stuff #
############################

###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

data_XCLASS   = fnunpickle(os.path.join(resultsdir, 'data_XCLASS.pickle'))
data_Gauss    = fnunpickle(os.path.join(resultsdir, 'data_Gauss.pickle'))
ratios_XCLASS = fnunpickle(os.path.join(resultsdir, 'ratios_XCLASS.pickle'))
ratios_Gauss  = fnunpickle(os.path.join(resultsdir, 'ratios_Gauss.pickle'))


###################################################################################################
# compare parameters
###################################################################################################

def plot_parameters(SSC):

    for q,ylab,ymin,ymax in [['column density','N [cm$^{-2}$]',1e11,1e21],['linewidth','$\sigma$ [km\,s$^{-1}$]',0,90],['velocity','v [km\,s$^{-1}$]',-50,50]]:

        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(8,6))
        ax.text(0.05, 0.9, 'SSC '+SSC['num']+': '+q, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

        # fitted_species = data_XCLASS[SSC['num']].keys()
        colors = [plt.cm.inferno(i/(len(fitable_species)+1)) for i,_ in enumerate(fitable_species)]

        for idx,(specie,c) in enumerate(zip(fitable_species, colors)):
            try:
                ncomps = len(data_XCLASS[SSC['num']][specie][q]['all'])
                for comp,q_all in enumerate(data_XCLASS[SSC['num']][specie][q]['all']):
                    spx_shift = 0.6/ncomps*(comp-(ncomps-1)/2.)
                    boxplot = ax.boxplot(q_all, patch_artist=True, positions=[idx+spx_shift], showfliers=False)
                    for box in boxplot['boxes']:
                        box.set(color=c, linewidth=1)
                        box.set(facecolor=mpl.colors.to_rgba(c, alpha=0.5))
                    for median in boxplot['medians']:
                        median.set(color=c, linewidth=2)
                    for whisker in boxplot['whiskers']:
                        whisker.set(color=c, linewidth=2)
                    for cap in boxplot['caps']:
                        cap.set(color=c, linewidth=2)
                    for flier in boxplot['fliers']:
                        flier.set(marker='o', color=c, alpha=0.5)
            except:
                pass

        # plot stripes for clarity
        for idx,spx in enumerate(fitable_species):
            if idx%2==0:
                ax.axvspan(idx-0.5, idx+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
            else:
                ax.axvspan(idx-0.5, idx+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

        ax.set_xlim(-0.5, len(fitable_species)-0.5)
        ax.set_ylim(ymin,ymax)
        ax.set_axisbelow(True)
        ax.grid(axis='y')
        ax.set_xticks(np.arange(len(fitable_species)))
        ax.set_xticklabels([specie_tex(f) for f in fitable_species])
        ax.tick_params(axis='x', rotation=90)
        if q=='column density':
            ax.set_yscale('log')
        ax.set_ylabel(ylab, fontsize=12)
        fig.tight_layout()

        savepath = os.path.join(plotdir, '10.results', 'parameters', 'SSC_'+SSC['num']+'.'+q+'.pdf')
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')


for SSC in tqdm(SSCs):
    plot_parameters(SSC)


###################################################################################################
# plot histogram of fit parameters
###################################################################################################

def plot_fit_hist(SSC, spx):

    def set_up_figure(SSC, spx):
        fig,axes = plt.subplots(nrows=5, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,12))
        fig.subplots_adjust(hspace=0.5)
        axes[0].text(0.05, 0.9, 'SSC '+SSC['num']+': '+specie_tex(spx), color='k', transform=axes[0].transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)
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

    def save_figure(fig, spx):
        savepath = escape_fname(os.path.join(plotdir, '10.results', 'histograms','SSC_'+SSC['num']+'.'+spx+'.histogram.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')

    temperatures     = data_XCLASS[SSC['num']][spx]['temperature']['all']
    column_densities = data_XCLASS[SSC['num']][spx]['column density']['all']
    linewidths       = data_XCLASS[SSC['num']][spx]['linewidth']['all']
    velocities       = data_XCLASS[SSC['num']][spx]['velocity']['all']
    opacities        = data_XCLASS[SSC['num']][spx]['integrated opacity']['all']
    fig,[ax_T, ax_N, ax_lw, ax_v, ax_o] = set_up_figure(SSC, spx)
    # fig,[ax_N, ax_lw, ax_v, ax_o] = set_up_figure(SSC, spx)
    plot_hist(ax_T, temperatures, r'T [K]')
    plot_hist(ax_N, column_densities, r'$\Sigma$ [cm$^{-2}$]')
    plot_hist(ax_lw, linewidths, r'$\sigma$ [km\,s$^{-1}$]')
    plot_hist(ax_v, velocities, r'v$_\mathrm{rel}$ [km\,s${-1}$]')
    plot_hist(ax_o, opacities, r'$\tau_\mathrm{int}$')
    format_figure([ax_T, ax_N, ax_lw, ax_v, ax_o])
    # format_figure([ax_N, ax_lw, ax_v, ax_o])
    save_figure(fig, spx)


for SSC in tqdm(SSCs):
    for spx in data_XCLASS[SSC['num']].keys():
        try:
            plot_fit_hist(SSC, spx)
        except:
            print(SSC['num'], spx)


###################################################################################################
# compare SO to CS opacity
###################################################################################################

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

for SSC in SSCs:
    c = colors[SSC['no']-1]
    try:
        tau_CS = data_XCLASS[SSC['num']]['CS;v=0']['integrated opacity']['all'][0]
        tau_SO = data_XCLASS[SSC['num']]['SO;v=0;#1']['integrated opacity']['all'][0]
        log_tau_ratio = np.log10(np.array(tau_CS)/np.array(tau_SO))
        boxplot = ax.boxplot(log_tau_ratio, patch_artist=True, positions=[SSC['no']], showfliers=False)
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
        ax.plot(SSC['no'], 0.85, marker='^', ms=8, lw=0, c=c)

ax.plot([0.5,14.5],[0,0], lw=1, ls='-', c='darkgrey', zorder=1)

# plot stripes for clarity
for idx,SSC in enumerate(SSCs):
    if idx%2==0:
        ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
    else:
        ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

ax.set_xlim(0.5, len(SSCs)+0.5)
ax.set_ylim(-0.9, 0.9)
ax.set_xticks(np.arange(1,len(SSCs)+1))
ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
ax.set_axisbelow(True)
ax.grid(axis='y')
ax.set_ylabel(r'$\log \frac{\tau\ (\mathrm{CS})}{\tau\ (\mathrm{SO})}$', fontsize=12)
ax.set_xlabel('SSC')
fig.tight_layout()

savepath = escape_fname(os.path.join(plotdir, '10.results', 'compare', 'opacity_CS_SO.pdf'))
os.system('mkdir -p '+os.path.dirname(savepath))
fig.savefig(savepath, dpi=300, bbox_inches='tight')


###################################################################################################
# component deviations
###################################################################################################

# line centroid
lps = {SSC['num']: {} for SSC in SSCs}
for s,data in data_XCLASS.items():
    for spx,dat in data.items():
        if len(dat['velocity']['median']) == 1:
            lps[s][spx] = dat['velocity']['median']


# linewidth
lws = {SSC['num']: {} for SSC in SSCs}
for s,data in data_XCLASS.items():
    for spx,dat in data.items():
        if len(dat['linewidth']['median']) == 1:
            lws[s][spx] = dat['linewidth']['median']

for SSC in SSCs:
    lw_mean = np.mean([lw for spx,lw in lws[SSC['num']].items()])
    print(SSC['num']+'\t'+'{:.2f}'.format(lw_mean))


###################################################################################################
# mean parameters
###################################################################################################

def plot_median_parameter(quantity):

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(8,6))
    colors = [plt.cm.inferno(i/(len(fitable_species)+1)) for i,_ in enumerate(fitable_species)]

    for idx,(specie,c) in enumerate(zip(fitable_species, colors)):
        values = []
        for SSC in SSCs:
            try:
                nearest_comp = np.argmin(np.abs(data_XCLASS[str(SSC['no'])][specie]['velocity']['median']))
                median = data_XCLASS[str(SSC['no'])][specie][quantity]['median'][nearest_comp]
                values.append(median)
            except:
                pass
        if not values==[]:
            values = np.array(values)
            lw16, lw, lw84 = np.percentile(values, (16,50,84))
            ax.errorbar(idx, lw, yerr=[[lw-lw16],[lw84-lw]], ls='', marker='o', ms=6, color=c, elinewidth=2, ecolor=c)
            ax.text((idx+0.5)/len(fitable_species), 0.95, len(values), color='k', transform=ax.transAxes, ha='center', va='center', weight='bold', fontsize=12) #, bbox=props)

    # plot stripes for clarity
    for idx,spx in enumerate(fitable_species):
        if idx%2==0:
            ax.axvspan(idx-0.5, idx+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(idx-0.5, idx+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(-0.5, len(fitable_species)-0.5)
    # ax.set_ylim(ymin,ymax)
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_xticks(np.arange(len(fitable_species)))
    ax.set_xticklabels([specie_tex(f) for f in fitable_species])
    ax.tick_params(axis='x', rotation=90)
    if quantity=='column density':
        ax.set_yscale('log')
    if quantity=='linewidth':
        ylab = r'$\sigma$ [km\,s$^{-1}$]'
    ax.set_ylabel(ylab, fontsize=12)
    fig.tight_layout()

    savepath = os.path.join(plotdir, '10.results', 'compare', 'median_'+quantity.replace(' ','_')+'.pdf')
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


plot_median_parameter('linewidth')


###################################################################################################

def compare_parameter(species, quantity):

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(8,6))
    # colors = [plt.cm.inferno(i/(len(species)+1)) for i,_ in enumerate(species)]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    def isplit(iterable,splitters):
        import itertools
        return [list(g) for k,g in itertools.groupby(iterable,lambda x:x in splitters) if not k]

    for s,SSC in enumerate(SSCs):
        c = colors[s]
        blocks = isplit(species,'')
        for b,block in enumerate(blocks):
            idxs   = []
            median = []
            p16    = []
            p84    = []
            for idx,specie in enumerate(block):
                idx = idx+sum([len(blocks[bb])+1 for bb in np.arange(b)])
                try:
                    nearest_comp = np.argmin(np.abs(data_XCLASS[str(SSC['no'])][specie]['velocity']['median']))
                    median.append( data_XCLASS[str(SSC['no'])][specie][quantity]['median'][nearest_comp] )
                    p16.append( data_XCLASS[str(SSC['no'])][specie][quantity]['16th'][nearest_comp] )
                    p84.append( data_XCLASS[str(SSC['no'])][specie][quantity]['84th'][nearest_comp] )
                    idxs.append( idx +0.6/len(SSCs)*(s-(len(SSCs)-1)/2.) )
                except:
                    pass
            if not median==[]:
                median = np.array(median)
                p16 = np.array(p16)
                p84 = np.array(p84)
                ax.errorbar(idxs, median, yerr=[median-p16,p84-median], ls='-', lw=0.5, marker='o', ms=3, color=c, elinewidth=1, ecolor=c, label=SSC['num'] if b==0 else '')
                # ax.plot(idxs, median, ls='-', marker='o', ms=6, color=c)

    # plot stripes for clarity
    for idx,spx in enumerate(fitable_species):
        if idx%2==0:
            ax.axvspan(idx-0.5, idx+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(idx-0.5, idx+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(-0.5, len(species)-0.5)
    # ax.set_ylim(ymin,ymax)
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_xticks(np.arange(len(species)))
    ax.set_xticklabels([specie_tex(f) for f in species])
    ax.tick_params(axis='x', rotation=90)
    if quantity=='column density':
        ax.set_yscale('log')
    if quantity=='linewidth':
        ylab = r'$\sigma$ [km\,s$^{-1}$]'
    ax.set_ylabel(ylab, fontsize=12)
    ax.legend(loc=3, bbox_to_anchor=(0.,1.1,1.,0.1), ncol=14, mode="expand", borderaxespad=0., fontsize=12)
    fig.tight_layout()

    savepath = os.path.join(plotdir, '10.results', 'compare', 'compare_linewidths.pdf')
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_parameter(species  = ['HCN;v=0','HC-13-N;v=0','HCN-15;v=0','HN-15-C;v=0','','HCN;v=0','HCN;v2=1','HCN;v2=2','','HCCCN;v=0','HCCCN;v6=1','HCCCN;v7=1','HCCCN;v7=2'],
                  quantity = 'linewidth'
                 )


###################################################################################################
# compare line centroids
###################################################################################################

def compare_velocity(species):
    """
    Compare line centroids of the given species and Leroy+18 values.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))

    colors = [plt.cm.inferno(i/(len(species)+0.2)) for i in np.arange(len(species)+1)]

    # my measurements
    for idx,specie in enumerate(species):
        spx_shift = -0.3+idx*0.6/len(species)

        SSC_no = [SSC['no'] for SSC in SSCs if specie in data_XCLASS[str(SSC['no'])]]

        med,low,upp = [],[],[]
        for SSC in SSCs:
            nearest_comp = np.argmin(np.abs(data_XCLASS[str(SSC['no'])][specie]['velocity']['median']))
            median = data_XCLASS[str(SSC['no'])][specie]['velocity']['median'][nearest_comp]
            p16    = data_XCLASS[str(SSC['no'])][specie]['velocity']['16th'][nearest_comp]
            p84    = data_XCLASS[str(SSC['no'])][specie]['velocity']['84th'][nearest_comp]
            med.append(median)
            low.append(median-p16)
            upp.append(p84-median)
        ax.errorbar([s+spx_shift for s in SSC_no], med, yerr=[low,upp], ls='', marker='o', ms=6, color=colors[idx], elinewidth=2, ecolor=colors[idx], label=specie_tex(specie))
        # ax.plot([s+spx_shift for s in SSC_no], med, ls='', marker='o', ms=6, color=colors[idx], label=specie_tex(specie))

    # Leroy measurements
    SSC_no    = [SSC['no'] for SSC in SSCs]
    SSC_velo  = [SSC['XXX'].value for SSC in SSCs]
    SSC_err   = [SSC['XXX_v'].value for SSC in SSCs]
    ax.errorbar([s+0.3 for s in SSC_no], SSC_velo, yerr=SSC_err, ls='', marker='_', ms=6, color=colors[-1], elinewidth=2, ecolor=colors[-1], label='Leroy+18')

    # plot stripes for clarity
    for idx,SSC in enumerate(SSCs):
        if idx%2==0:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, color='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(0.5, len(SSCs)+0.5)
    ax.set_ylim(-150,150)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'v [km\,s$^{-1}$]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.legend(loc=2, bbox_to_anchor=(0.02, 0.975), bbox_transform=ax.transAxes)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'compare', 'velocity_Leroy18.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_velocity(['CS;v=0','HC-13-N;v=0','CO;v=0','HCN;v=0','HCO+;v=0'])


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

        SSC_no = [SSC['no'] for SSC in SSCs if specie in data_XCLASS[str(SSC['no'])]]

        med,low,upp = [],[],[]
        for SSC in SSCs:
            nearest_comp = np.argmin(np.abs(data_XCLASS[str(SSC['no'])][specie]['velocity']['median']))
            median = data_XCLASS[str(SSC['no'])][specie]['linewidth']['median'][nearest_comp]
            p16    = data_XCLASS[str(SSC['no'])][specie]['linewidth']['16th'][nearest_comp]
            p84    = data_XCLASS[str(SSC['no'])][specie]['linewidth']['84th'][nearest_comp]
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

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'compare', 'linewidth_Leroy18.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_linewidths(['CS;v=0','HC-13-N;v=0','CO;v=0','HCN;v=0','HCO+;v=0'])


###################################################################################################
# opacity table
###################################################################################################

def opacity_table(species):
    # header
    header = r'\colhead{SSC} & '
    for specie in species:
        header += r'\colhead{'+specie_tex(specie)+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # column density ratios
    for SSC in SSCs:
        row = SSC['num'] +' & '
        for specie in species:
            try:
                nearest = np.argmin(np.abs(data_XCLASS[SSC['num']][specie]['velocity']['median']))
                median = data_XCLASS[SSC['num']][specie]['integrated opacity']['median'][nearest]
                p16    = data_XCLASS[SSC['num']][specie]['integrated opacity']['16th'][nearest]
                p84    = data_XCLASS[SSC['num']][specie]['integrated opacity']['84th'][nearest]
                upper  = p84-median
                lower  = median-p16
                # row += '{:#.2g}'.format(median)+r'$^{+'+'{:#.2g}'.format(upper)+'}_{-'+'{:#.2g}'.format(lower)+'}$ & '
                if median>10:
                    fmt = '{:2.1f}'
                else:
                    fmt = '{:1.2f}'
                row += fmt.format(median)+r'$^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
            except:
                row += ' & '
        row = row[:-3]
        row += r'\\'
        print(row)


opacity_table(['CO;v=0','HCO+;v=0','HCN;v=0','HC-13-N;v=0','HCN-15;v=0','HN-15-C;v=0'])


###################################################################################################
# plot opacity vs. opacity tracer ratio
###################################################################################################

def plot_opacity_check():
    """
    Plot XCLASS opacity vs. HCN/H13CN ratio as proxy for opacity.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    # ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for SSC,color in zip(SSCs,colors):
        nearest = np.argmin(np.abs(data_XCLASS[SSC['num']]['HCN;v=0']['velocity']['median']))
        opt     = data_XCLASS[SSC['num']]['HCN;v=0']['integrated opacity']['median'][nearest]
        opt16   = data_XCLASS[SSC['num']]['HCN;v=0']['integrated opacity']['16th'][nearest]
        opt84   = data_XCLASS[SSC['num']]['HCN;v=0']['integrated opacity']['84th'][nearest]
        opt_l   = opt-opt16
        opt_u   = opt84-opt
        ratio   = ratios_Gauss['HCN/HCNthin'][SSC['num']]['bestfit']
        ratio_e = ratios_Gauss['HCN/HCNthin'][SSC['num']]['error']
        ax.errorbar(opt,ratio,
                    xerr   = [[opt_l],[opt_u]],
                    yerr   = [[ratio_e],[ratio_e]],
                    marker = 'o', ms=5, color=color, elinewidth=2, ecolor=color, ls='',
                    label  = SSC['num']+r'\\',
                    zorder = 4
                   )

    ax.set_xlim(0.2,20)
    ax.set_ylim(0.7,200)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_ylabel(r'I (HCN) / I (H$^{13}$CN)', fontsize=12)
    ax.set_xlabel(r'$\tau$ (HCN)', fontsize=12)
    legend = ax.legend(loc=10, bbox_to_anchor=(0.,1.05,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12) #, frameon=False)
    for txt in legend.get_texts():
        txt.set_ha('center')
        txt.set_x(-22)
        txt.set_y(-15)
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'compare', 'opacity-HCN_ratio.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


plot_opacity_check()


###################################################################################################
# temperature
###################################################################################################

def compare_temperatures():

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    # colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    colors = ['orangered','darkorange']

    for SSC in SSCs:
        for a,(specie,color) in enumerate(zip(temperature_species,colors)):
            if specie in data_XCLASS[SSC['num']].keys():
                nearest = np.argmin(np.abs(data_XCLASS[SSC['num']][specie]['velocity']['median']))
                temp = data_XCLASS[SSC['num']][specie]['temperature']['all'][nearest]

                if np.nanmedian(temp)<250:
                    boxplot = ax.boxplot(temp, patch_artist=True, positions=[SSC['no']-0.2+a*0.4], showfliers=False, zorder=2)
                    for box in boxplot['boxes']:
                        box.set(color=color, linewidth=1)
                        box.set(facecolor=mpl.colors.to_rgba(color, alpha=0.5)) #+a*0.5))
                    for median in boxplot['medians']:
                        median.set(color=color, linewidth=3)
                    for whisker in boxplot['whiskers']:
                        whisker.set(color=color, linewidth=3)
                    for cap in boxplot['caps']:
                        cap.set(color=color, linewidth=3)
                    for flier in boxplot['fliers']:
                        flier.set(marker='.', size=8, color=color, alpha=0.5) #+a*0.5)
                else:
                    lowlim = np.nanpercentile(temp, 16)
                    ax.scatter([SSC['no']-0.2+a*0.4], [lowlim], s=40, c=color, marker='^')

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    tracers = [Line2D([0], [0],
                      marker='s', color=color, markerfacecolor=color, markersize=10,
                      label=specie_tex(specie))
               for a,(specie,color) in enumerate(zip(temperature_species,colors))]
    ax.legend(handles=tracers, loc=1, ncol=2) #loc=3, bbox_to_anchor=(0.,0.97,1.,0.05), ncol=14, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    # plot stripes for clarity
    for idx,SSC in enumerate(SSCs):
        if idx%2==0:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='lightgrey', lw=0, alpha=0.5, zorder=0)
        else:
            ax.axvspan(SSC['no']-0.5, SSC['no']+0.5, facecolor='darkgrey', lw=0, alpha=0.5, zorder=0)

    ax.set_xlim(0.5, len(SSCs)+0.5)
    ax.set_ylim(0,300)
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.set_axisbelow(True)
    ax.grid(axis='y', which='both')
    ax.set_ylabel(r'T$_\mathrm{rot}$ [K]', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'compare', 'compare_temperatures.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


compare_temperatures()


###################################################################################################
#
###################################################################################################
