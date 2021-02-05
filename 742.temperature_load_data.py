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


###################################################################################################
# load fitted data
###################################################################################################

def parse_molfit(SSC):

    import re

    def load_molfit(num, SSC):
        file = tempfiles[SSC['num']][str(num)]['model']['molfit']
        with open(file, 'r') as f:
            fdata = f.read()
            return fdata

    def filter_vals():
        # # consider fit good only if
        # #       temperature not within 2K of limits
        # #       column density not within 5% of limits
        # # ignore this selection for fixed fits
        # good_fit  = (t_val>t_llim+2 and t_val<t_ulim+2 and N_val>N_llim*1.05 and N_val<N_ulim*0.95 and N_val>2e12)
        # fixed_fit = (N_ulim/N_llim<5)
        # if good_fit or fixed_fit:
        data[spx]['temperature']['all'][ncomp].append(t_val)
        data[spx]['column density']['all'][ncomp].append(N_val)
        data[spx]['linewidth']['all'][ncomp].append(w_val)
        data[spx]['velocity']['all'][ncomp].append(v_val)
        # else:
        #     data[spx]['temperature']['all'][ncomp].append(np.nan)
        #     data[spx]['column density']['all'][ncomp].append(np.nan)
        #     data[spx]['linewidth']['all'][ncomp].append(np.nan)
        #     data[spx]['velocity']['all'][ncomp].append(np.nan)

    def split_at_molecule(mofit):
        return re.split(u'\n(?=[A-Z])', molfit)

    # load model molfit files
    molfits = [ load_molfit(num, SSC) for num in nums ]

    data = {}
    for molfit in molfits:

        # separate species
        species = split_at_molecule(molfit)
        for specie in species:
            spx = specie.split()[0]
            # spx = re.sub(u'(;$|;#1)', '', spx)
            spx = re.sub(u';$', '', spx)
            try:
                data[spx]
            except:
                data[spx] = {q: {'all': []} for q in ['temperature', 'column density', 'linewidth', 'velocity']}

            # separate components
            components = specie.split('\n')[1:]

            for ncomp,component in enumerate(filter(None,components)):
                c = component.split()

                # get values
                t_llim = float(c[5])
                t_ulim = float(c[6])
                t_val  = float(c[7])
                N_llim = float(c[9])
                N_ulim = float(c[10])
                N_val  = float(c[11])
                w_llim = float(c[13])
                w_ulim = float(c[14])
                w_val  = float(c[15])
                v_llim = float(c[17])
                v_ulim = float(c[18])
                v_val  = float(c[19])

                try:
                    data[spx]['temperature']['all'][ncomp]
                except:
                    for q in ['temperature', 'column density', 'linewidth', 'velocity']:
                        data[spx][q]['all'].append([])

                filter_vals()

    # get percentiles
    for spx,specie in data.items():
        for q,quantity in specie.items():
            for ncomp,component in enumerate(quantity['all']):
                p16,median,p84 = np.percentile(component, (16,50,84))
                try:
                    data[spx][q]['fit']
                except:
                    data[spx][q]['fit']    = []
                    data[spx][q]['16th']   = []
                    data[spx][q]['median'] = []
                    data[spx][q]['84th']   = []
                data[spx][q]['fit'].append( component[0] )
                data[spx][q]['16th'].append( p16 )
                data[spx][q]['median'].append( median )
                data[spx][q]['84th'].append( p84 )

    return data


def get_opacity(data_table):

    import re

    def load_opacity_files(num, SSC):
        taupath = escape_fname(os.path.join(tempdir,'SSC_'+str(SSC['no']),'run_'+str(num),'model','opacity.pickle'))
        with open(taupath, 'rb') as f:
            t_data = pickle.load(f, encoding="latin1")
            return t_data

    def extract_opacities(t_dat):
        int = t_dat[0]
        tau = t_dat[1]
        for specie in tau:
            spx = specie[0]
            # spx = re.sub(u'(;$|;#1)', '', spx)
            spx = re.sub(u';$', '', spx)
            c   = specie[1]-1
            opt = specie[2][:,1]
            int_opt  = np.sum(opt)
            peak_opt = np.max(opt)
            try:
                data_table[SSC['num']][spx]['peak opacity']['all'][c].append(peak_opt)
                data_table[SSC['num']][spx]['integrated opacity']['all'][c].append(int_opt)
            except:
                print(spx,c)

    def get_ncomps(SSC, spx):
        tau = load_opacity_files(0, SSC)[1]
        if spx+';' in [t[0] for t in tau]:
            ncomps = len([t for t in tau if spx+';'==t[0]])
        else:
            ncomps = len([t for t in tau if spx==t[0]])
        return ncomps

    def prepare_table(SSC):
        for spx, specie in data_table[SSC['num']].items():
            try:
                specie['peak opacity']
                specie['integrated opacity']
            except:
                ncomps = get_ncomps(SSC, spx)
                specie['peak opacity']       = {k: [[] for _ in np.arange(ncomps)] for k in ['all','16th','median','84th']}
                specie['integrated opacity'] = {k: [[] for _ in np.arange(ncomps)] for k in ['all','16th','median','84th']}

    def calc_percentiles(SSC):
        for spx in data_table[SSC['num']].keys():
            for c,comp in enumerate(data_table[SSC['num']][spx]['peak opacity']['all']):
                p16,median,p84 = np.percentile(comp, (16,50,84))
                data_table[SSC['num']][spx]['peak opacity']['16th'][c]   = p16
                data_table[SSC['num']][spx]['peak opacity']['median'][c] = median
                data_table[SSC['num']][spx]['peak opacity']['84th'][c]   = p84
            for c,comp in enumerate(data_table[SSC['num']][spx]['integrated opacity']['all']):
                p16,median,p84 = np.percentile(comp, (16,50,84))
                data_table[SSC['num']][spx]['integrated opacity']['16th'][c]   = p16
                data_table[SSC['num']][spx]['integrated opacity']['median'][c] = median
                data_table[SSC['num']][spx]['integrated opacity']['84th'][c]   = p84

    for SSC in tqdm(SSCs):
        prepare_table(SSC)
        t_data = [load_opacity_files(num,SSC) for num in nums]
        for t_dat in t_data:
            extract_opacities(t_dat)
        calc_percentiles(SSC)


temperature_data = {SSC['num']: parse_molfit(SSC) for SSC in SSCs}
get_opacity(temperature_data)
fnpickle(temperature_data, os.path.join(tempdir, 'temperature_data.pickle'))


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
        savepath = escape_fname(os.path.join(plotdir, '05.temperatures', 'histograms','SSC_'+SSC['num']+'.'+spx+'.histogram.pdf'))
        os.system('mkdir -p '+os.path.dirname(savepath))
        fig.savefig(savepath, dpi=300, bbox_inches='tight')

    temperatures     = temperature_data[SSC['num']][spx]['temperature']['all']
    column_densities = temperature_data[SSC['num']][spx]['column density']['all']
    linewidths       = temperature_data[SSC['num']][spx]['linewidth']['all']
    velocities       = temperature_data[SSC['num']][spx]['velocity']['all']
    opacities        = temperature_data[SSC['num']][spx]['integrated opacity']['all']
    fig,[ax_T, ax_N, ax_lw, ax_v, ax_o] = set_up_figure(SSC, spx)
    plot_hist(ax_T, temperatures, r'T [K]')
    plot_hist(ax_N, column_densities, r'$\Sigma$ [cm$^{-2}$]')
    plot_hist(ax_lw, linewidths, r'$\sigma$ [km\,s$^{-1}$]')
    plot_hist(ax_v, velocities, r'v$_\mathrm{rel}$ [km\,s${-1}$]')
    plot_hist(ax_o, opacities, r'$\tau_\mathrm{int}$')
    format_figure([ax_T, ax_N, ax_lw, ax_v, ax_o])
    save_figure(fig, spx)


for SSC in tqdm(SSCs):
    for spx in temperature_species:
        if spx in temperature_data[SSC['num']].keys():
            plot_fit_hist(SSC, spx)


###################################################################################################
#
###################################################################################################
