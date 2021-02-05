#######################
# GAS IN SSCS: XCLASS #
#######################

# use all_data to get column densities and intensities


###################################################################################################
# column density / intensity ratios
###################################################################################################

ratios = {
'CO/HCN':        {'a': {'specie': 'CO;v=0', 'transition': '3-2', 'vibration': 'v=0'},                'b': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'CO/HCO+':       {'a': {'specie': 'CO;v=0', 'transition': '3-2', 'vibration': 'v=0'},                'b': {'specie': 'HCO+;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'CO/CS':         {'a': {'specie': 'CO;v=0', 'transition': '3-2', 'vibration': 'v=0'},                'b': {'specie': 'CS;v=0', 'transition': '7-6', 'vibration': 'v=0'}},
'HCN/HCO+':      {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCO+;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'HCN/HNC':       {'a': {'specie': 'HCN-15;v=0', 'transition': '4-3', 'vibration': 'v=0'},            'b': {'specie': 'HN-15-C;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'HCN/H13CN':     {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HC-13-N;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'HCN/HC15N':     {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCN-15;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'CS/HCN':        {'a': {'specie': 'CS;v=0', 'transition': '7-6', 'vibration': 'v=0'},                'b': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'CS/HCO+':       {'a': {'specie': 'CS;v=0', 'transition': '7-6', 'vibration': 'v=0'},                'b': {'specie': 'HCO+;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'SO/S18O':       {'a': {'specie': 'SO;v=0;#2', 'transition': '8(8)-7(7)', 'vibration': '3Sum_v=0'},  'b': {'specie': 'SO-18;v=0', 'transition': '8(9)-7(8)', 'vibration': 'v=0'}},
'12C/13C':       {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HC-13-N;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'14N/15N':       {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCN-15;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'32S/33S':       {'a': {'specie': 'SO;v=0;#2', 'transition': '8(8)-7(7)', 'vibration': '3Sum_v=0'},  'b': {'specie': 'S-33-O;v=0', 'transition': '9(8)-8(7)', 'vibration': 'v=0'}},
'C/O':           {'a': {'specie': 'CS;v=0', 'transition': '7-6', 'vibration': 'v=0'},                'b': {'specie': 'SO;v=0;#2', 'transition': '8(8)-7(7)', 'vibration': '3Sum_v=0'}},
'SO/SO2':        {'a': {'specie': 'SO;v=0;#2', 'transition': '8(8)-7(7)', 'vibration': '3Sum_v=0'},  'b': {'specie': 'SO2;v=0', 'transition': '11(4,8)-11(3,9)', 'vibration': 'v=0'}},
'CS/SO2':        {'a': {'specie': 'CS;v=0', 'transition': '7-6', 'vibration': 'v=0'},                'b': {'specie': 'SO2;v=0', 'transition': '11(4,8)-11(3,9)', 'vibration': 'v=0'}},
'HCN/HC3N':      {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCCCN;v=0', 'transition': '38-37', 'vibration': 'v=0'}},
'HCN/HCNvib1':   {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCN;v2=1', 'transition': '4-3', 'vibration': 'v2=1,l=1f'}},
'HCN/HCNvib2':   {'a': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCN;v2=2', 'transition': '4-3', 'vibration': 'v2=2,l=2f'}},
'HC3N/HC3Nvib1': {'a': {'specie': 'HCCCN;v=0', 'transition': '39-38', 'vibration': 'v=0'},           'b': {'specie': 'HCCCN;v7=1', 'transition': '39-38', 'vibration': 'v7=1,l=1f'}},
'HC3N/HC3Nvib2': {'a': {'specie': 'HCCCN;v=0', 'transition': '39-38', 'vibration': 'v=0'},           'b': {'specie': 'HCCCN;v7=2', 'transition': '39-38', 'vibration': 'v7=2,l=2f'}},
# for XDR/PDR plot
'HNC/HCN':      {'a': {'specie': 'HN-15-C;v=0', 'transition': '4-3', 'vibration': 'v=0'},            'b': {'specie': 'HCN-15;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'HCO+/HCN':     {'a': {'specie': 'HCO+;v=0', 'transition': '4-3', 'vibration': 'v=0'},               'b': {'specie': 'HCN;v=0', 'transition': '4-3', 'vibration': 'v=0'}},
'H15NC/HCO+':   {'a': {'specie': 'HN-15-C;v=0', 'transition': '4-3', 'vibration': 'v=0'},            'b': {'specie': 'HCO+;v=0', 'transition': '4-3', 'vibration': 'v=0'}}
}


def get_nearest_component(SSC, specie):
    return np.argmin(np.abs(all_data[SSC['num']][specie]['velocity']['median']))

def get_strongest_component(SSC, specie):
    return np.argmax(np.abs(all_data[SSC['num']][specie]['column density']['median']))

def get_ratio_percentiles(R):
    med = np.nanpercentile(R, 50)
    bad_idx = np.append( np.where(R<med/5.)[0], np.where(R>med*5)[0] )
    p16,med,p84 = np.nanpercentile( np.delete(R,bad_idx), (16,50,84) )
    for bi in bad_idx:
        R[bi] = med
    return p16,med,p84

def add_ratio_to_dict(SSC,type,vals_a,vals_b):
    R = vals_a/vals_b
    p16,med,p84 = get_ratio_percentiles(R)

    rdict[type][SSC['num']]['median'] = med
    rdict[type][SSC['num']]['16th']   = p16
    rdict[type][SSC['num']]['84th']   = p84
    rdict[type][SSC['num']]['all']    = R

def add_dummy_ratio(SSC,type):
    rdict[type][SSC['num']]['median'] = np.nan
    rdict[type][SSC['num']]['16th']   = np.nan
    rdict[type][SSC['num']]['84th']   = np.nan
    rdict[type][SSC['num']]['all']    = [np.nan for _ in np.arange(101)]


# calculate ratios
####################################################################################################

for rname,rdict in tqdm(ratios.items()):
    spx_a = rdict['a']['specie']
    trn_a = rdict['a']['transition']
    vib_a = rdict['a']['vibration']
    spx_b = rdict['b']['specie']
    trn_b = rdict['b']['transition']
    vib_b = rdict['b']['vibration']
    rdict['column density']       = {}
    rdict['integrated intensity'] = {}

    for SSC in SSCs:
        rdict['column density'][SSC['num']]       = {}
        rdict['integrated intensity'][SSC['num']] = {}

        # test if both lines exist for this SSC
        try:

            # get nearest components
            near_a = get_nearest_component(SSC, spx_a)
            near_b = get_nearest_component(SSC, spx_b)
            # near_a = get_strongest_component(SSC, spx_a)
            # near_b = get_strongest_component(SSC, spx_b)

            # get column density
            N_a = np.array(all_data[SSC['num']][spx_a]['column density']['all'][near_a])
            N_b = np.array(all_data[SSC['num']][spx_b]['column density']['all'][near_b])

            # get intensity
            I_a = all_data[SSC['num']][spx_a]['integrated intensity'][trn_a][vib_a][near_a]['all'].value
            I_b = all_data[SSC['num']][spx_b]['integrated intensity'][trn_b][vib_b][near_b]['all'].value

            add_ratio_to_dict(SSC,'column density',N_a,N_b)
            add_ratio_to_dict(SSC,'integrated intensity',I_a,I_b)

        except:
            add_dummy_ratio(SSC, 'column density')
            add_dummy_ratio(SSC, 'integrated intensity')


fnpickle(ratios, os.path.join(refitdir, 'ratios.pickle'))


###################################################################################################
# plot ratios
###################################################################################################

def plot_ratios(rname, rdict, type):
    """
    Make a box-and-wiskers plot of a line ratio in all SSCs.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    sscs   = [SSC['no'] for SSC in SSCs]
    r_all  = [rdict[type][SSC['num']]['all'] for SSC in SSCs]
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
    if type=='column density':
        ax.set_ylabel(r'N ('+specie_tex(rdict['a']['specie'])+r') / N ('+specie_tex(rdict['b']['specie'])+r')', fontsize=12)
    elif type=='integrated intensity':
        ax.set_ylabel(r'I ('+specie_tex(rdict['a']['specie'])+' '+rdict['a']['transition']+' '+rdict['a']['vibration'].replace('_','\_')+r') / I ('+specie_tex(rdict['b']['specie'])+' '+rdict['b']['transition']+' '+rdict['b']['vibration'].replace('_','\_')+r')', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '13.redo_ratios', 'ratios', rname.replace('/','-')+'.'+type.replace(' ','_')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


for rname,rdict in tqdm(ratios.items()):
    plot_ratios(rname, rdict, type='column density')
    plot_ratios(rname, rdict, type='integrated intensity')


###################################################################################################
#
###################################################################################################
