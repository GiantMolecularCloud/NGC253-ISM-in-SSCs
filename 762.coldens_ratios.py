#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))

data_XCLASS = fnunpickle(os.path.join(resultsdir, 'data_XCLASS.pickle'))
data_Gauss  = fnunpickle(os.path.join(resultsdir, 'data_Gauss.pickle'))


###################################################################################################
# column density ratios
###################################################################################################

ratios_XCLASS = {
'CO/HCN':        {'a':'CO;v=0',     'b':'HCN;v=0'},
'CO/HCO+':       {'a':'CO;v=0',     'b':'HCO+;v=0'},
'CO/CS':         {'a':'CO;v=0',     'b':'CS;v=0'},
'HCN/HCO+':      {'a':'HCN;v=0',    'b':'HCO+;v=0'},
'HCN/HNC':       {'a':'HCN-15;v=0', 'b':'HN-15-C;v=0'},
'HCN/H13CN':     {'a':'HCN;v=0',    'b':'HC-13-N;v=0'},
'HCN/HC15N':     {'a':'HCN;v=0',    'b':'HCN-15;v=0'},
'CS/HCN':        {'a':'CS;v=0',     'b':'HCN;v=0'},
'CS/HCO+':       {'a':'CS;v=0',     'b':'HCO+;v=0'},
'SO/S18O':       {'a':'SO;v=0;#1',  'b':'SO-18;v=0'},
'12C/13C':       {'a':'HCN;v=0',    'b':'HC-13-N;v=0'},
'14N/15N':       {'a':'HCN;v=0',    'b':'HCN-15;v=0'},
'32S/33S':       {'a':'SO;v=0;#1',  'b':'S-33-O;v=0'},
'32S/34S':       {'a':'SO2;v=0',    'b':'S-34-O2;v=0'},
'C/O':           {'a':'CS;v=0',     'b':'SO;v=0;#1'},
'SO/SO2':        {'a':'SO;v=0;#1',  'b':'SO2;v=0'},
'CS/SO2':        {'a':'CS;v=0',     'b':'SO2;v=0'},
'HCN/HC3N':      {'a':'HCN;v=0',    'b':'HCCCN;v=0'},
'HCN/HCNvib1':   {'a':'HCN;v=0',    'b':'HCN;v2=1'},
'HCN/HCNvib2':   {'a':'HCN;v=0',    'b':'HCN;v2=2'},
'HC3N/HC3Nvib1': {'a':'HCCCN;v=0',   'b':'HCCCN;v7=1'},
'HC3N/HC3Nvib2': {'a':'HCCCN;v=0',   'b':'HCCCN;v7=2'},
# for XDR/PDR plot
'HNC/HCN':      {'a':'HN-15-C;v=0', 'b':'HCN-15;v=0'},
'HCO+/HCN':     {'a':'HCO+;v=0',    'b':'HCN;v=0'},
'H15NC/HCO+':   {'a':'HN-15-C;v=0', 'b':'HCO+;v=0'}
}

ratios_Gauss = {
'CO/HCN':       {'a':'CO 3-2',                  'b':'HCN 4-3'},
'CO/HCO+':      {'a':'CO 3-2',                  'b':'HCO+ 4-3'},
'CO/CS':        {'a':'CO 3-2',                  'b':'CS 7-6'},
'HCN/HCO+':     {'a':'HCN 4-3',                 'b':'HCO+ 4-3'},
'HCN/HNC':      {'a':'HC15N 4-3',               'b':'H15NC 4-3'},
'HCN/H13CN':    {'a':'HCN 4-3',                 'b':'H13CN 4-3'},
'HCN/HC15N':    {'a':'HCN 4-3',                 'b':'HC15N 4-3'},
'CS/HCN':       {'a':'CS 7-6',                  'b':'HCN 4-3'},
'CS/HCO+':      {'a':'CS 7-6',                  'b':'HCO+ 4-3'},
'SO/S18O':      {'a':'SO 8(8)-7(7) 3Sum_v=0',   'b':'S18O 8(9)-7(8)'},
'12C/13C':      {'a':'HCN 4-3',                 'b':'H13CN 4-3'},
'14N/15N':      {'a':'HCN 4-3',                 'b':'HC15N 4-3'},
'32S/33S':      {'a':'SO 8(8)-7(7) 3Sum_v=0',   'b':'33SO 9(8)-8(7)'},
'32S/34S':      {'a':'SO2 11(4,8)-11(3,9)',     'b':'34SO2 11(4,8)-11(3,9)'},
'C/O':          {'a':'CS 7-6',                  'b':'SO 8(8)-7(7) 3Sum_v=0'},
'SO/SO2':       {'a':'SO 8(8)-7(7) 3Sum_v=0',   'b':'SO2 11(4,8)-11(3,9)'},
'CS/SO2':       {'a':'CS 7-6',                  'b':'SO2 11(4,8)-11(3,9)'},
'HCN/HC3N':     {'a':'HCN 4-3',                 'b':'HC3N 38-37'},
'HCN/HCNvib1':  {'a':'HCN 4-3',                 'b':'HCN 4-3 v2=1,l=1f'},
'HCN/HCNvib2':  {'a':'HCN 4-3',                 'b':'HCN 4-3 v2=2,l=2f'},
'HC3N/HC3Nvib1': {'a':'HC3N 39-38',             'b':'HC3N 39-38 v7=1,l=1f'},
'HC3N/HC3Nvib2': {'a':'HC3N 39-38',             'b':'HC3N 39-38 v7=2,l=2f'},
# for XDR/PDR plot
'HNC/HCN':      {'a':'H15NC 4-3',               'b':'HC15N 4-3'},
'HCO+/HCN':     {'a':'HCO+ 4-3',                'b':'HCN 4-3'},
'H15NC/HCO+':   {'a':'H15NC 4-3',               'b':'HCO+ 4-3'}
}

def get_ratios_XCLASS(ratios):
    for rname,rdict in tqdm(ratios.items()):
        for SSC in SSCs:
            ratios[rname][SSC['num']] = {}
            try:
                nearest_comp_a = np.argmin(np.abs(data_XCLASS[SSC['num']][rdict['a']]['velocity']['median']))
                a_all = np.array(data_XCLASS[SSC['num']][rdict['a']]['column density']['all'][nearest_comp_a])
                nearest_comp_b = np.argmin(np.abs(data_XCLASS[SSC['num']][rdict['b']]['velocity']['median']))
                b_all = np.array(data_XCLASS[SSC['num']][rdict['b']]['column density']['all'][nearest_comp_b])

                r_all  = a_all/b_all
                rmed = np.percentile(r_all, 50)

                # remove outliers
                bad_idx = np.append( np.where(r_all<rmed/5.)[0], np.where(r_all>rmed*5)[0] )
                r16,rmed,r84 = np.percentile( np.delete(r_all,bad_idx), (16,50,84) )
                for bi in bad_idx:
                    r_all[bi] = rmed

            except KeyError:
                r_all = [np.nan]
                r16,rmed,r84 = np.nan,np.nan,np.nan

            ratios[rname][SSC['num']]['median'] = rmed
            ratios[rname][SSC['num']]['16th'] = r16
            ratios[rname][SSC['num']]['84th'] = r84
            ratios[rname][SSC['num']]['all'] = r_all

            try:
                nearest_comp_a = np.argmin(np.abs(data_XCLASS[SSC['num']][rdict['a']]['velocity']['median']))
                a_med = np.array(data_XCLASS[SSC['num']][rdict['a']]['column density']['median'][nearest_comp_a])
                a_p16 = np.array(data_XCLASS[SSC['num']][rdict['a']]['column density']['16th'][nearest_comp_a])
                a_p84 = np.array(data_XCLASS[SSC['num']][rdict['a']]['column density']['84th'][nearest_comp_a])
                nearest_comp_b = np.argmin(np.abs(data_XCLASS[SSC['num']][rdict['b']]['velocity']['median']))
                b_med = np.array(data_XCLASS[SSC['num']][rdict['b']]['column density']['median'][nearest_comp_a])
                b_p16 = np.array(data_XCLASS[SSC['num']][rdict['b']]['column density']['16th'][nearest_comp_a])
                b_p84 = np.array(data_XCLASS[SSC['num']][rdict['b']]['column density']['84th'][nearest_comp_a])
                rr   = a_med/b_med
                rlow = a_p16/b_p84
                rup  = a_p84/b_16
            except:
                rr,rlow,rup  = np.nan,np.nan,np.nan

            ratios[rname][SSC['num']]['ratio'] = rmed
            ratios[rname][SSC['num']]['error low'] = r16
            ratios[rname][SSC['num']]['error high'] = r84

    return ratios

def get_ratios_Gauss(ratios):
    for rname,rdict in tqdm(ratios.items()):
        for SSC in SSCs:
            ratios[rname][SSC['num']] = {}

            try:
                a     = data_Gauss[SSC['num']][rdict['a']]['integrated intensity']['bestfit'].value
                a_err = data_Gauss[SSC['num']][rdict['a']]['integrated intensity']['error'].value
                b     = data_Gauss[SSC['num']][rdict['b']]['integrated intensity']['bestfit'].value
                b_err = data_Gauss[SSC['num']][rdict['b']]['integrated intensity']['error'].value
            except:
                a, a_err, b, b_err = np.nan, np.nan, np.nan, np.nan

            r     = a/b
            r_err = r* np.sqrt( (a_err/a)**2 +(b_err/b)**2 )

            ratios[rname][SSC['num']]['bestfit'] = r
            ratios[rname][SSC['num']]['error']   = r_err

    return ratios


ratios_XCLASS = get_ratios_XCLASS(ratios_XCLASS)
fnpickle(ratios_XCLASS, os.path.join(resultsdir, 'ratios_XCLASS.pickle'))

# ratios_Gauss = fnunpickle(os.path.join(mandir, 'ratios_intensity.pickle'))
ratios_Gauss = get_ratios_Gauss(ratios_Gauss)
fnpickle(ratios_Gauss, os.path.join(resultsdir, 'ratios_Gauss.pickle'))


###################################################################################################
# plot ratios
###################################################################################################

def plot_ratios(rname, rdict):
    """
    Make a box-and-wiskers plot of a line ratio in all SSCs.
    """

    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    sscs   = [SSC['no'] for SSC in SSCs]
    r_all  = [ratios_XCLASS[rname][SSC['num']]['all'] for SSC in SSCs]
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

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'ratios', 'ratio_'+rname.replace('/','-')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


def ratios_table(rnames):
    # header
    header = r'\colhead{SSC} & '
    for rname in rnames:
        header += r'\colhead{'+rname+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # column density ratios
    for SSC in SSCs:
        row = SSC['num'] +' & '
        for rname in rnames:
            rdict = ratios_XCLASS[rname]
            median = rdict[SSC['num']]['median']
            if np.isnan(median):
                row += '... & '
            else:
                if median>100:
                    fmt = '{:.0f}'
                else:
                    fmt = '{:.1f}'
                p16    = rdict[SSC['num']]['16th']
                p84    = rdict[SSC['num']]['84th']
                upper  = p84-median
                lower  = median-p16
                row += fmt.format(median)+r'$^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
        row = row[:-3]
        row += r'\\'
        print(row)


def ratios_table2(rnames):
    # header
    header = r'\colhead{SSC} & '
    for rname in rnames:
        header += r'\colhead{'+rname+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # column density ratios
    for SSC in SSCs:
        row = SSC['num'] +'... & '
        for rname in rnames:
            rdict = ratios_XCLASS[rname]
            median = rdict[SSC['num']]['ratio']
            if np.isnan(median):
                row += ' & '
            else:
                if median>100:
                    fmt = '{:.0f}'
                else:
                    fmt = '{:.1f}'
                lower = rdict[SSC['num']]['error low']
                upper = rdict[SSC['num']]['error high']
                row += fmt.format(median)+r'$^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
        row = row[:-3]
        row += r'\\'
        print(row)


###################################################################################################
# plot ratios
###################################################################################################

for rname,rdict in tqdm(ratios_XCLASS.items()):
    plot_ratios(rname, rdict)


# dense gas ratio table
ratios_table(['CO/HCN','CO/CS','HCN/HCO+','CS/HCN'])

# optical depth ratio table
ratios_table(['HCN/HCNthin','HCN/HCNthin2','SO/S18O'])

# element ratio table
ratios_table(['C/O'])

# isotopologue ratio table
ratios_table(['12C/13C','14N/15N','32S/33S','32S/34S'])

# H,C,N chemistry table
ratios_table(['HCN/HNC'])

# HC3N super dense gas fraction
ratios_table(['HCN/HC3N'])

# sulfur chemistry
ratios_table(['SO/SO2','CS/SO2'])

# vibrational excitation table
ratios_table(['HCN/HCNvib1','HCN/HCNvib2','HC3N/HC3Nvib1','HC3N/HC3Nvib2'])


###################################################################################################
# XDR/PDR plots
###################################################################################################

def plot_XDR_PDR_XCLASS():
    """
    Line ratio plot to decide on XDR vs PDR and low density vs high density as in Baan+08
    """

    fig,axes = plt.subplots(nrows=2, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    fig.subplots_adjust(hspace=0, wspace=0) #, top=0.80, bottom=0.04, left=0.04, right=0.93)

    # get data
    sscs   = [SSC['no'] for SSC in SSCs]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    HCO_HCN, HNC_HCN, HNC_HCO = [],[],[]
    HCO_HCN_err, HNC_HCN_err, HNC_HCO_err = [],[],[]
    for SSC in SSCs:
        try:
            hco_hcn_med = ratios_XCLASS['HCO+/HCN'][SSC['num']]['median']
            hco_hcn_p16 = ratios_XCLASS['HCO+/HCN'][SSC['num']]['16th']
            hco_hcn_p84 = ratios_XCLASS['HCO+/HCN'][SSC['num']]['84th']
            hco_hcn_low = hco_hcn_med-hco_hcn_p16
            hco_hcn_hig = hco_hcn_p84-hco_hcn_med
            HCO_HCN.append(     np.log10(hco_hcn_med) )
            HCO_HCN_err.append( [0.434*hco_hcn_low/hco_hcn_med,0.434*hco_hcn_hig/hco_hcn_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )
        try:
            hnc_hcn_med = ratios_XCLASS['HNC/HCN'][SSC['num']]['median']
            hnc_hcn_p16 = ratios_XCLASS['HNC/HCN'][SSC['num']]['16th']
            hnc_hcn_p84 = ratios_XCLASS['HNC/HCN'][SSC['num']]['84th']
            hnc_hcn_low = hnc_hcn_med-hnc_hcn_p16
            hnc_hcn_hig = hnc_hcn_p84-hnc_hcn_med
            HNC_HCN.append(     np.log10(hnc_hcn_med) )
            HNC_HCN_err.append( [0.434*hnc_hcn_low/hco_hcn_med,0.434*hnc_hcn_hig/hco_hcn_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )
        try:
            hnc_hco_med = ratios_XCLASS['H15NC/HCO+'][SSC['num']]['median']*ratios_XCLASS['14N/15N'][SSC['num']]['median']
            hnc_hco_p16 = ratios_XCLASS['H15NC/HCO+'][SSC['num']]['16th']*ratios_XCLASS['14N/15N'][SSC['num']]['median']
            hnc_hco_p84 = ratios_XCLASS['H15NC/HCO+'][SSC['num']]['84th']*ratios_XCLASS['14N/15N'][SSC['num']]['median']
            hnc_hco_low = hnc_hco_med-hnc_hco_p16
            hnc_hco_hig = hnc_hco_p84=hnc_hco_med
            HNC_HCO.append(     np.log10(hnc_hco_med) )
            HNC_HCO_err.append( [0.434*hnc_hco_low/hnc_hco_med,0.434*hnc_hco_hig/hnc_hco_med] )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( [np.nan,np.nan] )

    # comparison from Baan+08
    B_hcn = [318.2, 14]
    B_hnc = [234.0, 7]
    B_hco = [276.1, 14]
    B_hco_hcn = [B_hco[0]/B_hcn[0], B_hco[0]/B_hcn[0]*np.sqrt((B_hco[1]/B_hco[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hcn = [B_hnc[0]/B_hcn[0], B_hnc[0]/B_hcn[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hco = [B_hnc[0]/B_hco[0], B_hnc[0]/B_hco[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hco[1]/B_hco[0])**2)]
    B_HCO_HCN = [np.log10(B_hco_hcn[0]), 0.434*B_hco_hcn[1]/B_hco_hcn[0]]
    B_HNC_HCN = [np.log10(B_hnc_hcn[0]), 0.434*B_hnc_hcn[1]/B_hnc_hcn[0]]
    B_HNC_HCO = [np.log10(B_hnc_hco[0]), 0.434*B_hnc_hco[1]/B_hnc_hco[0]]

    def format_panel(ax):
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_axisbelow(True)
        ax.grid(axis='both', which='both')

    def label_regions(ax):
        ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
        ax.text(0.05, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='left', va='bottom', weight='bold', fontsize=16)

    # panel 1: HCO+/HCN over HNC/HCO+
    ax = axes[0][0]
    ax.plot([-10,10],[10,-10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c,s in zip(HNC_HCO, HCO_HCN, HNC_HCO_err, HCO_HCN_err, colors, SSCs):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=[[a_err[0]],[a_err[1]]], yerr=[[b_err[0]],[b_err[1]]], marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, label='SSC '+str(s['no']), zorder=3)
    ax.errorbar(B_HCO_HCN[0],B_HNC_HCO[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', label=r'NGC 253 (Baan +08)', zorder=4)
    ax.set_xlim(-0.75,0.75)
    ax.set_ylim(-0.85,0.65)
    format_panel(ax)
    ax.set_ylabel(r'log N(HCO$^+$) / N(HCN)', fontsize=12)

    # panel 2: HNC/HCN over HCO/HCN
    ax = axes[0][1]
    ax.plot([0,0],[-10,10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([0,10],[-10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c in zip(HNC_HCN, HCO_HCN, HNC_HCN_err, HCO_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=[[a_err[0]],[a_err[1]]], yerr=[[b_err[0]],[b_err[1]]], marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCN[0],B_HCO_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-0.95,0.55)
    ax.set_ylim(-0.85,0.65)
    format_panel(ax)
    ax.tick_params(labelbottom=True)
    ax.set_xlabel(r'log N(HNC) / N(HCN)', fontsize=12)

    # panel 3: HNC/HCO over HNC/HCN
    ax = axes[1][0]
    ax.plot([-10,10],[0,0], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[0,0],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c in zip(HNC_HCO, HNC_HCN, HNC_HCO_err, HNC_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=[[a_err[0]],[a_err[1]]], yerr=[[b_err[0]],[b_err[1]]], marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCO[0],B_HNC_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-0.75,0.75)
    ax.set_ylim(-1.05,0.45)
    format_panel(ax)
    ax.set_xlabel(r'log N(HNC$^{**}$) / N(HCO$^+$)', fontsize=12)
    ax.set_ylabel(r'log N(HNC) / N(HCN)', fontsize=12)

    # panel 4: legend
    ax = axes[1][1]
    ax.set_axis_off()
    fig.legend(loc=3, bbox_to_anchor=(0.55,0.05,0.14,0.3), ncol=1, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'XDR-PDR_column_density.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


def plot_XDR_PDR_Gauss():
    """
    Line ratio plot to decide on XDR vs PDR and low density vs high density as in Baan+08
    """

    fig,axes = plt.subplots(nrows=2, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(6,6))
    fig.subplots_adjust(hspace=0, wspace=0) #, top=0.80, bottom=0.04, left=0.04, right=0.93)

    # get data
    sscs   = [SSC['no'] for SSC in SSCs]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    HCO_HCN, HNC_HCN, HNC_HCO = [],[],[]
    HCO_HCN_err, HNC_HCN_err, HNC_HCO_err = [],[],[]
    for SSC in SSCs:
        try:
            hco_hcn     = ratios_Gauss['HCO+/HCN'][SSC['num']]['bestfit']
            hco_hcn_err = ratios_Gauss['HCO+/HCN'][SSC['num']]['error']
            HCO_HCN.append(     np.log10(hco_hcn) )
            HCO_HCN_err.append( 0.434*hco_hcn_err/hco_hcn )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( np.nan )
        try:
            hnc_hcn     = ratios_Gauss['HNC/HCN'][SSC['num']]['bestfit']
            hnc_hcn_err = ratios_Gauss['HNC/HCN'][SSC['num']]['error']
            HNC_HCN.append(     np.log10(hnc_hcn) )
            HNC_HCN_err.append( 0.434*hnc_hcn_err/hnc_hcn )
        except:
            HNC_HCN.append(     np.nan )
            HNC_HCN_err.append( np.nan )
        try:
            hnc_hco     = ratios_Gauss['H15NC/HCO+'][SSC['num']]['bestfit']*ratios_Gauss['14N/15N'][SSC['num']]['bestfit']
            hnc_hco_err = np.sqrt( (ratios_Gauss['H15NC/HCO+'][SSC['num']]['error']/ratios_Gauss['H15NC/HCO+'][SSC['num']]['bestfit'])**2 +(ratios_Gauss['14N/15N'][SSC['num']]['error']/ratios_Gauss['14N/15N'][SSC['num']]['bestfit'])**2 )
            HNC_HCO.append(     np.log10(hnc_hco) )
            HNC_HCO_err.append( 0.434*hnc_hco_err/hnc_hco )
        except:
            HNC_HCO.append(     np.nan )
            HNC_HCO_err.append( np.nan )

    # comparison from Baan+08
    B_hcn = [318.2, 14]
    B_hnc = [234.0, 7]
    B_hco = [276.1, 14]
    B_hco_hcn = [B_hco[0]/B_hcn[0], B_hco[0]/B_hcn[0]*np.sqrt((B_hco[1]/B_hco[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hcn = [B_hnc[0]/B_hcn[0], B_hnc[0]/B_hcn[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hcn[1]/B_hcn[0])**2)]
    B_hnc_hco = [B_hnc[0]/B_hco[0], B_hnc[0]/B_hco[0]*np.sqrt((B_hnc[1]/B_hnc[0])**2+(B_hco[1]/B_hco[0])**2)]
    B_HCO_HCN = [np.log10(B_hco_hcn[0]), 0.434*B_hco_hcn[1]/B_hco_hcn[0]]
    B_HNC_HCN = [np.log10(B_hnc_hcn[0]), 0.434*B_hnc_hcn[1]/B_hnc_hcn[0]]
    B_HNC_HCO = [np.log10(B_hnc_hco[0]), 0.434*B_hnc_hco[1]/B_hnc_hco[0]]

    def format_panel(ax):
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_axisbelow(True)
        ax.grid(axis='both', which='both')

    def label_regions(ax):
        ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
        ax.text(0.05, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='left', va='bottom', weight='bold', fontsize=16)

    # panel 1: HCO+/HCN over HNC/HCO+
    ax = axes[0][0]
    ax.plot([-10,10],[10,-10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c,s in zip(HNC_HCO, HCO_HCN, HNC_HCO_err, HCO_HCN_err, colors, SSCs):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, label='SSC '+str(s['no']), zorder=3)
    ax.errorbar(B_HCO_HCN[0],B_HNC_HCO[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', label=r'NGC 253 (Baan +08)', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-0.80,0.80)
    format_panel(ax)
    ax.set_ylabel(r'log I(HCO$^+$) / I(HCN)', fontsize=12)

    # panel 2: HNC/HCN over HCO/HCN
    ax = axes[0][1]
    ax.plot([0,0],[-10,10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([0,10],[-10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c in zip(HNC_HCN, HCO_HCN, HNC_HCN_err, HCO_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCN[0],B_HCO_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-0.80,0.80)
    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    format_panel(ax)
    ax.set_xlabel(r'log I(HNC) / I(HCN)', fontsize=12)

    # panel 3: HNC/HCO over HNC/HCN
    ax = axes[1][0]
    ax.plot([-10,10],[0,0], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[0,0],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    label_regions(ax)
    for a,b,a_err,b_err,c in zip(HNC_HCO, HNC_HCN, HNC_HCO_err, HNC_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, lw=0, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCO[0],B_HNC_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, lw=0, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-1.00,0.60)
    format_panel(ax)
    ax.set_xlabel(r'log I(HNC$^{**}$) / I(HCO$^+$)', fontsize=12)
    ax.set_ylabel(r'log I(HNC) / I(HCN)', fontsize=12)

    # panel 4: legend
    ax = axes[1][1]
    ax.set_axis_off()
    fig.legend(loc=3, bbox_to_anchor=(0.55,0.05,0.14,0.3), ncol=1, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    savepath = escape_fname(os.path.join(plotdir, '10.results', 'XDR-PDR_line_ratio.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


plot_XDR_PDR_XCLASS()
plot_XDR_PDR_Gauss()


###################################################################################################
#
###################################################################################################
