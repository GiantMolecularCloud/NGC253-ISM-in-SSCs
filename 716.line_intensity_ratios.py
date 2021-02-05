#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Get physical quantities from fits, such as integrated intensity, mass, linewidth, line shift, ...


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
bandfit = fnunpickle(os.path.join(mandir, 'updated_band_fit_parameters.pickle'))
line_data_I = fnunpickle(os.path.join(mandir,'line_intensity_data.pickle'))


###################################################################################################
# line intensity ratios
###################################################################################################

ratios_I = {
'CO/HCN':       {'a':'CO 3-2 v=0',              'b':'HCN 4-3'},
'CO/HCO+':      {'a':'CO 3-2 v=0',              'b':'HCO+ 4-3'},
'CO/CS':        {'a':'CO 3-2 v=0',              'b':'CS 7-6'},
'HCN/HCO+':     {'a':'HCN 4-3',                 'b':'HCO+ 4-3'},
'HCN/HNC':      {'a':'HC15N 4-3',               'b':'H15NC 4-3'},
'HCN/HCNthin':  {'a':'HCN 4-3',                 'b':'H13CN 4-3'},
'HCN/HCNthin2': {'a':'HCN 4-3',                 'b':'HC15N 4-3'},
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

for rname,rdict in tqdm(ratios_I.items()):
    for SSC in SSCs:
        ratios_I[rname][str(SSC['no'])] = {}

        try:
            a     = line_data_I[str(SSC['no'])][rdict['a']]['integrated intensity']['bestfit'].value
            a_err = line_data_I[str(SSC['no'])][rdict['a']]['integrated intensity']['error'].value
            b     = line_data_I[str(SSC['no'])][rdict['b']]['integrated intensity']['bestfit'].value
            b_err = line_data_I[str(SSC['no'])][rdict['b']]['integrated intensity']['error'].value
        except:
            a, a_err, b, b_err = np.nan, np.nan, np.nan, np.nan

        r     = a/b
        r_err = r* np.sqrt( (a_err/a)**2 +(b_err/b)**2 )

        ratios_I[rname][str(SSC['no'])]['bestfit'] = r
        ratios_I[rname][str(SSC['no'])]['error']   = r_err

os.system('mkdir -p '+ratiodir)
fnpickle(ratios_I, os.path.join(ratiodir, 'ratios_intensity.pickle'))


###################################################################################################
# plot ratios
###################################################################################################

def plot_ratios_I(rname, rdict):
    fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    ax.text(0.05, 0.9, rname, color='k', transform=ax.transAxes, ha='left', va='top', weight='bold', fontsize=16, bbox=props)

    sscs   = [SSC['no'] for SSC in SSCs]
    ratios = [ratios_I[rname][str(SSC['no'])]['bestfit'] for SSC in SSCs]
    r_err  = [ratios_I[rname][str(SSC['no'])]['error'] for SSC in SSCs]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]

    for s,r,e,c in zip(sscs,ratios,r_err,colors):
        ax.errorbar(s, r, yerr=e, marker='o', ms=6, color=c, elinewidth=2, ecolor=c)

    ax.set_xlim(0, len(SSCs)+1)
    ax.set_yscale('log')
    ax.set_xticks(np.arange(1,len(SSCs)+1))
    ax.set_xticklabels([str(i) for i in np.arange(1, len(SSCs)+1)])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set_ylabel(r'I ('+lineID_tex(rdict['a'])+r') / I ('+lineID_tex(rdict['b'])+r')', fontsize=12)
    ax.set_xlabel('SSC')
    fig.tight_layout()

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'ratios_I', 'ratio_'+rname.replace('/','-')+'.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')


###################################################################################################
# plot ratios
###################################################################################################

for rname,rdict in tqdm(ratios_I.items()):
    plot_ratios_I(rname, rdict)


# dense gas ratio table
ratios_I_table(['CO/HCN','CO/CS','HCN/HCO+','CS/HCN'])

# optical depth ratio table
ratios_I_table(['HCN/HCNthin','HCN/HCNthin2','SO/S18O'])

# element ratio table
ratios_I_table(['C/O'])

# isotopologue ratio table
ratios_I_table(['12C/13C','14N/15N','32S/33S','32S/34S'])

# H,C,N chemistry table
ratios_I_table(['HCN/HNC'])

# HC3N super dense gas fraction
ratios_I_table(['HCN/HC3N'])

# sulfur chemistry
ratios_I_table(['SO/SO2','CS/SO2'])

# vibrational excitation table
ratios_I_table(['HCN/HCNvib1','HCN/HCNvib2','HC3N/HC3Nvib1','HC3N/HC3Nvib2'])


###################################################################################################
# XDR/PDR plots
###################################################################################################

def plot_XDR_PDR_I():
    """
    Line ratio plot to decide on XDR vs PDR and low density vs high density as in Baan+08
    """

    fig,axes = plt.subplots(nrows=2, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
    fig.subplots_adjust(hspace=0, wspace=0) #, top=0.80, bottom=0.04, left=0.04, right=0.93)

    # get data
    sscs   = [SSC['no'] for SSC in SSCs]
    colors = [plt.cm.inferno(i/(len(SSCs)+1)) for i in SSCs['no']]
    HCO_HCN, HNC_HCN, HNC_HCO = [],[],[]
    HCO_HCN_err, HNC_HCN_err, HNC_HCO_err = [],[],[]
    for SSC in SSCs:
        try:
            hco_hcn     = ratios_I['HCO+/HCN'][str(SSC['no'])]['bestfit']
            hco_hcn_err = ratios_I['HCO+/HCN'][str(SSC['no'])]['error']
            HCO_HCN.append(     np.log10(hco_hcn) )
            HCO_HCN_err.append( 0.434*hco_hcn_err/hco_hcn )
        except:
            HCO_HCN.append(     np.nan )
            HCO_HCN_err.append( np.nan )
        try:
            hnc_hcn     = ratios_I['HNC/HCN'][str(SSC['no'])]['bestfit']
            hnc_hcn_err = ratios_I['HNC/HCN'][str(SSC['no'])]['error']
            HNC_HCN.append(     np.log10(hnc_hcn) )
            HNC_HCN_err.append( 0.434*hnc_hcn_err/hnc_hcn )
        except:
            HNC_HCN.append(     np.nan )
            HNC_HCN_err.append( np.nan )
        try:
            hnc_hco     = ratios_I['H15NC/HCO+'][str(SSC['no'])]['bestfit']*ratios_I['14N/15N'][str(SSC['no'])]['bestfit']
            hnc_hco_err = np.sqrt( (ratios_I['H15NC/HCO+'][str(SSC['no'])]['error']/ratios_I['H15NC/HCO+'][str(SSC['no'])]['bestfit'])**2 +(ratios_I['14N/15N'][str(SSC['no'])]['error']/ratios_I['14N/15N'][str(SSC['no'])]['bestfit'])**2 )
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

    # panel 1: HCO+/HCN over HNC/HCO+
    ax = axes[0][0]
    ax.plot([-10,10],[10,-10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
    ax.text(0.95, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='right', va='bottom', weight='bold', fontsize=16)
    for a,b,a_err,b_err,c,s in zip(HNC_HCO, HCO_HCN, HNC_HCO_err, HCO_HCN_err, colors, SSCs):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, color=c, elinewidth=1, ecolor=c, label='SSC '+str(s['no']), zorder=3)
    ax.errorbar(B_HCO_HCN[0],B_HNC_HCO[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, color='lime', elinewidth=1, ecolor='lime', label=r'NGC 253 ($26^{\prime\prime}$, Baan +08)', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-0.80,0.80)
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_ylabel(r'log HCO$^+$ / HCN', fontsize=12)

    # panel 2: HNC/HCN over HCO/HCN
    ax = axes[0][1]
    ax.plot([0,0],[-10,10], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([0,10],[-10,-10],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
    ax.text(0.05, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='left', va='bottom', weight='bold', fontsize=16)
    for a,b,a_err,b_err,c in zip(HNC_HCN, HCO_HCN, HNC_HCN_err, HCO_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCN[0],B_HCO_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-0.80,0.80)
    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_xlabel(r'log HNC / HCN', fontsize=12)

    # panel 3: HNC/HCO over HNC/HCN
    ax = axes[1][0]
    ax.plot([-10,10],[0,0], ls='-', lw=1, c='grey', zorder=2)
    ax.fill_between([-10,10],[0,0],[10,10], color='lightgrey', alpha=0.5, zorder=1)
    ax.text(0.95, 0.9, 'XDR', color='k', transform=ax.transAxes, ha='right', va='top', weight='bold', fontsize=16)
    ax.text(0.95, 0.1, 'PDR', color='k', transform=ax.transAxes, ha='right', va='bottom', weight='bold', fontsize=16)
    for a,b,a_err,b_err,c in zip(HNC_HCO, HNC_HCN, HNC_HCO_err, HNC_HCN_err, colors):
        if np.isfinite(a) and np.isfinite(b):
            ax.errorbar(a,b, xerr=a_err, yerr=b_err, marker='o', ms=5, color=c, elinewidth=1, ecolor=c, zorder=3)
    ax.errorbar(B_HNC_HCO[0],B_HNC_HCN[0], xerr=B_HCO_HCN[1], yerr=B_HNC_HCO[1], marker='o', ms=5, color='lime', elinewidth=1, ecolor='lime', zorder=4)
    ax.set_xlim(-1.15,0.45)
    ax.set_ylim(-1.00,0.60)
    ax.set_axisbelow(True)
    ax.grid(axis='both')
    ax.set_xlabel(r'log HNC$^{**}$ / HCO$^+$', fontsize=12)
    ax.set_ylabel(r'log HNC / HCN)', fontsize=12)

    # panel 4: legend
    ax = axes[1][1]
    ax.set_axis_off()
    fig.legend(loc=3, bbox_to_anchor=(0.5,0.1,0.14,0.3), ncol=1, mode="expand", borderaxespad=0., fontsize=12, frameon=False)

    savepath = escape_fname(os.path.join(plotdir, '04.fit_results', 'XDR-PDR.pdf'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    fig.savefig(savepath, dpi=300, bbox_inches='tight')

plot_XDR_PDR_I()


###################################################################################################
#
###################################################################################################
