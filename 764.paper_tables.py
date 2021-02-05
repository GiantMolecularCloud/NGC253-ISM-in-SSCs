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
# paper tables
###################################################################################################

# table I: dense gas line ratio + column density ratio
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|cccc|cccc}
    \tablewidth{\linewidth}
    \tablecaption{Dense gas fractions specified by line intensity ratios and column density ratios.\label{table: ratios dense gas}}
    \tablehead{\colhead{SSC} & \multicolumn{4}{c}{line intensity ratio R$_\mathrm{I}$} & \multicolumn{4}{c}{column density ratio R$_\mathrm{N}$}\\
        \colhead{} & \colhead{CO/HCN} & \colhead{CO/CS} & \colhead{HCN/HCO$^+$} & \colhead{CS/HCN} & \colhead{CO/HCN} & \colhead{CO/CS} & \colhead{HCN/HCO$^+$} & \colhead{CS/HCN}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '

    # line ratios
    for rname in ['CO/HCN','CO/CS','HCN/HCO+','CS/HCN']:
        rdict = ratios_Gauss[rname]
        best = rdict[SSC['num']]['bestfit']
        if np.isnan(best):
            row += '... & '
        else:
            if best>100:
                fmt = '{:4.0f}'
            else:
                fmt = '{:4.1f}'
            error = rdict[SSC['num']]['error']
            row += '$'+fmt.format(best).replace(' ',r'\phantom{0}')+r'\pm'+fmt.format(error)+'$ & '

    # column density ratios
    for rname in ['CO/HCN','CO/CS','HCN/HCO+','CS/HCN']:
        rdict = ratios_XCLASS[rname]
        median = rdict[SSC['num']]['median']
        if np.isnan(median):
            row += '... & '
        else:
            if median>100:
                fmt = '{:5.0f}'
            else:
                fmt = '{:5.1f}'
            p16    = rdict[SSC['num']]['16th']
            p84    = rdict[SSC['num']]['84th']
            upper  = p84-median
            lower  = median-p16
            row += '$'+fmt.format(median).replace(' ',r'\phantom{0}')+r'^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
    \tablecomments{Errors of 0.0 are due to rounding. A few of the Gaussian fits (mainly for HCN and HCO$^+$) had to be fitted with fixed parameters (e.g. line centroid) to achieve convergence and thus do not have an error associated.\\
    \todo{Dividing the column densities such that the ratio is always $>1$ results in enormous errors. There is a certain spread in column densities that gets blown up by the division. Plotting these ratios the other way round does not cause the impression of giant or even wrong errors. A reasonable factor 1.5 uncertainty for $1.0 \times 10^{-3}$ appears much smaller than for 1000. If this should go in the paper, I will use log ratios for all column density ratios with errors in dex.}}
\end{deluxetable*}"""
print(footer)


# table II: isotopologue ratios
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|ccccccc}
    \tablewidth{\linewidth}
    \tablecaption{Line intensity ratios of selected species.\label{table: ratios other}}
    \tablehead{\colhead{SSC} & \colhead{HC$^{15}$N/H$^{15}$NC} & \colhead{HCN/H$^{13}$CN} & \colhead{HCN/HC$^{15}$N} & \colhead{SO/S$^{18}$O} & \colhead{HCN/HC$_3$N} & \colhead{SO/SO$_2$} & \colhead{CS/SO$_2$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '

    # line ratios
    for rname in ['HCN/HNC','HCN/H13CN','HCN/HC15N','SO/S18O','HCN/HC3N','SO/SO2','CS/SO2']:
        rdict = ratios_Gauss[rname]
        best = rdict[SSC['num']]['bestfit']
        if np.isnan(best):
            row += '... & '
        else:
            if best>100:
                fmt = '{:4.0f}'
            else:
                fmt = '{:4.1f}'
            error = rdict[SSC['num']]['error']
            row += '$'+fmt.format(best).replace(' ',r'\phantom{0}')+r'\pm'+fmt.format(error)+'$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
\end{deluxetable*}"""
print(footer)


# table III: optical depth
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|cccccc}
    \tablecaption{Total (integrated) optical depth $\tau$ for selected species obtained by \xclass fitting.\label{table: total optical depth}}
    \tablehead{\colhead{SSC} & \colhead{CO} & \colhead{HCO$^+$} & \colhead{HCN} & \colhead{H$^{13}$CN} & \colhead{HC$^{15}$N} & \colhead{H$^{15}$NC}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '
    for spx in ['CO;v=0','HCO+;v=0','HCN;v=0','HC-13-N;v=0','HCN-15;v=0','HN-15-C;v=0']:
        try:
            nearest_comp = np.argmin(np.abs(data_XCLASS[SSC['num']][spx]['velocity']['median']))
            median = data_XCLASS[SSC['num']][spx]['integrated opacity']['median'][nearest_comp]
        except:
            median = np.nan
        if np.isnan(median):
            row += ' ... & '
        else:
            p16    = data_XCLASS[SSC['num']][spx]['integrated opacity']['16th'][nearest_comp]
            p84    = data_XCLASS[SSC['num']][spx]['integrated opacity']['84th'][nearest_comp]
            upper  = p84-median
            lower  = median-p16
            row += '{:4.1f}'.format(median).replace(' ',r'\phantom{0}')+r'$^{+'+'{:4.1f}'.format(upper)+'}_{-'+'{:4.1f}'.format(lower)+'}$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
\end{deluxetable*}"""
print(footer)


# table IV: vibrational excitation ratios
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|cccc}
    \tablecaption{Vibrational excitation fraction for selected HCN and HC$_3$N lines.\label{table: ratios vibrational excitation}}
    \tablehead{\colhead{SSC} & \colhead{HCN/HCNvib} & \colhead{HCN/HCNvib} & \colhead{HC$_3$N/HC$_3$Nvib} & \colhead{HC$_3$N/HC$_3$N*}\\
    \colhead{} & \colhead{$\nu_2=1$, $l=1f$} & \colhead{$\nu_2=2$, $l=2f$} & \colhead{$\nu_7=1$, $l=1f$} & \colhead{$\nu_7=2$, $l=2f$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '

    # line ratios
    for rname in ['HCN/HCNvib1','HCN/HCNvib2','HC3N/HC3Nvib1','HC3N/HC3Nvib2']:
        rdict = ratios_Gauss[rname]
        best = rdict[SSC['num']]['bestfit']
        if np.isnan(best):
            row += r'... & '
        else:
            if best>100:
                fmt = '{:.0f}'
            else:
                fmt = '{:.1f}'
            error = rdict[SSC['num']]['error']
            row += '$'+fmt.format(best)+r'\pm'+fmt.format(error)+'$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
\end{deluxetable*}"""
print(footer)


# table V: temperature
####################################################################################################

header = r"""\floattable
\begin{deluxetable}{r|cccc}
    \tablecaption{Excitation temperatures obtained by \xclass fitting.\label{table: temperatures}}
    \tablehead{\colhead{SSC} & \colhead{H$_2$CS} & \colhead{SO$_2$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '

    # line ratios
    for spx in ['H2CS;v=0;#1','SO2;v=0']:
        try:
            nearest_comp = np.argmin(np.abs(data_XCLASS[SSC['num']][spx]['velocity']['median']))
            median = data_XCLASS[SSC['num']][spx]['temperature']['median'][nearest_comp]
        except:
            median = np.nan
        if np.isnan(median):
            row += ' ... & '
        else:
            p16    = data_XCLASS[SSC['num']][spx]['temperature']['16th'][nearest_comp]
            p84    = data_XCLASS[SSC['num']][spx]['temperature']['84th'][nearest_comp]
            upper  = p84-median
            lower  = median-p16
            if median<250:
                row += '{:.0f}'.format(median)+r'$^{+'+'{:.0f}'.format(upper)+'}_{-'+'{:.0f}'.format(lower)+'}$ & '
            else:
                row += '$>'+'{:.0f}'.format(p16)+'$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
    \tablecomments{We report lower limits in the cases where the temperature is weakly constrained.}
\end{deluxetable}"""
print(footer)


###################################################################################################
# Appendix
###################################################################################################

# intensities table
###################################################################################################

header = r"""\begin{longrotatetable}
\begin{deluxetable}{lLLCCCCCCCCCCCCCC}
    \tablecaption{Integrated intensities of the main component obtained by Gaussian fitting.\label{table: intensities}}
    \tablehead{\multirow{2}{*}{molecule} & \multicolumn{2}{c}{transition} & \multicolumn{14}{c}{SSC no.}\\
	& \colhead{rotational} & \colhead{vibrational} & \colhead{1}&\colhead{2}&\colhead{3}&\colhead{4}&\colhead{5}&\colhead{6}&\colhead{7}&\colhead{8}&\colhead{9}&\colhead{10}&\colhead{11}&\colhead{12}&\colhead{13}&\colhead{14}}
    \startdata"""
print(header)

for spx in unique_species:

    # get all lines of specie
    slines = [l for l in lines if l['XCLASS']==spx]
    slines = sorted(slines, key = lambda k: k['restfreq'])

    # get tex string
    for idx,line in enumerate(slines):
        m,t,v = tex_transition(line['molecule'],line['transition'],line['vibration'], return_all=True)
        row = ''
        if idx==0:
            row += m
            row += ' & '
            row += t
            row += ' & '
            row += v
            row += ' & '
        else:
            row += ' & '
            row += t
            row += ' & '
            row += v
            row += ' & '

        for SSC in SSCs:
            try:
                # nearest_comp = np.argmin(np.abs(data_Gauss[SSC['num']][line['ID']]['line shift']['bestfit']))
                # median = data_Gauss[SSC['num']][line['ID']]['integrated intensity']['bestfit'][nearest_comp]
                best = data_Gauss[SSC['num']][line['ID']]['integrated intensity']['bestfit'].value
            except:
                best = np.nan
            if np.isnan(best):
                row += ' ... & '
            else:
                # row += '{:4.0f}'.format(best).replace(' ',r'\phantom{0}')+' & '
                row += latex_float(best, fmt='{:4.0f}', phantom=True)
                row += ' & '
        row = row[:-3]
        row += r'\\'
        print(row)

footer = r"""    \enddata
\end{deluxetable}
\end{longrotatetable}"""
print(footer)


# column densities table
###################################################################################################

header = r"""\startlongtable
\begin{deluxetable*}{lLLCCCCCCC}
    \tablecaption{Column densities of the main component fitted with XCLASS.\label{table: column densities}}
    \tablehead{\multirow{2}{*}{molecule} & \multicolumn{2}{c}{transition} & \multicolumn{7}{c}{SSC no.}\\
	& \colhead{rotational} & \colhead{vibrational} & \colhead{1}&\colhead{2}&\colhead{3}&\colhead{4}&\colhead{5}&\colhead{6}&\colhead{7}}
    \startdata"""
print(header)

for spx in fitable_species:
    line = [l for l in lines if l['XCLASS']==spx][0]
    m,t,v = tex_transition(line['molecule'],line['transition'],line['vibration'], return_all=True)
    row = m+' & '+t+' & '+v+' & '

    for SSC in SSCs[:7]:
        try:
            nearest_comp = np.argmin(np.abs(data_XCLASS[SSC['num']][spx]['velocity']['median']))
            median = data_XCLASS[SSC['num']][spx]['column density']['median'][nearest_comp]
            row += latex_float(median, fmt='{:.1e}', phantom=True)
            row += ' & '
        except:
            row += ' ... & '
    row = row[:-3]
    row += r'\\'
    print(row)

sechead = r"""\tableline
\multirow{2}{*}{molecule} & \multicolumn{2}{c}{transition} & \multicolumn{7}{c}{SSC no.}\\
& \mathrm{rotational} & \mathrm{vibrational} & 8 & 9 & 10 & 11 & 12 & 13 & 14 \\
\tableline\\[2mm]"""
print(sechead)

for spx in fitable_species:
    line = [l for l in lines if l['XCLASS']==spx][0]
    m,t,v = tex_transition(line['molecule'],line['transition'],line['vibration'], return_all=True)
    row = m+' & '+t+' & '+v+' & '

    for SSC in SSCs[7:]:
        try:
            nearest_comp = np.argmin(np.abs(data_XCLASS[SSC['num']][spx]['velocity']['median']))
            median = data_XCLASS[SSC['num']][spx]['column density']['median'][nearest_comp]
            row += latex_float(median, fmt='{:.1e}', phantom=True)
            row += ' & '
        except:
            row += ' ... & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
\end{deluxetable*}"""
print(footer)


###################################################################################################
# relative errors for appendix tables
###################################################################################################

# integrated intensity
###################################################################################################

header = "relative errors in percent - Gaussian fit - integrated intensity\n----------------------------------------------------------------\n"
header += '{:<25}'.format('line')
for SSC in SSCs:
    header += '{:>10}'.format(SSC['num'])
header += '{:>10}'.format('mean')
header += '\n'
header += '-'*(25+len(SSCs)*10+10)
print(header)

all_I = []
all_e = []
means = []
for spx in unique_species:
    slines = [l for l in lines if l['XCLASS']==spx]
    slines = sorted(slines, key = lambda k: k['restfreq'])
    for idx,line in enumerate(slines):
        row = '{:<25}'.format(line['ID'])
        values = []
        for SSC in SSCs:
            try:
                best  = data_Gauss[SSC['num']][line['ID']]['integrated intensity']['bestfit'].value
                error = data_Gauss[SSC['num']][line['ID']]['integrated intensity']['error'].value
                rel_err = error/best*100
                if rel_err>1e4 or rel_err<1e-4:
                    rel_err = np.nan
                else:
                    all_I.append(best)
                    all_e.append(rel_err)
                row += '{:10.3f}'.format(rel_err)
                values.append(rel_err)
            except:
                row += '{:<10}'.format('')
        if not values==[] and not np.all(np.isnan(np.array(values))):
            mean = np.nanmean(values)
            row += '{:10.3f}'.format(mean)
            means.append(mean)
        print(row)
print(' '*(25+130)+'total mean:   '+'{:.3f}'.format(np.nanmean(means)))
print(' '*(25+130)+'total median: '+'{:.3f}'.format(np.nanmedian(means)))

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,12))
ax.scatter(all_I, all_e, s=20, marker='.')
ax.set_axisbelow(True)
ax.grid(ls=':', c='grey')
ax.set_xlabel(r'integrated intensity [K\,km\,s$^{-1}$]', fontsize=12)
ax.set_ylabel(r'relative error [\%]', fontsize=12)
ax.set_xscale('log')
ax.set_ylim([0,200])
fig.savefig(join(plotdir, 'errors_Gauss.pdf'), bbox_inches='tight', dpi=300)


# column density
###################################################################################################

header = "relative errors in percent - XCLASS fit - column density\n----------------------------------------------------------------\n"
header += '{:<25}'.format('line')
for SSC in SSCs:
    header += '{:>10}'.format(SSC['num'])
header += '{:>10}'.format('mean')
header += '{:>10}'.format('median')
header += '\n'
header += '-'*(25+len(SSCs)*10+10)
print(header)

all_I = []
all_e = []
means = []
for spx in unique_species:
    row = '{:<25}'.format(spx)
    values = []
    for SSC in SSCs:
        try:
            nearest_comp = np.argmin(np.abs(data_XCLASS[SSC['num']][spx]['velocity']['median']))
            median = data_XCLASS[SSC['num']][spx]['column density']['median'][nearest_comp]
            p16    = data_XCLASS[SSC['num']][spx]['column density']['16th'][nearest_comp]
            p84    = data_XCLASS[SSC['num']][spx]['column density']['84th'][nearest_comp]
            mean_error = np.nanmean([median-p16, p84-median])
            rel_err = mean_error/median*100
            if rel_err>1e4 or rel_err<1e-4:
                rel_err = np.nan
            else:
                all_I.append(median)
                all_e.append(rel_err)
            row += '{:10.3f}'.format(rel_err)
            values.append(rel_err)
        except:
            row += '{:<10}'.format('')
    if not values==[] and not np.all(np.isnan(np.array(values))):
        mean = np.nanmean(values)
        median = np.nanmedian(values)
        row += '{:10.3f}'.format(mean)
        row += '{:10.3f}'.format(median)
        means.append(mean)
    print(row)
print(' '*(25+140)+'total mean:   '+'{:.3f}'.format(np.nanmean(means)))
print(' '*(25+140)+'total median: '+'{:.3f}'.format(np.nanmedian(means)))

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,12))
ax.scatter(all_I, all_e, s=20, marker='.')
ax.set_axisbelow(True)
ax.grid(ls=':', c='grey')
ax.set_xlabel(r'column density [cm$^{-2}$]', fontsize=12)
ax.set_ylabel(r'relative error [\%]', fontsize=12)
ax.set_xscale('log')
ax.set_ylim([0,200])
fig.savefig(join(plotdir, 'errors_XCLASS.pdf'), bbox_inches='tight', dpi=300)


###################################################################################################
#
###################################################################################################
