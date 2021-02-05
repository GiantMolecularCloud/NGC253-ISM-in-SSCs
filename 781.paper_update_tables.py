#############################
# GAS IN SSCS: update paper #
#############################

# redo the paper tables


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
detected_species = fnunpickle(os.path.join(XCLASSdir, 'detected_species.pickle'))

all_data = fnunpickle(join(refitdir, 'all_data.pickle'))
ratios   = fnunpickle(join(refitdir, 'ratios.pickle'))


###################################################################################################
# table 1: integrated intensities
###################################################################################################

header = r"""% \begin{longrotatetable}
\floattable
\begin{deluxetable*}{lclCRRRRRRRRRRRRRR}
    \tabletypesize{\footnotesize}
    \tablecaption{Sample of the fitted integrated intensities in K\,\kms. The full table is available in electronic form including error intervals.
    \label{table: intensities}}
    \tablehead{\colhead{molecule} & \multicolumn{2}{c}{transition} & \colhead{comp.$^\dagger$} & \multicolumn{14}{c}{SSC no.}\\ \cline{2-3} \cline{5-18}
	& \colhead{rotational} & \colhead{vibrational} & & \colhead{1}&\colhead{2}&\colhead{3}&\colhead{4}&\colhead{5}&\colhead{6}&\colhead{7}&\colhead{8}&\colhead{9}&\colhead{10}&\colhead{11}&\colhead{12}&\colhead{13}&\colhead{14}}
    \startdata"""
print(header)

for row in data_transitions:
    line = ''
    line += '{:<14}'.format(molecule_tex(row['molecule'])) +' & '
    line += '{:<24}'.format(rotation_tex(row['rotation'])) +' & '
    line += '{:<14}'.format(vibration_tex(row['vibration'])) +' & '
    line += str(row['component']) +' & '
    for SSC in SSCs:
        intensity = row['SSC '+SSC['num']+' integrated intensity'][1].value
        if np.isnan(intensity):
            line += ' ... & '
        elif intensity<1.0:
            line += ' ... & '
        else:
            line += '{:4.0f}'.format(intensity) + ' & '
    line = line[:-3]
    line += r' \\'
    print(line)

footer = r"""    \enddata
    \tablenotetext{\dagger}{Components do not necessarily correspond to each other. For instance, the undetected marks (...) for component 5 of \co in SSC~1 merely indicate that five components fit the spectrum sufficiently well whereas in SSC~7 six components are required.}
    \tablecomments{The typical intensity error is 13.0\% and consists primarily of uncertainty due to line crowding. Additionally, the systematic flux uncertainty applies with $\lesssim 5\%$ for these observations according to the ALMA specifications.}
\end{deluxetable*}
% \end{longrotatetable}"""
print(footer)


###################################################################################################
# table 2: XCLASS fit results
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{rlccCRRRRR}
    \tablecaption{Sample of the parameters fitted with XCLASS. The full table is available in electronic form including error intervals for all quantities.
    \label{table: column densities}}
    \tablehead{\colhead{SSC} & \colhead{molecule} & \colhead{vibration} & \colhead{component} & \colhead{$\log \mathrm{N}$} & \colhead{$\sigma$} & \colhead{$\mu$} & \colhead{T$_\mathrm{kin}$} & \colhead{$\tau_\mathrm{int}$} \\
                & & & & [\mathrm{cm}^{-2}] & [\mathrm{km}\,\mathrm{s}^{-1}] & [\mathrm{km}\,\mathrm{s}^{-1}] & [\mathrm{K}] & \\
                & & & (1) & (2) & (3) & (4) & (5) & (6)}
    \startdata"""
print(header)

for row in data_species:
    if not row['SSC']=='14':
        line = '% '
    else:
        line = ''
    line += row['SSC'] +' & '
    line += '{:<12}'.format(molecule_tex(row['species'].split(';')[0])) +' & '
    line += '{:<12}'.format(vibration_tex(row['species'].split(';')[1])) +' & '
    line += str(row['component']) +' & '

    N = np.log10(row['column density'][1].value)
    s = row['linewidth'][1].value
    v = row['velocity'][1].value
    T = row['temperature'][1].value
    t = row['integrated opacity'][1]

    line += '... & ' if np.isnan(N) else ('{:5.2f}'.format(N) +' & ')
    line += '... & ' if np.isnan(s) else ('{:4.1f}'.format(s) +' & ')
    line += '... & ' if np.isnan(v) else ('{:6.1f}'.format(v) +' & ')
    line += '... & ' if np.isnan(T) else ('{:3.0f}'.format(T) +' & ')
    line += '... & ' if np.isnan(t) else ('{:5.2f}'.format(t) +' & ')

    line = line[:-3]
    line += r' \\'
    print(line)

footer = r"""    \enddata
    \tablecomments{(1) The numbers assigned to the fitted components do not necessarily correspond to each other.\\
    (2) Molecular column density.\\
    (3) Velocity dispersion.\\
    (4) Line centroid position with respect to the systemic velocity of each SSC (cf. \ref{section: spectra}).\\
    (5) Kinetic temperature is fitted for SO$_2$ and H$_2$CS only and fixed to 130\,K and 300\,K for rotational and ro-vibrational species, respectively.\\
    (6) Total line opacity, i.e. integrated over the line.\\
    The typical column density error is 14.3\% and consists primarily of uncertainty due to line crowding. Additionally, the systematic flux uncertainty applies with $\lesssim 5\%$ for these observations according to the ALMA specifications.}
\end{deluxetable*}"""
print(footer)


###################################################################################################
# table 3: dense gas line ratio + column density ratio
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
        median = ratios[rname]['integrated intensity'][SSC['num']]['median']
        if np.isnan(median):
            row += '... & '
        else:
            if median>100:
                fmt = '{:4.0f}'
            else:
                fmt = '{:4.1f}'
            p16 = ratios[rname]['integrated intensity'][SSC['num']]['16th']
            p84 = ratios[rname]['integrated intensity'][SSC['num']]['84th']
            upper = p84-median
            lower = median-p16
            row += '$'+fmt.format(median)+r'^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '

    # column density ratios
    for rname in ['CO/HCN','CO/CS','HCN/HCO+','CS/HCN']:
        median = ratios[rname]['column density'][SSC['num']]['median']
        if np.isnan(median):
            row += '... & '
        else:
            if median>100:
                fmt = '{:5.0f}'
            else:
                fmt = '{:5.1f}'
            p16 = ratios[rname]['column density'][SSC['num']]['16th']
            p84 = ratios[rname]['column density'][SSC['num']]['84th']
            upper = p84-median
            lower = median-p16
            median = round_significant(median, 2)
            upper  = round_significant(upper, 2)
            lower  = round_significant(lower, 2)
            row += '$'+fmt.format(median)+r'^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
    \tablecomments{The ratios above are derived from the following rotational transitions: \co, \hcn, \hco and \cs. Errors of 0.0 are due to rounding.}
\end{deluxetable*}"""
print(footer)


###################################################################################################
# table 4: selected ratios
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|ccccccc}
    \tablewidth{\linewidth}
    \tablecaption{Line intensity ratios of selected species.\label{table: ratios other}}
    \tablehead{\colhead{SSC} & \colhead{HC$^{15}$N/H$^{15}$NC} & \colhead{HCN/H$^{13}$CN} & \colhead{HCN/HC$^{15}$N} & \colhead{SO/S$^{18}$O} & \colhead{HCN/HC$_3$N} & \colhead{SO/SO$_2$} & \colhead{CS/SO$_2$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = '{:>2}'.format(SSC['num']) +' & '

    # line ratios
    for rname in ['HCN/HNC','HCN/H13CN','HCN/HC15N','SO/S18O','HCN/HC3N','SO/SO2','CS/SO2']:
        median = ratios[rname]['integrated intensity'][SSC['num']]['median']
        if np.isnan(median):
            row += '                   ... & '
        else:
            if median>100:
                fmt = '{:4.0f}'
            else:
                fmt = '{:4.1f}'
            p16 = ratios[rname]['integrated intensity'][SSC['num']]['16th']
            p84 = ratios[rname]['integrated intensity'][SSC['num']]['84th']
            upper = p84-median
            lower = median-p16
            row += '$'+fmt.format(median)+r'^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
    \tablecomments{The ratios above are derived from the following rotational transitions: \hcn, \mbox{H$^{13}$CN(4--3)}, HC$^{15}$N(4--3), H$^{15}$NC(4--3), HC$_{3}$N(38--37), \cs, SO($8_8$--$7_7$), S$^{18}$O ($8_9$--$7_8$) and SO$_2$($11_{4,8}$--$11_{3,9}$).}
\end{deluxetable*}"""
print(footer)


###################################################################################################
# table 5: vibrational excitation ratios
###################################################################################################

header = r"""\floattable
\begin{deluxetable*}{r|cccc}
    \tablecaption{Integrated intensity ro-vibrational over rotational excitation fraction for selected HCN and HC$_3$N lines. The ro-vibrational line used for the respective ratio is given in the second row.\label{table: ratios vibrational excitation}}
    \tablehead{\colhead{SSC} & \colhead{HCN/HCN*} & \colhead{HCN/HCN*} & \colhead{HC$_3$N/HC$_3$N*} & \colhead{HC$_3$N/HC$_3$N*}\\
    \colhead{} & \colhead{$\nu_2=1$, $l=1f$} & \colhead{$\nu_2=2$, $l=2f$} & \colhead{$\nu_7=1$, $l=1f$} & \colhead{$\nu_7=2$, $l=2f$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = SSC['num'] +' & '

    # line ratios
    for rname in ['HCN/HCNvib1','HCN/HCNvib2','HC3N/HC3Nvib1','HC3N/HC3Nvib2']:
        median = ratios[rname]['integrated intensity'][SSC['num']]['median']
        if np.isnan(median):
            row += r'... & '
        else:
            if median>100:
                fmt = '{:.0f}'
            else:
                fmt = '{:.1f}'
            p16 = ratios[rname]['integrated intensity'][SSC['num']]['16th']
            p84 = ratios[rname]['integrated intensity'][SSC['num']]['84th']
            upper = p84-median
            lower = median-p16
            row += '$'+fmt.format(median).replace(' ',r'\phantom{0}')+r'^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
\end{deluxetable*}"""
print(footer)


###################################################################################################
# table 7: temperatures
####################################################################################################

header = r"""\floattable
\begin{deluxetable}{r|cccc}
    \tablecaption{Excitation temperatures obtained by \xclass fitting.\label{table: temperatures}}
    \tablehead{\colhead{SSC} & \colhead{H$_2$CS} & \colhead{SO$_2$}}
    \startdata"""
print(header)

for SSC in SSCs:
    row = '{:>2}'.format(SSC['num']) +' & '

    # line ratios
    for spx in ['H2CS;v=0;#1','SO2;v=0']:
        try:
            median = all_data[SSC['num']][spx]['temperature']['median'][0]
        except:
            median = np.nan
        if np.isnan(median):
            row += '                ... & '
        else:
            p16    = all_data[SSC['num']][spx]['temperature']['16th'][0]
            p84    = all_data[SSC['num']][spx]['temperature']['84th'][0]
            upper  = p84-median
            lower  = median-p16
            if median<250:
                row += '{:3.0f}'.format(median)+r'$^{+'+'{:3.0f}'.format(upper)+'}_{-'+'{:3.0f}'.format(lower)+'}$ & '
            else:
                row += '             $>'+'{:3.0f}'.format(p16)+'$ & '
    row = row[:-3]
    row += r'\\'
    print(row)

footer = r"""    \enddata
    \tablecomments{We report lower limits in the cases where the temperature is weakly constrained.}
\end{deluxetable}"""
print(footer)



###################################################################################################
# estimate errors
###################################################################################################

# column density
###################################################################################################

dat = np.array([row['column density'].value for row in data_species])
medians = dat[:,1]
p16s    = dat[:,0]
p84s    = dat[:,2]
uppers  = p84s-medians
lowers  = medians-p16s
rel_ups = uppers/medians*100
rel_los = lowers/medians*100
med_rel_ups = np.median(rel_ups)
med_rel_los = np.median(rel_los)
print("\ncolumn density\nmedian relavtive error\nupper: "+str(med_rel_ups)+"   lower:"+str(med_rel_los))

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,12))
ax.scatter(medians, rel_ups, c='blue', s=20, marker='.', label='lower')
ax.scatter(medians, rel_los, c='red', s=20, marker='.', label='upper')
ax.axhline(y=med_rel_ups, c='grey')
ax.axhline(y=med_rel_los, c='grey')
ax.text(1e20, med_rel_ups, 'upper', color='grey', ha='right', va='bottom')
ax.text(1e20, med_rel_los, 'lower', color='grey', ha='right', va='top')
ax.set_axisbelow(True)
ax.grid(ls=':', c='grey')
ax.set_xlabel(r'column density [cm$^{-2}$]', fontsize=12)
ax.set_ylabel(r'relative error [\%]', fontsize=12)
ax.set_xscale('log')
ax.set_ylim([0,200])
fig.legend(loc='top right')
fig.savefig(join(plotdir, 'errors_species.pdf'), bbox_inches='tight', dpi=300)


# integrated intensity
###################################################################################################

dat = np.array([row['SSC '+SSC['num']+' integrated intensity'].value for row in data_transitions for SSC in SSCs])
no_nans = [False if np.isnan(x[1]) else True for x in dat]
dat = dat[no_nans]
medians = dat[:,1]
p16s    = dat[:,0]
p84s    = dat[:,2]
uppers  = p84s-medians
lowers  = medians-p16s
rel_ups = uppers/medians*100
rel_los = lowers/medians*100
med_rel_ups = np.median(rel_ups)
med_rel_los = np.median(rel_los)
print("\nintegrated intensity\nmedian relavtive error\nupper: "+str(med_rel_ups)+"   lower:"+str(med_rel_los))

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(8,12))
ax.scatter(medians, rel_ups, c='blue', s=20, marker='.', label='lower')
ax.scatter(medians, rel_los, c='red', s=20, marker='.', label='upper')
ax.axhline(y=med_rel_ups, c='grey')
ax.axhline(y=med_rel_los, c='grey')
ax.text(1e5, med_rel_ups, 'upper', color='grey', ha='right', va='bottom')
ax.text(1e5, med_rel_los, 'lower', color='grey', ha='right', va='top')
ax.set_axisbelow(True)
ax.grid(ls=':', c='grey')
ax.set_xlabel(r'integrated intensity [K\,km\,s$^{-1}$]', fontsize=12)
ax.set_ylabel(r'relative error [\%]', fontsize=12)
ax.set_xscale('log')
ax.set_xlim([1e-2,1e5])
ax.set_ylim([0,200])
fig.legend(loc='top right')
fig.savefig(join(plotdir, 'errors_transitions.pdf'), bbox_inches='tight', dpi=300)


# line width
###################################################################################################

dat = np.array([row['linewidth'].value for row in data_species])
medians = dat[:,1]
p16s    = dat[:,0]
p84s    = dat[:,2]
uppers  = p84s-medians
lowers  = medians-p16s
rel_ups = uppers/medians*100
rel_los = lowers/medians*100
med_rel_ups = np.median(rel_ups)
med_rel_los = np.median(rel_los)
print("\nlinewidth\nmedian relavtive error\nupper: "+str(med_rel_ups)+"   lower:"+str(med_rel_los))


# line position
###################################################################################################

dat = np.array([row['velocity'].value for row in data_species])
medians = dat[:,1]
p16s    = dat[:,0]
p84s    = dat[:,2]
uppers  = p84s-medians
lowers  = medians-p16s
rel_ups = uppers/medians*100
rel_los = lowers/medians*100
med_rel_ups = np.median(rel_ups)
med_rel_los = np.median(rel_los)
print("\nvelocity\nmedian relavtive error\nupper: "+str(med_rel_ups)+"   lower:"+str(med_rel_los))


# integrated opacity
###################################################################################################

dat = np.array([row['integrated opacity'] for row in data_species])
medians = dat[:,1]
p16s    = dat[:,0]
p84s    = dat[:,2]
uppers  = p84s-medians
lowers  = medians-p16s
rel_ups = uppers/medians*100
rel_los = lowers/medians*100
med_rel_ups = np.median(rel_ups)
med_rel_los = np.median(rel_los)
print("\nintegrated opacity\nmedian relavtive error\nupper: "+str(med_rel_ups)+"   lower:"+str(med_rel_los))


###################################################################################################
#
###################################################################################################
