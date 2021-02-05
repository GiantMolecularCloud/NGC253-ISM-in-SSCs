###################################################################################################
# table 1: integrated intensities
###################################################################################################

tab1 = open('data_transitions.MRF.dat','w+')


# generate format
###################################################################################################

tab1.write('FORMAT\n======\n\n')

fmt = 'A14,X3,A24,X3,A14,X3,I9,X3'
for SSC in SSCs:
    for i in range(3):
        fmt += ',I4,X3'
fmt = fmt[:-3]
tab1.write(fmt)
tab1.write('\n\n\n')



# generate units
###################################################################################################

tab1.write('UNITS\n=====\n\n')
units = '---\n'*4
for SSC in SSCs:
    for i in range(3):
        units += 'K.km/s\n'
tab1.write(units)
tab1.write('\n\n\n')


# generate labels
###################################################################################################

tab1.write('LABELS\n======\n\n')
labels = 'molecule\nrotation\nvibration\ncomponent\n'
for SSC in SSCs:
        labels += 'I-16th(SSC'+SSC['num']+')\n'
        labels += 'I-50th(SSC'+SSC['num']+')\n'
        labels += 'I-84th(SSC'+SSC['num']+')\n'
tab1.write(labels)
tab1.write('\n\n\n')


# generate explanations
###################################################################################################

tab1.write('EXPLANATIONS\n======\n\n')
explanations = ''
explanations += 'molecule\n'
explanations += 'rotational state\n'
explanations += 'vibrational state\n'
explanations += 'spectral component\n'
for SSC in SSCs:
    explanations += 'SSC'+SSC['num']+' 16th percentile lower bound (lower 1sigma error) on integrated intensity for the intensity distribution across 100 fits\n'
    explanations += 'SSC'+SSC['num']+' 50th percentile median of integrated intensity for the intensity distribution across 100 fits\n'
    explanations += 'SSC'+SSC['num']+' 84th percentile upper bound (upper 1sigma error) on integrated intensity for the intensity distribution across 100 fits\n'
tab1.write(explanations)
tab1.write('\n\n\n')


# generate table data
###################################################################################################

delimiter = '   '

header = ''
header += '{:<14}'.format('molecule') +delimiter
header += '{:<24}'.format('rotation') +delimiter
header += '{:<14}'.format('transition') +delimiter
header += '{:<9}'.format('component')+delimiter
for SSC in SSCs:
    header += '{:^21}'.format('SSC '+SSC['num'])
header += '\n'+' '*(14+24+14+9)+delimiter*4
for SSC in SSCs:
    header += '16th'+delimiter+'50th'+delimiter+'84th'+delimiter
header += '\n'+'='*(14+24+14+9+len(delimiter)*3 +len(SSCs)*3*(4+len(delimiter)))
tab1.write(header)

for row in data_transitions:
    line = '\n'
    line += '{:<14}'.format(row['molecule']) +delimiter
    line += '{:<24}'.format(row['rotation']) +delimiter
    line += '{:<14}'.format(row['vibration']) +delimiter
    line += '{:9d}'.format(row['component']) +delimiter
    for SSC in SSCs:
        for val in row['SSC '+SSC['num']+' integrated intensity']:
            if np.isnan(val.value):
                line += ' nan'
            elif val.value<1.0:
                line += ' nan'
            else:
                line += '{:4.0f}'.format(val.value)
            line += delimiter
    line = line[:-3]
    tab1.write(line)

tab1.close()



###################################################################################################
# table 2: species table
###################################################################################################

tab2 = open('data_species.MRF.dat','w+')


# generate format
###################################################################################################

tab2.write('FORMAT\n======\n\n')
fmt = 'I3,X3,A12,X3,A12,X3,I9,X3'
for i in range(3):
    fmt += ',F5.2,X3'
for i in range(3):
    fmt += ',F4.1,X3'
for i in range(3):
    fmt += ',F6.1,X3'
for i in range(3):
    fmt += ',F3.0,X3'
for i in range(3):
    fmt += ',F5.2'
tab2.write(fmt)
tab2.write('\n\n\n')


# generate units
###################################################################################################

tab2.write('UNITS\n=====\n\n')
units = '---\n'*4
for i in range(3):
    units += '[cm-2]\n'
for i in range(3):
    units += 'km/s\n'
for i in range(3):
    units += 'km/s\n'
for i in range(3):
    units += 'K\n'
for i in range(3):
    units += '---\n'
tab2.write(units)
tab2.write('\n\n\n')


# generate labels
###################################################################################################

tab2.write('LABELS\n======\n\n')
labels = 'SSC\nmolecule\nvibration\ncomponent\n'
for l in ['logN','lw','v','T','tau']:
    labels += l+'-16th\n'
    labels += l+'-50th\n'
    labels += l+'-84th\n'
tab2.write(labels)
tab2.write('\n\n\n')


# generate explanations
###################################################################################################

tab2.write('EXPLANATIONS\n======\n\n')
explanations = ''
explanations += 'SSC\n'
explanations += 'molecule\n'
explanations += 'vibrational state\n'
explanations += 'spectral component\n'
explanations += 'log10(column density) 16th percentile lower bound (lower 1sigma error) for the distribution across 100 fits\n'
explanations += 'log10(column density) 50th percentile median for the distribution across 100 fits\n'
explanations += 'log10(column density) 84th percentile upper bound (upper 1sigma error) for the distribution across 100 fits\n'
explanations += 'linewidth 16th percentile lower bound (lower 1sigma error) for the distribution across 100 fits\n'
explanations += 'linewidth 50th percentile median for the distribution across 100 fits\n'
explanations += 'linewidth 84th percentile upper bound (upper 1sigma) error for the distribution across 100 fits\n'
explanations += 'velocity 16th percentile lower bound (lower 1sigma error) for the distribution across 100 fits\n'
explanations += 'velocity 50th percentile median for the distribution across 100 fits\n'
explanations += 'velocity 84th percentile upper bound (upper 1sigma error) for the distribution across 100 fits\n'
explanations += 'kinetic temperature 16th percentile lower bound (lower 1sigma error) for the distribution across 100 fits\n'
explanations += 'kinetic temperature 50th percentile median for the distribution across 100 fits\n'
explanations += 'kinetic temperature 84th percentile upper bound (upper 1sigma error) for the distribution across 100 fits\n'
explanations += 'integrated opacity 16th percentile lower bound (lower 1sigma error) for the distribution across 100 fits\n'
explanations += 'integrated opacity 50th percentile median for the distribution across 100 fits\n'
explanations += 'integrated opacity 84th percentile upper bound (upper 1sigma error) for the distribution across 100 fits\n'
tab2.write(explanations)
tab2.write('\n\n\n')


# generate table data
###################################################################################################

delimiter = '   '

header = '{:<3}'.format('SSC')
header += '{:<12}'.format('molecule') +delimiter
header += '{:<12}'.format('vibration') +delimiter
header += '{:<9}'.format('component') +delimiter
header += '{:>5}'.format('N') +delimiter
header += '{:>4}'.format('lw') +delimiter
header += '{:>6}'.format('v') +delimiter
header += '{:>3}'.format('T') +delimiter
header += '{:>5}'.format('tau') +delimiter
header += '\n'+'='*(14+24+14+9+len(delimiter)*3 +len(SSCs)*3*(4+len(delimiter)))
tab2.write(header)

for row in data_species:
    line = '\n'
    line += '{:<3}'.format(row['SSC']) +delimiter
    line += '{:<12}'.format(row['species'].split(';')[0]) +delimiter
    line += '{:<12}'.format(row['species'].split(';')[1]) +delimiter
    line += '{:9d}'.format(row['component']) +delimiter

    # column density
    for N in row['column density']:
        if np.isnan(N.value):
            line += '  nan'
        else:
            line += '{:5.2f}'.format(np.log10(N.value))
        line += delimiter

    # linewidth
    for s in row['linewidth']:
        if np.isnan(s.value):
            line += '  nan'
        else:
            line += '{:4.1f}'.format(s.value)
        line += delimiter

    # velocity
    for v in row['velocity']:
        if np.isnan(v.value):
            line += '  nan'
        else:
            line += '{:6.1f}'.format(v.value)
        line += delimiter

    # temperature
    for T in row['temperature']:
        if np.isnan(T.value):
            line += '  nan'
        else:
            line += '{:3.0f}'.format(T.value)
        line += delimiter

    # integrated opacity
    for t in row['integrated opacity']:
        if np.isnan(t):
            line += '  nan'
        else:
            line += '{:5.2f}'.format(t)
        line += delimiter

    line = line[:-3]
    tab2.write(line)

tab2.close()
