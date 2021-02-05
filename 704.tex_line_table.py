def tex_transition(molecule, transition=None, vibration=None, return_all=False):
    m_replace = [['O2','O$_2$'],['13C','$^{13}$C'],['O+','O$^+$'],['15N','$^{15}$N'],['C3','C$_3$'],['33S','$^{33}$S'],['34S','$^{34}$S'],['18O','$^{18}$O'],['H2','H$_2$']]
    t_replace = [['(','_{'],[')','}']]
    v_replace = [[',',', '],['v2',r'\nu_2'],['v4',r'\nu_4'],['v6',r'\nu_6'],['v7',r'\nu_7'],['3Sum_v',r'\sum_{3} \nu'],['v=',r'\nu=']]

    for m_r in m_replace:
        molecule = molecule.replace(m_r[0],m_r[1])
    for t_r in t_replace:
        transition = transition.replace(t_r[0],t_r[1])
    for v_r in v_replace:
        vibration = vibration.replace(v_r[0],v_r[1])
    if return_all:
        return molecule,transition,vibration
    else:
        if not transition and not vibration:
            return molecule
        if transition and not vibration:
            return molecule,transition
        if transition and vibration:
            return molecule,transition,vibration


def tex_line_string(str):
    x_replace = [[';',' '],['#1','']]
    m_replace = [['O2','O$_2$'],['13C','$^{13}$C'],['O+','O$^+$'],['15N','$^{15}$N'],['C3','C$_3$'],['33S','$^{33}$S'],['34S','$^{34}$S'],['18O','$^{18}$O'],['H2','H$_2$']]
    t_replace = [['(','_{'],[')','}']]
    v_replace = [[',',', '],['v2',r'\nu_2'],['v4',r'\nu_4'],['v6',r'\nu_6'],['v7',r'\nu_7'],['3Sum_v',r'\sum_{3} \nu'],['v=',r'\nu=']]

    for r in x_replace:
        str = str.replace(r[0],r[1])
    for r in m_replace:
        str = str.replace(r[0],r[1])
    for r in t_replace:
        str = str.replace(r[0],r[1])
    for r in v_replace:
        str = str.replace(r[0],r[1])
    return str


def tex_line_table(lines, species):

    for specie in species:

        # get all lines of specie
        slines = [l for l in lines if l['XCLASS']==specie]
        slines = sorted(slines, key = lambda k: k['restfreq'])

        # get tex string
        for idx,line in enumerate(slines):
            m,t,v = tex_transition(line['molecule'],line['transition'],line['vibration'])
            tex = ''
            if idx==0:
                tex += m
                tex += ' & '
                tex += t
                tex += ' & '
                tex += v
                tex += ' & '
                tex += '{:.5f}'.format(line['restfreq'].to(u.GHz).value)
                tex += r'\\'
            else:
                tex += ' & '
                tex += t
                tex += ' & '
                tex += v
                tex += ' & '
                tex += '{:.5f}'.format(line['restfreq'].to(u.GHz).value)
                tex += r'\\'
            print(tex)


def detected_lines_table(detected_species):

    # get list of all detected_species
    all_species = []
    for SSC in SSCs:
        for specie in detected_species[str(SSC['no'])]:
            if not specie in all_species:
                all_species.append(specie)

    for specie in all_species:

        # get all lines of specie
        slines = [l for l in lines if l['XCLASS']==specie]
        slines = sorted(slines, key = lambda k: k['restfreq'])

        # get tex string
        for idx,line in enumerate(slines):
            m,t,v = tex_transition(line['molecule'],line['transition'],line['vibration'])
            # write line, transition, frequency
            tex = ''
            if idx==0:
                tex += m
                tex += ' & '
                tex += t
                tex += ' & '
                tex += v
                tex += ' & '
                tex += '{:.5f}'.format(line['restfreq'].to(u.GHz).value)
            else:
                tex += ' & '
                tex += t
                tex += ' & '
                tex += v
                tex += ' & '
                tex += '{:.5f}'.format(line['restfreq'].to(u.GHz).value)

            tex += ' && '

            # mark if detected
            for SSC in SSCs:
                tex += ' & '
                if specie in detected_species[str(SSC['no'])]:
                    tex += r'\checkmark'

            tex += r'\\'

            # print table row
            print(tex)


def ratios_N_table(ratios):
    # header
    header = r'\colhead{SSC} & '
    for rname in ratios:
        header += r'\colhead{'+rname+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # column density ratios
    for SSC in SSCs:
        row = SSC['num'] +' & '
        for rname in ratios:
            rdict = ratios_N[rname]
            median = rdict['median'][SSC['no']-1]
            if np.isnan(median):
                row += ' & '
            else:
                if median>100:
                    fmt = '{:.0f}'
                else:
                    fmt = '{:.1f}'
                p16    = rdict['16th'][SSC['no']-1]
                p84    = rdict['84th'][SSC['no']-1]
                upper  = p84-median
                lower  = median-p16
                row += fmt.format(median)+r'$^{+'+fmt.format(upper)+'}_{-'+fmt.format(lower)+'}$ & '
        row = row[:-3]
        row += r'\\'
        print(row)


def ratios_I_table(ratios='all'):
    if ratios=='all':
        # header
        header = r'\colhead{SSC} & '
        for SSC in SSCs:
            header += r'\colhead{'+str(SSC['no'])+r'} & '
        header = header[:-3]
        header += r'\\'
        print(header)
        # line ratios
        for rname,rdict in ratios_I.items():
            row = tex_line_string(rname) +' & '
            for SSC in SSCs:
                bestfit = rdict[str(SSC['no'])]['bestfit']
                if np.isnan(bestfit):
                    row += ' & '
                else:
                    if bestfit>100:
                        fmt = '{:.0f}'
                    else:
                        fmt = '{:.1f}'
                    error = rdict[str(SSC['no'])]['error']
                    row += '$'+fmt.format(bestfit)+r'\pm'+fmt.format(error)+'$ & '
            row = row[:-3]
            row += r'\\'
            print(row)
    else:
        # header
        header = r'\colhead{SSC} & '
        for ratio in ratios:
            header += r'\colhead{'+tex_line_string(ratio)+r'} & '
        header = header[:-3]
        header += r'\\'
        print(header)
        for SSC in SSCs:
            row = str(SSC['no']) +' & '
            for ratio in ratios:
                bestfit = ratios_I[ratio][str(SSC['no'])]['bestfit']
                if np.isnan(bestfit):
                    row += ' & '
                else:
                    if bestfit>100:
                        fmt = '{:.0f}'
                    else:
                        fmt = '{:.1f}'
                    error = ratios_I[ratio][str(SSC['no'])]['error']
                    row += '$'+fmt.format(bestfit)
                    if not np.isnan(error):
                        row += r'\pm'+fmt.format(error)
                    row += '$ & '
            row = row[:-3]
            row += r'\\'
            print(row)



def temperature_table(t_species):
    header = r'\colhead{SSC} & '
    for specie in t_species:
        header += r'\colhead{'+tex_line_string(specie)+r'} & '
    header = header[:-3]
    header += r'\\'
    print(header)
    # temperatures
    for SSC in SSCs:
        row = str(SSC['no'])+' & '
        for specie in t_species:
            try:
                nearest_comp = np.argmin(np.abs(line_data_N[str(SSC['no'])][specie]['velocity']['median']))
                median = line_data_N[str(SSC['no'])][specie]['temperature']['median'][nearest_comp]
            except:
                median = np.nan
            if np.isnan(median):
                row += ' & '
            else:
                p16    = line_data_N[str(SSC['no'])][specie]['temperature']['16th'][nearest_comp]
                p84    = line_data_N[str(SSC['no'])][specie]['temperature']['84th'][nearest_comp]
                upper  = p84-median
                lower  = median-p16
                row += '{:.0f}'.format(median)+r'$^{+'+'{:.0f}'.format(upper)+'}_{-'+'{:.0f}'.format(lower)+'}$ & '
        row = row[:-3]
        row += r'\\'
        print(row)
