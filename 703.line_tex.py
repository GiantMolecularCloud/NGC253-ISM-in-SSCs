def line_tex(line):
    """Convert line dictionary to LaTeX string.

    Parameters
    ----------
    line : dict
        Line dictionary

    Returns
    -------
    string
        Formatted LaTeX compatible string including line and transition.
    """

    m_replace = [['O2','O$_2$'],['13C','$^{13}$C'],['O+','O$^+$'],['15N','$^{15}$N'],['C3','C$_3$'],['33S','$^{33}$S'],['34S','$^{34}$S'],['18O','$^{18}$O'],['H2','H$_2$']]
    t_replace = [['(','_{'],[')','}']]
    v_replace = [[',',', '],['v2',r'\nu_2'],['v4',r'\nu_4'],['v6',r'\nu_6'],['v7',r'\nu_7'],['3Sum_v',r'\sum_{3} \nu'],['v=',r'\nu=']]

    m = line['molecule']
    t = line['transition']
    v = line['vibration']

    for m_r in m_replace:
        m = m.replace(m_r[0],m_r[1])
    for t_r in t_replace:
        t = t.replace(t_r[0],t_r[1])
    for v_r in v_replace:
        v = v.replace(v_r[0],v_r[1])

    tex = m+' $('+t
    if not ( v == '' ):
        tex += '\ '+v
    tex += ')$'

    return tex


def lineID_tex(lineID):
    """Convert line ID to LaTeX string.

    Parameters
    ----------
    line : str
        Line ID

    Returns
    -------
    string
        Formatted LaTeX compatible string including line and transition.
    """

    x_replace = [['O2','O$_2$'],['13C','$^{13}$C'],['O+','O$^+$'],['15N','$^{15}$N'],['C3','C$_3$'],['33S','$^{33}$S'],['34S','$^{34}$S'],['18O','$^{18}$O'],
                 ['(','_{'],[')','}'],[',l',', l'],
                 ['v2',r'\nu_2'],['v4',r'\nu_4'],['v6',r'\nu_6'],['v7',r'\nu_7'],['3Sum_v',r'\sum_{3} \nu'],['v=',r'\nu=']]

    for x in x_replace:
        lineID = lineID.replace(x[0],x[1])

    lineID = lineID.replace(' ',' $',1)
    lineID = lineID.replace('$ $',' ')
    lineID += '$'

    return lineID
