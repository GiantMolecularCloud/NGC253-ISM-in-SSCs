def specie_tex(specie):
    """Convert line dictionary to LaTeX string.

    Parameters
    ----------
    specie : string
        Specie label as used in XCLASS

    Returns
    -------
    string
        Formatted LaTeX compatible string.
    """

    replace = [['O2','O$_2$'],['C-13-','$^{13}$C'],['O+','O$^+$'],['N-15-','$^{15}$N'],['N-15','$^{15}$N'],['C3','C$_3$'],['S-33-','$^{33}$S'],['S-34-','$^{34}$S'],['O-18-','$^{18}$O'],['O-18','$^{18}$O'],['H2','H$_2$'],['CCC','C$_3$'],
               ['O2','O$_2$'],['13C','$^{13}$C'],  ['O+','O$^+$'],['15N','$^{15}$N'], ['32S','$^{32}$S'],  ['C3','C$_3$'],['33S','$^{33}$S'],  ['34S','$^{34}$S'],  ['18O','$^{18}$O'], ['12C','$^{12}$C'],  ['H2','H$_2$'],['C3','C$_3$'], ['14N','$^{14}$N'],
               [';',' '],['#1',''],['#2',''],
               ['v=0',''],['v2=1',r'$\nu_2=1$'],['v2=2',r'$\nu_2=2$'],['v4=1',r'$\nu_4=1$'],['v6=1',r'$\nu_6=1$'],['v6=2',r'$\nu_6=2$'],['v7=1',r'$\nu_7=1$'],['v7=2',r'$\nu_7=2$']]

    for r in replace:
        specie = specie.replace(r[0],r[1])

    return specie


def molecule_tex(molecule):
    replace = [['O2','O$_2$'],['C-13-','$^{13}$C'],['O+','O$^+$'],['N-15-','$^{15}$N'],['N-15','$^{15}$N'],['C3','C$_3$'],['S-33-','$^{33}$S'],['S-34-','$^{34}$S'],['O-18-','$^{18}$O'],['O-18','$^{18}$O'],['H2','H$_2$'],['CCC','C$_3$'],
               ['O2','O$_2$'],['13C','$^{13}$C'],  ['O+','O$^+$'],['15N','$^{15}$N'], ['32S','$^{32}$S'],  ['C3','C$_3$'],['33S','$^{33}$S'],  ['34S','$^{34}$S'],  ['18O','$^{18}$O'], ['12C','$^{12}$C'],  ['H2','H$_2$'],['C3','C$_3$'], ['14N','$^{14}$N']]
    for r in replace:
        molecule = molecule.replace(r[0],r[1])
    return molecule


def rotation_tex(rotation):
    replace = [['-','--']]
    for r in replace:
        rotation = rotation.replace(r[0],r[1])
    rotation = '}$'.join(('$_{'.join(rotation.split('('))).split(')'))
    return rotation


def vibration_tex(vibration):
    replace = [['3Sum_v=0',r'$\nu=0$'],['v=0',r'$\nu=0$'],['v2=1',r'$\nu_2=1$'],['v2=2',r'$\nu_2=2$'],['v4=1',r'$\nu_4=1$'],['v6=1',r'$\nu_6=1$'],['v6=2',r'$\nu_6=2$'],['v7=1',r'$\nu_7=1$'],['v7=2',r'$\nu_7=2$']]
    for r in replace:
        vibration = vibration.replace(r[0],r[1])
    return vibration
