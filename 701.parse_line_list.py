def parse_line_list(line_file, vsys=250*u.km/u.s):
    """Parse a line list.

    Parameters
    ----------
    line_file : string
        Txt file containing the lines and their parameters.
    vsys      : astropy velocity
        Systemic velocity. Optional. Needed to convert observed frequency.

    Returns
    -------
    list of dictionaries
        Each dict contains the keys "molecule" (name), "transistion", "restfreq", "obsfreq" (based on
        given vsys), "confidence" (how certain is the detection).
    """

    # load line file
    lines = []
    line_table = np.genfromtxt(line_file,
                 delimiter=(11,)+(25,)+(13,)+(15,)+(14,)+(29,)+(5,)+(20,),
                 dtype=None,
                 comments='%'
                )

    # clean up und convert
    for row in line_table:
        if not ( row[0] == b'\n' ):
            molecule   = row[0].strip().decode('UTF-8')
            transition = row[1].strip().decode('UTF-8')
            vibration  = row[2].strip().decode('UTF-8')
            restfreq   = float(row[3])*u.GHz
            Eupper     = float(row[4])*u.K
            confidence = row[5].strip().decode('UTF-8')
            seenby     = row[6].strip().decode('UTF-8')
            XCLASS     = row[7].strip().decode('UTF-8').replace('\n','')
            ID         = molecule+' '+transition
            if not vibration=='':
                ID += ' '+vibration

            if not ( vsys == None ):
                obsfreq = restfreq/(1.+vsys/c.c)                 # optical convention
                lines.append({'ID': ID, 'molecule': molecule, 'transition': transition, 'vibration': vibration, 'restfreq': restfreq, 'obsfreq': obsfreq, 'Eupper': Eupper, 'confidence': confidence, 'seen_by': seenby, 'XCLASS': XCLASS})
            else:
                lines.append({'ID': ID, 'molecule': molecule, 'transition': transition, 'vibration': vibration, 'restfreq': restfreq,                     'Eupper': Eupper, 'confidence': confidence, 'seen_by': seenby, 'XCLASS': XCLASS})

    return lines
