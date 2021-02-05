###################################################################################################
# create molfit file
###################################################################################################

def create_molfile_file(SSC, specie, components=1, fit_with=None):

    def lpos_deviation(specie):
        if specie=='CO;v=0':
            return 150.0
        elif specie=='HCN;v=0':
            return 50.0
        elif specie=='HCO+;v=0':
            return 50.0
        elif specie=='CS;v=0':
            return 50.0
        else:
            return 10.0

    def width(specie):
        if specie=='CO;v=0':
            return 80.0
        elif specie=='HCN;v=0':
            return 50.0
        elif specie=='HCO+;v=0':
            return 50.0
        elif specie=='CS;v=0':
            return 50.0
        else:
            return 50.0

    def coldens(specie):
        if specie=='CO;v=0':
            return ['1.0E+15','1.0E+23','1.0E+19']
        elif specie=='HCN;v=0':
            return ['1.0E+12','1.0E+20','1.0E+18']
        elif specie=='HCO+;v=0':
            return ['1.0E+12','1.0E+20','1.0E+18']
        elif specie=='CS;v=0':
            return ['1.0E+12','1.0E+20','1.0E+17']
        elif ('=1' in specie) or ('=2' in specie):
            return ['1.0E+12','1.0E+18','1.0E+14']
        else:
            return ['1.0E+12','1.0E+19','1.0E+16']

    def fix_temp(specie):
        if specie in reliable_Ts:
            return "y  25.0  1000.0  130.0"
        else:
            return "n  25.0  1000.0  130.0"

    if '#1' in specie:
        spx = specie
    else:
        spx = specie+';'

    molfit = f"""% fit parameter?(y/n)   lower_limit upper_limit starting_value
% source_size rotation_temperature column_density line_width velocity_offset
{spx}   {components}"""

    # write multiple components if components>1
    lpos = lpos_deviation(specie)
    w    = width(specie)
    col1,col2,col3 = coldens(specie)
    temp = fix_temp(specie)
    for c in np.arange(components):
        molfit += f"""
n  0.01  10.0  1.0   {temp}   y  {col1}  {col2}  {col3}   y  20.0  {w}  25.0   y  {-lpos}  {lpos}  0.0   c"""

    # add additional lines if blended
    if fit_with is not None:
        for fw in fit_with:
            lp = lpos_deviation(fw)
            w    = width(fw)
            col1,col2,col3 = coldens(fw)
            temp = fix_temp(fw)
            molfit += f"""
{fw};   1
n  0.01  10.0  1.0   {temp}   y  {col1}  {col2}  {col3}   y  20.0  {w}  25.0   y  {-lp}  {lp}  0.0   c"""

    # save to disk
    savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules.molfit'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(molfit)
    return savepath


###################################################################################################
#
###################################################################################################
