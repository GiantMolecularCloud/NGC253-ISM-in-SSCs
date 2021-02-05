###################################################################################################
# create observations file
###################################################################################################

def create_observations_file(num, SSC, specie, fit_range, fstep=2.5*u.MHz):

    if np.array(fit_range).shape == (2,):
        fit_range = [fit_range]

    observation = f"""<?xml version="1.0" encoding="UTF-8"?>
<ExpFiles>
<NumberExpFiles>{len(fit_range)}</NumberExpFiles>"""

    for frange in fit_range:
        fmin = frange[0]
        fmax = frange[1]

        if fmin=='min':
            if fmax<350:
                dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'LSB.dat'))
            elif fmax>350:
                dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'USB.dat'))
            frequency = np.genfromtxt(dat_file, usecols=(0))
            fmin = np.min(frequency)
        else:
            fmin = fmin*1000
        if fmax=='max':
            if fmin<350000:
                dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'LSB.dat'))
            elif fmin>350000:
                dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'USB.dat'))
            frequency = np.genfromtxt(dat_file, usecols=(0))
            fmax = np.max(frequency)
        else:
            fmax = fmax*1000

        if fmin<350000 and fmax<350000:
            dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'LSB.dat'))
        elif fmin>350000 and fmax>350000:
            dat_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'USB.dat'))
        else:
            raise ValueError("SSC "+str(SSC['no'])+", "+specie+": fmin and fmax are not on the same band: ["+str(fmin)+','+str(fmax)+']')

        observation += f"""
<file>
    <FileNamesExpFiles>{dat_file}</FileNamesExpFiles>
    <ImportFilter>xclassASCII</ImportFilter>
    <NumberExpRanges>1</NumberExpRanges>
    <FrequencyRange>
        <MinExpRange>{fmin}</MinExpRange>
        <MaxExpRange>{fmax}</MaxExpRange>
        <StepFrequency>{fstep.to(u.MHz).value}</StepFrequency>
        <t_back_flag>False</t_back_flag>
        <BackgroundTemperature>0.0</BackgroundTemperature>
        <TemperatureSlope>0.0</TemperatureSlope>
        <HydrogenColumnDensity>0.e+0</HydrogenColumnDensity>
        <DustBeta>0.0</DustBeta>
        <Kappa>0.0</Kappa>
    </FrequencyRange>
    <GlobalvLSR>0.0</GlobalvLSR>
    <TelescopeSize>0.15</TelescopeSize>
	<Inter_Flag>True</Inter_Flag>
	<ErrorY>no</ErrorY>
	<NumberHeaderLines>0</NumberHeaderLines>
	<SeparatorColumns> </SeparatorColumns>
</file>"""

    observation += f"""
<iso_flag>False</iso_flag>
</ExpFiles>"""

    # save to disk
    savepath = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'observation.xml'))
    os.system('mkdir -p '+os.path.dirname(savepath))
    with open(savepath, 'w') as f:
        f.write(observation)
    return savepath


###################################################################################################
#
###################################################################################################
