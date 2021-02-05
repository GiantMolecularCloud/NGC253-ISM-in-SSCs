###################################################################################################
# helpers
###################################################################################################

def escape_fname(filename):
    """
    escape the ; symbol in file names
    """
    return filename.replace(';','_').replace('_#1','')

def get_samples(flux_uncertainty=0.10, size=100):
    """
    bootstrap error distribution by scaling the spectra
    sample 0 is always the actual data without scaling (scale factor 1.0)
    """
    samples = np.random.normal(loc=1.0, scale=flux_uncertainty, size=size+1)
    samples[0] = 1.0
    nums    = np.arange(size+1)
    joint_list = []
    for n,s in zip(nums, samples):
        joint_list.append([int(n),s])
    return joint_list

def bands_from_specie(specie):
    """
    return the band or bands that contain the given specie
    """
    # get lines for specie
    specie_lines = []
    for line in lines:
        if line['XCLASS']==specie:
            specie_lines.append(line)
    # get bands
    bands = []
    for line in specie_lines:
        if line['restfreq']<350*u.GHz:
            bands.append('LSB')
        elif line['restfreq']>350*u.GHz:
            bands.append('USB')
    # make unique
    bands = list(dict.fromkeys(bands))
    return bands

def get_final_files(num, SSC, specie):
    mol_file  = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules.molfit'))
    obs_file  = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_0','observation.xml'))
    man_mol   = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules_manual.molfit'))
    man_obs   = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_0','observation_manual.xml'))
    final_mol = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules_final.molfit'))
    final_obs = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'observation_final.xml'))

    os.system('mkdir -p '+escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num))))

    if os.path.exists(man_mol):
        os.system('cp -r '+man_mol+' '+final_mol)
    else:
        os.system('cp -r '+mol_file+' '+final_mol)

    if os.path.exists(man_obs):
        os.system('cp -r '+man_obs+' '+final_obs)
    else:
        os.system('cp -r '+obs_file+' '+final_obs)

    # replace run number in observations file
    os.system("sed -i 's/run_0/run_"+str(num)+"/g' "+final_obs)


###################################################################################################
#
###################################################################################################
