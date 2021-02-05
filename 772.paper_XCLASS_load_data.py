#######################
# GAS IN SSCS: XCLASS #
#######################


###################################################################################################
# load fitted data
###################################################################################################

def parse_molfit(SSCnum):

    import re

    def load_molfit(num, SSCnum):
        file = refitfiles[SSCnum][str(num)]['model']['molfit']
        with open(file, 'r') as f:
            fdata = f.read()
            return fdata

    def filter_vals():
        # # consider fit good only if
        # #       temperature not within 2K of limits
        # #       column density not within 5% of limits
        # # ignore this selection for fixed fits
        # good_fit  = (t_val>t_llim+2 and t_val<t_ulim+2 and N_val>N_llim*1.05 and N_val<N_ulim*0.95 and N_val>2e12)
        # fixed_fit = (N_ulim/N_llim<5)
        # if good_fit or fixed_fit:
        data[spx]['temperature']['all'][ncomp].append(t_val)
        data[spx]['column density']['all'][ncomp].append(N_val)
        data[spx]['linewidth']['all'][ncomp].append(w_val)
        data[spx]['velocity']['all'][ncomp].append(v_val)
        # else:
        #     data[spx]['temperature']['all'][ncomp].append(np.nan)
        #     data[spx]['column density']['all'][ncomp].append(np.nan)
        #     data[spx]['linewidth']['all'][ncomp].append(np.nan)
        #     data[spx]['velocity']['all'][ncomp].append(np.nan)

    def split_at_molecule(mofit):
        return re.split(u'\n(?=[A-Z])', molfit)

    # load model molfit files
    molfits = [ load_molfit(num, SSCnum) for num in nums ]

    data = {}
    for molfit in molfits:

        # separate species
        species = split_at_molecule(molfit)
        for specie in species:
            spx = specie.split()[0]
            # spx = re.sub(u'(;$|;#1)', '', spx)
            spx = re.sub(u';$', '', spx)

            if spx==refitspecies:
                try:
                    data[spx]
                except:
                    data[spx] = {q: {'all': []} for q in ['temperature', 'column density', 'linewidth', 'velocity']}

                # separate components
                components = specie.split('\n')[1:]

                for ncomp,component in enumerate(filter(None,components)):
                    c = component.split()

                    # get values
                    t_llim = float(c[5])
                    t_ulim = float(c[6])
                    t_val  = float(c[7])
                    N_llim = float(c[9])
                    N_ulim = float(c[10])
                    N_val  = float(c[11])
                    w_llim = float(c[13])
                    w_ulim = float(c[14])
                    w_val  = float(c[15])
                    v_llim = float(c[17])
                    v_ulim = float(c[18])
                    v_val  = float(c[19])

                    try:
                        data[spx]['temperature']['all'][ncomp]
                    except:
                        for q in ['temperature', 'column density', 'linewidth', 'velocity']:
                            data[spx][q]['all'].append([])

                    filter_vals()

    # get percentiles
    for spx,specie in data.items():
        for q,quantity in specie.items():
            for ncomp,component in enumerate(quantity['all']):
                p16,median,p84 = np.percentile(component, (16,50,84))
                try:
                    data[spx][q]['fit']
                except:
                    data[spx][q]['fit']    = []
                    data[spx][q]['16th']   = []
                    data[spx][q]['median'] = []
                    data[spx][q]['84th']   = []
                data[spx][q]['fit'].append( component[0] )
                data[spx][q]['16th'].append( p16 )
                data[spx][q]['median'].append( median )
                data[spx][q]['84th'].append( p84 )

    return data


def get_opacity(data_table):

    import re

    def load_opacity_files(num, SSCnum):
        taupath = escape_fname(join(refitdir,refitspecies,'SSC_'+SSCnum,'run_'+str(num),'model','opacity.pickle'))
        with open(taupath, 'rb') as f:
            t_data = pickle.load(f, encoding="latin1")
            return t_data

    def extract_opacities(t_dat):
        int = t_dat[0]
        tau = t_dat[1]
        for specie in tau:
            spx = specie[0]
            # spx = re.sub(u'(;$|;#1)', '', spx)
            spx = re.sub(u';$', '', spx)
            c   = specie[1]-1
            opt = specie[2][:,1]
            int_opt  = np.sum(opt)
            peak_opt = np.max(opt)
            try:
                data_table[SSCnum][spx]['peak opacity']['all'][c].append(peak_opt)
                data_table[SSCnum][spx]['integrated opacity']['all'][c].append(int_opt)
            except:
                print(spx,c)

    def get_ncomps(SSCnum, spx):
        tau = load_opacity_files(0, SSCnum)[1]
        if spx+';' in [t[0] for t in tau]:
            ncomps = len([t for t in tau if spx+';'==t[0]])
        else:
            ncomps = len([t for t in tau if spx==t[0]])
        return ncomps

    def prepare_table(SSCnum):
        for spx, specie in data_table[SSCnum].items():
            try:
                specie['peak opacity']
                specie['integrated opacity']
            except:
                ncomps = get_ncomps(SSCnum, spx)
                specie['peak opacity']       = {k: [[] for _ in np.arange(ncomps)] for k in ['all','16th','median','84th']}
                specie['integrated opacity'] = {k: [[] for _ in np.arange(ncomps)] for k in ['all','16th','median','84th']}

    def calc_percentiles(SSCnum):
        for spx in data_table[SSCnum].keys():
            for c,comp in enumerate(data_table[SSCnum][spx]['peak opacity']['all']):
                p16,median,p84 = np.percentile(comp, (16,50,84))
                data_table[SSCnum][spx]['peak opacity']['16th'][c]   = p16
                data_table[SSCnum][spx]['peak opacity']['median'][c] = median
                data_table[SSCnum][spx]['peak opacity']['84th'][c]   = p84
            for c,comp in enumerate(data_table[SSCnum][spx]['integrated opacity']['all']):
                p16,median,p84 = np.percentile(comp, (16,50,84))
                data_table[SSCnum][spx]['integrated opacity']['16th'][c]   = p16
                data_table[SSCnum][spx]['integrated opacity']['median'][c] = median
                data_table[SSCnum][spx]['integrated opacity']['84th'][c]   = p84

    for SSCnum in tqdm(refitfiles.keys()):
        prepare_table(SSCnum)
        t_data = [load_opacity_files(num,SSCnum) for num in nums]
        for t_dat in t_data:
            extract_opacities(t_dat)
        calc_percentiles(SSCnum)


refit_data = {SSCnum: parse_molfit(SSCnum) for SSCnum in refitfiles.keys()}
get_opacity(refit_data)         # ignore non-existent species
fnpickle(refit_data, escape_fname(join(refitdir, refitspecies, 'refit_data.pickle')))


###################################################################################################
# merge refit data into older data structure
###################################################################################################

data_XCLASS_refit = copy.deepcopy(data_XCLASS)

for SSCnum in refit_data.keys():
    for spx in refit_data[SSCnum].keys():
        data_XCLASS_refit[SSCnum][spx] = copy.deepcopy(refit_data[SSCnum][spx])

# fix the renamed species: SO;v=0;#1 --> SO;v=0;#2
for _,data in data_XCLASS_refit.items():
    if 'SO;v=0;#1' in data.keys():
        data['SO;v=0;#2'] = data.pop('SO;v=0;#1')

fnpickle(data_XCLASS_refit, escape_fname(join(refitdir, refitspecies, 'data_XCLASS_refit.pickle')))


###################################################################################################
# FOR RE_RUN ONLY: merge the two H13CN runs
###################################################################################################

part_1 = fnunpickle(escape_fname(join(refitdir, 'HC-13-N_v=0', 'data_XCLASS_refit_part1.pickle')))
part_2 = fnunpickle(escape_fname(join(refitdir, 'HC-13-N_v=0', 'data_XCLASS_refit_part2.pickle')))

merged = copy.deepcopy(part_1)
for SSCnum in ['1','2','3','4','6','10','12','13']:
    merged[SSCnum]['HC-13-N;v=0'] = copy.deepcopy(part_2[SSCnum]['HC-13-N;v=0'])
for SSCnum in ['1','2','3','4','10','12','13']:
    merged[SSCnum]['SO2;v=0'] = copy.deepcopy(part_2[SSCnum]['SO2;v=0'])

fnpickle(merged, escape_fname(join(refitdir, 'HC-13-N;v=0', 'data_XCLASS_merged.pickle')))


###################################################################################################
# FOR RE_RUN ONLY: merge all re-fitted lines
###################################################################################################

data_XCLASS_refit = copy.deepcopy(data_XCLASS)

# fix the renamed species: SO;v=0;#1 --> SO;v=0;#2
for _,data in data_XCLASS_refit.items():
    if 'SO;v=0;#1' in data.keys():
        data['SO;v=0;#2'] = data.pop('SO;v=0;#1')


updated = {'CO;v=0':      {'SSCs': ['3','6','7','9','14'],                'data': fnunpickle(escape_fname(join(refitdir, 'CO_v=0', 'data_XCLASS_refit.pickle')))},
           'HC-13-N;v=0': {'SSCs': ['1','2','3','4','5','6','8','9','10','11','12','13','14'], 'data': fnunpickle(escape_fname(join(refitdir, 'HC-13-N_v=0', 'data_XCLASS_merged.pickle')))},
           'HCN;v2=1':    {'SSCs': ['11'],                                'data': fnunpickle(escape_fname(join(refitdir, 'HCN_v2=1', 'data_XCLASS_refit.pickle')))},
           'SO2;v=0':     {'SSCs': ['1','2','3','4','10','11','12','13'], 'data': fnunpickle(escape_fname(join(refitdir, 'HC-13-N_v=0', 'data_XCLASS_merged.pickle')))}
          }

for spx,dat in updated.items():
    for SSCnum in dat['SSCs']:
        data_XCLASS_refit[SSCnum][spx] = dat['data'][SSCnum][spx]

fnpickle(data_XCLASS_refit, escape_fname(join(refitdir, 'data_XCLASS_refit.pickle')))


###################################################################################################
# FOR RE_RUN ONLY: remove already fitted species/SSCs
###################################################################################################

data_only_refit = {s: {} for s in list(dict.fromkeys([i for _,x in updated.items() for i in x['SSCs']]))}

for spx,dat in updated.items():
    for SSCnum in dat['SSCs']:
        data_only_refit[SSCnum][spx] = copy.deepcopy(dat['data'][SSCnum][spx])

fnpickle(data_only_refit, escape_fname(join(refitdir, 'data_only_refit.pickle')))


###################################################################################################
#
###################################################################################################
