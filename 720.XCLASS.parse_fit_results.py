###################################################################################################
# parse fit results
###################################################################################################

def parse_molfit(SSC, specie, return_all=False, mode='final'):
    """
    Load bestfit molfit files and calculate bestfit of original data, 16th percentile, median and
    84th percentile for temperature, column density, linewidth and velocity. If multiple components
    are present list all of that for all present components.
    If return_all=True return all fit parameters and percentiles instead of only percentiles.
    """

    if mode=='final':
        modestr = '_final'
    elif mode=='manual':
        modestr = '_manual'
    else:
        modestr = ''

    # load bestfit molfit file
    molfit_files = glob.glob(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_*','results','molecules'+modestr+'__LM__call_1.out.molfit')))
    molfit_files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    bestfits = {'temperature': [], 'column density': [], 'linewidth': [], 'velocity': []}
    for f in molfit_files:

        # handle specie names correctly
        if '#1' in specie:
            spx = specie
        else:
            spx = specie+';'

        # load molfit file
        molfit = open(f, 'r').readlines()
        for idx,ll in enumerate(molfit):
            ll = ll.replace('\n','')
            ll = [l for l in ll.split(' ') if not l=='']
            molfit[idx] = ll

        # get which lines contain the components of this specie
        for idx,ll in enumerate(molfit):
            if ll[0]==spx:
                idx_start = idx+1
            if ll[0]!='n' and idx>1:
                idx_stop  = idx
                break
            if idx==len(molfit)-1:
                idx_stop  = idx+1

        try:
            idx_start, idx_stop
        except:
            print(SSC, specie)

        # get values
        for idx in np.arange(idx_start,idx_stop):
            for name in ['temperature','column density','linewidth','velocity']:
                try:
                    bestfits[name][idx-1]
                except:
                    bestfits[name].append([])

            t_llim = float(molfit[idx][5])
            t_ulim = float(molfit[idx][6])
            t_val  = float(molfit[idx][7])
            N_llim = float(molfit[idx][9])
            N_ulim = float(molfit[idx][10])
            N_val  = float(molfit[idx][11])
            w_llim = float(molfit[idx][13])
            w_ulim = float(molfit[idx][14])
            w_val  = float(molfit[idx][15])
            v_llim = float(molfit[idx][17])
            v_ulim = float(molfit[idx][18])
            v_val  = float(molfit[idx][19])

            # consider fit good only if
            #       temperature not within 2K of limits
            #       column density not within 5% of limits
            # good_fit  = (t_val>t_llim+2 and t_val<t_ulim+2 and N_val>N_llim*1.05 and N_val<N_ulim*0.95 and N_val>2e12)
            # good_fit  = (N_val>N_llim*1.05 and N_val<N_ulim*0.95)
            # fixed_fit = (N_ulim/N_llim<5)
            # an_exception = True if (SSC['no']==3 and specie=='HC-13-N;v=0') else False
            # if good_fit or fixed_fit or an_exception:
            bestfits['temperature'][idx-1].append(t_val)
            bestfits['column density'][idx-1].append(N_val)
            bestfits['linewidth'][idx-1].append(w_val)
            bestfits['velocity'][idx-1].append(v_val)

    # get percentiles
    percentiles = {'temperature': [], 'column density': [], 'linewidth': [], 'velocity': []}
    for name in percentiles.keys():
        for idx,values in enumerate(bestfits[name]):
            lower,median,upper = np.percentile(values, (16,50,84))
            percentiles[name].append([values[0], lower, median, upper])

    if return_all:
        return bestfits,percentiles
    else:
        return percentiles


def parse_log(SSC, specie):
    """
    Parse log file and report chi2 evolution for the Levenberg-Marqard iterations. Returns chi2 for
    fit to original data, 16th percentile, median and 84th percentile for all iterations.
    The chi2 file is incomplete randomly for no reason! Have to use the log file instead.
    """

    # load LM log file
    log_files = glob.glob(escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_*','results','fit__LM__call_1.log')))
    log_files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    chi2s = []
    for f in log_files:
        log = np.genfromtxt(f,
                            skip_header = 17,
                            skip_footer = 3
                           )
        chi2 = log[:,1]
        chi2s.append(chi2)
    chi2s = np.transpose(chi2s)

    # get percentiles
    percentiles = []
    for it,values in enumerate(chi2s):
        lower, median, upper = np.percentile(values, (16,50,84))
        percentiles.append([values[0], lower, median, upper])

    return np.array(percentiles)


def get_physical(SSCs, species, mode='final'):
    """
    Get a dictionary of the physical parameters for the specified SSCs and species.
    """

    line_data = {}

    for SSC in tqdm(SSCs):
        line_data[str(SSC['no'])] = {}

        for specie in species:
            line_data[str(SSC['no'])][specie] = {}

            try:
                bestfits,percentiles = parse_molfit(SSC, specie, return_all=True, mode=mode)

                for Q in ['temperature','column density','linewidth','velocity']:
                    line_data[str(SSC['no'])][specie][str(Q)] = {'fit':    [p[0] for p in  percentiles[str(Q)]],
                                                                 '16th':   [p[1] for p in  percentiles[str(Q)]],
                                                                 'median': [p[2] for p in  percentiles[str(Q)]],
                                                                 '84th':   [p[3] for p in  percentiles[str(Q)]],
                                                                 'all':    bestfits[str(Q)]}
            except ValueError:
                print("SSC "+str(SSC['no'])+" "+specie+"\tCould not load due to varying columns in molfit file.")
            except TypeError:
                print("SSC "+str(SSC['no'])+" "+specie+"\tIteration over 0-d array")
            except:
                print("Problem in SSC "+str(SSC['no'])+" "+specie)
    return line_data


###################################################################################################
#
###################################################################################################
