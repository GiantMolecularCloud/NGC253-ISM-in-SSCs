#######################
# GAS IN SSCS: XCLASS #
#######################

# Merge independent fitting, temperature line fitting and ro-vib line fitting.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))

data_I   = fnunpickle(os.path.join(mandir,'line_intensity_data.pickle'))
data_N   = fnunpickle(os.path.join(Xfinaldir, 'data.pickle'))
data_T   = fnunpickle(os.path.join(tempdir, 'temperature_data.pickle'))
data_vib = fnunpickle(os.path.join(vibdir, 'vibrational_data.pickle'))

temperature_species = fnunpickle(os.path.join(tempdir, 'temperature_species.pickle'))
vibrational_species = fnunpickle(os.path.join(vibdir, 'vibrational_species.pickle'))


###################################################################################################
# merge data
###################################################################################################

data_XCLASS = copy.deepcopy(data_N)
data_Gauss  = copy.deepcopy(data_I)

for SSC in tqdm(SSCs):
    for specie in data_XCLASS[SSC['num']].keys():

        if specie in temperature_species:
            # for temperature specie fits:
            #       T   fitted
            #       N   fitted
            #       w   (semi-)fixed
            #       v   fixed
            #       tau fitted implicitely
            for q in ['temperature','column density','peak opacity','integrated opacity']:
                data_XCLASS[SSC['num']][specie][q] = data_T[SSC['num']][specie][q]

        if specie in vibrational_species:
            # for temperature specie fits:
            #       T   fixed at 300K
            #       N   fitted
            #       w   (semi-)fixed
            #       v   fixed
            #       tau fitted implicitely
            for q in ['temperature','column density','peak opacity','integrated opacity']:
                data_XCLASS[SSC['num']][specie][q] = data_vib[SSC['num']][specie][q]

os.system('mkdir -p '+resultsdir)
fnpickle(data_XCLASS, os.path.join(resultsdir, 'data_XCLASS.pickle'))
fnpickle(data_Gauss, os.path.join(resultsdir, 'data_Gauss.pickle'))


###################################################################################################
# filter bad fits
###################################################################################################

# bad fits are most easily identified by way too low column densities
# remove a fit if its column density is lower than 10% of the median
# replace these values with the sample median to keep the data structure equal in size, update in place

for s, ss in tqdm(data_XCLASS.items()):
    for spx, data in ss.items():
        for comp in np.arange(len(data['velocity']['median'])):

            # get indeces of bad fits: column density
            all    = data['column density']['all'][comp]
            median = data['column density']['median'][comp]
            bad_cd = np.append( np.where(all < 0.1*median)[0], np.where(all > 10*median)[0] )

            if not spx=='CO;v=0' or (s in ['1','2'] and spx=='CS;v=0'):
                # get indeces of bad fits: linewidth
                all    = data['linewidth']['all'][comp]
                median = data['linewidth']['median'][comp]
                bad_lw = np.append( np.where(all > 0.*median+59.)[0], np.where(all < 0.*median+11.)[0] )
                if len(bad_lw) > 50:
                    bad_lw = np.array([], dtype='int64')
            else:
                bad_lw = np.array([], dtype='int64')

            bad_idx = np.append(bad_cd, bad_lw)

            # for all quantities
            for q,quantity in data.items():
                # get good statistics (without the bad values)
                good_p16, good_median, good_p84 = np.percentile( np.delete(quantity['all'][comp],bad_idx), (16,50,84) )
                quantity['median'][comp] = good_median
                quantity['16th'][comp]   = good_p16
                quantity['84th'][comp]   = good_p84
                # replace bad values with good median
                for bi in bad_idx:
                    # quantity['all'][comp][bi] = good_median
                    data_XCLASS[s][spx][q]['all'][comp][bi] = good_median


fnpickle(data_XCLASS, os.path.join(resultsdir, 'data_XCLASS.pickle'))


###################################################################################################
#
###################################################################################################
