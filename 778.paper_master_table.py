###############################
# GAS IN SSCS: data structure #
###############################

# build a simple data structure using astropy tables


###################################################################################################
# get species data
###################################################################################################

# get column data
###################################################################################################

columns = {'SSC':                {'values': [], 'unit': None},
           'species':            {'values': [], 'unit': None},
           'component':          {'values': [], 'unit': None},
           'column density':     {'values': [], 'unit': 'cm^-2'},
           'linewidth':          {'values': [], 'unit': 'km/s'},
           'velocity':           {'values': [], 'unit': 'km/s'},
           'temperature':        {'values': [], 'unit': 'K'},
           'integrated opacity': {'values': [], 'unit': None},
           'peak opacity':       {'values': [], 'unit': None}
           }

for ssc in all_data.keys():
    for specie in all_data[ssc].keys():
        components = all_data[ssc][specie]['components']
        for component in np.arange(components):

            # append base data
            columns['SSC']['values'].append(       ssc       )
            columns['species']['values'].append(   specie   )
            columns['component']['values'].append( component )

            # append fit data
            for quantity in ['column density', 'linewidth', 'velocity', 'temperature', 'integrated opacity', 'peak opacity']:
                p16    = all_data[ssc][specie][quantity]['16th'][component]
                median = all_data[ssc][specie][quantity]['median'][component]
                p84    = all_data[ssc][specie][quantity]['84th'][component]

                columns[quantity]['values'].append( [p16,median,p84] )


# construct table
###################################################################################################

data_species = QTable()
for colname,coldata in columns.items():
    col = Column(data = coldata['values'],
                 name = colname,
                 unit = coldata['unit']
                )
    data_species.add_column(col)

data_species.write(join(subprojectdir,'data_species.fits'), overwrite=True)
data_species.write(join(plotdir,'data_species.txt'), format='ascii.aastex', overwrite=True)


###################################################################################################
# get transitions data
###################################################################################################

# get column data
###################################################################################################

columns = {'molecule':                    {'values': [], 'unit': None},
           'rotation':                    {'values': [], 'unit': None},
           'vibration':                   {'values': [], 'unit': None},
           'component':                   {'values': [], 'unit': None},
           'SSC 1 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 2 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 3 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 4 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 5 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 6 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 7 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 8 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 9 integrated intensity':  {'values': [], 'unit': 'K km/s'},
           'SSC 10 integrated intensity': {'values': [], 'unit': 'K km/s'},
           'SSC 11 integrated intensity': {'values': [], 'unit': 'K km/s'},
           'SSC 12 integrated intensity': {'values': [], 'unit': 'K km/s'},
           'SSC 13 integrated intensity': {'values': [], 'unit': 'K km/s'},
           'SSC 14 integrated intensity': {'values': [], 'unit': 'K km/s'},
           }

def get_max_components(specie):
    comps = []
    for s in all_data.keys():
        try:
            comps.append( all_data[s][specie]['components'] )
        except:
            pass
    return np.max(comps)

for specie in all_data['14'].keys():
    molecule = specie.split(';')[0]
    for rotation in all_data['14'][specie]['integrated intensity'].keys():
        for vibration in all_data['14'][specie]['integrated intensity'][rotation].keys():
            max_comps = get_max_components(specie)
            for component in np.arange(max_comps):

                # append base data
                columns['molecule']['values'].append( molecule )
                columns['rotation']['values'].append( rotation )
                columns['vibration']['values'].append( vibration )
                columns['component']['values'].append( component )

                for ssc in all_data.keys():

                    # append fit data
                    try:
                        p16    = all_data[ssc][specie]['integrated intensity'][rotation][vibration][component]['16th'].value
                        median = all_data[ssc][specie]['integrated intensity'][rotation][vibration][component]['median'].value
                        p84    = all_data[ssc][specie]['integrated intensity'][rotation][vibration][component]['84th'].value
                    except:
                        p16, median, p84 = np.nan, np.nan, np.nan

                    columns['SSC '+ssc+' integrated intensity']['values'].append( [p16,median,p84] )


# construct table
###################################################################################################

data_transitions = QTable()
for colname,coldata in columns.items():
    col = Column(data = coldata['values'],
                 name = colname,
                 unit = coldata['unit']
                )
    data_transitions.add_column(col)

data_transitions.write(join(subprojectdir,'data_transitions.fits'), overwrite=True)
data_transitions.write(join(plotdir,'data_transitions.txt'), format='ascii.aastex', overwrite=True)


###################################################################################################
#
###################################################################################################
