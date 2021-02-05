#######################
# GAS IN SSCS: XCLASS #
#######################

# Prepare the data in python3 in a way that keeps the XCLASS code simple and manageable.


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(scriptdir, '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(mandir, 'spectra.pickle'))
bandfit = fnunpickle(os.path.join(mandir, 'band_fit_parameters.pickle'))
linedat = fnunpickle(os.path.join(mandir,'linedat.pickle'))

fitable_species = ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v7=2', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'S-34-O;v=0', 'SO-18;v=0', 'SO2;v=0', 'S-34-O2;v=0', 'H2CS;v=0;#1']
reliable_Ts     = ['SO2;v=0','S-34-O2;v=0','H2CS;v=0;#1']
fnpickle(fitable_species, os.path.join(XCLASSdir, 'fitable_species.pickle'))
fnpickle(reliable_Ts, os.path.join(XCLASSdir, 'reliable_Ts.pickle'))


###################################################################################################
# load functions
###################################################################################################

execfile(os.path.join(scriptdir, '720.XCLASS.create_CASA_script.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.create_model_spectrum.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.create_molfit_file.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.create_observations_file.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.create_XCLASS_script.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.export_data.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.helpers.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.parse_fit_results.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_combined_spectrum.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_convergence.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_parameter_histogram.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_parameters_N.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_ratios_N.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_restfreq_spectrum.py'))
execfile(os.path.join(scriptdir, '720.XCLASS.plot_spectrum_fit.py'))


###################################################################################################
# run single XCLASS fit to get parameters
###################################################################################################

# do this automatic approach only for the lines that can be fit independently (well separated lines
# and lines that can be separated by appropriate cuts in the fit_ranges table)

components = {
'1':  {'CO;v=0': 6, 'CS;v=0': 2, 'HCN;v=0': 2, 'HCO+;v=0': 2, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'2':  {'CO;v=0': 5, 'CS;v=0': 2, 'HCN;v=0': 2, 'HCO+;v=0': 2, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'3':  {'CO;v=0': 6, 'CS;v=0': 3, 'HCN;v=0': 2, 'HCO+;v=0': 2, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'4':  {'CO;v=0': 5, 'CS;v=0': 1, 'HCN;v=0': 2, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'5':  {'CO;v=0': 5, 'CS;v=0': 1, 'HCN;v=0': 1, 'HCO+;v=0': 1, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'6':  {'CO;v=0': 5, 'CS;v=0': 1, 'HCN;v=0': 1, 'HCO+;v=0': 2, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'7':  {'CO;v=0': 5, 'CS;v=0': 2, 'HCN;v=0': 3, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'8':  {'CO;v=0': 6, 'CS;v=0': 1, 'HCN;v=0': 4, 'HCO+;v=0': 4, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'9':  {'CO;v=0': 6, 'CS;v=0': 1, 'HCN;v=0': 3, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'10': {'CO;v=0': 5, 'CS;v=0': 2, 'HCN;v=0': 3, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'11': {'CO;v=0': 7, 'CS;v=0': 1, 'HCN;v=0': 3, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'12': {'CO;v=0': 6, 'CS;v=0': 1, 'HCN;v=0': 3, 'HCO+;v=0': 4, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'13': {'CO;v=0': 6, 'CS;v=0': 3, 'HCN;v=0': 4, 'HCO+;v=0': 3, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1},
'14': {'CO;v=0': 5, 'CS;v=0': 2, 'HCN;v=0': 2, 'HCO+;v=0': 2, 'H2CS;v=0;#1': 1, 'HC-13-N;v=0': 1, 'HCN-15;v=0': 1, 'HN-15-C;v=0': 1, 'SO;v=0;#1': 1, 'S-33-O;v=0': 1, 'SO2;v=0': 1, 'HCN;v2=2': 1, 'HCN;v2=1': 1, 'HC-13-N;v=0': 1, 'HCCCN;v=0': 1, 'HCCCN;v7=2': 1, 'HCCCN;v6=1': 1, 'HCCCN;v7=1': 1, 'HCCCN;v6=2': 1, 'HN-15-C;v=0': 1, 'S-34-O;v=0': 1, 'SO-18;v=0': 1, 'S-34-O2;v=0': 1}
}
fnpickle(components, os.path.join(XCLASSdir, 'components.pickle'))

fit_ranges = {
'1':  {'CO;v=0': [345.62,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.90], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'2':  {'CO;v=0': [345.62,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.90], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'3':  {'CO;v=0': [345.63,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.90], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'4':  {'CO;v=0': [345.50,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'5':  {'CO;v=0': [345.50,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'6':  {'CO;v=0': [345.50,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'7':  {'CO;v=0': [345.50,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'8':  {'CO;v=0': [345.68,346.10], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.68], 'HCO+;v=0': [356.50,356.91], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [[345.58,345.70],[354.66,354.80]], 'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'9':  {'CO;v=0': [345.50,345.95], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [[345.58,345.70],[354.66,354.80]], 'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'10': {'CO;v=0': [345.65,346.05], 'CS;v=0': [242.80,343.05], 'HCN;v=0': ['min',354.75], 'HCO+;v=0': [356.50,357.00], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'11': {'CO;v=0': [345.50,345.95], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.85], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'12': {'CO;v=0': [345.50,345.95], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.90], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [354.66,354.80],                   'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'13': {'CO;v=0': [345.68,346.00], 'CS;v=0': [242.80,343.00], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.88], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [[345.58,345.70],[354.66,354.80]], 'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]},
'14': {'CO;v=0': [345.70,345.90], 'CS;v=0': [242.80,343.05], 'HCN;v=0': ['min',354.62], 'HCO+;v=0': [356.50,356.88], 'H2CS;v=0;#1': [343.12,343.6], 'HC-13-N;v=0': [345.30,345.50], 'HCN-15;v=0': [344.05,344.26], 'HN-15-C;v=0': [355.35,355.52], 'SO;v=0;#1': [344.25, 344.50], 'S-33-O;v=0': [343.02,343.15], 'SO2;v=0': [[342.70,342.80],[355.00,355.23],[357.12,'max']], 'HCN;v2=2': [356.10,356.40], 'HCN;v2=1': [356.20,356.40], 'HCCCN;v=0': [[345.58,345.70],[354.66,354.80]], 'HCCCN;v6=1': [[355.22,355.32],[355.50,355.70]], 'HCCCN;v6=2': [[355.78,355.88],[356.00,356.20]], 'HCCCN;v7=1': [[355.50,355.70],[356.00,356.20]], 'HCCCN;v7=2': [[343.70,343.90],[356.00,356.20],[356.88,357.12]], 'HN-15-C;v=0': [355.38,355.50], 'S-34-O;v=0': [343.70,343.95], 'SO-18;v=0': [355.50,355.70], 'S-34-O2;v=0': [[344.50,345.30],[345.90,'max'],[357.02,357.15]]}
}
fnpickle(fit_ranges, os.path.join(XCLASSdir, 'fit_ranges.pickle'))

fit_with = {
'CO;v=0':       None,
'CS;v=0':       None,
'HCN;v=0':      None,
'HCO+;v=0':     None,
'H2CS;v=0;#1':  None,
'HC-13-N;v=0':  None,
'HCN-15;v=0':   None,
'HN-15-C;v=0':  None,
'SO;v=0;#1':    None,
'S-33-O;v=0':   None,
'SO2;v=0':      None,
'HCN;v2=2':     ['HCN;v2=1','HCCCN;v7=2'],
'HCN;v2=1':     None,
'HCCCN;v=0':    ['S-34-O2;v=0'],
'HCCCN;v6=1':   ['SO-18;v=0','HCCCN;v7=1'],
'HCCCN;v6=2':   ['HCCCN;v7=1','HCN;v2=2'],
'HCCCN;v7=1':   ['SO-18;v=0','HCCCN;v6=1','HCCCN;v7=2','HCN;v2=2'],
'HCCCN;v7=2':   ['S-34-O;v=0','HCCCN;v7=1','HCN;v2=1','HCN;v2=2','S-34-O2;v=0'],
'S-34-O;v=0':   ['HCCCN;v7=2'],
'SO-18;v=0':    ['HCCCN;v6=1','HCCCN;v7=1'],
'S-34-O2;v=0':  ['SO2;v=0','HCCCN;v7=2']
}
fnpickle(fit_with, os.path.join(XCLASSdir, 'fit_with.pickle'))

detected_species = {
'1':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0',               'HCCCN;v=0', 'HCCCN;v6=1',                             'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1',               'SO-18;v=0', 'SO2;v=0'                              ],
'2':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'3':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'4':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2',                'CS;v=0', 'SO;v=0;#1',               'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'5':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1',               'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'6':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',                         'HC-13-N;v=0',                                                                                     'CS;v=0'                                                                                 ],
'7':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',                         'HC-13-N;v=0',                                                                                     'CS;v=0'                                                                                 ],
'8':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1',               'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'9':  ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',                         'HC-13-N;v=0',                                                                                     'CS;v=0', 'SO;v=0;#1'                                                                    ],
'10': ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',                         'HC-13-N;v=0',               'HCCCN;v=0',                                                          'CS;v=0', 'SO;v=0;#1',                            'SO2;v=0'                              ],
'11': ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',             'HCN;v2=1', 'HC-13-N;v=0',               'HCCCN;v=0',               'HCCCN;v7=1',                              'CS;v=0', 'SO;v=0;#1',               'SO-18;v=0', 'SO2;v=0'                              ],
'12': ['CO;v=0', 'HCO+;v=0', 'HCN;v=0',                         'HC-13-N;v=0',                                                                                     'CS;v=0', 'SO;v=0;#1',                            'SO2;v=0'                              ],
'13': ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'SO-18;v=0', 'SO2;v=0',                'H2CS;v=0;#1'],
'14': ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HCCCN;v7=2', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#1', 'S-33-O;v=0', 'SO-18;v=0', 'SO2;v=0', 'S-34-O2;v=0', 'H2CS;v=0;#1']
}
fnpickle(detected_species, os.path.join(XCLASSdir, 'detected_species.pickle'))


XCL_files = []
for SSC in tqdm(SSCs):
    for specie in detected_species[SSC['num']]:

        # export spectrum
        dat_files = export_spectrum(0, SSC, specie,
                                    rms = 0.46*u.K
                                  )

        # build observation file
            obs_file = create_observations_file(0, SSC, specie,
                                                fit_range = fit_ranges[str(SSC['no'])][specie]
                                               )

        # build molfit file
        mol_file = create_molfile_file(SSC, specie,
                                       components = components[str(SSC['no'])][specie],
                                       fit_with   = fit_with[specie]
                                      )

        # build XCLASS script
        XCL_file = create_XCLASS_script(0, SSC, specie,
                                        alg_file = os.path.join(basescriptdir,'NGC253','paper_20a','algorithm.xml'),
                                        obs_file = obs_file,
                                        mol_file = mol_file
                                       )
        XCL_files.append(XCL_file)

# create_CASA_script(XCL_files, savepath=os.path.join(XCLASSdir,'all_XCLASS.py'))
# run XCLASS script in CASA

def run_parallel_XCLASS(SSCs):

    def CASA_command(SSC):
        these_XCL_files = [f for f in XCL_files if 'SSC_'+SSC['num']+'/' in f]
        Xcommand = f"""
import time,datetime

start = time.time()
for f in {these_XCL_files}:
\texecfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
"""
        return Xcommand

    commands = [CASA_command(SSC) for SSC in SSCs]
    pool = Pool(len(SSCs))
    pool.map(run_in_casa, commands)
    pool.close()
    pool.join()


run_parallel_XCLASS(SSCs)


# check fit
for SSC in tqdm(SSCs):
    for specie in detected_species[SSC['num']]:
        bands = bands_from_specie(specie)
        for band in bands:
            try:
                plot_spectrum_fit(SSC, specie, band)
            except:
                print("error in SSC "+str(SSC['no'])+" "+specie)


###################################################################################################
# manually re-fit overlapping lines and failed fits
###################################################################################################

SSC    = get_SSC('14')
specie = 'HCN;v2=2'
band   = 'USB'
XCL_file = create_XCLASS_script(0, SSC, specie,
                                alg_file = os.path.join(basescriptdir,'NGC253','paper_20a','algorithm.xml'),
                                obs_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_0','observation.xml')),
                                mol_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules_manual.molfit')),
                                mode     = 'manual'
                               )
run_in_casa(f"execfile('{XCL_file}')")
plot_spectrum_fit(SSC, specie, band)


###################################################################################################
# run 100 times to bootstrap errors using final setup files
###################################################################################################

nums = np.arange(0,101)

def final_100_run(SSC):
    XCL_files = []
    for specie in detected_species[str(SSC['no'])]:
        print('SSC '+SSC['num']+': '+specie)
        for num in nums:

            get_final_files(num, SSC, specie)

            # export spectrum
            dat_files = export_spectrum(num, SSC, specie,
                                        rms = 0.46*u.K
                                      )

            # build XCLASS script
            XCL_file = create_XCLASS_script(num, SSC, specie,
                                            alg_file = os.path.join(basescriptdir,'NGC253','paper_20a','algorithm.xml'),
                                            obs_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'run_'+str(num),'observation_final.xml')),
                                            mol_file = escape_fname(os.path.join(XCLASSdir,'SSC_'+str(SSC['no']),specie,'molecules_final.molfit')),
                                            mode     = 'final'
                                           )
            XCL_files.append(XCL_file)
    fnpickle(XCL_files, 'temp.XCL_files.SSC_'+SSC['num']+'.pickle')


pool = Pool(len(SSCs))
pool.map(final_100_run, SSCs)
pool.close()
pool.join()
XCL_files = [fnunpickle('temp.XCL_files.SSC_'+SSC['num']+'.pickle') for SSC in SSCs]
XCL_files = [a for b in XCL_files for a in b]


def run_parallel_XCLASS(SSCs):

    def CASA_command(SSC):
        these_XCL_files = [f for f in XCL_files if 'SSC_'+SSC['num']+'/' in f]
        Xcommand = f"""
import time,datetime

start = time.time()
for f in {these_XCL_files}:
\texecfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
"""
        return Xcommand

    commands = [CASA_command(SSC) for SSC in SSCs]
    pool = Pool(len(SSCs))
    pool.map(run_in_casa, commands)
    pool.close()
    pool.join()


run_parallel_XCLASS(SSCs)


###################################################################################################
# examine fit
###################################################################################################

# plot fit
for SSC in tqdm(SSCs):
    for specie in detected_species[str(SSC['no'])]:
        bands = bands_from_specie(specie)
        for band in bands:
            try:
                plot_spectrum_fit(SSC, specie, band)
            except:
                pass

# check convergence
for SSC in tqdm(SSCs):
    for specie in detected_species[str(SSC['no'])]:
        plot_chi2_evolution(SSC, specie)

# plot histogram of fit parameters
for SSC in tqdm(SSCs):
    for specie in detected_species[str(SSC['no'])]:
        try:
            plot_fit_hist(SSC, specie)
        except:
            print("Could not plot histogram: SSC "+str(SSC['no'])+" "+specie)


###################################################################################################
# plot bestfit model spectrum
###################################################################################################

# XCL_files = []
# for SSC in tqdm(SSCs):
#     model_molfit = create_model_spectrum(0, SSC, detected_species[str(SSC['no'])], mode='final')
#     XCL_files.append(model_spectrum_XCLASS(0, SSC, model_molfit))
#
# create_CASA_script(XCL_files, savepath=os.path.join(XCLASSdir,'model_spectrum_XCLASS.py'))
#
# # run XCLASS script in CASA
#
# for SSC in tqdm(SSCs):
#     for band in ['LSB','USB']:
#         plot_combined_spectrum(SSC, band)


###################################################################################################
# plot model/data variation
###################################################################################################

XCL_files = []
for SSC in tqdm(SSCs):
    for num in nums:
        model_molfit = create_model_spectrum(num, SSC, detected_species[str(SSC['no'])], mode='final')
        XCL_files.append(model_spectrum_XCLASS(num, SSC, model_molfit))

def run_parallel_XCLASS(SSCs):

    def CASA_command(SSC):
        these_XCL_files = [f for f in XCL_files if 'SSC_'+SSC['num']+'/' in f]
        Xcommand = f"""
import time,datetime

start = time.time()
for f in {these_XCL_files}:
\texecfile(f)

stop = time.time()
exec_time = np.round(stop-start, 1)
casalog.post("\\nFinished XCLASS runs \\nExecution took "+str(datetime.timedelta(seconds=exec_time))+"hours.\\n")
"""
        return Xcommand

    commands = [CASA_command(SSC) for SSC in SSCs]
    pool = Pool(len(SSCs))
    pool.map(run_in_casa, commands)
    pool.close()
    pool.join()


run_parallel_XCLASS(SSCs)

for SSC in tqdm(SSCs):
    for band in ['LSB','USB']:
        plot_combined_variation(nums, SSC, band, rms=0.46*u.K)


###################################################################################################
# load fitted data
###################################################################################################

line_data_N = get_physical(SSCs, fitable_species)

# remove undetected lines from line_data
# not really necessary anymore after running the fits only for detected lines
ld2 = copy.deepcopy(line_data_N)
for SSC in SSCs:
    for k in line_data_N[str(SSC['no'])].keys():
        if not k in detected_species[str(SSC['no'])]:
            del ld2[str(SSC['no'])][k]
line_data_N = copy.deepcopy(ld2)

# cross-check that all detected lines where fitted
for SSC in SSCs:
    for k in detected_species[str(SSC['no'])]:
        if not k in line_data_N[str(SSC['no'])].keys():
            print("Line detected but not fitted: SSC "+str(SSC['no'])+" "+k)

os.system('mkdir -p '+XCLASSdir)
fnpickle(line_data_N, os.path.join(XCLASSdir, 'line_column_density_data.pickle'))


###################################################################################################
# plot fit parameters
###################################################################################################

for SSC in tqdm(SSCs):
    plot_parameters_N(SSC)


###################################################################################################
# column density ratios
###################################################################################################

ratios_N = {
'CO/HCN':       {'a':'CO;v=0',     'b':'HCN;v=0'},
'CO/HCO+':      {'a':'CO;v=0',     'b':'HCO+;v=0'},
'CO/CS':        {'a':'CO;v=0',     'b':'CS;v=0'},
'HCN/HCO+':     {'a':'HCN;v=0',    'b':'HCO+;v=0'},
'HCN/HNC':      {'a':'HCN-15;v=0', 'b':'HN-15-C;v=0'},
'CS/HCN':       {'a':'CS;v=0',     'b':'HCN;v=0'},
'CS/HCO+':      {'a':'CS;v=0',     'b':'HCO+;v=0'},
'SO/S18O':      {'a':'SO;v=0;#1',  'b':'S-18-O;v=0'},
'12C/13C':      {'a':'HCN;v=0',    'b':'HC-13-N;v=0'},
'14N/15N':      {'a':'HCN;v=0',    'b':'HCN-15;v=0'},
'32S/33S':      {'a':'SO;v=0;#1',  'b':'S-33-O;v=0'},
'32S/34S':      {'a':'SO2;v=0',    'b':'S-34-O2;v=0'},
'C/O':          {'a':'CS;v=0',     'b':'SO;v=0;#1'},
'SO/SO2':       {'a':'SO;v=0;#1',  'b':'SO2;v=0'},
'HCN/HC3N':     {'a':'HCN;v=0',    'b':'HCCCN;v=0'},
'HCNvib1/HCN':  {'a':'HCN;v2=1',   'b':'HCN;v=0'},
'HCNvib2/HCN':  {'a':'HCN;v2=2',   'b':'HCN;v=0'},
'HC3Nvib1/HC3N': {'a':'HCCCN;v7=1', 'b':'HCCCN;v=0'},
'HC3Nvib2/HC3N': {'a':'HCCCN;v7=2', 'b':'HCCCN;v=0'}
}

for rname,rdict in tqdm(ratios_N.items()):
    ratios_N[rname]['median'] = []
    ratios_N[rname]['16th']   = []
    ratios_N[rname]['84th']   = []
    ratios_N[rname]['all']    = []

    for SSC in SSCs:
        try:
            nearest_comp_a = np.argmin(np.abs(line_data_N[str(SSC['no'])][rdict['a']]['velocity']['median']))
            a_all = np.array(line_data_N[str(SSC['no'])][rdict['a']]['column density']['all'][nearest_comp_a])
            nearest_comp_b = np.argmin(np.abs(line_data_N[str(SSC['no'])][rdict['b']]['velocity']['median']))
            b_all = np.array(line_data_N[str(SSC['no'])][rdict['b']]['column density']['all'][nearest_comp_b])

            # removing bad fits results in varying length arrays
            # fill in the missing values with the median
            if len(a_all)<len(b_all):
                for i in np.arange(np.abs(len(a_all)-len(b_all))):
                    a_all = np.append(a_all, np.median(a_all))
            if len(a_all)>len(b_all):
                for i in np.arange(np.abs(len(a_all)-len(b_all))):
                    b_all = np.append(b_all, np.median(b_all))

            r_all  = a_all/b_all
            r16,rmed,r84 = np.percentile(r_all, (16,50,84))
        except KeyError:
            r_all = [np.nan]
            r16   = np.nan
            rmed  = np.nan
            r84   = np.nan

        ratios_N[rname]['median'].append(rmed)
        ratios_N[rname]['16th'].append(r16)
        ratios_N[rname]['84th'].append(r84)
        ratios_N[rname]['all'].append(r_all)

os.system('mkdir -p '+ratiodir)
fnpickle(ratios_N, os.path.join(ratiodir, 'ratios_column_density.pickle'))


###################################################################################################
# plot ratios
###################################################################################################

for rname,rdict in tqdm(ratios_N.items()):
    plot_ratios_N(rname, rdict)

ratios_N_table(['12C/13C','14N/15N','32S/33S','32S/34S'])


###################################################################################################
# plot fitted spectrum oerview
###################################################################################################

execfile(os.path.join(scriptdir, '725.spectrum_overview.py'))
spectrum_overview()


###################################################################################################
# compare quantities across SSCs
###################################################################################################

execfile(os.path.join(scriptdir, '720.XCLASS.plot_SSC_quantity_N.py'))

# temperature
for T_specie in tqdm(reliable_Ts):
    plot_SSC_quantity_N(T_specie, 'temperature')

# column density
for N_specie in tqdm(['CO;v=0','HCN;v=0','CS;v=0']):
    plot_SSC_quantity_N(N_specie, 'column density')

# linewidth
plot_SSC_quantity_N('CS;v=0', 'linewidth')
compare_linewidths_N(['CS;v=0','HC-13-N;v=0','CO;v=0','HCN;v=0','HCO+;v=0'])

# velocity
plot_SSC_quantity_N('CS;v=0', 'velocity')

# temperature table
temperature_table(['SO2;v=0','S-34-O2;v=0','H2CS;v=0;#1'])

###################################################################################################
#
###################################################################################################
