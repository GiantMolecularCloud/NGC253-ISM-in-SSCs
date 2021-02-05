##################################
# SSC CHEMISTRY SUB-PROJECT INFO #
##################################

####################################################################################################
# directories
####################################################################################################

subprojectdir = 'XXXXX'
datadir   = os.path.join(subprojectdir, '01.data/')
mandir    = os.path.join(subprojectdir, '02.manual_fit/')
XCLASSdir = os.path.join(subprojectdir, '03.XCLASS_fit/')
Xfinaldir = os.path.join(subprojectdir, '04.XCLASS_final/')
tempdir   = os.path.join(subprojectdir, '05.temperatures/')
vibdir    = os.path.join(subprojectdir, '06.vibrations/')
resultsdir = os.path.join(subprojectdir, '10.results/')
refitdir  = os.path.join(subprojectdir, '11.refit/')
intensitydir = os.path.join(subprojectdir, '12.intensities/')
paperdir  = os.path.join(subprojectdir, '13.paper/')

scriptdir = os.path.join(basescriptdir,'NGC253/paper_20a/')
plotdir   = os.path.join(baseplotdir,'NGC253/paper_20a/')


####################################################################################################
# additional imports
####################################################################################################

from astropy.table import QTable
from astropy.table.column import Column
import astropy.constants as c
execfile(os.path.join(scriptdir,'701.parse_line_list.py'))
execfile(os.path.join(scriptdir,'703.line_tex.py'))
execfile(os.path.join(scriptdir,'704.tex_line_table.py'))
execfile(os.path.join(scriptdir,'705.specie_tex.py'))


####################################################################################################
# data
####################################################################################################

vsys     = 250*u.km/u.s
distance = 3.5*u.Mpc

lines = parse_line_list(os.path.join(scriptdir,'702.lines.txt'))
SSCs  = QTable.read(os.path.join(subprojectdir,'SSCs_Leroy+18.fits'))
LSB   = SpectralCube.read(os.path.join(datadir, 'NGC253.band7.TP+12m-mid+12m-high.LSB.K.image.pbcor.fits'))
USB   = SpectralCube.read(os.path.join(datadir, 'NGC253.band7.TP+12m-mid+12m-high.USB.K.image.pbcor.fits'))

unique_species = []
for line in lines:
    if not line['XCLASS'] in unique_species:
        unique_species.append(line['XCLASS'])

fitable_species = ['CO;v=0', 'HCO+;v=0', 'HCN;v=0', 'HCN;v2=2', 'HCN;v2=1', 'HC-13-N;v=0', 'HCN-15;v=0', 'HCCCN;v=0', 'HCCCN;v7=2', 'HCCCN;v6=1', 'HCCCN;v7=1', 'HN-15-C;v=0', 'CS;v=0', 'SO;v=0;#2', 'S-33-O;v=0', 'S-34-O;v=0', 'SO-18;v=0', 'SO2;v=0', 'S-34-O2;v=0', 'H2CS;v=0;#1']


####################################################################################################
