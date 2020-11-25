import sys
import numpy as np
import netCDF4 as nc

################################################################################
'''
DESCRIPTION:
-----------
Build an ocean heat content cost function for the NEMOTAM "TAM_PRED" configuration
for a given basin over a given depth range.

PARAMETERS:
----------
BASIN: str 
   takes one of the following values:
   'arc'   - Arctic Ocean           [ 70, 90]N
   'natl'  - North Atlantic         [ 35, 70]N
   'itatl' - Intertropical Atlantic [-35, 35]N
   'npac'  - North Pacific          [ 35, 70]N
   'itpac' - Intertropical Pacific  [-35, 35]N
   'ind'   - Indian Ocean          ~[ 27,-70]N
   'sou'   - Southern Ocean         [-90,-70]N
DEPTH_LEVEL: int 
    the depth range over which to evaluate the heat content, determined by the
    index (0-30)
filename_out: str
    The output filename, supplied to NEMOTAM using the cn_tam_input in the 
    adjoint namelist.

'''
################################################################################
### PARAMETERS:
DEPTH_LEVEL = 30                   #Index of cost fn. max depth level [0-30]
BASIN  ='natl'                     #Basin key (see above)
filename_out='your/output/file.nc' #Output filename

# Constants
cp=  3850 # Specific heat capacity (JK^-1kg^-1)
rho0=1026 # Reference density

################################################################################

# Load grid geometry
O2G  =nc.Dataset('mesh_mask.nc')
latt =O2G.variables['gphit'][:]
lont =O2G.variables['glamt'][:]
e1t  =O2G.variables['e1t'  ][:]
e2t  =O2G.variables['e2t'  ][:]
e3t  =O2G.variables['e3t'  ][:]
tmask=O2G.variables['tmask'][:]
gdept=O2G.variables['gdept'][:]

# Load basin masks
O2M  =nc.Dataset('ORCA2_subbasins.nc')
AT   =O2M.variables['atlmsk_nomed'][:]
PA   =O2M.variables['pacmsk'][:]
IN   =O2M.variables['indmsk'][:]

depmask=np.ones(np.shape(gdept))
depmask[:,DEPTH_LEVEL+1:,:,:]=0

if   (BASIN == 'arc'  ): #Arctic Ocean
    latmask=np.zeros(np.shape(latt));latmask[latt>70]=1
    MASK=latmask*tmask

elif (BASIN == 'natl' ): #North Atlantic
    latmask=np.zeros(np.shape(latt));latmask[(latt>35) & (latt<70)]=1
    MASK=latmask*AT*tmask

elif (BASIN == 'itatl'): #Intertropical Atlantic
    latmask=np.zeros(np.shape(latt));latmask[(latt>-35) & (latt<35)]=1
    MASK=latmask*AT*tmask

elif (BASIN == 'npac' ): #North Pacific
    latmask=np.zeros(np.shape(latt));latmask[(latt>35) & (latt<70)]=1
    MASK=latmask*PA*tmask

elif (BASIN == 'itpac'): #Intertropical Pacific
    latmask=np.zeros(np.shape(latt));latmask[(latt>-35) & (latt<35)]=1
    MASK=latmask*PA*tmask

elif (BASIN == 'ind'  ): #Indian ocean    
    latmask=np.ones(np.shape(latt));latmask[(latt<-35)]=0
    MASK=IN*tmask*latmask

elif (BASIN == 'sou'  ): #Southern ocean
    latmask=np.zeros(np.shape(latt));latmask[(latt<-35)]=1
    MASK=latmask*tmask
else:
    raise Exception("No valid ocean basin key provided")

#Tinit=(e1t*e2t*e3t*MASK*depmask)/np.sum(e1t*e2t*e3t*MASK*depmask) #Heat content in K
Tinit=(e1t*e2t*e3t*MASK*depmask*cp*rho0) #Heat content in J

# SAVE TO NETCDF
with nc.Dataset(filename_out,'w') as NC:
    NC.createDimension('z',31)
    NC.createDimension('y',149)
    NC.createDimension('x',182)
    Tinit_nc=NC.createVariable('t0_ad',np.float64(),('z','y','x'))
    vinit_nc=NC.createVariable('v0_ad',np.float64(),('z','y','x'))
    Tinit_nc[:]=Tinit[:]
    vinit_nc[:]=np.zeros(np.shape(Tinit))
