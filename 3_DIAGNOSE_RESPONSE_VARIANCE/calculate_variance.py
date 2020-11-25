import sys
import numpy as np
import netCDF4 as nc
import scipy.sparse as sp
###############################################################################
'''
DESCRIPTION: This script takes "nested" adjoint outputs along with stochastic
covariance matrices and calculates the variance accumulated due to white noise
forcing at the surface.

"Nesting" is used as metrics are likely most sensitive to change at short lags.
As the variance is calculated using numerical integration (see Eq. 3 of our
manuscript), this means large time-steps will exaggerate this early stage
sensitivity. The model is therefore run three times: once for 5 days with
output every timestep, once for 1 year with output every 5 days, and once for
60 years with output every year. This allows the integral over the first period
of each run to be replaced by a much finer numerical integral from the shorter
run with higher frequency output.

PARAMETERS:
----------
filename_05d: str
    The 5 day adjoint output file
filename_01y: str
    The 1  y. adjoint output file
filename_60y: str
    The 60 y. adjoint output file
cov_mat_dir:  str
    The location of the files COV_SSX_SSY.npz containing the white noise 
    covariance matrix for variables X,Y (e.g. COV_SST_SSS.npz)
filename_out: str
    The desired output netCDF file name.

OUTPUT:
-------
<filename_out>.nc netCDF file containing buoyancy-forced and momentum-forced
variance accumulated during each of the 60 years of the adjoint run.

'''
### PARAMETERS

var_calcs=True #Calculate variance or not
save_vars=True #Save calculated variances to netCDF file or not

filename_05d='/your/05d/adjoint/output.nc'
filename_01y='/your/01y/adjoint/output.nc'
filename_60y='/your/60y/adjoint/output.nc'

filename_out='/your/output/filename.nc'

cov_mat_dir='/location/of/your/covariance/matrix/npz/files'

################################################################################

def variance_sources(var1_ad,var2_ad,cov_12,dt):
    '''
    Calculate the covariance between two adjoint sensitivity fields (var1_ad,
    var2_ad) in response to a stochastic forcing covariance matrix (cov_12).

    Parameters:
    -----------
    var1_ad,var2_ad: ndarray [N,149,182]
           Adjoint sensitivity fields (M*|F> for the cost function |F>)
    cov_12: sparse matrix [27118,27118]
           Covariance matrix of white noise forcing (units of flux^2*time)
    
    Returns:
    -------
    var_map: 3d array [N,149,182]
           The contribution to the total variance of the adjoint metric at each
           output snapshot.
    '''
                    
    NT=np.shape(       var1_ad)[0]             #number of output snapshots
    var1_ad=np.reshape(var1_ad,(NT,149*182)) # first adjoint variable
    var2_ad=np.reshape(var1_ad,(NT,149*182)) # second adjoint variable
    # Calculate <F|M S M*|F>dt :
    var_map=np.cumsum((\
                     ((var2_ad.T)*(cov_12.dot(var1_ad.T))).T\
                                                         ).reshape(NT,149,182)\
                                                                     *dt,axis=0)
    return var_map
################################################################################

if var_calcs:
    #################
    # Calculate the integral for each of the 3 adjoint outputs and add result
    # to a dictionary
    ##################

    # Open the adjoint outputs:
    NC05d=nc.Dataset(filename_05d)
    NC01y=nc.Dataset(filename_01y)
    NC60y=nc.Dataset(filename_60y)
    
    #Load the ice cover mask:
    ice05d=np.flip(NC05d.variables['ice_fraction'][:],axis=0)
    ice01y=np.flip(NC01y.variables['ice_fraction'][:],axis=0)
    ice60y=np.flip(NC60y.variables['ice_fraction'][:],axis=0)
    V05d={};V01y={};V60y={}

    #Loop through all variable combinations, skip ones which throw a loading
    # error (e.g. covariance between T and u) as these do not exist
    print('Calculating variance in response to covariance between')
    for A in ['t','s','u','v']:
        for B in ['t','s','u','v']:
            print(A+' and '+B+':')
            try:
                # Load covariance matrix:
                C=sp.load_npz(cov_mat_dir+\
                              '/COV_SS'+A.upper()+'_SS'+B.upper()+'.npz')
            except:
                print('Either '+cov_mat_dir+\
                              '/COV_SS'+A.upper()+'_SS'+B.upper()+'.npz'+\
                      ' does not exist (no problem) or an error occured')
                continue
                
            # Load adjoint variable 1:
            X05d=np.flip(NC05d.variables[A+'_ad'][:,0,:],axis=0)*(1-ice05d)
            X01y=np.flip(NC01y.variables[A+'_ad'][:,0,:],axis=0)*(1-ice01y)
            X60y=np.flip(NC60y.variables[A+'_ad'][:,0,:],axis=0)*(1-ice60y)
            # Load adjoint variable 2:
            Y05d=np.flip(NC05d.variables[B+'_ad'][:,0,:],axis=0)*(1-ice05d)
            Y01y=np.flip(NC01y.variables[B+'_ad'][:,0,:],axis=0)*(1-ice01y)
            Y60y=np.flip(NC60y.variables[B+'_ad'][:,0,:],axis=0)*(1-ice60y)
            # Calculate integral
            V05d[A+'_'+B]=variance_sources(X05d,Y05d,C,5760)
            V01y[A+'_'+B]=variance_sources(X01y,Y01y,C,5760*75)
            V60y[A+'_'+B]=variance_sources(X60y,Y60y,C,5760*5475)
            print('success')


# Combine the three different runs into a single nested integral:
VAR={}
for i in V60y.keys():
    VAR[i]=np.zeros(np.shape(V60y[i]))    
    VAR[i][1 ,:]=(V01y[i][-1,:]- V01y[i][1,:])+(V05d[i][-1,:]-V05d[i][1,:])
    VAR[i][1:,:]=VAR[i][1,:]+(V60y[i][1:,:]- V60y[i][1,:])

# Calculate momentum and buoyancy components:
BUO=VAR['t_t']+VAR['t_s']*2+VAR['s_s']
MOM=VAR['u_u']+VAR['u_v']*2+VAR['v_v']
################################################################################
if save_vars:
    #Write output to netCDF file
    print('Saving variance calculations to netCDF')
    NC=nc.Dataset(filename_out,'w')
    NC.createDimension('t',None)
    NC.createDimension('y',149)
    NC.createDimension('x',182)
    BO=NC.createVariable('buoyancy_var',np.float64(),('t','y','x'))
    MO=NC.createVariable('momentum_var',np.float64(),('t','y','x'))
    ti=NC.createVariable('time',np.float64(),('t'))

    BO[:]=BUO*1e-42
    MO[:]=MOM*1e-42
    ti[:]=np.arange(61)

    ti.setncattr('Units','yr')
    BO.setncattr('Units','(ZJ)^2')
    MO.setncattr('Units','(ZJ)^2')

    BO.setncattr('Description','Cumulative heat content variance due to'+\
                 ' stochastic surface buoyancy forcing')
    BO.setncattr('Description','Cumulative heat content variance due to'+\
                 'stochastic surface momentum forcing')
    NC.close()

