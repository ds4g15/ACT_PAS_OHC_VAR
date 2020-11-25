import numpy as np
import netCDF4 as nc
import scipy.sparse as sp
import matplotlib.pyplot as plt

'''
DESCRIPTION:
------------
Combine flux covariance matrices with flux decorrelation time to make
stochastic white noise matrix representation. 

REQUIREMENTS:
-------------
Sparse flux covariance matrices (produced by 1_calculate_flux_covariance.py)
flux_cov_SST_SST.npz 
flux_cov_SST_SSS.npz
flux_cov_SSS_SSS.npz
flux_cov_SSU_SSU.npz
flux_cov_SSU_SSV.npz
flux_cov_SSV_SSV.npz

e-folding decorrelation times (produced by 2_calculate_flux_e-folding_time.py)
flux_e_folding_time_SST.npy
flux_e_folding_time_SSS.npy
flux_e_folding_time_SSU.npy
flux_e_folding_time_SSV.npy

OUTPUT:
-------
Sparse white noise forcing covariance matrices of shape (27118,27118):
COV_SST_SST.npz [Units K^2/s]
COV_SST_SSS.npz [Units (K psu)/s]
COV_SSS_SSS.npz [Units psu^2/s]
COV_SSU_SSU.npz [Units m^4/s^3]
COV_SSU_SSV.npz [Units m^4/s^3]
COV_SSV_SSV.npz [Units m^4/s^3]
'''

################################################################################
#Loop through all variable combinations, skip ones which throw a loading
# error (e.g. covariance between T and u) as these do not exist
for A in ['SST','SSS','SSU','SSV']:
    for B in ['SST','SSS','SSU','SSV']:
        print('COV_'+A+'_'+B)
        try:
            print('loading flux covariance matrix')
            X=sp.load_npz('flux_cov_'+A+'_'+B+'.npz')
        except:
            print('flux_cov_'+A+'_'+B+'.npz does not exist or an error occurred')
            continue
        print('loading e-folding time arrays')
        t1=np.load('flux_e_folding_time_'+A+'.npy')
        t2=np.load('flux_e_folding_time_'+B+'.npy')
        #broadcast t1+t2 in sparse format [27118,27118]
        T=sp.csr_matrix(t1.reshape(-1,1)+t2.reshape(1,-1))
        T.eliminate_zeros()
        print('Combine with flux covariance matrix')
        Z=X.multiply(T)
        print('saving white noise matrix to npz file')
        sp.save_npz(('COV_'+A+'_'+B+'.npz'),Z)
        print('done')
