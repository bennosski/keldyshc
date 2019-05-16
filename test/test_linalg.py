import numpy as np
import h5py

f = h5py.File("results/test_linalg.h5")
arr1 = f['/arr1'][...]
arr1 = arr1['real'] + 1j*arr1['imag']
arr2 = f['/arr2'][...]
arr2 = arr2['real'] + 1j*arr2['imag']
ret  = f['/ret'][...]
ret  = ret['real'] + 1j*ret['imag']

n = int(np.sqrt(len(ret)))

arr1 = np.reshape(arr1, [n, n], order='F')
arr2 = np.reshape(arr2, [n, n], order='F') 
ret  = np.reshape(ret, [n, n], order='F')

print('diff : {:1.5e}'.format(np.mean(np.abs(ret-arr1@arr2))))
