import time
import numpy as np
import h5py
from read_params import read_params
from util import *

"""
import asyncio

async def run(cmd):
    proc = await asyncio.create_subprocess_shell( \
        cmd, \
        stdout=asyncio.subprocess.PIPE, \
        stderr=asyncio.subprocess.PIPE)

    stdout, stderr = await proc.communicate()

    print('Running keldysh C code')
    print(f'[stdout]\n{stdout.decode()}')

    '''
    print(f'[{cmd!r} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')
    '''

loop = asyncio.get_event_loop()
loop.run_until_complete(asyncio.wait([run('./keldysh_test')]))
loop.close()
"""

class A:
    pass

A.tmax, A.nt, A.beta, A.ntau, A.norb, A.order = read_params()

A.sig   = -1
A.dtau  = A.beta/(A.ntau-1)

A.M = zload('results/GM.h5', '/M')

rcorr = zload('results/integ.h5', '/rcorr')
gregory_matrix_M = zload('results/integ.h5', '/gmM')
gregory_matrix_R = zload('results/integ.h5', '/gmR')

rcorr = np.reshape(rcorr, [A.order, A.order, A.order], order='F')

A.M = np.reshape(A.M, [A.ntau, A.norb, A.norb], order='F')

gregory_matrix_M = np.reshape(gregory_matrix_M, [A.ntau, A.ntau], order='F')

gregory_matrix_R = np.reshape(gregory_matrix_R, [A.nt, A.nt], order='F')

'''
print('rcorr')
print(rcorr)

print('A.M')
print(A.M)

print('gmM')
print(gregory_matrix_M)

print('gmR')
print(gregory_matrix_R)
'''

#A.M = np.random.randn(A.ntau, A.norb, A.norb)
#rcorr = np.random.randn(A.order, A.order, A.order)
#gregory_matrix_M = np.random.randn(A.ntau, A.ntau)

def prep_MxM(A):
    '''
    prepare A for higher-order (Gregory and boundary correction) convolution using matrix multiplication
    
    A : a matsubara or langreth object with a Matsubara matrix of size (ntau x norb, norb)

    output : a matrix of size (ntau x norb x ntau x norb) which can be multiplied by a second vector to produce the convolution
    '''
        
    ntau  = A.ntau
    dtau  = A.dtau
    norb  = A.norb
    order = A.order
    
    Cmk = np.zeros((ntau,norb,ntau,norb),dtype=np.complex128)        
        
    for iorb in range(0,norb):
        for korb in range(0,norb):

            for m in range(ntau-order,ntau-1):
                for k in range(ntau-order,ntau):
                    for l in range(0,order):                            
                        Cmk[m,iorb,k,korb] += -1j*dtau * rcorr[ntau-1-m,l,ntau-1-k] * A.sig * A.M[ntau-1-l,iorb,korb]

            for m in range(0,ntau-order):
                Cmk[m,iorb,m:ntau,korb] += -1j*dtau * gregory_matrix_M[ntau-1-m,:ntau-m] * A.sig * A.M[np.arange(ntau-1,m-1,-1),iorb,korb]

            for m in range(1,order):
                for k in range(0,order):
                    for l in range(0,order):
                        Cmk[m,iorb,k,korb] += -1j*dtau * rcorr[m,l,k] * A.M[l,iorb,korb]

            for m in range(order,ntau):
                Cmk[m,iorb,:m+1,korb] += -1j*dtau * gregory_matrix_M[m,:m+1] * A.M[np.arange(m,-1,-1),iorb,korb]                        
    return np.reshape(Cmk, [ntau*norb, ntau*norb], order='F')

t = time.time()
Cmk = prep_MxM(A)
print('took {:.4e}'.format(time.time()-t))

print('Cmk')
print(Cmk)

f = h5py.File("results/Cmk.h5")
arr = f['/Cmk'][...]
Cmkc = arr[:]['real'] + 1j*arr[:]['imag']
f.close()

Cmkc = np.reshape(Cmkc, [A.ntau*A.norb, A.ntau*A.norb], order='F')
print('Cmkc')
print(Cmkc)

print('difference')
print(np.mean(np.abs(Cmk-Cmkc)))
