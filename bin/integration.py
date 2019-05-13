import time
import numpy as np
import h5py

class A:
    pass

A.tmax  = 10.0; 
A.nt    = 20;
A.beta  = 4.0
A.ntau  = 40
A.norb  = 3
A.order = 6

A.sig   = -1
A.dtau  = A.beta/(A.ntau-1)

f = h5py.File("results/GM.h5")
arr = f['/M'][...]
A.M = arr[:]['real'] + 1j*arr[:]['imag']
f.close()

f = h5py.File("results/integ.h5")
arr = f['/rcorr'][...]
rcorr = arr[:]
arr = f['/gmM'][...]
gregory_matrix_M = arr[:]
arr = f['/gmR'][...]
gregory_matrix_R = arr[:]
f.close()

rcorr = np.reshape(rcorr, [A.order, A.order, A.order], order='F')

A.M = np.reshape(A.M, [A.ntau, A.norb, A.norb], order='F')

gregory_matrix_M = np.reshape(gregory_matrix_M, [A.ntau, A.ntau], order='F')

gregory_matrix_R = np.reshape(gregory_matrix_R, [A.nt, A.nt], order='F')

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
