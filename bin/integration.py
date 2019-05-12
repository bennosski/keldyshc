import time
import numpy as np

class A:
    pass

A.ntau  = 1000
A.dtau  = 0.1
A.norb  = 4
A.order = 6
A.M = np.random.randn(A.ntau, A.norb, A.norb)
A.sig = -1

rcorr = np.random.randn(A.order, A.order, A.order)
gregory_matrix_M = np.random.randn(A.ntau, A.ntau)

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
    return np.reshape(Cmk, [ntau*norb, ntau*norb])


t = time.time()
Cmk = prep_MxM(A)
print('took {:.4e}'.format(time.time()-t))

