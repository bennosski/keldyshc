import numpy as np

def _lagrange(tpts,j,t):
    Lj = 1.0
    for m in range(0,len(tpts)):
        if not m==j:
            Lj *= (t-tpts[m])/(tpts[j]-tpts[m])
    return Lj

def compute_rcorr(p):
    from scipy import integrate

    #Rcorr = np.zeros((p,p,p), dtype=np.complex128)
    Rcorr = np.zeros((p,p,p))

    tpts = np.linspace(0.0,(p-1),p)

    for i in range(0,p):
        for j in range(0,p):
            for m in range(1,p):
                def kern(t):
                    return _lagrange(tpts,i,m-t)*_lagrange(tpts,j,t)

                Rcorr[m,i,j],err = integrate.quad(kern, 0.0, m)
    return Rcorr

def construct_rcorr():    
    rcorr = {}
    for p in range(7):
        rcorr[p] = compute_rcorr(p)
    np.save('rcorr.npy', rcorr)

def h5pysave(dr):
    import h5py

    f = h5py.File('weights.h5', 'w')
    for p in range(1, 7):
        arr = dr[p]
        arr = np.reshape(arr, [-1]) # flatten first
        rcorr = np.reshape(arr, np.shape(arr), order='F')
        f.create_dataset(f'/rcorr_p{p}', data=rcorr)

    #f = h5py.File('gregory.h5', 'w')

    #for p in range(1, 7):
    
        wstart, omega = weights(p)
        
        print('shapes ', np.shape(wstart), np.shape(omega), p)
        assert np.shape(wstart) == (2*p-1, p)
        assert np.shape(omega) == (p, )
        
        arr = np.reshape(wstart, [-1])
        arr = np.reshape(arr, np.shape(wstart), order='F')
        f.create_dataset(f'/wstart_p{p}', data=arr)
        
        f.create_dataset(f'/omega_p{p}', data=omega)

    f.close()

def h5pyload():
    import h5py

    f = h5py.File('rcorr.h5', 'r')
    for p in range(7):
        arr = f[f'p{p}'][...]
        print(arr)
    f.close()

def npload():    
    dr = np.load('rcorr.npy').item()
    for p in range(7):
        arr = dr[p]
        print(arr)
    return dr
    


def StartPolyW(p):
    w = np.zeros((p,p))
    if p==2:
        w[1,:] = [0.5,0.5]
    elif p==3:
        w[1,:] = [5.0/12.0,2.0/3.0,-1.0/12.0]
        w[2,:] = [1.0/3.0,4.0/3.0,1.0/3.0]
    elif p==4:
        w[1,:] = [3.0/8.0,19.0/24.0,-5.0/24.0,1.0/24.0]
        w[2,:] = [1.0/3.0,4.0/3.0,1.0/3.0,0.0]
        w[3,:] = [3.0/8.0,9.0/8.0,9.0/8.0,3.0/8.0]
    elif p==5:
        w[1,:] = [251.0/720.0,323.0/360.0,-11.0/30.0,53.0/360.0,-19.0/720.0]
        w[2,:] = [29.0/90.0,62.0/45.0,4.0/15.0,2.0/45.0,-1.0/90.0]
        w[3,:] = [27.0/80.0, 51.0/40.0, 9.0/10.0, 21.0/40.0, -(3.0/80.0)]
        w[4,:] = [14.0/45.0, 64.0/45.0, 8.0/15.0, 64.0/45.0, 14.0/45.0]
    elif p==6:
        w[1,:] = [95.0/288.0, 1427.0/1440.0, -(133.0/240.0), 241.0/720.0, -(173.0/1440.0), 3.0/160.0]
        w[2,:] = [14.0/45.0, 43.0/30.0, 7.0/45.0, 7.0/45.0, -(1.0/15.0), 1.0/90.0]
        w[3,:] = [51.0/160.0, 219.0/160.0, 57.0/80.0, 57.0/80.0, -(21.0/160.0), 3.0/160.0]
        w[4,:] = [14.0/45.0, 64.0/45.0, 8.0/15.0, 64.0/45.0, 14.0/45.0, 0.0]
        w[5,:] = [95.0/288.0, 125.0/96.0, 125.0/144.0, 125.0/144.0, 125.0/96.0, 95.0/288.0]
    return w
#============================================================
def AdamsMoultonW(p):
    if p==1:
        w = [1.0]
        return w
    elif p==2:
        w = [0.5,0.5]
        return w
    elif p==3:
        w = [-(1.0/12.0), 2.0/3.0, 5.0/12.0]
        return w
    elif p==4:
        w = [1.0/24.0, -(5.0/24.0), 19.0/24.0, 3.0/8.0]
        return w
    elif p==5:
        w = [-(19.0/720.0), 53.0/360.0, -(11.0/30.0), 323.0/360.0, 251.0/720.0]
        return w
    elif p==6:
        w = [3.0/160.0, -(173.0/1440.0), 241.0/720.0, -(133.0/240.0), 1427.0/1440.0, 95.0/288.0]
        return w
#============================================================
def OmegaW(p):
    if p==1:
        w = [1.0]
    elif p==2:
        w = [0.5,1.0]
    elif p==3:
        w = [5.0/12.0, 13.0/12.0, 1.0]
    elif p==4:
        w = [3.0/8.0, 7.0/6.0, 23.0/24.0, 1.0]
    elif p==5:
        w = [251.0/720.0, 299.0/240.0, 211.0/240.0, 739.0/720.0, 1.0]
    elif p==6:
        w = [95.0/288.0, 317.0/240.0, 23.0/30.0, 793.0/720.0, 157.0/160.0, 1.0]
    return w
#============================================================
def weights(p):
    wpoly = StartPolyW(p)

    wam = AdamsMoultonW(p)

    omega = OmegaW(p)

    wstart = np.zeros((2*p-1,p))

    wstart[0:p,0:p] = wpoly
    for k in range(0,p-1):
        for j in range(0,k+1):
            wstart[p+k,j] = wstart[p+k-1,j]
        for j in range(k+1,p):
            wstart[p+k,j] = wstart[p+k-1,j] + wam[j-k-1]

    return wstart,omega

'''
def construct_and_save_gregory_weights():
    import h5py

    f = h5py.File('gregory.h5', 'w')

    for p in range(1, 7):
        wstart, omega = weights(p)
        
        print('shapes ', np.shape(wstart), np.shape(omega), p)
        assert np.shape(wstart) == (2*p-1, p)
        assert np.shape(omega) == (p, )
        
        arr = np.reshape(wstart, [-1])
        arr = np.reshape(arr, np.shape(wstart), order='F')
        f.create_dataset(f'/wstart_p{p}', data=arr)
        
        f.create_dataset(f'/omega_p{p}', data=omega)
        
    f.close()

    #np.save('gregory.npy', gregory)
'''

def test():
    #arr = np.random.randn(5,5)

    arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])


    arr1 = np.reshape(arr, [-1])
    arr1 = np.reshape(arr1, np.shape(arr), order='F')
    
    arr2 = np.reshape(arr, np.shape(arr), order='F')

    print('diff ', np.mean(np.abs(arr1-arr2)))
    
    print(arr1)
    print(' ')
    print(arr2)


if __name__=='__main__':
    # for the rcorr weights: #

    try:
        from scipy import integrate

        print('scipy exists. constructing rcorr')

        construct_rcorr()

    except:
        print('no scipy. saving weights')

        dr = npload()
        h5pysave(dr)
        h5pyload()


    
    


