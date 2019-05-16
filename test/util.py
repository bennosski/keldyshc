
import h5py


def zload(filename, dset):
    f = h5py.File(filename)
    arr = f[dset][...]
    arr = arr[:]['real'] + 1j*arr[:]['imag']
    f.close()
    return arr

def dload(filename, dset):
    f = h5py.File(filename)
    arr = f[dset][...]
    f.close()
    return arr
