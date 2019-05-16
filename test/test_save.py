from util import *

arr = dload('results/test_dsave.h5', '/dset')

print(arr)

arr1 = dload('results/test_zsave.h5', '/dset1')
arr2 = dload('results/test_zsave.h5', '/dset2')

print(arr1)
print(arr2)
