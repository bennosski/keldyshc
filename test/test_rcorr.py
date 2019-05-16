import time
import numpy as np
import h5py
from read_params import read_params
from util import *

dr = np.load("../bin/rcorr.npy").item();
rcorr = dr[6]
order = np.shape(rcorr)[0]
rcorr = np.reshape(rcorr, [-1])
rcorr = np.reshape(rcorr, [order, order, order], order='F')

crcorr = zload("../test/results/test_rcorr.h5", "/rcorr")
crcorr = np.reshape(crcorr, [order, order, order], order='F')

print('diff = {}'.format(np.mean(np.abs(crcorr - rcorr))))



