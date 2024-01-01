import numpy as np
import multiprocessing as mp
from functools import partial

zarrObj = None
data = None


def getZarrRegion(zarrObj, regionStart, regionEnd):
    # TODO - extend to ND
    if len(regionStart) == 3:
        data = zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]]   
    elif len(regionStart) == 2:
        data = zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1]]   
    
    # add support for big-endian byte order
    if data.dtype.byteorder == '>':
        # data = data.newbyteorder().byteswap() 
        data = data.astype(data.dtype.str.replace('>', '<'));
    return data
    

def setZarrRegion(zarrObj, regionStart, regionEnd, data):
    # TODO - extend to ND
    if data.ndim == 3:
        zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]] = data
    elif data.ndim == 2:
        if len(zarrObj.shape) == 3:
            zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]] = np.expand_dims(data, axis=-1)
        elif len(zarrObj.shape) == 2:
            zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1]] = data



def setZarrData(zarrObj_in, data_in):
    # write whole data for zarr with multiprocessing

    global zarrObj
    global data

    zarrObj = zarrObj_in
    data = data_in

    blockSize = zarrObj.chunks
    sz = zarrObj.shape
    nbx = np.ceil(sz[0] / blockSize[0])
    nby = np.ceil(sz[1] / blockSize[1])
    nbz = np.ceil(sz[2] / blockSize[2])

    [X, Y, Z] = np.meshgrid(np.arange(nbx), np.arange(nby), np.arange(nbz))
    # print(X, Y, Z)
    X = X.ravel()
    Y = Y.ravel()
    Z = Z.ravel()

    len(X)

    blockStarts = [None] * len(X)
    blockEnds = [None] * len(X)
    for i in range(len(X)):
        blockStarts[i] = [X[i], Y[i], Z[i]]
        blockEnds[i] = [X[i] + 1, Y[i] + 1, Z[i] + 1]
    
    # f = partial(setZarrBlock, zarrObj, data)
    with mp.Pool(mp.cpu_count()) as p:
        # p.starmap(f, zip(blockStarts, blockEnds))
        p.starmap(setZarrBlock, zip(blockStarts, blockEnds))
        
    return 0


def setZarrBlock(blockStart, blockEnd):
    # write one or more isolated blocks
    blockSize = np.array(zarrObj.chunks)

    regionStart = np.array(blockStart) * blockSize
    regionEnd = np.minimum(np.array(blockEnd) * blockSize, np.array(data.shape))
    regionStart = regionStart.astype('int32')
    regionEnd = regionEnd.astype('int32')

    # print(regionStart, regionEnd)
    if data.ndim == 3:
        zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]] = data[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]]
    elif data.ndim == 2:
        if len(zarrObj.shape) == 3:
            zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1], regionStart[2]:regionEnd[2]] = np.expand_dims(data, axis=-1)
        elif len(zarrObj.shape) == 2:
            zarrObj[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1]] = data[regionStart[0]:regionEnd[0], regionStart[1]:regionEnd[1]]

