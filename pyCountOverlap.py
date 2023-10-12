# -*- coding: utf-8 -*-
""" __author__ =  'Ron Hoebe <r.a.hoebe@amsterdamumc.nl>'
    __version__ = '1.01'
    __website__ = 'https://github.com/orgs/Cellular-Imaging-Amsterdam-UMC/repositories'
    __license__ = 'GNU 3.0 license'        
"""
from skimage import measure
import numpy as np
import cv2
from cv2 import imread
import pandas as pd


def count_overlap(bigMask: str = 'images/CellsNucleiLabels.tif',
                  smallMask: str = 'images/CellsGranulesLabels.tif',
                  columnName: str = 'Count'
                  ) -> pd.DataFrame:
    # -*- coding: utf-8 -*-
    """Count overlap of 2 masks

        input:
            bigMask:    str     file to read
            smallMask:  str     file to read
            columnName: str     column name
        output:
            result:     DataFrame
    """

    # Count overlapping
    bigLabels = imread(bigMask, cv2.IMREAD_ANYDEPTH)
    smallLabels = imread(smallMask, cv2.IMREAD_ANYDEPTH)
    numBig = np.max(bigLabels)
    BigNumSmall = np.zeros([numBig, 1], dtype=np.int16)
    smallStats = pd.DataFrame(measure.regionprops_table(
        smallLabels, properties=('centroid',)))
    smallStatsnp = np.ndarray.astype(
        np.round(smallStats.to_numpy()), dtype=np.uint16)
    smallStatsInBigLabel = bigLabels[smallStatsnp[:, 0],
                                     smallStatsnp[:, 1]]
    for i in range(1, numBig+1):
        BigNumSmall[i-1, 0] = np.count_nonzero(smallStatsInBigLabel == i)
    result = pd.DataFrame(BigNumSmall, columns=[columnName])

    return result
