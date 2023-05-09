"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import math
import sys
import time
import os

import MDAnalysis
import argparse
import numpy as np

import scipy.interpolate as si
import scipy.stats
import warnings
from collections import defaultdict
from mpl_toolkits.axes_grid1 import make_axes_locatable
from multiprocessing.dummy import Pool as ThreadPool

import CDPL.Base as Base
import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Grid as Grid
import CDPL.Math as Math
import CDPL.Pharm as Pharm

import pyevtk.hl as hl


if __name__ == '__main__':
    print('> All dependencies are satisfied.')
