#!/usr/bin/env python

'''
This library file contains functions for generating and manipulating Biopython
graphics objects.
'''

import colorsys
import itertools
import math
import numpy

def RGB_to_hex(RGBlist):
    '''
    Convert an RGB color into a HEX string (required for some display functions)
    '''
    n = lambda x: int(x*255)
    RGB256 = [(n(r),n(g),n(b)) for r,g,b in RGBlist]
    colors = ['#%02x%02x%02x' % (r, g, b) for r, g, b in RGB256]
    return colors

def colormap(valuelist):
    '''
    Generate a list of divergent colors for use with labeling SeqFeature objects
    '''
    values = numpy.unique(valuelist)
    N = len(values)
    #we will vary in 2 dimensions, so this is how many steps in each
    perm = int(math.ceil(math.sqrt(N)))
    #need offset, as humans can't tell colors that are unsaturated apart
    H = [(x*1.0/perm) for x in range(perm)]
    S = [(x*1.0/perm)+0.2 for x in range(perm)]
    #we will use this to truncate at correct length
    V = [0.7]*N
    # Create all combinations of our colors.                                                                                                                                                           
    HS = itertools.product(H, S)
    H, S = zip(*HS)
    HSV = zip(H,S,V)
    RGB = [colorsys.hsv_to_rgb(h,s,v) for h, s, v in HSV]
    colorlookup = dict(zip(values, RGB[:N]))
    return colorlookup
