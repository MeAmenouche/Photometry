#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 09:24:45 2019

@author: melissa

Soft to create the standalone catalogue for the SSP Model
"""


import numpy as np
import os
import pandas as pd
import Cat_Extraction.photometry_test as pht
from astropy.io import fits
from astropy import units as u


#We will import cat_extraction.py as a module to use all its functions to build up
#a standalone catalogue

class Star_list(object):
    """
    """
    
    def __init__(self, Ra, Dec):
        """
        """
        self.Ra = Ra
        self.Dec = Dec
    
    def get_data(self, Ra, Dec):
        """"We will use the ztf_method from photometry_test to get the data around this region
        """
        star_obj = pht.Star(Ra, Dec)
        star_obj.ztf_data(Ra, Dec,  ['sciimg.fits'])
        
        #should return alist of all the sciimg, downloaded via ztfqeury or present...
    def select_seeing(self, listImg):
            """Look for a "refrence" image with the best seeing (need to check its unity)
            ############## Check for the listImg..... should contain the path
            """
            fits_data_ = fits.open(listImg[0])
            seeing_max= fits_data_[0].header['SEEING']
                
            for i in range(1, len(listImg)):
                fits_data = fits.open(listImg[i])
                seeing_img[i] = fits_data[0].header['SEEING']
            
                if seeing_img[i] > seeing_max:
                    seeing_max = seeing_img[i]
                    img = listImg[i]
            
            return seeing_max, img
        
    def source_extraction(self, img, bk_gr, cat, radius, mag, mypath='/Users/melissa/Phd/'):
        """!!!!!!!!! attention si img (apth de sciimg) démarre de Users....
        """
        source_cat = pht.gaia_cat(img, bk_gr, cat, radius)
        source_cut = #source_cat mettre le critère en magnitude pour éviter les étoiles qui saturent
        
        fits_img = fits.open(img)
        
        
        
        
        obs_src = SkyCoord(ra=source_cat['ra']*u.degree, dec=source_cat['dec']*u.degree)
        wcs_img = wcs.WCS(header=fits_img[0].header)

        xp, yp = obs_src.to_pixel(wcs=wcs_img)
        xp,yp