#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 13:49:47 2019

@author: melissa
"""

import astropy.units as u
import numpy as np
import os.path as op
from astropy.io import fits
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import coordinates, units, wcs, time
import os
import astrobject

def get_gaia(path, source_cat, path_cat ='/Users/melissa/Phd/ZTF/catalogs/gaia_astrobject.csv'):
    """To load GAIA cat using another  ZTF SCI Image. Otherwise, use the path_cat
    """
    i = astrobject.get_image(path, background=0)
    i.download_catalogue(source = source_cat)
    data_cat = pd.DataFrame(data=i.catalogue.data)
    data_cat.to_csv(path_cat)
    return data_cat

def get_files(path_rep):
    """Construct a list with the names of all the ZTF files in sci rep
    """
    listImg = []
    for (repertoire, sousRepertoires, fichiers) in os.walk(path_rep):
        listImg.extend(fichiers)
    listImg.pop(0)
    return listImg

def get_path(list_files, path_rep):
    """Get the paths of all the files from get_files
    """
    path = []
    for i in range(len(list_files)):
        for (repertoire, sousRepertoires, fichiers) in os.walk(path_rep):
            if  list_files[i] in fichiers:
                path.append(os.path.join(repertoire, list_files[i]))
    return path

def extract_sources(file, gaia_cat):
    """Special psfcat
    """
    obsGAIA = coordinates.SkyCoord(ra=gaia_cat['RA_ICRS'], dec=gaia_cat['DE_ICRS'], unit=units.deg)
    fImg = fits.open(file)
    data_img = fImg[1].data
    raImg, decImg, fluxImg,  = data_img["ra"], data_img["dec"], data_img["flux"]
        
    errfluxImg, magImg, errmagImg =  data_img["sigflux"], data_img['mag'], data_img['sigmag']
    
    obsImg = SkyCoord(ra=raImg*u.degree, dec=decImg*u.degree)
    idx, d2d, d3d = obsImg.match_to_catalog_sky(obsGAIA)
    d2darcsec = d2d.to("arcsec").value

    flagin = (d2darcsec<1.0)
    d2dflag = d2darcsec[flagin]          
    
    idxImg = np.arange(len(idx))[flagin]
    idxGAIA = idx[flagin]
    
    return idxImg, idxGAIA

def get_cood(id_gaia):
    """To get the GAIA cood of a ZTF source using gaia dataframe and the indexes got from 
    the match to cat sky. Here, id_gaia is one element from idxGAIA got from extract_sources
    """
    cood = tab.to_string().split()
    ra = float(cood[0])
    dec = float(cood[1])
    return ra, dec
    
class Source(object):
    
    def __init__(self, ra_gaia, dec_gaia, listImg):
        self.ra_gaia = ra_gaia
        self.dec_gaia = dec_gaia
        #self.path = path
        self.cood_gaia = SkyCoord(ra=self.ra_gaia*u.degree, dec=self.dec_gaia*u.degree)
        self.listImg = listImg
        self.ra = []
        self.dec = []
        self.x = []
        self.y = []
        self.f = []
        self.sig_f = []
        self.m = []
        self.sig_m = []
        #self.date = [] à rajouter si on veut
    
    def get_data(self, list_path):
        """
        """
        data_type = "psfcat.fits"
        fail_file= []
        
        for i in range(len(self.listImg)):
            print(i)
            
            if (self.listImg[i][38] == 'p' or i == 211):
                try:
                    fImg = fits.open(list_path[i])
                    dataImg = fImg[1].data
                    dataImgH = fImg[0].header
                    obsmjdImg = dataImgH["OBSMJD"]
                    fieldImg = dataImgH["FIELDID"] 
                    filterImg = dataImgH["FILTERID"]
                    
                    raImg, decImg = dataImg["ra"], dataImg["dec"]
                    X, Y = dataImg["xpos"], dataImg["ypos"]
                    fluxImg, errfluxImg = dataImg["flux"], dataImg["sigflux"]
                    magImg, errmagImg = dataImg['mag']+dataImgH["MAGZP"], dataImg['sigmag'] 
                    obsImg = SkyCoord(ra=raImg*u.degree, dec=decImg*u.degree)
                    idx, d2dSN, d3dSN = self.cood_gaia.match_to_catalog_sky(obsImg)
                    self.ra.append(raImg[idx])
                    self.dec.append(decImg[idx])
                    self.x.append(X[idx])
                    self.y.append(Y[idx])
                    self.f.append(fluxImg[idx])
                    self.sig_f.append(errfluxImg[idx])
                    self.m.append(magImg[idx])
                    self.sig_m.append(errmagImg[idx])
                    
                except:
                    #continue
                    fail_file.append(list_path[i])
                
        source_ = {"ra": self.ra, "dec": self.dec, "X": self.x, 
                   "Y": self.y, "flux": self.f, "err_f": self.sig_f,
                   "mag" : self.m, "err_m" : self.sig_m}
        self.data_source = pd.DataFrame(data=source_)
            
        return self.data_source, fail_file

if __name__ == "__main__":
    gaia_cat = pd.read_csv('/Users/melissa/Phd/ZTF/catalogs/gaia_astrobject.csv')
    list_img = get_files('/Users/melissa/Phd/ZTF/2018/')
    paths = get_path(list_img, '/Users/melissa/Phd/ZTF/2018/')
    ########## à régelr cette histoire de pop#####
    paths.pop(0)
    paths.pop(0)
    ########## à régler cette histoire de pop #########
    index_img, index_gai = extract_sources(
            '/Users/melissa/Phd/ZTF/2018/0609/351076/ztf_20180609351076_000762_zg_c06_o_q2_psfcat.fits',
            gaia_cat)
    
      