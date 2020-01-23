#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:13:19 2019

@author: melissa
"""

import astrobject
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from astropy.io import fits
from astroquery import vizier
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astrobject.catalogue import catalogues as cat
import Cat_Extraction.photometry_test as pht


class ZTF_field(object):
    """Creation of an object :
        
        *ra dec : attributes  of the object
        *We will extract a gaia catalog using the fetch_gaia_catalogue 
        *We will match those gaia sources with one on the several psfcats 
    """
    def __init__(self, ra, dec):
        """
        """
        self.ra = ra
        self.dec = dec
    
    ######################################################################
    
    def fetch_gaia_catalogue(self, radius, r_unit="deg", extracolumns=["Var"],column_filters={"Gmag":"0..25"},**kwargs):
        """ query online gaia-catalogue in Vizier (I/345, DR2) using astroquery.
        This function requieres an internet connection.
        
        Parameters
        ----------
        center: [string] 'ra dec'
        position of the center of the catalogue to query.
        
        radius: [string] 'value unit'
        radius of the region to query. For instance '1d' means a
        1 degree raduis

        extracolumns: [list-of-string] -optional-
        Add extra column from the V/139 catalogue that will be added to
        the basic query (default: position, ID, object-type, magnitudes)
        column_filters: [dict] -optional-
        Selection criterium for the queried catalogue.

        **kwargs goes to astroquery.vizier.Vizier

        Returns
        -------
        GAIA Catalogue (child of Catalogue)
        """
        columns = ["Source","RA_ICRS","e_RA_ICRS","DE_ICRS","e_ED_ICRS","Var"]
        #for band in SDSS_INFO["bands"]:

        #try:
        coord = SkyCoord(ra=self.ra,dec=self.dec, unit=(u.deg,u.deg))
        angle = Angle(radius,r_unit)
        v = vizier.Vizier(columns, column_filters=column_filters)
        v.ROW_LIMIT = -1
        t = v.query_region(coord, radius=angle,catalog="I/345/gaia2").values()[0]
        
        return t

    ######################################################################
    
    def cat_extraction(self, radius, catalogue='GAIA', **kwargs):
        """Extract a catalogue using the ra dec coordinates.  
        We start by calling the fetch_gaia_catalogue to apply it to our object.
        We create an astrobject catalogue via (astrobject.catalogue 
        import catalogues as cat) from what we got with fetch_gaia_catalogue.
        We apply the fovmask on the catalogue and use the define_around with 
        the mask of isolated stars. This returns a gaia catalogue.
        """
        t_cat = self.fetch_gaia_catalogue(radius)
        astro_cat = cat.get_catalogue(t_cat, catalogue)
        cat_ztf = astro_cat.data[astro_cat.fovmask]
        astro_cat.define_around(radius*u.arcsec)
        mask_stars = astro_cat.get_mask(isolated_only=True)
        gaia_cat = cat_ztf[mask_stars]
        
        return gaia_cat
    
    ######################################################################
    
    def match_data(self, obsGAIA, Img, path_img, dist, gaia_sources):
        """
        Match one ZTF psfcat with the gaia catalogue.
        obsGAIA : the gaia sources skycoord (stpe on the notebook)
        Img : the name of the psfcat image we want to match to
        path_img : the path  to the repository where the image is
        dist : the min distance required between one source from the psfcat and another
        from the gaia cat
        gaia_sources : the gaia  catalogue
        
        It makes the identification and write on a file (saved on the repository
        of the psfcat) the id of the gaia source and the line where the ztf source
        had been found on the psfcat
        """
        #print(Img, path_img)
        fImg = fits.open(os.path.join(path_img, Img))
        #print(os.path.join(path_img, Img))
        dataImg = fImg[1].data
        dataImgH = fImg[0].header
        obsmjdImg = dataImgH["OBSMJD"]
        #print(obsmjdImg)
        fieldImg = dataImgH["FIELDID"] 
        filterImg = dataImgH["FILTERID"]
                    
        #raImg, decImg,X, Y, fluxImg, errfluxImg, magImg, errmagImg = dataImg["ra"], dataImg["dec"], dataImg["xpos"], dataImg["ypos"], dataImg["flux"], 
        #dataImg["sigflux"], dataImg['mag']+dataImgH["MAGZP"], dataImg['sigmag']
                
        raImg, decImg = dataImg["ra"], dataImg["dec"]
        X, Y = dataImg["xpos"], dataImg["ypos"]
        fluxImg, errfluxImg = dataImg["flux"], dataImg["sigflux"]
        magImg, errmagImg = dataImg['mag']+dataImgH["MAGZP"], dataImg['sigmag']
        obsImg = SkyCoord(ra=raImg*u.degree, dec=decImg*u.degree)
        #idx, d2dSN, d3dSN = obsImg.match_to_catalog_sky(obsGAIA)
        idx, d2dSN, d3dSN = obsImg.match_to_catalog_sky(obsGAIA)
        d2dSNarcsec = d2dSN.to("arcsec").value
        flagin = (d2dSNarcsec<dist)
        d2dflag = d2dSNarcsec[flagin]
        idx_Img = np.arange(len(idx))[flagin]
        idx_source = idx[flagin]
        new_cat = gaia_sources.iloc[idx_source]
        data_ = {"id_source": new_cat['Source'], "idpsf": idx_Img}
        data__ = pd.DataFrame(data=data_)
        #data__.to_csv(os.path.join(path_img, 'gaia_sources.csv'))
        data__.to_csv(os.path.join(path_img, Img[0:44]+'_gaia.csv'))
        
        #We need to remove a hidden file ".DS_Store" in each directory
        #ds_file = os.path.exists(os.path.join(path_img, '.DS_Store'))
        
    
    

class LC(ZTF_field):
    
    def __init__(self, ra, dec):
        
        ZTF_field.__init__(self, ra, dec)
        
    def parcours(self, name_csv, path_csv, name_img, dicmag, filter_ = 1):
        """
        """
        value = bool(dicmag)
        #Load the csv file that contain the number of line 
        data = pd.read_csv(path_csv+'/'+name_csv)
        
        #Read the psfcat from the same directory than the csv file
        file_img = fits.open(os.path.join(path_csv, name_img))
        data_img = file_img[1].data
        head_img = file_img[0].header
        filterImg = head_img["FILTERID"]
        
        #Studies done filter by filter
        if filterImg == filter_: 
            imgmag = data_img['mag']+head_img['MAGZP']
            mag_match = imgmag[data['idpsf'].values] 
            
            #Creation ou remplissage of a dictionnary with KEY :id_source, VALUE : 
            #list of magnitudes
        
            """if value == True:
                for r, l in zip(data['id_source'], mag_match):
                    dicmag[r].append(l) if (r in dicmag.keys()) else dicmag.update(r = l)
        
            else:
                dicmag = {data['id_source'][i]: [mag_match[i]] for i in range(len(mag_match))}
                
            #if 'r' in dicmag.keys():
                #print(name_csv, path_csv, name_img)"""
            
            for k, l in zip(data['id_source'], mag_match):
                dicmag[k].append(l)
    
        return dicmag
        
        
        
        
        
        
        
        
        
        
        
        
        
 
