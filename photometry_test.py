#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 17:27:27 2019

@author: melissa

Test of teh ztf photometry
"""
import astrobject
import os
import ztfquery
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u



def gaia_cat(sci_img, bk_gr, cat, radius, mypath='/Users/melissa/Phd/'):
    """Extraction of catalogue sources in a science ZTF image using astrobject.
        sci_img : ZTF sciimg from which we want to extract sources     
        bk_gr : 0 or 1; When 0 astrobject doesn't estimate the back_ground of teh image
        (becaise ZTF sciimage are not substracted from back-ground)
        cat :  corresponds to the catalogur from which astrobject exratct the sources
            (can be gaia, sdss...)
        radius : 
            
        ******** Attention to path and sci_img. Depending on the path taht shows
        sci_img, need to configure this point after tests.
    """
    i = astrobject.get_image(os.path.join(mypath, sci_img), background= bk_gr)
    i.download_catalogue(source= cat)
    #if doesn't work successive several times 
    first_data = i.catalogue.data
    #We need to search for soyrces only in the field of view of the camera
    second_data = i.catalogue.data[i.catalogue.fovmask]
    
    i.catalogue.define_around(radius*u.arcsec)
    #We want only solated stars (for the psf forward)
    mask_stars = i.catalogue.get_mask(isolated_only=True)
    gaia_stars_ = second_data[mask_stars]
    
    return gaia_stars_



class Star(object):
    
    def __init__(self, RA, DEC):
        """
        """
        self.RA = RA
        self.DEC = DEC
    
    def ztf_data(self, RA, DEC, data_type):
        """RA DEC from GAIA catalogs, data_type can be : sciimg, psfcat, sexcat...
            it returns the odered list of all the data files
            Need to test (because first use of ztfquery)
        """

        zquery = ztfquery.query.ZTFQuery()
        zquery.load_metadata(radec=[RA, DEC], size=0.01, 
                     sql_query="obsjd BETWEEN 2458270 and 2458370 AND fid=1")
        zquery.download_data(data_type, show_progress=False, nprocess=4, overwrite=True)
        #Use download or get_local, don't know the difference yet
        zquery.get_local_data(data_type)
        listImg = sorted(zquery.get_local_data(data_type))
        
        return listImg
    
    def data_loaded(self,my_path = '/Users/melissa/Phd/'):
        """To use when the data are availible on memory (instead of ztf_data).
        But need to delete the first file and even the second 'ds_store'
        """
        list_Img  = []
       
        for (repertoire, sousRepertoires, fichiers) in os.walk(os.path.join(my_path, 'ZTF/2018/')):
           list_Img.extend(fichiers)
       
       #Time to get the paths of all those files
       
        path = []
       
        for i in range(len(list_Img)):
            for (repertoire, sousRepertoires, fichiers) in os.walk(os.path.join(my_path, 'ZTF/2018/')):
                if  list_Img[i] in fichiers:
                    #path.append(os.path.join(repertoire, list_Img[i]))
                    path.append(repertoire)
      
        return list_Img, path
        
   
    
    def matching_data(self, obsGAIA, listImg, path, dist, gaia_sources, my_path = '/Users/melissa/Phd/', 
                      data_type = "psfcat.fits", ffcq = "F_000762_zg_c06_q2"):
        """Use it if data_loaded used before
        
        osbImg : corresponds to the skycoord of the Gaia sources
        listImg: the list/path of images (need to specify the type of image) where 
        to search for the GAIA sources
        path : paths of the images
        dist : (in arsec) corresponds to the min 2d_dist (more details on match_to_catalog_sky
        from astropy) between the gaia source and the matched ztf source
        gaia_sources : the Gaia catalog exyracted from a ZTF sciimg from which 
        we do the matching 
        ffcq : field filter ccd number quadrant number of the data (we process per
        quadrant per ccd per filter)
        Broowse the images and extract data from the images corresponding 
        to the gaia sources. WE DO IT ONE QUADRANT PER ONE CCD
        
        """
        
        count = 0
        fail = []
        
        for i in range(len(listImg)):
            try:
                
                if listImg[i][38:len(listImg[i])] == data_type :
                
                    
                    count += 1
                    #fImg = fits.open(path[i])
                    fImg = fits.open(os.path.join(path[i], list_Img[i]))
                
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
                
                    ra_sources = [gaia_sources['RA_ICRS'][i] for i in idx_source]
                    dec_sources = [gaia_sources['DE_ICRS'][i] for i in idx_source]
                
                    for k,j,l in zip(idx_Img, ra_sources, dec_sources):
                        #print(j)
                        dat_2 = {"ra": [raImg[k]], "dec": [decImg[k]], "X": [X[k]], 
                                 "Y": [Y[k]], "flux": [fluxImg[k]], "err_f": [errfluxImg[k]], 
                                 "mag": [magImg[k]], "err_m": [errmagImg[k]], "filter_id": filterImg, "src": path[i]}
                    
                        #dat_2 = {'mag': [magImg[k]], "ra": [raImg[k]], "src": path[i]}
                        data_2 = pd.DataFrame(data=dat_2)
                    
                        #data_2 = pd.DataFrame(data=dat_2)
        
                        #path_source = my_path"/ZTF/stars_info/"+ffcq+"/source_"+str(j)+"_"+str(l)+".csv"
                        path_source = os.path.join(my_path, "ZTF/stars_info/", ffcq, "source_"+str(j)+"_"+str(l)+".csv")
                    
                        test_file = os.path.exists(path_source)
                        if test_file == True:
                            data_1 = pd.read_csv(path_source)
                            data_1.drop(data_1.columns[data_1.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
                            new = data_1.append(data_2)
                            new.to_csv(path_source)
                            print('existed', i, path_source)
                        if test_file == False:
                            print('new', i, path_source)
                            data_new = data_2
                            data_new.to_csv(path_source)
            
            except:
                fail.append(path[i])
                
        
        return fail
        
        
    def test_match(self,  obsGAIA, listImg, path, dist, gaia_sources, my_path = '/Users/melissa/Phd/', 
                      data_type = "psfcat.fits"):
        
        """Same function as matching_data but for one source (cad my gaia cat is one source)
        """
        list_idx = []
        count = 0
        for i in range(len(listImg)):
            #try:
                
            if listImg[i][38:len(listImg[i])] == data_type :
                print(i)
                count += 1
                fImg = fits.open(path[i])
                
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
                idx, d2dSN, d3dSN = obsGAIA.match_to_catalog_sky(obsImg)
                
                d2dSNarcsec = d2dSN.to("arcsec").value
                
                flagin = (d2dSNarcsec<dist)
                
                #idx_Img = np.arange(len(idx))[flagin]
                idx_ = np.array([idx])
                idx_Img = idx_[flagin]
                list_idx.append(len(idx_Img))
                #ra_sources = [gaia_sources['RA_ICRS'][i] for i in idx_source]
                #dec_sources = [gaia_sources['DE_ICRS'][i] for i in idx_source]
                
        return list_idx, count
    
    def matching_ztf(self, obsGAIA, listImg, dist, gaia_sources, my_path = '/Users/melissa/Phd/',
                     ffcq = "F_000762_zg_c06_q2"):
        """To use if we got the data using ztfquery. No need to specify the data_type,
        we do it when downloading data
        """
        
        count = 0
        fail = []
        
        for i in range(len(listImg)):
            #try:
            count += 1
            fImg = fits.open(path[i])
            dataImg = fImg[1].data
            dataImgH = fImg[0].header
            obsmjdImg = dataImgH["OBSMJD"]
            #print(obsmjdImg)
            fieldImg = dataImgH["FIELDID"] 
            filterImg = dataImgH["FILTERID"]
                
            raImg, decImg = dataImg["ra"], dataImg["dec"]
            X, Y = dataImg["xpos"], dataImg["ypos"]
            fluxImg, errfluxImg = dataImg["flux"], dataImg["sigflux"]
            magImg, errmagImg = dataImg['mag']+dataImgH["MAGZP"], dataImg['sigmag']
                    
            obsImg = SkyCoord(ra=raImg*u.degree, dec=decImg*u.degree)
            idx, d2dSN, d3dSN = obsImg.match_to_catalog_sky(obsGAIA)
            d2dSNarcsec = d2dSN.to("arcsec").value
            flagin = (d2dSNarcsec<dist)
            d2dflag = d2dSNarcsec[flagin]
            idx_Img = np.arange(len(idx))[flagin]
            idx_source = idx[flagin]
        
            ra_sources = [gaia_sources['RA_ICRS'][i] for i in idx_source]
            dec_sources = [gaia_sources['DE_ICRS'][i] for i in idx_source]
                
            for k,j,l in zip(idx_Img, ra_sources, dec_sources):
                
                dat_2 = {"ra": [raImg[k]], "dec": [decImg[k]], "X": [X[k]], 
                         "Y": [Y[k]], "flux": [fluxImg[k]], "err_f": [errfluxImg[k]], 
                         "mag": [magImg[k]], "err_m": [errmagImg[k]], "filter_id": filterImg, "src": path[i]}
                    
                data_2 = pd.DataFrame(data=dat_2)
                    
                path_source = my_path/"ZTF/stars_info/"+ffcq+"/source_"+str(j)+"_"+str(l)+".csv"
                    
                test_file = pathi.exists(path_source)
                if test_file == True:
                    data_1 = pd.read_csv(path_source)
                    data_1.drop(data_1.columns[data_1.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
                    new = data_1.append(data_2)
                    new.to_csv(path_source)
                    print('existed', i, path_source)
                if test_file == False:
                    print('new', i, path_source)
                    data_new = data_2
                    data_new.to_csv(path_source)
            
            #except:
            #    fail.append(path[i])
                
        
    def stat_sources(self, gaia_sources, ffcq = "F_000762_zg_c06_q2",
                     my_path = "/Users/melissa/Phd/", filter_ = 1 ):
        """Browse all the sources file got from the matching to get : mean,ùedian, 
        std_median, number of data, ra, dec per source in ONE filter (here 1 : G)
        """
        sources_ra = gaia_sources['RA_ICRS']
        sources_dec = gaia_sources['DE_ICRS']
        ra = []
        dec =[]
        mean_mag = []
        median_mag = []
        std_mag = []
        file_unfound = []
        nb_data = []
        path_= []
        count = 0
        for w,v in zip(sources_ra, sources_dec):
            count += 1
            try:
                path_ = os.path.join(my_path, "ZTF/stars_info/", ffcq, "source_"+str(w)+"_"+str(v)+".csv")
                #path_ = my_path/"ZTF/stars_info/"+ffcq+"/source_"+str(w)+"_"+str(v)+".csv"
                files = pd.read_csv(path_)
                data_g = files[files['filter_id'] == filter_]
                nb_data_, __ = data_g.shape
                nb_data.append(nb_data_)
                mean_source = np.mean(data_g['mag'])
                mean_mag.append(mean_source)
                median_source = np.median(data_g['mag'])
                median_mag.append(median_source)
                std_source = np.std(data_g['mag'])
                std_mag.append(std_source)
                ra.append(str(w))
                dec.append(str(v))
            
            except:
                file_unfound.append(path_)
                
        #build up a pd dataframe with the collected data
        donnees_ = {"ra": ra, "dec": dec, "median_m": median_mag, "std_m": std_mag,
             "mean_m": mean_mag ,"npoints": nb_data}
        
        data_sources = pd.DataFrame(data=donnees_)
        #data_sources.to_csv(my_path/'ZTF/stars_info/ztf_id_sources/data'+ffcq+'.csv')
        
        return data_sources, path_
    

