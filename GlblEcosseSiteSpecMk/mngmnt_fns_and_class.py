#-------------------------------------------------------------------------------
# Name:        mngmnt_fns_and_class.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'mngmnt_fns_and_class.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import join, normpath, isfile
from netCDF4 import Dataset
from glob import glob
from numpy import seterr, ma, arange

ngranularity = 120
month_names_short = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
crop_name_dict = {'Winter Wheat': 'Wheat_Winter', 'Spring Wheat': 'Wheat_Spring', 'Grain Maize': 'Maize',
                                                 'Spring Oats': 'Oats_Spring', 'Spring Barley': 'Barley_Spring'}
var_name_dict = {'Winter Wheat': 'wheat', 'Spring Wheat': 'wheat', 'Grain Maize': None,
                                                 'Spring Oats': None, 'Spring Barley': 'barley'}

def check_mask_location(mask_defn, site_rec, varname, resol_deg):
    '''
    resol_deg is size of cell in degrees
    masked elements indicate zero, i.e. not pasture, therefore if all the elements are zero then count as not pasture
    '''
    gran_lat, gran_lon, lat, lon, dummy, dummy = site_rec

    res_d2 = resol_deg/2

    lat_ur = lat + res_d2
    lon_ur = lon + res_d2
    lat_ur_indx, lon_ur_indx, ret_code = mask_defn.get_nc_coords(lat_ur, lon_ur)

    lat_ll = lat - res_d2
    lon_ll = lon - res_d2
    lat_ll_indx, lon_ll_indx, ret_code = mask_defn.get_nc_coords(lat_ll, lon_ll)

    vals = mask_defn.nc_dset.variables[varname][lat_ll_indx:lat_ur_indx + 1, lon_ll_indx:lon_ur_indx + 1]
    if ma.is_masked(vals):
        if vals.count() == 0:  # Count the non-masked elements of the array along the given axis.
            return False
        else:
            return True
    else:
        val_mean = vals.mean()
        if val_mean > 0:
            return True
        else:
            return False

def fetch_vals_from_dset(metric_defn, lat, lon, nyears, varname = None):
    '''
    read data from NC file depending on resource

    NB  cnvrsn_fctr = metric_defn.cnvrsn_fctr
        vals = [round(val_raw.item() * cnvrsn_fctr, 3) for val_raw in vals_raw]
    '''
    lat_indx, lon_indx, ret_code = metric_defn.get_nc_coords(lat, lon)

    if varname is None:
        varname = metric_defn.varname

    rsrce = metric_defn.rsrce
    if rsrce == 'fertiliser':
        val_raw = metric_defn.nc_dset.variables[varname][:, :, lat_indx, lon_indx]
        vals = nyears*[float(val_raw[0][0].item())]
    elif rsrce == 'yields':
        vals_raw = metric_defn.nc_dset.variables[varname][:, lat_indx, lon_indx]
        vals = [round(val_raw.item(), 3) for val_raw in vals_raw]
    else:       # sowing and harvest
        vals_raw = metric_defn.nc_dset.variables[varname][:, lat_indx, lon_indx]
        vals = [round(val_raw.item()) for val_raw in vals_raw]

    # check length of vals and extend if short
    # ========================================
    if rsrce != 'fertiliser':
        nvals = len(vals)
        if nvals < nyears:
            nyrs_add = nyears - nvals
            vals = nyrs_add*[vals[0]] + vals

    return vals

def identify_datasets(project_path, mask_fname):
    '''
    called at startup 
    '''
    descriptor = '\tNC files - mask file: '

    if isfile(mask_fname):
        descriptor += 'OK\t'
        print('\tmask file: ' + mask_fname)
    else:
        descriptor += 'none\t'
        print('\tmask file: none')

    # yields  NB yield is a Python reserved word
    # ==========================================
    resource = 'yields'
    descriptor += resource + ': '

    fns_path = join(project_path, 'Yield data')
    yield_fnames = glob(fns_path + '\\*.nc')
    if len(yield_fnames) == 0:
        descriptor += 'none\t'
        print('\t' + resource + ' file: none')
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' file: ' + yield_fnames[0])

    # fertiliser
    # ==========
    resource = 'fertiliser'
    descriptor += resource + ': '

    fert_path = join(project_path, 'Fertiliser')
    fert_fnames = glob(fert_path + '\\Napprate*.nc')
    if len(fert_fnames) == 0:
        descriptor += 'none\t'
        print('\t' + resource + ' file: none')
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' file: ' + fert_fnames[0])

    # sowing_harvest
    # ==============
    resource = 'phenology'
    descriptor += resource + ': '

    pheno_path = join(project_path, 'Phenology')
    pheno_fnames = glob(pheno_path + '\\' + '*.nc')
    nfiles = len(pheno_fnames)
    if nfiles == 0:
        descriptor += 'none\t'
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' first of {} files: {}'.format(nfiles, pheno_fnames[0]))

    return descriptor

def create_proj_data_defns(project_path, mask_fname, crop_name, req_resol_deg):
    '''

    '''
    if crop_name not in crop_name_dict:
        print('*** Error *** ' + crop_name + ' not mappable to phenology NC name')
        return None
    crop_mapped = crop_name_dict[crop_name]

    crop = crop_name.lower()
    resol = str(req_resol_deg)

    resource = 'cropmasks'
    mask_defn = ManagementSet(mask_fname, resource)

    # for now we just use the mean of the yields from 2000 to 2014;  NB yield is a Python reserved word
    # =================================================================================================
    resource = 'yields'
    fns_path = join(project_path, 'Yield data')
    yield_fname = glob(fns_path + '\\*' + crop_mapped + '*.nc')
    if len(yield_fname) == 0:
        print('*** Error *** ' + fns_path + ' no yield file found for ' + crop_mapped)
        return None

    yield_defn = ManagementSet(yield_fname[0], resource, 'YIELD', units = 't/ha')

    # fertiliser
    # ==========
    resource = 'fertiliser'
    fert_path = join(project_path, 'Fertiliser')
    fert_fname = glob(fert_path + '\\*' + crop_mapped + '*.nc')
    if len(fert_fname) == 0:
        print('*** Error *** ' + fert_path + ' no ' + resource + ' files found')
        return None

    fert_defn = ManagementSet(fert_fname[0], resource, var_name_dict[crop_name] + 'Napprate', 'N kg/ha')

    # sowing_harvest
    # ==============
    resource = 'phenology'
    pheno_path = join(project_path, 'Phenology')
    pheno_fname = glob(pheno_path + '\\' + crop_mapped + '*.nc')
    if len(pheno_fname) == 0:
        print('*** Error *** ' + pheno_path + ' no phenology file found')
        return None

    pheno_defn = ManagementSet(pheno_fname[0], resource)

    return mask_defn, yield_defn, pheno_defn, fert_defn

def open_proj_NC_sets(mask_defn, yield_defn, pheno_defn, fert_defn):
    '''
    '''
    mask_defn.nc_dset  = Dataset(mask_defn.nc_fname, mode='r')
    yield_defn.nc_dset = Dataset(yield_defn.nc_fname, mode='r')
    pheno_defn.nc_dset = Dataset(pheno_defn.nc_fname, mode='r')
    fert_defn.nc_dset = Dataset(fert_defn.nc_fname, mode='r')

    return

def close_proj_NC_sets(mask_defn, yield_defn, pheno_defn, fert_defn):

    mask_defn.nc_dset.close()
    yield_defn.nc_dset.close()
    pheno_defn.nc_dset.close()
    fert_defn.nc_dset.close()

    return

class ManagementSet(object, ):
    '''

    '''
    def __init__(self, nc_fname, resource, rqrd_varname = None, units = None, cnvrsn_fctr = None):
        '''

        '''
        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname, mode='r')
        if 'lat' in nc_dset.variables:
            lat_var = 'lat'
            lon_var = 'lon'
        else:
            lat_var = 'latitude'
            lon_var = 'longitude'
        lats = nc_dset.variables[lat_var][:]
        lons = nc_dset.variables[lon_var][:]

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        start_year = None
        end_year   = None
        varnames  = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                varnames.append(var)

            if var == 'time':
                time_var = nc_dset.variables[var]
                if hasattr(time_var, 'units'):
                    time_units = time_var.units
                else:
                    time_units = 'year'

                if resource == 'yields' or resource == 'fertiliser':
                    start_year = time_var[0].item()
                    end_year = time_var[-1].item()
                else:
                    '''
                    since_time = time_units.split('since')[1]
                    start_year = int(since_time.split('-')[0]) + time_var[0]    # messy way to get to 1961
                    '''
                    start_year = 1981
                    end_year = start_year + len(time_var) - 1

        nc_dset.close()

        lat_frst = float(lats[0])
        lon_frst = float(lons[0])
        lat_last = float(lats[-1])
        lon_last = float(lons[-1])

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.lat_frst = float(lats[0])
        self.lon_frst = float(lons[0])
        self.lat_last = float(lats[-1])
        self.lon_last = float(lons[-1])

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.varnames = varnames
        self.varname = rqrd_varname
        self.units = units
        self.cnvrsn_fctr = cnvrsn_fctr
        self.nc_dset = None

        # resolutions
        # ===========
        self.resol_lon = (lons[-1] - lons[0])/(len(lons) - 1)
        self.resol_lat = (lats[-1] - lats[0])/(len(lats) - 1)
        self.max_lat_indx = len(lats) - 1
        self.max_lon_indx = len(lons) - 1

        #
        self.lats = list(lats)
        self.lons = list(lons)
        self.start_year = start_year
        self.end_year   = end_year
        self.rsrce = resource

    def get_nc_coords(self, latitude, longitude):

        ret_code = 'OK'

        lat_indx = int(round((latitude -  self.lat_frst)/self.resol_lat))
        lon_indx = int(round((longitude - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            ret_code = '*** Warning *** latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                                                                                round(latitude, 4), self.max_lat_indx)
        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            ret_code = '*** Warning *** longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                                round(longitude, 4), self.max_lon_indx)
        return lat_indx, lon_indx, ret_code
