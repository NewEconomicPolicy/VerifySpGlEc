#-------------------------------------------------------------------------------
# Name:        validate.py
# Purpose:     Class to organise run parameters
# Author:      Mike Martin
# Created:     28/09/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'validate.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
import os
import json

def check_harv_doys_adjust(sow_days, harv_days):
    '''
    add 365 days if sowing date later than harvest
    '''
    if sow_days[0] > harv_days[0]:
        harv_days_adust = [doy + 365 for doy in harv_days]
    else:
        harv_days_adust = harv_days

    return harv_days_adust

def check_phenology_fnames(form):
    '''
    fertiliser and crop input files should look include the crop type in their names
    '''

    sowing_crop_type = 'unknown'
    fert_crop_type = 'unknown';
    fert_field = ''
    form.crop_type_fert = {'type': fert_crop_type, 'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    fert_fname = form.w_lbl14.text()
    if not os.path.isfile(fert_fname):
        return 'Fertiliser input file ' + fert_fname + ' does not exist'
    try:
        with open(fert_fname, 'r') as fobj:
            print('Fertiliser file ' + fert_fname + ' is readable')

    except (OSError, IOError) as err:
        print(err)
        return 'Could not read fertiliser input file'

    fert_crop_type = 'spring wheat'
    # fert_crop_type = 'maize'

    # sowing dates
    # ============
    sowing_fname = form.w_lbl15.text()
    if not os.path.isfile(sowing_fname):
        return 'Sowing_dates input file ' + sowing_fname + ' does not exist'
    try:
        with open(sowing_fname, 'r') as fobj:
            print('Sowing dates input file ' + sowing_fname + ' is readable')

    except (OSError, IOError) as err:
        print(err)
        return 'Could not read sowing dates input file'

    sowing_crop_type = 'spring wheat'
    #    sowing_crop_type = 'maize'

    form.w_create_files.setEnabled(False)
    # form.w_create_files.setEnabled(True)

    form.crop_type_fert = {'type': fert_crop_type, 'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    # disable simulation years
    # ========================
    if form.years_from_flag == 'sowing':
        form.combo11s.setEnabled(False)
        form.combo11e.setEnabled(False)

    return 'tout est bon'

def check_fert_sowing_fnames(form):
    '''
    fertiliser and crop input files should look include the crop type in their names
    '''

    sowing_crop_type = 'unknown'
    fert_crop_type = 'unknown'; fert_field = ''
    form.crop_type_fert   = {'type': fert_crop_type,   'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    fert_fname = form.w_lbl14.text()
    if not os.path.isfile(fert_fname):
        return 'Fertiliser input file ' + fert_fname + ' does not exist'
    try:
        with open(fert_fname, 'r') as fobj:
            print('Fertiliser file ' + fert_fname + ' is readable')

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read fertiliser input file'

    if fert_fname.lower().rfind('swhe') > 0:
        fert_crop_type = 'spring wheat'
        fert_field = 'SWHE'
    elif fert_fname.lower().rfind('maiz') > 0:
        fert_crop_type = 'maize'
        fert_field = 'MAIZ'

    # sowing dates
    # ============
    sowing_fname = form.w_lbl15.text()
    if not os.path.isfile(sowing_fname):
        return 'Sowing_dates input file ' + sowing_fname + ' does not exist'
    try:
        with open(sowing_fname, 'r') as fobj:
            print('Sowing dates input file ' + sowing_fname + ' is readable - simulation years will be disabled')

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read sowing dates input file'
    
    if sowing_fname.lower().rfind('swheat') > 0:
        sowing_crop_type = 'spring wheat'
    elif sowing_fname.lower().rfind('maize') > 0:
        sowing_crop_type = 'maize'

    # check file contents
    # ===================
    mess = ''
    if fert_crop_type == 'unknown' or sowing_crop_type == 'unknown':
        mess = '*** Error *** unknown crop type or sowing crop type'
        if fert_crop_type == 'unknown':
            print('*** Error *** fertiliser crop type is unknown - cannot enable simulation file creation')
        if sowing_crop_type == 'unknown':
            print('*** Error *** sowing crop type is unknown - cannot enable simulation file creation')
        form.w_create_files.setEnabled(False)
    else:
        mess = 'fertiliser for crop ' + fert_crop_type + ' and sowing dates for crop ' + sowing_crop_type + \
                                                                                            ' input files are valid'
        form.w_create_files.setEnabled(True)

    form.crop_type_fert   = {'type': fert_crop_type,   'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    # disable simulation years
    # ========================
    if form.years_from_flag == 'sowing':
        form.combo11s.setEnabled(False)
        form.combo11e.setEnabled(False)

    return mess