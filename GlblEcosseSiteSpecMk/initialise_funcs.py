"""
#-------------------------------------------------------------------------------
# Name:        initialise_funcs.py
# Purpose:     script to read read and write the setup and configuration files
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'initialise_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os import getcwd, remove, makedirs, mkdir, name as os_name
from os.path import exists, normpath, getmtime, splitext, isfile, isdir, join, lexists

import json
import logging
from time import sleep, strftime
from sys import exit
from glob import glob
from hwsd_mu_globals_fns import HWSD_mu_globals_csv
from shape_funcs import format_bbox, calculate_area
from weather_datasets import read_weather_dsets_detail
from mngmnt_fns_and_class import identify_datasets

bbox_default = [116.90045, 28.2294, 117.0, 29.0] # bounding box default - somewhere in SE Europe
sleepTime = 5

GLBL_ECSSE_STR = 'global_ecosse_config_site_specific_'
APPLIC_STR = 'glbl_ecss_site_spec_mk'
ERROR_STR = '*** Error *** '

def initiation(form):
    '''
    this function is called to initiate the programme to process non-GUI settings.
    '''
    # retrieve settings
    _read_setup_file(form)
    if len(form.weather_sets) == 0:
        print(ERROR_STR + 'No valid weather sets found')
        sleep(sleepTime)
        exit(0)

    form.parms_settings = _read_site_specific_parms()

    # set years_from_flag to gui or 'sowing' to take simulation period from sowing dates file
    # =======================================================================================
    form.years_from_flag = 'gui'
    form.req_resol_deg = None

    # facilitate multiple config file choices
    # =======================================
    config_files = glob(form.config_dir + '/' + GLBL_ECSSE_STR + '*.txt')
    studies = []
    for fname in config_files:
        dummy, remainder = fname.split(GLBL_ECSSE_STR)
        study, dummy = splitext(remainder)
        if study != '':
            studies.append(study)
    form.studies = studies

    if len(config_files) > 0:
        form.config_file = config_files[0]
    else:
        form.config_file = form.config_dir + '/' + GLBL_ECSSE_STR + 'dummy.txt'

    # for log file
    date_stamp = strftime('_%Y_%m_%d_%I_%M_%S')
    log_file_name = 'global_ecosse_min' + date_stamp + '.log'

    # crop_pars, fnames.dat and model_switches must be present
    # =========================================================
    fname_model_switches = 'Model_Switches.dat'
    cwDir = getcwd()
    dflt_model_switches = join(cwDir, fname_model_switches)
    if isfile(dflt_model_switches):
        form.dflt_model_switches = dflt_model_switches
    else:
        print('{} file does not exist in directory {}'.format(fname_model_switches, cwDir))
        sleep(sleepTime)
        exit(0)

    fnames_dat = 'fnames.dat'
    dflt_fnames_dat = join(cwDir, fnames_dat)
    if isfile(dflt_fnames_dat):
        form.dflt_fnames_dat = dflt_fnames_dat
    else:
        print('{} file does not exist in directory {}'.format(fnames_dat, cwDir))
        sleep(sleepTime)
        exit(0)

    dflt_crop_sun = 'CROP_SUN.DAT'
    crop_pars_fname = join(cwDir, dflt_crop_sun)
    if isfile(crop_pars_fname):
        print('Will use {} ECOSSE crop definition file'.format(crop_pars_fname))
        form.dflt_crop_sun = crop_pars_fname
    else:
        print('{} file does not exist in directory {}'.format(crop_pars_fname, cwDir))
        sleep(sleepTime)
        exit(0)

    _read_crop_codes(form, crop_pars_fname)

    # set up logging
    # ==============
    log_dir = form.log_dir
    log_fname = join(log_dir, log_file_name)

    # Setup up initial logger to handle logging prior to setting up the
    # full logger using config options
    form.lgr = logging.getLogger('glbl_ecsse')
    form.lgr.setLevel(logging.INFO)
    fh = logging.FileHandler(log_fname)    # send log recs (created by loggers) to appropriate destination
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    # formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)  # specify layout of records
    form.lgr.addHandler(fh)

    if not check_sims_dir(form):
        sleep(sleepTime)
        exit(0)

    # if more than 15 or so files then delete oldest
    # ==============================================
    log_flist = glob(normpath(log_dir + '/global_ecosse_*.log'))
    list_len = len(log_flist)
    max_log_files = 15
    num_to_delete = list_len - max_log_files
    if num_to_delete > 0:
        # date order list
        log_flist.sort(key=getmtime)
        for ifile in range(num_to_delete):
            try:
                remove(log_flist[ifile])
                form.lgr.info('removed log file: ' + log_flist[ifile])
            except (OSError, IOError) as e:
                print('Failed to delete log file: {0}'.format(e))

    # create dump files for grid point with mu_global 0
    form.fobjs = {}
    output_fnames = list(['nodata_muglobal_cells.csv'])
    if form.zeros_file:
        output_fnames.append('zero_muglobal_cells.csv')
    for file_name in output_fnames:
        long_fname = join(log_dir, file_name)
        key = file_name.split('_')[0]
        if exists(long_fname):
            try:
                remove(long_fname)
            except (PermissionError) as e:
                mess = 'Failed to delete mu global zeros dump file: {}\n\t{} '.format(long_fname, e)
                print(mess + '\n\t- check that there are no other instances of GlblEcosse'.format(long_fname, e))
                sleep(sleepTime)
                exit(0)

        form.fobjs[key] = open(long_fname,'w')

    return

def change_config_file(form, new_study = None):
    '''
       identify and read the new configuration file
       '''
    if new_study is None:
        new_study = form.combo00s.currentText()

    new_config = 'global_ecosse_config_site_specific_' + new_study
    config_file = form.config_dir + '/' + new_config + '.txt'

    if isfile(config_file):
        form.config_file = config_file
        read_config_file(form)
        form.study = new_study
        form.w_study.setText(new_study)
    else:
        print('Could not locate ' + config_file)

def build_and_display_studies(form):
    '''
    is called at start up and when user creates a new project
    '''

    config_files = glob(form.config_dir + '/' + GLBL_ECSSE_STR + '*.txt')
    studies = []
    for fname in config_files:
        dummy, remainder = fname.split(GLBL_ECSSE_STR)
        study, dummy = splitext(remainder)
        if study != '':
            studies.append(study)
    form.studies = studies

    # display studies list
    # ====================
    if hasattr(form, 'combo00s'):
        form.combo00s.clear()
        for study in studies:
            form.combo00s.addItem(study)

    return config_files

def _read_crop_codes(form, crop_pars_fname):
    '''
    read crop names and their corresponding codes
    '''
    required_crops = ['Winter Wheat', 'Spring Wheat', 'Grain Maize', 'Spring Oats', 'Spring Barley']

    with open(crop_pars_fname) as fhand:
        lines = fhand.readlines()

    crop_codes = {}
    for iline in range(0,len(lines), 9):
        crop_name = lines[iline].strip()
        code = int(lines[iline + 1].strip())
        if crop_name in required_crops:
            crop_codes[crop_name] = code

    form.crop_codes = crop_codes

    return

def _default_parms_settings ():

    _default_parms = {
          'site' : {
            'iom_c': 0.0,
            'toc': 2.7,
            'drain_class': 2,
            'depth_imperm_lyr': 3,
            'wtr_tbl_dpth': 300,
            'prev_lu': 1,
            'prev_crop_code': 1,
            'yield_prev_crop': 8.0,
            'prev_crop_harvest_doy': 0,
            'equil_mode': 2
          },
          'cultivation': {
            'cult_type': 3,
            'cult_vigor': 0.5
          }
        }
    return _default_parms

def _read_site_specific_parms():
    """
    read programme run settings from the parameters file, if it exists
    """
    func_name =  __prog__ +  '  _read_site_specific_parms'

    # look for setup file here...
    parms_setup_file = join(getcwd(), 'additional_setup', 'site_specific_parms.json')

    if exists(parms_setup_file):
            with open(parms_setup_file, 'r') as fsetup:
                try:
                    parms_settings = json.load(fsetup)
                except (json.decoder.JSONDecodeError, OSError, IOError) as err:
                    print(err)
    else:
        parms_settings = _default_parms_settings ()

    return parms_settings

def _read_setup_file(form):
    """
    # read settings used for programme from the setup file, if it exists,
    # or create setup file using default values if file does not exist
    """
    func_name =  __prog__ +  ' _read_setup_file'

    # validate setup file
    # ===================
    fname_setup = APPLIC_STR + '_setup.json'
    setup_file = join(getcwd(), fname_setup)

    if exists(setup_file):
        with open(setup_file, 'r') as fsetup:
            try:
                settings = json.load(fsetup)
            except (json.decoder.JSONDecodeError, OSError, IOError) as err:
                print(err)
                sleep(sleepTime)
                exit(0)
    else:
        settings = _write_default_setup_file(setup_file)
        print('Read setup file ' + setup_file)

    # validate setup file
    # ===================
    grp = 'setup'
    settings_list = ['config_dir', 'fname_png', 'hwsd_dir', 'log_dir', 'mask_fn', 'proj_path',
                                                'python_exe', 'runsites_py', 'shp_dir', 'sims_dir', 'weather_dir']
    for key in settings_list:
        if key not in settings[grp]:
            print('*** Error *** setting {} is required in setup file {} '.format(key, setup_file))
            sleep(sleepTime)
            exit(0)

    # initialise vars
    # ===============
    config_dir      = settings[grp]['config_dir']
    form.fname_png  = settings[grp]['fname_png']
    form.hwsd_dir   = settings[grp]['hwsd_dir']
    log_dir         = settings[grp]['log_dir']
    form.mask_fn    = settings[grp]['mask_fn']
    form.proj_path  = settings[grp]['proj_path']
    form.python_exe = settings[grp]['python_exe']
    runsites_py     = settings[grp]['runsites_py']
    form.shp_dir    = settings[grp]['shp_dir']
    form.sims_dir   = settings[grp]['sims_dir']
    weather_dir     = settings[grp]['weather_dir']

    # make sure directories exist for configuration and log files
    # ===========================================================
    if not lexists(log_dir):
        makedirs(log_dir)
    form.log_dir = log_dir

    if not lexists(config_dir):
        makedirs(config_dir)
    form.config_dir = config_dir

    # file to run Ecosse
    # ==================
    if isfile(runsites_py):
        print('Will use ECOSSE run script : ' + runsites_py)
    else:
        print('Error reading {}\tECOSSE run script {} must exist'.format(setup_file, runsites_py))
        sleep(sleepTime)
        exit(0)

    form.runsites_py = runsites_py

    # check that the runsites configuration file exists
    # =================================================
    runsites_config_file = join(config_dir, 'global_ecosse_site_spec_runsites_config.txt')

    mess = 'Run sites configuration file ' + runsites_config_file
    if exists(runsites_config_file):
        mess += ' exists'
        form.runsites_config_file = runsites_config_file
    else:
        mess += ' does not exist - cannot run ECOSSE'
        form.runsites_config_file = None
    print(mess)

    # weather is crucial
    # ==================
    if lexists(weather_dir):
        form.weather_dir = weather_dir
    else:
        print('Error reading {}\tClimate directory {} must exist'.format(setup_file, weather_dir))
        sleep(sleepTime)
        exit(0)

    # check weather data
    # ==================
    read_weather_dsets_detail(form)

    # TODO: most of these are not used
    # ================================
    grp = 'run_settings'
    try:
        check_space_every = settings[grp]['check_space_every']
        form.completed_max = settings[grp]['completed_max']
        form.kml_flag      = settings[grp]['kml_flag']
        form.soilTestFlag = settings[grp]['soil_test_flag']
        form.space_remaining_limit = settings[grp]['space_remaining_limit']
        form.zeros_file   = settings[grp]['zeros_file']
    except KeyError:
        print('{0}\tError in group: {1}'.format(func_name,grp))
        sleep(sleepTime)
        exit(0)
    from glbl_ecss_cmmn_cmpntsGUI import print_resource_locations

    lta_nc_fname = None
    print_resource_locations(setup_file, config_dir, form.hwsd_dir, weather_dir, lta_nc_fname, form.sims_dir, log_dir)

    return True

def _write_default_setup_file(setup_file):
    """
    #  stanza if setup_file needs to be created
    """
    # Windows only for now
    # =====================
    if os_name != 'nt':
        print('Operating system is ' + os_name + 'should be nt - cannot proceed with writing default setup file')
        sleep(sleepTime)
        exit(0)

    # return list of drives
    # =====================
    import win32api

    drives = win32api.GetLogicalDriveStrings()
    drives = drives.split('\000')[:-1]
    if 'S:\\' in drives:
        root_dir_app  = 'S:\\tools\\'     # Global Ecosse installed here
        root_dir_user = 'H:\\'            # user files reside here
    else:
        root_dir_app  = 'E:\\'
        root_dir_user = 'C:\\AbUniv\\'

    suite_path = root_dir_app + 'GlobalEcosseSuite\\'
    data_path  = root_dir_app + 'GlobalEcosseData\\'
    outputs_path  = root_dir_app + 'GlobalEcosseOutputs\\'
    root_dir_user += 'GlobalEcosseSuite\\'
    runsites_py = ''

    _default_setup = {
        'setup': {
            'root_dir_user' : root_dir_user,
            'fname_png'     : join(suite_path + 'Images', 'Tree_of_life.PNG'),
            'shp_dir'       : data_path + 'CountryShapefiles',
            'images_dir'    : outputs_path + 'images',
            'python_exe'   : 'E:\\Python36\\python.exe',
            'runsites_py'  : runsites_py,
            'sims_dir'      : outputs_path + 'EcosseSims',
            'log_dir'       : root_dir_user + 'logs',
            'config_dir'    : root_dir_user + 'config',
            'hwsd_csv_fname': '',
            'hwsd_dir'   : data_path + 'HWSD_NEW',
            'weather_dir': 'E:\\GlobalEcosseData'
        },
        'run_settings': {
            'completed_max'         : 5000000000,
            'check_space_every'     : 10,
            'space_remaining_limit' : 1270,
            'kml_flag'      : True,
            'soil_test_flag': False,
            'zeros_file'    : False
        }
    }
    # create setup file
    # =================
    with open(setup_file, 'w') as fsetup:
        json.dump(_default_setup, fsetup, indent=2, sort_keys=True)
        fsetup.close()
        return _default_setup

def _write_default_config_file(config_file):
    """
    #        ll_lon,    ll_lat  ur_lon,ur_lat
    # stanza if config_file needs to be created
    """
    _default_config = {
        'minGUI': {
            'bbox': bbox_default,
            'cordexFlag': False,
            'aveWthrFlag': False,
            'hwsdCsvFname': '',
            'fertiliserFname': '',
            'daily_mode': True
        },
        'cmnGUI': {
            'climScnr' : 0,
            'cropIndx' : 0,
            'cruStrtYr': 0,
            'cruEndYr' : 0,
            'eqilMode' : 9.5,
            'futStrtYr': 0,
            'futEndYr' : 0,
            'gridResol': 0,
            'study'    : ''
        }
    }
    # if config file does not exist then create it...
    with open(config_file, 'w') as fconfig:
        json.dump(_default_config, fconfig, indent=2, sort_keys=True)
        fconfig.close()
        return _default_config

def read_config_file(form):
    """
    # read widget settings used in the previous programme session from the config file, if it exists,
    # or create config file using default settings if config file does not exist
    """
    func_name =  __prog__ +  ' read_config_file'
    config_file = form.config_file
    if exists(config_file):
        try:
            with open(config_file, 'r') as fconfig:
                config = json.load(fconfig)
                fconfig.close()
                print('Read config file ' + config_file)
        except (OSError, IOError) as err:
                print(err)
                return False
    else:
        config = _write_default_config_file(config_file)
        print('Wrote configuration file ' + config_file)

    grp = 'minGUI'
    try:
        wthr_rsrce_indx = config[grp]['cordexFlag']
        ave_weather = config[grp]['aveWthrFlag']
        form.bbox = config[grp]['bbox']
        fertiliser_fname = config[grp]['fertiliserFname']
        sowing_dates_fname = config[grp]['sowingDatesFname']
        daily_mode  =  config[grp]['dailyMode']
        manure_flag  =  config[grp]['manureFlag']
        hwsd_csv_fname = config[grp]['hwsdCsvFname']
    except KeyError:
        print('{}\tError in group: {}'.format(func_name,grp))
        return False

    # patch to permit backwards compatibility
    if not isinstance(wthr_rsrce_indx, int):
        wthr_rsrce_indx = 0   # sets to CRU, the default

    # make sure index is within the permissable range of entries
    nitems = form.combo10w.count()
    if wthr_rsrce_indx >= 0 and wthr_rsrce_indx < nitems:
        form.combo10w.setCurrentIndex(wthr_rsrce_indx)

    # redundant - fertiliser and sowing dates file feedback
    # ======================================================
    descr = identify_datasets(form.proj_path, form.mask_fn)  # displays file info
    print(descr)

    # common area
    # ===========
    grp = 'cmnGUI'
    try:
        form.w_study.setText(str(config[grp]['study']))
        form.combo09s.setCurrentIndex(config[grp]['cruStrtYr'])
        form.combo09e.setCurrentIndex(config[grp]['cruEndYr'])
        form.combo10.setCurrentIndex(config[grp]['climScnr'])
        form.combo11s.setCurrentIndex(config[grp]['futStrtYr'])
        form.combo11e.setCurrentIndex(config[grp]['futEndYr'])
        form.combo12.setCurrentIndex(config[grp]['cropIndx'])
        form.combo16.setCurrentIndex(config[grp]['gridResol'])
        form.w_equimode.setText(str(config[grp]['eqilMode']))
    except KeyError:
        print('{}\tError in group: {}'.format(func_name,grp))
        return False

    # settings for single point 0peration and use of shape file

    # bounding box set up
    area = calculate_area(form.bbox)
    ll_lon, ll_lat, ur_lon, ur_lat = form.bbox
    form.w_ll_lon.setText(str(ll_lon))
    form.w_ll_lat.setText(str(ll_lat))
    form.w_ur_lon.setText(str(ur_lon))
    form.w_ur_lat.setText(str(ur_lat))
    form.lbl03.setText(format_bbox(form.bbox,area))

    # reset widgets associated with the HWSD file
    # ===========================================
    form.fstudy = ''
    if hwsd_csv_fname != '':
        form.w_lbl06.setText(hwsd_csv_fname)
        if isfile(hwsd_csv_fname):
            # read CSV file using pandas and create obj
            form.hwsd_mu_globals = HWSD_mu_globals_csv(form, hwsd_csv_fname)
            form.w_lbl07.setText(form.hwsd_mu_globals.aoi_label)
        else:
            print('HWSD csv file ' + hwsd_csv_fname + ' does not exist')
            hwsd_csv_fname = ''

    if hwsd_csv_fname == '':
        form.hwsd_mu_globals = None
        form.w_lbl06.setText('')
        form.w_lbl07.setText('')

    form.hwsd_csv_fname = hwsd_csv_fname

    # set check boxes
    # ===============
    if ave_weather:
        form.w_ave_wthr.setCheckState(2)
    else:
        form.w_ave_wthr.setCheckState(0)

    if daily_mode:
        form.w_daily.setChecked(True)
    else:
        form.w_mnthly.setChecked(True)

    if form.runsites_config_file == None:
        form.w_run_ecosse.setEnabled(False)
        form.w_process_files.setEnabled(False)

    # avoids errors when exiting
    # ==========================
    form.req_resol_deg = None
    form.req_resol_granul = None

    form.w_use_dom_soil.setChecked(True)
    form.w_use_high_cover.setChecked(True)
    form.w_use_mask.setChecked(True)

    DFLT_MAX_BAND  = 999

    start_at_band = 0
    stop_after_band = DFLT_MAX_BAND
    # stop_after_band = 8
    if start_at_band != 0 or stop_after_band < DFLT_MAX_BAND:
        print('***Warning*** will start simulations at band {} and stop after {}'.format(start_at_band, stop_after_band))
    form.start_at_band = start_at_band
    form.stop_after_band = stop_after_band

    return True

def write_config_file(form):
    """
    # write current selections to config file
    """
    study = form.w_study.text()
    if study == '':
        study = form.combo00s.currentText()

    # facilitate multiple config file choices
    config_file = join(form.config_dir, GLBL_ECSSE_STR + study + '.txt')

    # prepare the bounding box
    try:
        ll_lon = float(form.w_ll_lon.text())
        ll_lat = float(form.w_ll_lat.text())
        ur_lon = float(form.w_ur_lon.text())
        ur_lat = float(form.w_ur_lat.text())
    except ValueError:
        ll_lon = 0.0
        ll_lat = 0.0
        ur_lon = 0.0
        ur_lat = 0.0
    form.bbox =  list([ll_lon,ll_lat,ur_lon,ur_lat])

    # TODO:
    # print('Weather choice Id: {}'.format(form.w_weather_choice.checkedId()))
    config = {
        'minGUI': {
            'bbox'         : form.bbox,
            'cordexFlag'      : form.combo10w.currentIndex(),
            'aveWthrFlag'     : form.w_ave_wthr.isChecked(),
            'hwsdCsvFname'    : normpath(form.w_lbl06.text()),
            'dailyMode'       : form.w_daily.isChecked(),
            'manureFlag'      : True,
            'fertiliserFname' : None,
            'sowingDatesFname' : None
        },
        'cmnGUI': {
            'study'    : form.w_study.text(),
            'cruStrtYr': form.combo09s.currentIndex(),
            'cruEndYr' : form.combo09e.currentIndex(),
            'climScnr' : form.combo10.currentIndex(),
            'futStrtYr': form.combo11s.currentIndex(),
            'futEndYr' : form.combo11e.currentIndex(),
            'cropIndx' : form.combo12.currentIndex(),
            'eqilMode' : form.w_equimode.text(),
            'gridResol': form.combo16.currentIndex()
            }
        }
    if isfile(config_file):
        descriptor = 'Overwrote existing'
    else:
        descriptor = 'Wrote new'
    if study != '':
        with open(config_file, 'w') as fconfig:
            json.dump(config, fconfig, indent=2, sort_keys=True)
            fconfig.close()
            print('\n' + descriptor + ' configuration file ' + config_file)

def write_study_definition_file(form):
    """
    # write current selections to config file
    """
    # prepare the bounding box
    try:
        ll_lon = float(form.w_ll_lon.text())
        ll_lat = float(form.w_ll_lat.text())
        ur_lon = float(form.w_ur_lon.text())
        ur_lat = float(form.w_ur_lat.text())
    except ValueError:
        ll_lon = 0.0
        ll_lat = 0.0
        ur_lon = 0.0
        ur_lat = 0.0
    bbox =  list([ll_lon,ll_lat,ur_lon,ur_lat])
    study = form.w_study.text()

    wthr_rsrce = form.combo10w.currentText()
    if wthr_rsrce == 'CRU':
        fut_clim_scen = form.combo10.currentText()
    else:
        fut_clim_scen = wthr_rsrce

    # construct land_use change - not elegant but adequate
    # =========================
    land_use = ''
    if hasattr(form, 'lu_pi_content'):
        for indx in form.lu_pi_content['LandusePI']:
            lu, pi = form.lu_pi_content['LandusePI'][indx]
            land_use += form.lu_type_abbrevs[lu] + '2'
    else:
        land_use = 'unk2unk'

    land_use = land_use.rstrip('2')

    # TODO: replace "luCsvFname": "" with 'luPiJsonFname': form.fertiliser_fname
    if hasattr(form, 'combo12'):
        crop_name = form.combo12.currentText()
    else:
        crop_name = 'Unknown'

    study_defn = {
        'studyDefn': {
            'bbox'     : bbox,
            "luCsvFname": "",
            'hwsdCsvFname' : form.w_lbl06.text(),
            'study'    : study,
            'land_use' : land_use,
            'histStrtYr': form.combo09s.currentText(),
            'histEndYr' : form.combo09e.currentText(),
            'climScnr' : fut_clim_scen,
            'futStrtYr': form.combo11s.currentText(),
            'futEndYr' : form.combo11e.currentText(),
            'cropName' : crop_name,
            'province' : 'xxxx',
            'resolution' : form.req_resol_deg,
            'shpe_file': 'xxxx',
            'version'  : form.version
            }
        }

    # copy to sims area
    if study != '':
        study_defn_file = join(form.sims_dir, study + '_study_definition.txt')
        with open(study_defn_file, 'w') as fstudy:
            json.dump(study_defn, fstudy, indent=2, sort_keys=True)
            fstudy.close()
            print('\nWrote study definition file ' + study_defn_file)

    return

def check_sims_dir(form):
    '''
    called from _initiation func.
    '''
    retFlag = False
    sims_dir = form.sims_dir
    form.sims_dir = '' # default in case of failure

    # make sure directory has write permissions and it exists
    if not lexists(sims_dir):
        mkdir(sims_dir)
        form.sims_dir = sims_dir
        retFlag = True
    else:
        if isdir(sims_dir):
            form.lgr.info('Directory {0} already exists'.format(sims_dir))
            form.sims_dir = sims_dir
            retFlag = True
        elif isfile(sims_dir):
            print('{0} already exists but as a file' + sims_dir)
        else:
            print('{0} is not a directory or a file' + sims_dir)

    return retFlag

def write_runsites_config_file(form, identifer = None):

    func_name =  __prog__ +  ' write_runsites_config_file'

    # read the runsites config file and edit one line
    # ======================================
    runsites_config_file = form.runsites_config_file
    try:
        with open(runsites_config_file, 'r') as fconfig:
            config = json.load(fconfig)
            print('Read config file ' + runsites_config_file)
    except (OSError, IOError) as err:
            print(err)
            return False

    # overwrite config file
    # =====================
    study = form.w_study.text()
    if identifer is not None:
        study += '_' + identifer

    sims_dir = normpath(join(form.sims_dir, study))
    if not isdir(sims_dir):
        print('Simulations directory ' + sims_dir + ' must exist')
        return False

    config['Simulations']['simdir'] = sims_dir
    # config['Simulations']['resume_frm_prev'] = form.w_skip_sites.isChecked()
    # config['General']['cropName'] = 'unknown'
    config['General']['cropName'] = form.combo12.currentText()

    with open(runsites_config_file, 'w') as fconfig:
        json.dump(config, fconfig, indent=2, sort_keys=True)
        print('Edited ' + runsites_config_file + ' with simulation location: ' + sims_dir)
        return True
