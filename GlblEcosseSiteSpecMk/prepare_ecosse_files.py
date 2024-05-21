'''
#-------------------------------------------------------------------------------
# Name:
# Purpose:
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
'''
__version__ = '1.0.00'
__prog__ = 'prepare_ecosse_files.py'

import csv
from os import makedirs, fsync, remove
from os.path import lexists, join, normpath, basename, exists, split as split_dir
import time
import sys
from calendar import isleap
from shutil import copyfile
import json
from validate import check_harv_doys_adjust
from glob import glob

from thornthwaite import thornthwaite
from make_site_specific_files import MakeSiteFiles
from crop_soil_fert_funcs import site_data_modify_daily, site_data_modify_mnthly, create_site_soil_layers
from mngmnt_fns_and_class import fetch_vals_from_dset

_MONTHDAYS = [31,28,31,30,31,30,31,31,30,31,30,31]
_LEAP_MONTHDAYS = [31,29,31,30,31,30,31,31,30,31,30,31]
SPACER_LEN = 12

def make_ecosse_files(form, climgen, kml_output, site_rec, province, fert_defn, pheno_defn, yield_defn,
                                                                        pettmp_hist_grid_cell, pettmp_fut_grid_cell):
    '''
    generate sets of Ecosse files for each site
    where each site has one or more soils and each soil can have one or more dominant soils
    pettmp_grid_cell is climate data for this soil grid point
    '''
    func_name = 'make_ecosse_files'

    # simulation start and end year
    # =============================
    start_year = climgen.fut_start_year
    end_year = climgen.fut_end_year
    if climgen.start_from_1901:
        ntot_years = end_year - 1901 + 1
    else:
        ntot_years = end_year - start_year + 1

    crop_name = form.combo12.currentText()
    crop_code = form.crop_codes[crop_name]
    wthr_rsrce = climgen.wthr_rsrce
    if wthr_rsrce not in list(['NASA','EWEMBI','EObs','HARMONIE']):
        print('Weather resource {} not recognised in function {}'.format(wthr_rsrce, func_name))
        return

    if not isinstance(pettmp_hist_grid_cell, dict) or not isinstance(pettmp_hist_grid_cell, dict):
        print('Bad input to ' + func_name)
        return

    gran_lat, gran_lon, latitude, longitude, area, mu_globals_props = site_rec
    sims_dir = form.sims_dir
    fut_clim_scen = climgen.fut_clim_scen

    # get sowing dates nearest to this lat/lon
    # ========================================
    fert_amnts = fetch_vals_from_dset(fert_defn, latitude, longitude, ntot_years)
    yields = fetch_vals_from_dset(yield_defn, latitude, longitude, ntot_years)     # NB yield is Python inbuilt function
    sow_days = fetch_vals_from_dset(pheno_defn, latitude, longitude, ntot_years, 'sowing')
    harv_days = fetch_vals_from_dset(pheno_defn, latitude, longitude, ntot_years, 'harvest')
    harv_days = check_harv_doys_adjust(sow_days, harv_days)

    # calculate historic average weather
    # ==================================
    dset_start_year = form.weather_sets[wthr_rsrce + '_Mnth']['start_year']

    hist_start_year = climgen.hist_start_year
    indx_start = 12*(hist_start_year - dset_start_year)

    hist_end_year = climgen.hist_end_year
    indx_end   = 12*(hist_end_year - dset_start_year + 1) # end year includes all 12 months - TODO: check

    # use dict-comprehension to initialise precip. and temperature dictionaries
    # =========================================================================
    hist_precip = {mnth: 0.0 for mnth in climgen.months}
    hist_tmean  = {mnth: 0.0 for mnth in climgen.months}

    for indx in range(indx_start, indx_end, 12):

        for imnth, month in enumerate(climgen.months):
            try:
                hist_precip[month] += pettmp_hist_grid_cell['precipitation'][indx + imnth]
                hist_tmean[month]  += pettmp_hist_grid_cell['temperature'][indx + imnth]
            except (IndexError, TypeError) as err:
                print('\n' + str(err) + 'indx: {}\timnth: {}\tat lat/lon: {} {}'.format(indx, imnth, latitude, longitude))
                return

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    hist_weather_recs = []
    num_hist_years = hist_end_year - hist_start_year + 1
    hist_lta_precip = []
    for month in climgen.months:
        ave_precip = hist_precip[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_precip,1)), \
                                            '{} long term average monthly precipitation [mm]'.format(month)))
        hist_lta_precip.append(ave_precip)

    hist_lta_tmean = []
    for month in climgen.months:
        ave_tmean = hist_tmean[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_tmean,2)), \
                                            '{} long term average monthly temperature [degC]'.format(month)))
        hist_lta_tmean.append(ave_tmean)

    # write a single set of met files for all simulations for this grid cell
    # ======================================================================
    study = form.study
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    met_rel_path = '..\\..\\' + climgen.rgn_wthr_dir + '\\' + gran_coord + '\\'
    clim_dir = normpath(join(sims_dir, climgen.rgn_wthr_dir, gran_coord))
    met_fnames = make_met_files(clim_dir, latitude, climgen, pettmp_fut_grid_cell) # future weather
    nyears = len(met_fnames)

    site = MakeSiteFiles(form, climgen, comments = True)
    # TODO: these need sorted
    site.equil_mode = form.w_equimode.text()
    site.soil_code = 1
    site.prev_crop_harvest_doy = 0

    site.nyears = nyears
    site.ncrops = nyears # one crop per year
    site.ncultivations = nyears # one cultivation per year

    # create additional weather related files from already existing met files
    irc = climgen.create_FutureAverages(clim_dir, latitude, site, hist_lta_precip, hist_lta_tmean)
    site.start_year = start_year
    site.end_year = end_year

    #------------------------------------------------------------------
    # Create a set of simulation input files for each dominant
    # soil-land use type combination
    #------------------------------------------------------------------
    # construct directory name with all dominant soils
    manure_flag = False
    for pair in mu_globals_props.items():
        mu_global, proportion = pair
        area_for_soil = area*proportion
        soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

        for soil_num, soil in enumerate(soil_list):
            identifer = 'lat{0:0=7d}_lon{1:0=7d}_mu{2:0=5d}_s{3:0=2d}'.format(gran_lat, gran_lon, mu_global, soil_num + 1)

            sim_dir = join(sims_dir, study, identifer)

            if not lexists(sim_dir):
                makedirs(sim_dir)

            create_site_soil_layers(site, soil, form.depths)
            if site.timestep == 1:
                site = site_data_modify_daily(site, latitude, longitude, crop_code, crop_name, manure_flag,
                                                            fert_amnts, sow_days, harv_days, yields, met_fnames)
            else:
                site = site_data_modify_mnthly(site, latitude, longitude, crop_code, crop_name, manure_flag,
                                                                                        fert_amnts, yields, met_fnames)

            site.write_sim_files(sim_dir, soil, latitude, hist_weather_recs, met_rel_path)

            # write kml file if requested and signature file
            # ==============================================
            if form.kml_flag and soil_num == 0:
                write_kml_file(form.lgr, kml_output, sim_dir,  str(mu_global), mu_global, latitude, longitude)

            write_signature_file(form.lgr, sim_dir, mu_global, soil, latitude, longitude, province)

            # copy across Model_Switches.dat file
            # ===================================
            out_model_switches = join(sim_dir, basename(form.dflt_model_switches))
            copyfile(form.dflt_model_switches, out_model_switches)
            out_fnames_dat = join(sim_dir, basename(form.dflt_fnames_dat))
            copyfile(form.dflt_fnames_dat, out_fnames_dat)

        # manifest file is essential for subsequent processing
        # ====================================================
        write_manifest_file(form, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area_for_soil)

    # end of Soil loop
    # ================

    return

def make_met_files(clim_dir, latitude, climgen, pettmp_fut_grid_cell, skip_fut_wthr_flag = True):
    '''
    feed annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    '''
    start_year = climgen.fut_start_year
    end_year = climgen.fut_end_year

    # do not repeat met file creation if already in existence
    # =======================================================
    if lexists(clim_dir):
        if skip_fut_wthr_flag:
            met_fnames = [split_dir(met_fname)[1] for met_fname in glob(clim_dir + '/met*s.txt')]
            nyears = end_year - start_year + 1
            if len(met_fnames) >= nyears:
                return met_fnames
    else:
        makedirs(clim_dir)

    precip = pettmp_fut_grid_cell['precipitation']
    temp   = pettmp_fut_grid_cell['temperature']

    indx1 = 0
    met_fnames = []
    for year in range(start_year, end_year + 1):
        fname = 'met{}s.txt'.format(year)
        met_fnames.append(fname)
        met_path = join(clim_dir, fname)

        if climgen.fut_wthr_mnthly_flag:
            ntime_incrs = 12
        else:
            if isleap(year):
                ntime_incrs = 366
            else:
                ntime_incrs = 365

        indx2 = indx1 + ntime_incrs

        # precipitation and temperature
        precipitation = precip[indx1:indx2]
        temp_mean     = temp[indx1:indx2]

        # pet
        # ===
        if climgen.fut_wthr_mnthly_flag:
            pet = thornthwaite(temp_mean, latitude, year)
        else:
            pet = _get_thornthwaite(temp_mean, latitude, year)

        # TODO: do something about occasional runtime warning...
        pot_evapotrans = [round(p, 2) for p in pet]
        precip_out     = [round(p, 2) for p in precipitation]
        tmean_out      = [round(t, 2) for t in temp_mean]

        # write file
        output = []
        for tstep, mean_temp in enumerate(tmean_out):
            output.append([tstep+1, precip_out[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_path, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        indx1 += ntime_incrs

    return  met_fnames

def _get_thornthwaite(temp_mean, latitude, year):
    '''
    feed daily annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    '''
    func_name =  __prog__ + ' _get_thornthwaite'

    ntime_steps = len(temp_mean)

    if ntime_steps == 365:
        month_days = _MONTHDAYS
    else:
        month_days =  _LEAP_MONTHDAYS

    indx1 = 0
    monthly_t = []
    for ndays in month_days:
        indx2 = indx1 + ndays
        accum_temp = 0.0
        for temp in temp_mean[indx1:indx2]:
            accum_temp +=temp
        monthly_t.append(accum_temp/ndays)
        indx1 = indx2

    pet_mnthly = thornthwaite(monthly_t, latitude, year)

    # now inflate pet to daily
    # ========================
    pet_daily = []
    for pet_month, ndays in zip(pet_mnthly, month_days):
        pet_day = pet_month/ndays
        for ic in range(ndays):
            pet_daily.append(pet_day)

    # print('{} {} {} {}'.format(year, ntime_steps, len(pet_daily),func_name))
    return pet_daily

def input_txt_line_layout(data, comment):

        spacer_len = max(SPACER_LEN - len(data), 2)
        spacer = ' ' * spacer_len
        return '{}{}# {}\n'.format(data, spacer, comment)

def write_line_summary(form, coord_frst, coord_last, max_cells_in_line, max_cells_in_cluster):
    '''
    write line summary; function not used
    '''
    gran_lat_last, latitude_last, gran_lon_last, longitude_last = coord_last
    gran_lat_frst, latitude_frst, gran_lon_frst, longitude_frst = coord_frst
    form.fstudy[1].write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.
                      format(gran_lat_frst, gran_lon_frst, round(latitude_frst,6), round(longitude_frst,6),
                                gran_lat_last, gran_lon_last, round(latitude_last,6), round(longitude_last,6),
                                                    max_cells_in_line, max_cells_in_cluster))
    form.fstudy[1].flush()
    fsync(form.fstudy[1].fileno())
    return

def write_study_manifest_files(form, lon_lats):
    '''
    create study manifest file <study>_manifest.csv and summary <study>_summary_manifest.csv
    write first line and leave file object open for subsequent records
    '''
    ngranularity = 120
    NoData = -999

    # set up mainfest file for the study
    form.fstudy = []
    for frag_name in list(['_','_summary_']):
        study_fname = join(form.sims_dir, form.study + frag_name + 'manifest.csv')
        if exists(study_fname):
            remove(study_fname)
        form.fstudy.append(open(study_fname,'w',100))

    # write first two lines so they are distinguisable for sorting purposes
    for lon_lat_pair in lon_lats:
        longitude, latitude = lon_lat_pair
        gran_lon = round((180.0 + longitude)*ngranularity)
        gran_lat = round((90.0 - latitude)*ngranularity)
        form.fstudy[1].write('{}\t{}\t{}\t{}\t'.format(gran_lat, gran_lon, latitude, longitude))
        form.fstudy[0].write('{}\t{}\t{}\t{}\t{}\n'.format(gran_lat, gran_lon, latitude, longitude, NoData))

    form.fstudy[1].write('{}\t{}\n'.format(form.req_resol_deg, form.req_resol_granul))

    return

def write_manifest_file(form, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area):
    '''
    write json consisting of mu_global and soil shares
    '''

    # location etc.
    # =============
    manifest = {
        'location': {
            'longitude' : round(longitude,6),
            'latitude'  : round(latitude,6),
            'area' : round(area,8),
            'area_description' : form.study,
            'province' : 'province',
            'scenario' : fut_clim_scen
        }
    }
    # soil shares
    # ===========
    smu_global = str(mu_global)
    manifest[smu_global] =  {}
    for soil_num, soil in enumerate(soil_list):
        manifest[smu_global][soil_num + 1] =  soil[-1]

    # deprecated
    # ==========
    manifest['longitudes'] =  {}
    manifest['granular_longs'] = {}

    # construct file name and write
    # =============================
    manif_dir, fname_part2 = split_dir(sim_dir)
    manifest_fname = join(manif_dir, 'manifest_' + fname_part2[:-4] + '.txt')
    with open(manifest_fname, 'w') as fmanif:
        json.dump(manifest, fmanif, indent=2, sort_keys=True)
        fmanif.close()

    return

def write_signature_file(lgr, sim_dir, mu_global, soil, latitude, longitude, province = '', bad_val = 0):
    '''
    write json consisting of mu_global and soil details
    '''

    config = {
        'location': {
            'province' : province,
            'longitude' : round(longitude,6),
            'latitude'  : round(latitude,6),
            'share' : soil[-1]
        },
        'soil_lyr1': {
            'C_content': soil[0],
            'Bulk_dens': soil[1],
            'pH'       : soil[2],
            '%_clay': soil[3],
            '%_silt': soil[4],
            '%_sand': soil[5]
        }
    }
    # add subsoil layer, if it exists
    if len(soil) >= 12:
       config['soil_lyr2'] =  {
                'C_content': soil[6],
                'Bulk_dens': soil[7],
                'pH'       : soil[8],
                '%_clay': soil[9],
                '%_silt': soil[10],
                '%_sand': soil[11]
            }

    signature_fname = join(sim_dir, str(mu_global) + '.txt')
    with open(signature_fname, 'w') as fsig:
        json.dump(config, fsig, indent=2, sort_keys=True)
        fsig.close()

    return

def write_kml_file(lgr, output, sim_dir, fname_short, mu_global, latitude, longitude):
    '''
    write kml consisting of mu_global and soil details
    '''
    # set Icon size
    scaleVal = 0.5

    output[0] = '<?xml version="1.0" encoding="utf-8"?>\n'
    output[1] = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    output[2] = '\t<Document>\n'
    output[3] = '\t<Style id="defaultStyle">\n'
    output[4] = '\t\t<LabelStyle><scale>' + str(scaleVal) + '</scale></LabelStyle>\n'
    output[5] = '\t\t<IconStyle id="hoverIcon"><scale>' + str(scaleVal) + '</scale></IconStyle>\n'
    output[6] = '\t</Style>\n'

    output[7] = '\t\t<Placemark>\n'
    output[8] = '\t\t\t<name>' + str(mu_global) + '</name>\n'
    output[9] = '\t\t\t<description>\n' + 'Long: ' + str(round(longitude,4)) + \
            ',   Lat:' + str(round(latitude,4)) +'\nMu_global:' + str(mu_global) + '\n\t\t\t</description>\n'

    output[10] = '\t\t\t<Point><coordinates>' + str(longitude) + ',' + str(latitude) + '</coordinates></Point>\n'

    output[11] = '\t\t\t<styleUrl>#defaultStyle</styleUrl>\n'
    output[12] = '\t\t</Placemark>\n'
    output[13] = '\t</Document>\n'
    output[14] = '</kml>\n'

    kml_fname = join(sim_dir, fname_short + '.kml')
    fpout = open(kml_fname, 'w')
    fpout.writelines(output)
    fpout.close()

    # tmp_list = kml_fname.split('\\')
    # lgr.info('Wrote kml file {} {}'.format(tmp_list[-2],tmp_list[-1]))

    return

def update_progress(last_time, start_time, completed, est_num_sims, no_climate, no_cropland):
    '''
    Update progress bar             num_meta_cells, no_climate, no_cropland)
    '''
    new_time = time.time()
    if new_time - last_time > 5:
        # used to be: Estimated remaining
        mess = '\rCompleted: {:<5d}\tRemaining: {:<5d}\tNo climate: {:<4d}\tNot cropland: {:<5d}'\
                                    .format(completed, est_num_sims - completed, no_climate, no_cropland)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time