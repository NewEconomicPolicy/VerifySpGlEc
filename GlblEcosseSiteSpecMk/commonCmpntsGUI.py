"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'commonCmpntsGUI.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import os

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QLabel, QComboBox, QRadioButton, QButtonGroup, QPushButton, QCheckBox

from shape_funcs import calculate_area, format_bbox
from initialise_funcs import write_study_definition_file, read_config_file, write_config_file

STD_FLD_SIZE = 60
STD_BTN_SIZE = 100
RESOLUTIONS = {120:'30"', 30:'2\'', 20:'3\'', 10:'6\'', 6:'10\'', 4:'15\'', 3:'20\'', 2:'30\''}
reverse_resols = {}

for key in RESOLUTIONS:
    reverse_resols[RESOLUTIONS[key]] = key

def calculate_grid_cell(form, granularity = 120):

    latitude = 52.0
    # use current lower left latitude for reference
    latText = form.w_ll_lat.text()
    try:
        latitude = float(latText)
    except ValueError:
        print(latText)

    resol = form.combo16.currentText()
    try:
        granul = reverse_resols[resol]
    except KeyError as e:
        print(str(e))
        return

    resol_deg = 1.0/float(granul)  # units of decimal degrees
    bbox = list([0.0, latitude, resol_deg, latitude + resol_deg])
    area = calculate_area(bbox)
    form.lbl16a.setText(format_bbox(bbox,area,2))

    form.req_resol_upscale = int(granularity/granul)    # number of base granuals making up one side of a cell
    form.req_resol_granul = granul                            # number of cells per degree
    form.req_resol_deg = resol_deg                            # size of each trapizoid cell in decimal degrees

    return

def commonSection(form, grid, irow):
    
    # TODO: some of these could be transferred to the settings file

    # start_year = form.weather_sets['EObs_Mnth']['start_year']
    # end_year   = form.weather_sets['EObs_Mnth']['end_year']
    start_year = form.weather_sets['HARMONIE_Day']['start_year']
    end_year   = form.weather_sets['HARMONIE_Day']['end_year']
    hist_syears = list(range(start_year, end_year))
    hist_eyears = list(range(start_year + 1, end_year + 1))
    fut_syears = range(start_year, end_year)
    fut_eyears = list(range(start_year + 1, end_year + 1))
    scenarios = list(['A1B_MG1','A2_MG1','B1_MG1','B2_MG1'])  #
    form.depths = list([30,100]) # soil depths

    luTypes = {}; lu_type_abbrevs = {}
    for lu_type, abbrev, ilu in zip(
                    ['Arable','Forestry','Miscanthus','Grassland','Semi-natural', 'SRC', 'Rapeseed', 'Sugar cane'],
                    ['ara',   'for',      'mis',      'gra',      'nat',          'src', 'rps',      'sgc'],
                    [1,        3,          5,          2,          4,              6,     7,          7]):
        luTypes[lu_type] = ilu
        lu_type_abbrevs[lu_type] = abbrev

    form.land_use_types = luTypes
    form.lu_type_abbrevs = lu_type_abbrevs

    # line 9: resources
    # ===================
    irow += 1
    lbl10w = QLabel('Weather resource')
    lbl10w.setAlignment(Qt.AlignRight)
    helpText = 'permissable weather dataset resources are limited to EObs only'
    lbl10w.setToolTip(helpText)
    grid.addWidget(lbl10w, irow, 0)

    combo10w = QComboBox()
    for wthr_rsrce in form.valid_wthr_dset_rsrces:
        combo10w.addItem(wthr_rsrce)
    form.combo10w = combo10w
    grid.addWidget(combo10w, irow, 1)

    # line 9: scenarios
    # =================
    lbl10 = QLabel('Climate Scenario')
    lbl10.setAlignment(Qt.AlignRight)
    helpText = 'Ecosse requires future average monthly precipitation and temperature derived from climate models.\n' \
        + 'The data used here is ClimGen v1.02 created on 16.10.08 developed by the Climatic Research Unit\n' \
        + ' and the Tyndall Centre. See: http://www.cru.uea.ac.uk/~timo/climgen/'

    lbl10.setToolTip(helpText)
    grid.addWidget(lbl10, irow, 2)

    combo10 = QComboBox()
    for scen in scenarios:
        combo10.addItem(str(scen))
    combo10.setEnabled(False)
    combo10.setFixedWidth(STD_FLD_SIZE)
    form.combo10 = combo10
    grid.addWidget(combo10, irow, 3)

    # historic
    # ========
    irow += 1
    lbl09s = QLabel('Historic start year')
    lbl09s.setAlignment(Qt.AlignRight)
    helpText = 'Ecosse requires long term average monthly precipitation and temperature\n' \
            + 'which is derived from datasets managed by Climatic Research Unit (CRU).\n' \
            + ' See: http://www.cru.uea.ac.uk/about-cru'
    lbl09s.setToolTip(helpText)

    combo09s = QComboBox()
    for year in hist_syears:
        combo09s.addItem(str(year))
    combo09s.setFixedWidth(STD_FLD_SIZE)
    form.combo09s = combo09s

    grid.addWidget(lbl09s, irow,  0)
    # row, column, rowSpan, columnSpan
    grid.addWidget(combo09s, irow,  1)

    lbl09e = QLabel('End year')
    lbl09e.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl09e, irow, 2)

    combo09e = QComboBox()
    for year in hist_eyears:
        combo09e.addItem(str(year))
    combo09e.setFixedWidth(STD_FLD_SIZE)
    form.combo09e = combo09e
    grid.addWidget(combo09e, irow,  3)

    # years into future
    # =================
    irow += 1
    lbl11s = QLabel('Simulation start year')
    helpText = 'Simulation start and end years determine the number of growing seasons to simulate\n' \
            + 'CRU and CORDEX resources run to 2100 whereas EObs resource runs to 2017'
    lbl11s.setToolTip(helpText)
    lbl11s.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl11s, irow, 0)

    combo11s = QComboBox()
    for year in fut_syears:
        combo11s.addItem(str(year))
    combo11s.setFixedWidth(STD_FLD_SIZE)
    form.combo11s = combo11s
    grid.addWidget(combo11s, irow,  1)

    lbl11e = QLabel('End year')
    lbl11e.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl11e, irow, 2)

    combo11e = QComboBox()
    for year in fut_eyears:
        combo11e.addItem(str(year))
    combo11e.setFixedWidth(STD_FLD_SIZE)
    form.combo11e = combo11e
    grid.addWidget(combo11e, irow,  3)
    
    w_ave_wthr = QCheckBox('Start simulation from 1901')
    helpText = 'Select this option to use average weather, from the CRU year range, for\n' \
               ' the climate file for each of the simulation years'
    w_ave_wthr.setToolTip(helpText)
    grid.addWidget(w_ave_wthr, irow,  4, 1, 2)
    form.w_ave_wthr = w_ave_wthr

    # crop choice
    # ===========
    irow += 1
    w_lbl12 = QLabel('Crop name')
    w_lbl12.setAlignment(Qt.AlignRight)
    helpText = 'list of crop names from CROP_SUN.DAT file. This name is mapped to \n' \
                                                            'its corresponding crop code in the management file'
    w_lbl12.setToolTip(helpText)
    grid.addWidget(w_lbl12, irow, 0)

    combo12 = QComboBox()
    for crop_name in form.crop_codes:
        combo12.addItem(crop_name)
    form.combo12 = combo12
    combo12.setToolTip(helpText)
    grid.addWidget(combo12, irow, 1)

    # spacer
    # ======
    irow += 3
    grid.addWidget(QLabel(), irow, 0)

    return irow

def save_clicked(form):

        # write last GUI selections
        write_config_file(form)
        calculate_grid_cell(form)
        write_study_definition_file(form)

def grid_coarseness(form, grid, irow):
    '''
   function to lay out grid resolution dropdown and reporting
    '''

    # row 16 - grid coarseness
    # ========================
    irow += 1
    lbl16 = QLabel('Grid resolution')
    lbl16.setAlignment(Qt.AlignRight)
    helpText = 'The size of each grid cell is described in arc minutes and arc seconds. The smallest cell resolution \n' \
        + 'corresponds to that of the HWSD database (30 arc seconds) and the largest to that used by the climate data ' \
        + '(30 arc minutes)'
    lbl16.setToolTip(helpText)
    grid.addWidget(lbl16, irow,  0)

    combo16 = QComboBox()
    for resol in sorted(RESOLUTIONS, reverse = True):
        combo16.addItem(str(RESOLUTIONS[resol]))
    combo16.setToolTip(helpText)
    combo16.setFixedWidth(STD_FLD_SIZE)
    combo16.currentIndexChanged[str].connect(form.resolutionChanged)
    grid.addWidget(combo16, irow,  1)
    form.combo16 = combo16

    form.lbl16a = QLabel('')
    form.lbl16a.setToolTip(helpText)
    grid.addWidget(form.lbl16a, irow,  2, 1, 2)

    return irow

def exit_clicked(form, write_config_flag = True):

    # write last GUI selections
    if write_config_flag:
        write_config_file(form)
        calculate_grid_cell(form)
        write_study_definition_file(form)

    # close various files
    if hasattr(form, 'fobjs'):
        for key in form.fobjs:
            form.fobjs[key].close()

    # close logging
    try:
        form.lgr.handlers[0].close()
        # print('Log file written to: ' + form.log_dir)
    except AttributeError:
        pass

    form.close()
