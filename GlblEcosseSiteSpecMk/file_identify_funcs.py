#-------------------------------------------------------------------------------
# Name:        file_identify_funcs.py
# Purpose:
# Author:      Mike Martin
# Created:     25/10/2019
# Licence:     <your licence>
# defs:
#   _retrieve_abbreviated_data_range
#   _analyse_line
#   _infer_column_meanings  guesses xcol, ycol and zcol; work out number of header lines
#   identify_file           reads the first n (e.g. five) thousand lines and identifies file
#
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'file_identify_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
import os
from locale import LC_ALL, setlocale, format_string

setlocale(LC_ALL, '')

coord_max_val = 999.0
ecosse_large_number = 999999999.0 # TODO: very crude!
max_lines_to_read = 200

def _retrieve_abbreviated_data_ranges(line_recs, separator, num_header_recs, field_dscrptrs):
    '''
    identifies fields
    '''
    # examine first 10 columns of each line
    field_types = field_dscrptrs[:10]
    field_vals = {}
    for icol, field_type in enumerate(field_types):
        if field_type == 'Number':
            field_vals[icol] = []

    for line in line_recs[num_header_recs:]:

        rec = line.rstrip('\n').split(separator)
        for icol, field_type in enumerate(field_types):
            if field_type == 'Number':
                try:
                    field_vals[icol].append(float(rec[icol]))
                except ValueError as e:
                    field_vals[icol].append(-999)

    # step through each field to identify columns for lat/long
    field_extents = {}
    for icol in field_vals:
        field_extents[icol] = list([min(field_vals[icol]), max(field_vals[icol])])
        # print('{} {} {}'.format(icol, field_extents[icol][0], field_extents[icol][1]))

    return field_extents

def _analyse_line(line_rec, separator):

    field_dscrptrs = []

    num_nums = 0
    num_nums_first9 = 0
    num_strs = 0
    for ic, el in enumerate(line_rec.split(separator)):

        try:
            val = float(el)
            field_dscrptrs.append('Number')
            num_nums += 1
            if ic < 10:
                num_nums_first9 += 1
        except (ValueError) as e:
            field_dscrptrs.append('String')
            num_strs += 1

    # if the data record flag is false then it is probably header line
    if num_nums_first9 >= 3:
        data_record_flag = True
    else:
        data_record_flag = False

    return num_nums, num_strs, data_record_flag, field_dscrptrs

def _infer_column_meanings(line_recs, separator):
    '''
    a) make reasonable guesses for xcol, ycol and zcol
    b) work out number of header lines
    c) identify non-number columns
    '''
    num_lines = len(line_recs)

    # is first line a header?
    # =======================
    num_nums, num_strs, data_record_flag, field_dscrptrs = _analyse_line(line_recs[0], separator)

    # data_record_flag indicates at least 3 numbers are present
    # =========================================================
    if data_record_flag:
        num_header_recs = 0
    else:
        num_header_recs = 1

    # choose 4th line to make sure we are past the headers
    # ====================================================
    num_nums, num_strs, data_record_flag, field_dscrptrs = _analyse_line(line_recs[3], separator)
    num_fields = num_nums + num_strs

    field_extents = _retrieve_abbreviated_data_ranges(line_recs, separator, num_header_recs, field_dscrptrs)

    if data_record_flag:

        # first look for latitude, then longitude
        # =======================================
        for indx in field_extents:
            v_min, v_max = field_extents[indx]
            if (v_min > -90.0 and v_max < 90.0):
                ycol = indx
                del field_extents[indx]
                break

        for indx in field_extents:
            v_min, v_max = field_extents[indx]
            if (v_min > -180.0 and v_max < 180.0):
                xcol = indx
                del field_extents[indx]
                break

        for indx in field_extents:
            zcol = indx

    else:
        # file is not valid
        # =================
        xcol = None; ycol = None; zcol = None
        num_header_recs = num_lines

    return xcol, ycol, zcol, num_fields, num_header_recs, field_dscrptrs

def identify_csv_file(csv_file, max_lines = max_lines_to_read):
    '''
    reads the specified number of lines and identifes number of columns, their ranges and types
    '''

    if not os.path.isfile(csv_file):
        print('File {} does not exist'.format(csv_file))
        return None

    # read first n (eg 200) lines to estimate number of lines and fields
    # ==================================================================
    fobj = open(csv_file, 'r')
    file_size_read = 0
    line_recs = []
    for nline in range(0, max_lines):
        raw_line = fobj.readline()

        # empty string indicates an end of line
        # =====================================
        if raw_line == '':
            break
        file_size_read += len(raw_line)
        line = raw_line.rstrip('\n')
        line_recs.append(line)
    fobj.close()

    if file_size_read == 0 or nline == 0:
        print('File {} is empty'.format(csv_file))
        return None

    ave_line_len = file_size_read/nline

    # identify which type of file
    # ===========================
    num_spaces = line.count(' ')
    num_tabs = line.count('\t')
    num_commas = line.count(',')

    if num_tabs == 0 and num_commas == 0 and num_spaces == 0:
        print('No tabs or commas or spaces on line {} of file {} - file not identified'.format(nline, csv_file))
        return None

    if num_spaces > num_tabs and num_spaces > num_commas:
        separator = ' '
        separator_name = 'space'
    else:
        if num_tabs > num_commas:

            # tabs dominate - probably Ecosse output
            # ======================================
            separator = '\t'
            separator_name = 'tab'
        else:
            # commas dominate - probably CSV of HWSD grid cells
            # =================================================
            separator = ','
            separator_name = 'comma'

    xcol, ycol, zcol, num_fields, num_header_recs, field_dscrptrs = _infer_column_meanings(line_recs, separator)

    # construct descriptor
    # ====================
    fname_info = os.stat(csv_file)
    setlocale(LC_ALL, '')
    nlines_est = format_string('%d', fname_info.st_size/ave_line_len, grouping=True)
    descriptor = '   estimated # lines: {}'.format(nlines_est)
    descriptor += '   # columns: {}'.format(str(num_fields))

    # construct extended file description
    # ===================================
    mess = 'File: ' + os.path.split(csv_file)[1] + '\tseparator: ' + separator_name
    mess += '\txcol: {}\tycol: {}\tzcol: {}\t# fields: {}'.format(xcol, ycol, zcol, num_fields)
    mess += '\t# num_header_recs: {}\t# tabs: {}\t# commas: {}'.format(num_header_recs, num_tabs, num_commas)
    for field_dscrptr in field_dscrptrs:
        mess += '\n\t' + field_dscrptr

    return list([nlines_est, num_fields, separator_name, descriptor])
