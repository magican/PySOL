#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################################
#Author: Denis Spiridonov
#Copyright 2012, SOLab, http://solab.rshu.ru/. All rights reserved.
#
#License:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License (version 3) 
# as published by the Free Software Foundation. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details:
#   http://www.gnu.org/licenses/gpl-3.0.html
########################################################

# slbtiles --hirlam-dir temp/ -o /mnt/d/output/ -z 3-9
# slbtiles -i temp/ -o /mnt/d/output/ -z 3-9
# slbtiles -i temp/ -o /mnt/d/output/ --move-red --red1 7-9 --red2 5-6 --red4 3-4

import subprocess, sys, os
import xml.dom.minidom
from PIL import Image
import ConfigParser, argparse
import datetime
import logging, logging.handlers
import shutil
 
def copytree(src, dst, symlinks=0):
    print ("copy tree " + src)
    names = os.listdir(src)
    if not os.path.exists(dst):
        os.mkdir(dst)
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks)
            else:
                shutil.copy2(srcname, dstname)
        except (IOError, os.error) as why:
            print "Can't copy %s to %s: %s" %(srcname, dstname, str(why))

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-i', type=str, dest='input', help="input directory")
    parser.add_argument(
            '--image', type=str, dest='input_image', help="input image file")
    parser.add_argument(
            '--kml', type=str, dest='input_kml', help="input kml file")
    parser.add_argument(
            '--hirlam-dir', type=str, dest='red_dir',
            help="input hirlam product directory")
    parser.add_argument(
             '-o', type=str, dest='output', help='output directory')
    parser.add_argument(
             '-z', '--zoom', type=str, dest='zoom', default = '3-10',
             help="zoom levels to render (format:'2-5' or '10')")
    parser.add_argument(
             '--log', type=str, dest='log_file', help="path to log file")
    parser.add_argument(
             '--notiles', action='store_true', dest='notiles',
             default = False, help="will be not created tiles")
    parser.add_argument(
             '-t_srs', type=str, dest='t_srs',
             help="target spatial reference set")

    parser.add_argument(
             '--move-red', action='store_true', dest='move_red',
             default = False, help="selective movement tiles from " +\
             "red1, red2 and red4 directories to output directory")
    parser.add_argument(
             '--red1', type=str, dest='red1', help="range for red1")
    parser.add_argument(
             '--red2', type=str, dest='red2', help="range for red2")
    parser.add_argument(
             '--red4', type=str, dest='red4', help="range for red4")
    return parser

def create_log(args):
    file_logger = logging.getLogger('slbtilesLogger')
    file_logger.setLevel(logging.DEBUG)
    
    if not args.log_file:
        log_filename = os.path.join (os.path.expanduser('~'), 'slbtiles.log')
    else:
        log_filename = args.log_file
    
    # Add the log message handler to the logger
    handler = logging.handlers.RotatingFileHandler(
                        log_filename, maxBytes=1000000, backupCount=1)
    
    file_logger.addHandler(handler)
    return file_logger

class SLBTiles():
    def __init__(self, args, kml_file, png_file, file_logger):
        self.input_dir = args.input
        self.output_dir = args.output
        self.zoom = args.zoom
        self.kml_file = kml_file
        self.png_file = png_file
        self.file_logger = file_logger
        self.t_srs = args.t_srs if args.t_srs else None

        now_time = datetime.datetime.now().isoformat()
        start_time = now_time[:10] +' '+ now_time[11:19]
        str_to_log = 'slbtiles is running %s\n'%start_time +\
                     'Command-line options:\n' +\
                     'input directory: %s\n'%self.input_dir +\
                     'output directory: %s\n'%self.output_dir +\
                     'notiles: %s\n'%str(args.notiles)
        if not args.notiles:
            str_to_log += 'zoom levels: %s\n'%self.zoom
        self.file_logger.debug(str_to_log)

    def getValue(self, nodelist):
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return ''.join(rc)

    def gdal_translate(self):
        conf = ''
        for line in open(self.kml_file, 'r').readlines():
            conf += line.replace('gx:','')

        self.kml_file_temp = os.path.join(self.output_dir,
                                     os.path.basename(self.kml_file+'_temp'))
        fd = open(self.kml_file_temp, 'w')
        fd.write(conf)
        fd.close()

        self.dom = xml.dom.minidom.parseString(\
                                        open(self.kml_file_temp, 'r').read())

        north = self.dom.getElementsByTagName('north')[0]
        south = self.dom.getElementsByTagName('south')[0]
        east = self.dom.getElementsByTagName('east')[0]
        west = self.dom.getElementsByTagName('west')[0]

        self.n_val = self.getValue(north.childNodes).strip()
        self.s_val = self.getValue(south.childNodes).strip()
        self.w_val = self.getValue(west.childNodes).strip()
        self.e_val = self.getValue(east.childNodes).strip()

        png_image = Image.open(self.png_file)

        width, height = png_image.size

        self.vrt_file = os.path.join(self.output_dir,
                      os.path.basename(self.png_file.split('.')[0]+'.vrt'))

        of = 'VRT'
        a_srs = 'EPSG:4326'
        cmd = 'gdal_translate -of '+ of +' -a_srs '+a_srs+\
        ' -gcp 0 0 %s %s' %(self.w_val, self.n_val) +\
        ' -gcp %s 0 %s %s' %(width, self.e_val, self.n_val)+\
        ' -gcp %s %s %s %s' %(width, height, self.e_val, self.s_val)+\
        ' -b 1 -b 2 -b 3 -b 4 %s %s' %(self.png_file, self.vrt_file)

        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        proc.wait()

        self.file_logger.debug('Step 1/3: gdal_translate')
        error = False
        for line in proc.stderr.readlines():
            if line:
                if not error:
                    self.file_logger.debug('Error: ')
                self.file_logger.debug(line)
                error = True
            print line
        if not error:
            self.file_logger.debug('success\n')

        os.remove(self.kml_file_temp)

    def gdalwrap(self):
        of = 'VRT'
        # 'EPSG:3857' - Spherical Mercator
        # 'EPSG:4326'
        # 'EPSG:900913'
        t_srs = self.t_srs if self.t_srs else 'EPSG:3857' #'EPSG:3857'
        self.vrt_file2 = self.vrt_file.rsplit('.',1)[0] + '.2.vrt'
        if os.path.isfile(self.vrt_file2):
            os.remove(self.vrt_file2)
        cmd = 'gdalwarp -of %s -t_srs %s %s %s' \
              %(of, t_srs, self.vrt_file, self.vrt_file2)
        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        proc.wait()

        self.file_logger.debug('Step 2/3: gdalwarp')
        error = False
        for line in proc.stderr.readlines():
            if line:
                if not error:
                    self.file_logger.debug('Error: ')
                self.file_logger.debug(line)
                error = True
            print line
        if not error:
            self.file_logger.debug('success\n')

    def gdal2tiles(self):
        zoom_text = '-z %s' %self.zoom if self.zoom else ''
        if args.output:
            self.tiles_dir = os.path.join(args.output, 'tiles')
        else:
            self.tiles_dir = 'tiles'

        cmd = 'gdal2tiles.py -w none --no-kml %s %s %s' \
                                %(zoom_text, self.vrt_file2, self.tiles_dir)
        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        proc.wait()

        self.file_logger.debug('Step 3/3: gdal2tiles.py')
        error = False
        for line in proc.stderr.readlines():
            if line:
                if not error:
                    self.file_logger.debug('Error: ')
                self.file_logger.debug(line)
                error = True
            print line
        if not error:
            self.file_logger.debug('success\n')

    def create_ini(self, notiles):
        self.ini_file = os.path.join(self.output_dir, 'desc.ini')
        now_time = datetime.datetime.now().isoformat()
        time_to_ini = now_time[:10] +' '+ now_time[11:19]
        about_ini = '; Generated by slbtiles\n' +\
                    '; %s\n' %time_to_ini +\
                    '; version 0.1\n\n'

        config = ConfigParser.RawConfigParser()

        config.add_section('general')
        dataset_name = os.path.split(os.path.abspath(self.input_dir))[1]
        config.set('general', 'dataset_name', dataset_name)

        if not notiles:
            config.add_section('zoom_level')
            tiles_files = os.listdir(self.tiles_dir)
            temp_zoom_list = filter (lambda x: x.isdigit(), tiles_files)
            zoom_levels = []
            for zoom in temp_zoom_list:
                try:
                    zoom_levels.append(int(zoom))
                except:
                    pass
            min_zoom_level = min(zoom_levels)
            max_zoom_level = max(zoom_levels)
            config.set('zoom_level', 'min_zoom_level', str(min_zoom_level))
            config.set('zoom_level', 'max_zoom_level', str(max_zoom_level))

        config.add_section('time_coverage')
        time_from_kml = self.dom.getElementsByTagName('when')
        if len(time_from_kml):
            time_from_kml = time_from_kml[0]
            t_c = self.getValue(time_from_kml.childNodes).strip()
            if t_c[4] == '-':
                time_coverage = t_c[:10] + ' ' + t_c[11:19]
            else: time_coverage = t_c[:8] + ' ' + t_c[9:19]
            config.set('time_coverage', 'begin_datetime', time_coverage)
            config.set('time_coverage', 'end_datetime', time_coverage)

        config.add_section('spatial_coverage')
        config.set('spatial_coverage', 'type', 'polygon')
        coord_ini = self.dom.getElementsByTagName('coordinates')
        if len(coord_ini):
            coord_ini = coord_ini[0]
            coordinates = self.getValue(coord_ini.childNodes).strip()
            coordinates_list = coordinates.split()
            config.set('spatial_coverage', 'coordinates_length',
                       str(len(coordinates_list)))
            for i in range(len(coordinates_list)):
                config.set('spatial_coverage', 'coordinates[%d]' %i,
                           coordinates_list[i])

        config.add_section('bbox')
        config.set('bbox', 'north', self.n_val)
        config.set('bbox', 'south', self.s_val)
        config.set('bbox', 'west', self.w_val)
        config.set('bbox', 'east', self.e_val)

        with open(self.ini_file, 'wb') as configfile:
            configfile.write(about_ini)
            config.write(configfile)

def create_output_dir(args):
    if not args.output:
        print "No output directory specified. Use -o"
        sys.exit(1)

    if os.path.isdir(args.output):
        print "Output dataset %s exists" %args.output
#        sys.exit(1)
    elif not os.path.isdir(args.output):
        try:
            os.makedirs(args.output)
#            uid = os.getuid()
#            gid = os.getgid()
#            os.chown(args.output, uid, gid)
        except OSError as e:
            print "Could not create directory %s:" %args.output, e
            sys.exit(1)

def red_dir(args):
    if not os.path.isdir(args.red_dir):
        print "Input directory %s does not exist" %args.red_dir
        sys.exit(1)

    create_output_dir(args)
    files_list = os.listdir(args.red_dir)
    file_logger = create_log(args)
    args.input = args.red_dir
    output_dir = args.output

    for _file in files_list:
        if _file.endswith('.png'):
            png_file = os.path.join(args.red_dir, _file)
            kml_file = os.path.join(args.red_dir, _file[:-3]+'kml')
            if os.path.isfile(png_file) and os.path.isfile(kml_file):

                _name, _extension = png_file.split('.')
#                _name_split = _name.split('_')
#                if len(_name_split) > 5:
#                    _satellite, _place, _date, _daytime, _hour, _red = \
#                                                            _name_split[-6:]

                args.output = os.path.join(output_dir, os.path.basename(_name))
                if not os.path.isdir(args.output):
                    try:
                        os.makedirs(args.output)
                    except OSError as e:
                        print "Could not create directory %s:" %args.output, e
                        sys.exit(1)

                try:
                    slb = SLBTiles(args, kml_file, png_file, file_logger)

                    slb.gdal_translate()

                    slb.gdalwrap()

                    if not args.notiles:
                        slb.gdal2tiles()
                    else:
                        notiles_log = 'Step 3/3: gdal2tiles.py\n'+\
                                      'using key --notiles, not running'
                        file_logger.debug(notiles_log)

                    slb.create_ini(args.notiles)
                    del slb
                except KeyboardInterrupt:
                    error_code = '\nInterrupted by the user'
                    print error_code
                    file_logger.debug(error_code+ "\n\n")
                    sys.exit(1)

                file_logger.debug("\n")
    sys.exit(0)

def get_range(red):
    # get list range from string, e.g. 5-9
    range_red = []
    if '-' in red:
        red_split = red.split('-')
        if not red_split[0].isdigit() or not red_split[1].isdigit():
            print "range for red1 is incorrect (format:'2-5' or '7')"
        for i in range(int(red_split[0]), int(red_split[1])+1):
            range_red.append(i)
    elif not red.isdigit():
        print "range for red1 is incorrect (format:'2-5' or '7')"
    else:
        range_red.append(red)
    return range_red

def red_copy(red, red_range, input_dir, output_dir):
    for zoom in red_range:
        src = os.path.join(input_dir, red, 'tiles', str(zoom))
        dst = os.path.join(output_dir, red[:-5], 'tiles', str(zoom))
        print 'FROM: ', src
        print 'TO: ', dst
        red_dir = os.path.join(output_dir, red[:-5])
        tiles_dir = os.path.join(output_dir, red[:-5], 'tiles')
        for path in [red_dir, tiles_dir, dst]:
            if not os.path.isdir(path):
                try:
                    os.mkdir(path)
#                    uid = os.getuid()
#                    gid = os.getgid()
#                    os.chown(args.output, uid, gid)
                except OSError as e:
                    print "Could not create directory %s:" %path, e
                    sys.exit(1)
        copytree(src, dst)

def move_red(args):
    # for moved tiles from red1, red2 and red4
    list_dir = os.listdir(args.input)
    for _dir in list_dir:
        if _dir.endswith('red1'):
            red1 = _dir
            red2 = _dir.rstrip('1') + '2'
            red4 = _dir.rstrip('1') + '4' if args.red4 else None
            if not os.path.isdir(os.path.join(args.input, red2)):
                print 'for %s not exists %s' %(red1, red2)
            if red4:
                if not os.path.isdir(os.path.join(args.input, red4)):
                    print 'for %s not exists %s' %(red1, red4)

            range_red1 = get_range(args.red1)
            range_red2 = get_range(args.red2)
            range_red4 = get_range(args.red4)

            red_copy(red1, range_red1, args.input, args.output)
            red_copy(red2, range_red2, args.input, args.output)
            red_copy(red4, range_red4, args.input, args.output)

            # copy ini file
            red1_dir = os.path.join(args.input, red1)
            for _file in os.listdir(red1_dir):
                if _file.endswith('.ini'): # or _file.endswith('.vrt'):
                    src_ini = os.path.join(args.input, red1, _file)
                    dst_ini = os.path.join(args.output, red1[:-5], _file)
                    shutil.copy2(src_ini, dst_ini)

    sys.exit(1)

if __name__=='__main__':
    parser = parse()
    args = parser.parse_args()
    kml_file, png_file = False, False

    if args.red_dir:
        red_dir(args)
        sys.exit(1)

    if not args.input:
        if not args.input_image or not args.input_kml:
            print "No input directory specified. Use -i or keys --image and "\
                  "--kml"
            sys.exit(1)
        else:
            args.input = os.path.split(args.input_image)[0]
            png_file = args.input_image
            kml_file = args.input_kml

    if not os.path.isdir(args.input):
        print "Input directory %s does not exist" %args.input
        sys.exit(1)

    create_output_dir(args)

    if args.move_red:
        move_red(args)

#    if not args.output:
#        print "No output directory specified. Use -o"
#        sys.exit(1)
#
#    if os.path.isdir(args.output):
#        print "Output dataset %s exists" %args.output
##        sys.exit(1)
#    elif not os.path.isdir(args.output):
#        try:
#            os.makedirs(args.output)
#        except OSError:
#            print "Could not create directory %s" %args.output
#            sys.exit(1)

    files_list = os.listdir(args.input)
    if not kml_file or not png_file:
        for files in files_list:
            if files.endswith('.kml'):
                kml_file = os.path.join (args.input, files)
            elif files.endswith('.png'):
                png_file = os.path.join (args.input, files)

    file_logger = create_log(args)

    if not png_file or not os.path.isfile(png_file):
        print "png file not found"
        file_logger.debug("png file not found")
        sys.exit(1)
    elif not kml_file or not os.path.isfile(kml_file):
        print "kml file not found"
        file_logger.debug("kml file not found")
        sys.exit(1)

    try:
        slb = SLBTiles(args, kml_file, png_file, file_logger)

        slb.gdal_translate()

        slb.gdalwrap()

        if not args.notiles:
            slb.gdal2tiles()
        else:
            notiles_log = 'Step 3/3: gdal2tiles.py\n'+\
                          'using key --notiles, not running'
            file_logger.debug(notiles_log)

        slb.create_ini(args.notiles)
    except KeyboardInterrupt:
        error_code = '\nInterrupted by the user'
        print error_code
        file_logger.debug(error_code+ "\n\n")
        sys.exit(1)

    file_logger.debug("\n")
    sys.exit(0)