#!/usr/bin/python
import re, sys
from os import path
from urllib import urlopen, urlretrieve
from datetime import datetime
from argparse import ArgumentParser

today = datetime.now().strftime('%d.%m.%Y')
date_pattern='^\d{2}\.\d{2}\.\d{4}$'
uri_pattern='a_date1={day_start}&a_date2={day_end}&f_ed3=4&f_ed4=4&f_ed5=3&f_pe=1&f_pe1=2&lng_id=2&wmo_id={station_id}'
href_pattern='<a href=([^>]+)>.*</a>'

parser = ArgumentParser(description='Download weather information from rp5.ru')
parser.add_argument('--station_id', dest='station_id', default=26063, help='Meteostation WMO ID', type=int)
parser.add_argument('--day', dest='day_start', default=today, help='Day of interest as d.m.Y (or start of interval)')
parser.add_argument('--day_end', dest='day_end', help='End day on interval(if any)')
parser.add_argument('--dest', dest='dest', default='/tmp', help='Destination directory')

args = parser.parse_args()

if not args.station_id:
    print 'station_id argument is mandatory'
    sys.exit(1)
if not args.day_end:
    args.day_end = args.day_start
if (not re.match(date_pattern, args.day_start)) or (not re.match(date_pattern, args.day_end)):
    print 'Date(s) are not valid, please check again'
    sys.exit(1)

dict_args = vars(args)
uri = uri_pattern.format(**dict_args)
h = urlopen('http://rp5.ru/inc/f_archive.php', uri)
data = h.read()

matches = re.search(href_pattern, data, re.MULTILINE)
if not matches:
    print 'Can\'t parse link from response: \n%s' % data
    sys.exit(1)

url = matches.group(1)

filename = url.split('/')[-1]
dest = path.join(args.dest, filename)
urlretrieve(url, dest)

print "File has to be located as path %s" % dest
