import os
import glob
import datetime
import logging
import numpy
from dateutil import relativedelta
import pdb
import re

class URLSeries(object):
    """
    implement a chronologically ordered series of URLs (files,...)
    """
    def __init__(self, urlpattern=None, timepattern=None, periodicity=None, time_reference=None, **kwargs):
        """
        Create a chronologically ordered series of URLs (files,...) from which a feature (PointTimeSeries, GridTimeSeries,Trajectory) overlapping several files can be instantiated from.

        A time series can be defined in different ways:
        * providing a regular expression - or pattern

        SUPPRIMER TIME PATTERN ?

        :param urlpattern: pattern of the files, including the time information (use python datetime directives) for the feature start or reference time.
        :param time_increment: 
        """
        self.urlpattern = urlpattern
        self.timepattern = timepattern
        self.intervals = None
        if periodicity:
            if periodicity[-1] == 'S':
                self.periodicity = periodicity
            elif periodicity[-1] == 'M':
                self.periodicity = periodicity
            else:
                raise NotImplementedError
            if time_reference == None:
                pass
                #raise Exception("time_reference must be provided")
            else:
                self.timereference = time_reference
        else:
            self.periodicity = None
            self.timereference = None

    def get_example_url(self):
        """
        Return the reference URL of a URL series. This reference is used to
        extract the feature information (dimensions, variables, etc...).
        It can be any relevant file from the file series.

        If the reference time was not provided at the URL series creation,
        the latest time in the series is used
        """
        if not self.timereference is None:
            refURL = self.urlpattern.replace(
                        '%TIME',
                        self.timereference.strftime(self.timepattern['%TIME'])
                        )
            logging.debug("reference url : %s", refURL)
            return refURL
        else:
            raise NotImplementedError


    def get_url_for_time(self,time, proximity='exact',forbid_overlap=True, delta=relativedelta.relativedelta(day=1)):
        """
        :keyword forbid_overlap:forbids having several overlapping files for the same time step (if True raise exception). In such case, a liste of files is sent rather than a single file.
        proximity:
        -before: echeance before given time 
        -after: echeance after the given time
        -exact: given time correspond to an echeance
        -closest: look for the closest
        @return: datetime
        """
        # CASE 1 : time of URL is unique
#         print self.urlpattern
        if not '#S' in self.urlpattern:
            f = self.urlpattern.replace('%TIME', time.strftime(self.timepattern['%TIME']))
            fpattern = time.strftime( f )
#             print fpattern
            if os.path.exists(fpattern) and proximity=='exact':
                return fpattern, time
            if '*' in fpattern or '?' in fpattern:
                files = glob.glob(fpattern)
                if len(files) > 1:
                    if forbid_overlap:
                        raise Exception('Several files existing for the same time')
                    else:
                        return files, time
                if len(files) == 0:
                    return None,time
                else:
                    f = files[0]
                return f,time
            if proximity == 'exact':
                return None, time
            if self.timereference == None:
                raise Exception("timereference must be provided")
            if self.periodicity[-1] == 'S' or self.periodicity[-1] == 'M':
                periodicity = int(self.periodicity[:-1])
            else:
                raise Exception("Only periodicity expressed in seconds (S) are implemented")
            if self.periodicity[-1] == 'S':
#                 nbStep = int((time - self.timereference).total_seconds() / periodicity)
                nbStep = numpy.floor((time - self.timereference).total_seconds() / periodicity)
#                 logging.debug( 'number of step: %s %s',nbStep,nbStep*periodicity)
                time1 = (self.timereference + datetime.timedelta(seconds=nbStep*periodicity))
                time2 = self.timereference + datetime.timedelta(seconds=(nbStep+1)*periodicity)
#                 logging.debug('time1: %s time2: %s',time1,time2)
            elif self.periodicity[-1] == 'M':
                difftime = relativedelta.relativedelta(self.timereference, time)
                difftime_in_month = abs(difftime.years) * 12 + abs(difftime.months)
                nbstep = abs(difftime_in_month / periodicity)
                time1 = self.timereference\
                            + relativedelta.relativedelta(
                                    months=nbstep % 12,
                                    years=nbstep / 12
                                    )
                time2 = self.timereference\
                            + relativedelta.relativedelta(
                                    months=(nbstep + 1) % 12,
                                    years=(nbstep + 1) / 12
                                    )
            f1 = self.urlpattern.replace('%TIME',time1.strftime(self.timepattern['%TIME']))
#             print "FOUND : ", f1, time, time1,time2
            if time1 == time or proximity=='before':
                return f1,time1
            f2 = self.urlpattern.replace('%TIME',time2.strftime(self.timepattern['%TIME']))
            if proximity=='after':
                return f2,time2
            if proximity=='closest':
                if abs((time1-time).total_seconds()) <= abs((time2-time).total_seconds()):
                    return f1,time1
                else:
                    return f2,time2
    
        # CASE 2 : time of URL is expressed as an interval
        # ------------------------------------------------
        elif '#S' in self.urlpattern and '#E' in self.urlpattern:
            # search all possible files if no time interval was provided and build the list of time intervals
            if self.intervals is None:
                startDirectives = self.urlpattern.split('#S')
                url = ''
                for s in startDirectives:
                    ss = s
                    y = ss.rfind('%Y')
                    m = ss.rfind('%m')
                    e = ss.rfind("#E")
                    if not (y != -1 and y > e):
                        ss = ss.replace('%Y','????')
                    if not(m != -1 and m>e):
                        ss = ss.replace('%m','??')
                    ss = ss.replace('%d','??')
                    ss = ss.replace('%H','??')
                    ss = ss.replace('%M','??')
                    ss = ss.replace('%S','??')
                    url += ss
                url = url.replace('#E','')            
                urllist = {}
                start = time - delta
                end = time + delta
                curdate = start
                while True:
                    flist = glob.glob(curdate.strftime(url))
                    for f in flist:
                        if not urllist.has_key(f):
                            urllist[f] = ''
                    if (curdate.year == end.year and curdate.month == end.month):
                        break
                    curdate = curdate + relativedelta.relativedelta(month=1)
#                 print "URL ", url
#                 for f in urllist.keys():
#                     print f
#                     print URLSeries._get_interval_from_pattern( self.urlpattern, f )
#                 print urllist
            if self.intervals != None:
                # search the interval framing the time value (or the closest interval if not found)
                k = 0
                minDT = sys.maxint
                best = None
                for start,end in self.intervals:
                    if start <= time <= end:
                        best = k
                        return self._fill_url_pattern(self.urlpattern, self.intervals[k])
                    else:
                        tmpMin=min([abs(start-time), minDT, abs(end-time)])
                        if tmpMin <= minDT:
                            minDT = tmpMin
                            best = k
                #fill in URL pattern
                if proximity == 'exact':
                    return None
                elif proximity == 'closest':
                    return self._fill_url_pattern(self.urlpattern, self.intervals[k])
                else:
                    raise Exception("To be Implemented")
            else:
                #TBD
                pass

    def _get_time_directives(self, pattern_struct, code):
        '''
        :param type: type of time directive ('start' or 'end')
        '''
        strdirectives = []
        for s in pattern_struct:
            if isinstance(s,tuple):
                directive,c = s
                if c == code:
                    strdirectives.append(directive)
        return strdirectives

    def _decode_pattern(cls, pattern):
        '''
        return the structure of the pattern
        '''
        pattern_struct = []
        i = pattern.find('%')
        pattern_struct.append(pattern[0:i])
        k = 0
        while i != -1:
            print i,k
            k = pattern.find('#',i+1)
            code = pattern[k+1]
            directive = pattern[i:k]
            pattern_struct.append( (directive,'#'+code) )
            i = pattern.find('%',k+2)
            if i != -1:
                pattern_struct.append( pattern[k+2:i] )
        if k != 0:
            pattern_struct.append( pattern[k+2:] )
        return pattern_struct
    _decode_pattern=classmethod(_decode_pattern)


    def _get_time_from_directives(cls, url, pattern_struct, code):
        '''
        return the start or end time from an url following the provided pattern
        
        :param code: type of time directive (#S or #E)
        '''
        directives = URLSeries._get_time_directives (pattern_struct)
##        for s in pattern_struct:
##            if isinstance(s,tuple):
        
        subpattern = ''
        if code == '#S':
            anticode = '#E'
        else:
            anticode = '#S'
        for s in directives:
            ss = s
            y = ss.rfind('%Y')
            m = ss.rfind('%m')
            e = ss.rfind(anticode)
            if not (y != -1 and y > e):
                ss = ss.replace('%Y','????')
            if not(m != -1 and m>e):
                ss = ss.replace('%m','??')
            d = ss.rfind('%d')
            if not (d != -1 and d > e):
                ss = ss.replace('%d','??')
            H = ss.rfind('%H')
            if not (H != -1 and H > e):
                ss = ss.replace('%H','??')
            M = ss.rfind('%M')
            if not (M != -1 and M > e):
                ss = ss.replace('%M','??')
            S = ss.rfind('%S')
            if not (S != -1 and S > e):
                ss = ss.replace('%S','??')
            subpattern += ss
        subpattern = subpattern.replace(anticode,'')
        print "SUB : ", url, pattern, subpattern
        date = datetime.datetime.strptime( url,subpattern )
        return date
    _get_time_from_directives=classmethod(_get_time_from_directives)


    def _get_interval_from_pattern(cls, pattern, url):
        '''
        return the corresponding time interval of a file (using the file pattern to extract this information)
        '''
        start = URLSeries._get_time_from_directives(url,pattern,'#S')
        end = URLSeries._get_time_from_directives(url,pattern,'#E')
        return (start,end)
    _get_interval_from_pattern=classmethod(_get_interval_from_pattern)
    

    def _fill_url_pattern(cls, pattern, interval ):
        '''
        reconstruct an url from a pattern with time start/end fields
        '''
        start,end = interval
        startDirectives = pattern.split('#S')
        url = ''
        for s in startDirectives:
            i=-2
            while s[i-2] == '%':
                i = i-2
            start.strftime(s[i:])
            url = url + s
        endDirectives = url.split('#E')
        url = ''
        for s in endDirectives:
            i=-2
            while s[i-2] == '%':
                i = i-2
            end.strftime(s[i:])
            url = url + s
        return url    
    _fill_url_pattern=classmethod(_fill_url_pattern)


    def _get_urllist_and_timesteps_in_interval(self,start, end,forbid_overlap=True,unique_result=False):
        '''
        return the list of URLs of files intersecting or contained within a time interval
        '''
        res = []
#         logging.debug('start %s, end %s',start,end)
        # CASE 1 : file time pattern is regular and periodic (gridded files)
        if self.periodicity:
            if unique_result==True:
                middle = start + (end-start)/2
                logging.debug('get_url_for_time %s',middle)
                unique_file,t = self.get_url_for_time( middle, 'closest' ,forbid_overlap=forbid_overlap )
                res = (unique_file,t)
            else:
                periodicity = int(self.periodicity[:-1])
                #find the number of files between start and end
                fstart,tstart = self.get_url_for_time( start, 'before' ,forbid_overlap=forbid_overlap )
                fend,tend = self.get_url_for_time( end, 'after',forbid_overlap=forbid_overlap )
                tcurrent = tstart
    #             print  'interval file start file end',tstart,fstart,fend,tend
                if self.periodicity[-1] == 'S':
                    steprange = xrange( int(abs((tend-tstart).total_seconds()) / periodicity) + 1)
    #                 print 'steprange',steprange
                elif self.periodicity[-1] == 'M':
                    difftime = relativedelta.relativedelta(tend, tstart)
                    difftime_in_month = difftime.years*12 + difftime.months
                    steprange = xrange(abs(difftime_in_month) / periodicity)
                for i in steprange:
    #                 print 'tcurrent',tcurrent
                    f,t = self.get_url_for_time(tcurrent, 'exact', forbid_overlap=forbid_overlap )
                    if f is not None:
                        if type(f) == 'list':
                            for ff in f:
                                res.append((ff,tcurrent))
                        else:
                            res.append((f,tcurrent))
                    if self.periodicity[-1] == 'S':
                        tcurrent += datetime.timedelta(seconds=periodicity)
                    elif self.periodicity[-1] == 'M':
                        tcurrent += relativedelta.relativedelta(months=periodicity)
        # CASE 2 : file time pattern is random (orbit files)
        # repositories have to be scanned
        else:
            day1 = datetime.datetime(start.year,start.month,start.day)
            day2 = datetime.datetime(end.year,end.month,end.day)
            d = day1
#             logging.debug('[urlseries.py] day1:%s, day2:%s',day1,day2)
            while d <= day2:
                pattern = self.urlpattern.replace( '%TIME', d.strftime(self.timepattern['%TIME']) )
                logging.debug("URL series pattern : %s" % os.path.join( os.path.dirname( pattern ), '*' ))
                urls = glob.glob(os.path.join( os.path.dirname( pattern ), '*' ))
                urls.sort()
                patdep,patindep = self.urlpattern.split('%TIME')
                patdep = patdep + self.timepattern['%TIME']
                if self.timepattern['%TIME'][-1]=='/':
                    flag_filename_time_precision = False
                else:
                    flag_filename_time_precision = True
                logging.debug('urlseries| the filename contains a time pattern: %s',flag_filename_time_precision)
                # filter redundant time information in time pattern
#                 logging.debug('patdep %s',patdep)
                idy = patdep.find('%')
                timeInfo = {}
                while idy != -1:
#                         logging.debug('inot loop whgile %s %s',idx,patdep)
                    #print idx, patdep
                    field = patdep[idy+1]
                    if timeInfo.has_key(field):
                        #print '%'+field,d.strftime('%'+field)
                        replstr = d.strftime('%'+field)
                        patdep = patdep.replace( '%'+field,replstr, 1 )
                        idy = patdep.find('%', idy+len(replstr))
                    else:
                        timeInfo[field] = 0
                        idy = patdep.find('%', idy+2)
#                 logging.debug('patdep trasnformed %s',patdep)
#                 pdb.set_trace()
                for u in urls:
#                     logging.debug('url to treat %s',u)
                    #print u, self.urlpattern.replace( '%TIME', self.timepattern['%TIME'] )
                    # remove the time independant part of url and pattern to extract url reference time
                    patdepname, patdepext = os.path.splitext(patdep)
                    uname, uext = os.path.splitext(u)
#                     logging.debug('patdepext %s uext %s',patdepext,uext)
#                     if uext == patdepext: 
#                     idx = pattern.find('*')
#                     pdb.set_trace()
                    idx = [i for i, letter in enumerate(pattern) if letter == '*'][-1]
#                     for 
                    
#                     idx = [m.start() for m in re.finditer(pattern,'*')][-1]
                    udep = u[:idx]
#                     logging.debug('udep %s patdep %s', udep, patdep)
                    if '?' in patdep:
                        #remove variable character not depending on time, to help find the time of the file
                        dir_u = os.path.dirname(udep)
                        dir_p = os.path.dirname(patdep)
                        base_u = os.path.basename(udep)
                        base_p = os.path.basename(patdep)
#                         print 'basep',base_p
                        idx_question_mark = [i for i, letter in enumerate(base_p) if letter == '?']
#                         print 'idx',idx_question_mark
                        base_u = base_u[0:idx_question_mark[0]]+base_u[idx_question_mark[-1]+1:]
                        base_p = base_p[0:idx_question_mark[0]]+base_p[idx_question_mark[-1]+1:]
#                         print 'base_u',base_u
                        udep = os.path.join(dir_u,base_u)
                        patdep_simplified = os.path.join(dir_p,base_p)
                    else:
                        patdep_simplified = patdep
#                         logging.debug('udep %s patdep %s', udep, patdep_simplified)
#                     utime = datetime.datetime.strptime( udep, patdep )
#                     try:
                    utime = datetime.datetime.strptime( udep, patdep_simplified )
#                     logging.debug('utime %s start %s end %s',utime,start,end)
                    if flag_filename_time_precision == False:
                        #collect all the files of the directory
                        res.append( (u,utime) )
                    else:
                        #enlarge the 
#                         start.replace(hour=0,minute=0,second=0,microsecond=0)
#                         end.replace(hour=0,minute=0,second=0,microsecond=0)
                        if utime >= start and utime < end:
                            if unique_result == True:
                                middle = start + (end-start)/2
                                
                                if res != []:
    #                                 pdb.set_trace()
                                    if abs(middle-utime)<abs(middle-res[0][1]):
                                        res.append( (u,utime) )
                                else:
                                    res.append( (u,utime) )
                            else:
                                if not (u,utime) in res:
                                    res.append( (u,utime) )
#                                 logging.debug('[urlseries.py] url %s', u)
#                     except:
#                         pass
                d = d + datetime.timedelta(days=1)
        return res

    def get_urllist_in_interval(self, start, end, times=None, forbid_overlap=True,unique_result=False):
        '''
        return the list of files (or URLs) found in the *start/end* interval.

        If times is not None, it will contains the corresponding list of
        reference time associated with each found URL.
        '''
        tmp = self._get_urllist_and_timesteps_in_interval(start, end, forbid_overlap=forbid_overlap ,unique_result=unique_result)
        logging.debug('get_urllist_in_interval propose: %s',tmp)
#         if times:
#             del times[:]
#             times.extend([t for u, t in tmp])
#         logging.debug("Interval : %s", times)
        if isinstance(tmp,list):
            res = [u for u, t in tmp]
        else:
            res = [tmp[0]]
        final_res = []
        for uu in res:
            if os.path.exists(uu):
                final_res.append(uu)
        return final_res

    def get_timesteps_in_interval(self,start, end,forbid_overlap=True):
        return [t for u,t in self._get_urllist_and_timesteps_in_interval(start,end,forbid_overlap=forbid_overlap)]


if __name__ == '__main__' :
    
#    series = URLSeries(urlpattern="/home/cercache/project/ww3/public/HINDCAST/GLOBAL/2002_CFSR/hs/ww3.200209_hs.nc", \
#                        timepattern={"%TIME":'%Y_CFSR/hs/ww3.%Y%m'},periodicity="1M", \
#                        timereference=datetime.datetime(2002,09,01))
#    t = series.get_url_for_time(datetime.datetime(2003,09,01))
    
    t = datetime.datetime(2012,5,23,23,0,0)
    series = URLSeries(urlpattern="/home/cerdata/provider/ghrsst/satellite/l4/glob/odyssea-nrt/data/%TIME-IFR-L4_GHRSST-SSTfnd-ODYSSEA-GLOB_010-v2.0-fv1.0.nc", \
                        timePattern={"%TIME":'%Y/%j/%Y%m%d'},periodicity="86400S", \
                        timeReference=datetime.datetime(2011,9,1,0) )
    f = series.get_url_for_time(t, proximity='after')    
    print t, f 
    
    t1 = datetime.datetime(2012,5,23,23,0,0)
    t2 = datetime.datetime(2012,5,30,23,0,0)
    files = series.getURLInTimeInterval( t1, t2 )
    for f in files:
        print f


    
    
