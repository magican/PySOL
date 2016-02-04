#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .ncfile import NCFile


class SaralNCFile(NCFile):

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
        '''
        fieldnames = NCFile.get_fieldnames(self)
        fieldnames_1hz = []
        for f in fieldnames:
            # filter out 40 Hz measurements
            #print f, self.get_dimensions(f)
            vdims = self.get_handler().variables[f].dimensions
            if not 'meas_ind' in list(vdims):
                fieldnames_1hz.append(f)
        return fieldnames_1hz

    def get_dimensions(self, fieldname=None):
        """
        """
        dims = NCFile.get_dimensions(self, fieldname)
        dims_1hz = []
        for d in dims:
            if not d == 'meas_ind':
                dims_1hz.append(d)
        return tuple(dims_1hz)
