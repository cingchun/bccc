#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   ugvbbcccdf.py
#            Des:   ugvbbcccdf
#           Mail:   qingchun720@foxmail.com
#   Created Time:   21:00 八月-25/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


from bccc import gvbbcccdf

class UGVBBCCCDF(gvbbcccdf.GVBBCCCDF):
    def __init__(self, mf, na=0, nt=4, nc=0, df=None):
        gvbbcccdf.GVBBCCCDF.__init__(mf, nt, nc)
        self.thd,self.ncyc,self.xE = 1.0e-7,500,0.015
        self.na,self.df = na,df


