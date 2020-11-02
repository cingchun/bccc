#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   gvbbcccdf.py
#            Des:   gvbbcccdf
#           Mail:   qingchun720@foxmail.com
#   Created Time:   21:00 八月-25/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


from bccc import gvbbccc


class GVBBCCCDF(gvbbccc.GVBBCCC):
    '''
    Density fitting GVB-BCCC

    Attributes:


    Saved results：
    '''

    def __init__(self, mf, nt=4, nc=0, df=None):
        gvbbccc.GVBBCCC.__init__(mf, nt, nc)
        self.thd,self.ncyc,self.xE = 1.0e-7,500,0.015
        self.df = df




