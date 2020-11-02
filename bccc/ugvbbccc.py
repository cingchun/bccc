#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   ugvbbccc.py
#            Des:   ugvbbccc
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:18 八月-24/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


from bccc import gvbbccc


class UGVBBCCC(gvbbccc.GVBBCCC):
    '''
    High-spin GVB-BCCC

    Attributes:

    Saved results：
    '''

    def __init__(self, mf, na=0, nt=4, nc=0):
        gvbbccc.GVBBCCC.__init__(mf, nt, nc)
        self.thd,self.ncyc,self.xE = 1.0e-7,500,0.015
        self.na = na

