#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
# 
# Author: Qingchun Wang 
# E-mail: qingchun720@foxmail.com 
# 


import numpy


nbs = 16
nrf = 312


def rep(ci):
    r = {}
    inv = numpy.linalg.inv(ci)
    
    # O0:  <Bp|  \hat{0}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = ci[0,1]
    m[12,1] = ci[1,1]
    m[11,2] = ci[2,3]
    m[11,3] = ci[3,3]
    m[13,5] = 1
    m[7,6] = 1
    m[4,8] = 1
    m[0,9] = inv[0,0]
    m[1,9] = inv[0,1]
    m[2,10] = inv[2,2]
    m[3,10] = inv[2,3]
    m[15,14] = 1
    r[0] = m
    
    # O1:  <Bp|  \hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,0] = ci[0,0]
    m[9,1] = ci[1,0]
    m[10,2] = ci[2,2]
    m[10,3] = ci[3,2]
    m[8,4] = 1
    m[6,7] = 1
    m[2,11] = inv[3,2]
    m[3,11] = inv[3,3]
    m[0,12] = inv[1,0]
    m[1,12] = inv[1,1]
    m[5,13] = 1
    m[14,15] = 1
    r[1] = m
    
    # O2:  <Bp|  \hat{0}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = ci[0,1]
    m[14,1] = ci[1,1]
    m[13,2] = -ci[2,2]*1
    m[13,3] = -ci[3,2]*1
    m[11,4] = -1.0
    m[9,6] = 1
    m[0,7] = -1.0*inv[0,0]
    m[1,7] = -1.0*inv[0,1]
    m[2,8] = inv[3,2]
    m[3,8] = inv[3,3]
    m[5,10] = 1
    m[15,12] = -1.0
    r[2] = m
    
    # O3:  <Bp|  \hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,0] = -ci[0,0]*1
    m[7,1] = -ci[1,0]*1
    m[8,2] = ci[2,3]
    m[8,3] = ci[3,3]
    m[10,5] = 1
    m[6,9] = 1
    m[4,11] = -1.0
    m[2,13] = -1.0*inv[2,2]
    m[3,13] = -1.0*inv[2,3]
    m[0,14] = inv[1,0]
    m[1,14] = inv[1,1]
    m[12,15] = -1.0
    r[3] = m
    
    # O4:  <Bp|  \hat{1}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = ci[0,0]
    m[11,1] = ci[1,0]
    m[12,2] = -ci[2,2]*1
    m[12,3] = -ci[3,2]*1
    m[14,5] = -1.0
    m[8,6] = 1
    m[4,7] = -1.0
    m[2,9] = -1.0*inv[3,2]
    m[3,9] = -1.0*inv[3,3]
    m[0,10] = inv[1,0]
    m[1,10] = inv[1,1]
    m[15,13] = 1
    r[4] = m
    
    # O5:  <Bp|  \hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,0] = ci[0,1]
    m[10,1] = ci[1,1]
    m[9,2] = -ci[2,3]*1
    m[9,3] = -ci[3,3]*1
    m[7,4] = -1.0
    m[6,8] = 1
    m[0,11] = inv[0,0]
    m[1,11] = inv[0,1]
    m[2,12] = -1.0*inv[2,2]
    m[3,12] = -1.0*inv[2,3]
    m[5,14] = -1.0
    m[13,15] = 1
    r[5] = m
    
    # O6:  <Bp|  \hat{1}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = ci[0,0]
    m[13,1] = ci[1,0]
    m[14,2] = ci[2,3]
    m[14,3] = ci[3,3]
    m[12,4] = 1
    m[10,6] = 1
    m[2,7] = -1.0*inv[2,2]
    m[3,7] = -1.0*inv[2,3]
    m[0,8] = -1.0*inv[1,0]
    m[1,8] = -1.0*inv[1,1]
    m[5,9] = -1.0
    m[15,11] = -1.0
    r[6] = m
    
    # O7:  <Bp|  \hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,0] = -ci[0,1]*1
    m[8,1] = -ci[1,1]*1
    m[7,2] = -ci[2,2]*1
    m[7,3] = -ci[3,2]*1
    m[9,5] = -1.0
    m[6,10] = 1
    m[4,12] = 1
    m[0,13] = inv[0,0]
    m[1,13] = inv[0,1]
    m[2,14] = inv[3,2]
    m[3,14] = inv[3,3]
    m[11,15] = -1.0
    r[7] = m
    
    # O8:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}  |Bq> = 
    
    # O9:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[0,0]
    m[1,0] = ci[0,0]*inv[0,1]
    m[0,1] = ci[1,0]*inv[0,0]
    m[1,1] = ci[1,0]*inv[0,1]
    m[2,2] = ci[2,2]*inv[2,2]
    m[3,2] = ci[2,2]*inv[2,3]
    m[2,3] = ci[3,2]*inv[2,2]
    m[3,3] = ci[3,2]*inv[2,3]
    m[4,4] = 1
    m[7,7] = 1
    m[11,11] = 1
    m[12,12] = 1
    m[13,13] = 1
    m[15,15] = 1
    r[9] = m
    
    # O10:  <Bp|  \hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O11:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,0] = ci[0,1]
    m[15,1] = ci[1,1]
    m[0,6] = inv[0,0]
    m[1,6] = inv[0,1]
    m[11,8] = 1
    m[13,10] = 1
    r[11] = m
    
    # O12:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,2] = ci[2,3]
    m[4,3] = ci[3,3]
    m[2,5] = inv[2,2]
    m[3,5] = inv[2,3]
    m[7,9] = 1
    m[12,14] = 1
    r[12] = m
    
    # O13:  <Bp|  \hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,0] = -ci[0,0]*1
    m[6,1] = -ci[1,0]*1
    m[8,11] = -1.0
    m[10,13] = -1.0
    m[0,15] = -1.0*inv[1,0]
    m[1,15] = -1.0*inv[1,1]
    r[13] = m
    
    # O14:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,5] = -1.0
    m[4,6] = 1
    m[11,9] = -1.0
    m[12,10] = 1
    r[14] = m
    
    # O15:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,1]*inv[2,2]
    m[3,0] = ci[0,1]*inv[2,3]
    m[2,1] = ci[1,1]*inv[2,2]
    m[3,1] = ci[1,1]*inv[2,3]
    m[0,2] = -ci[2,3]*inv[0,0]
    m[1,2] = -ci[2,3]*inv[0,1]
    m[0,3] = -ci[3,3]*inv[0,0]
    m[1,3] = -ci[3,3]*inv[0,1]
    m[7,8] = 1
    m[13,14] = -1.0
    r[15] = m
    
    # O16:  <Bp|  \hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,4] = -1.0
    m[9,11] = 1
    m[10,12] = -1.0
    m[5,15] = 1
    r[16] = m
    
    # O17:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,2] = ci[2,3]
    m[15,3] = ci[3,3]
    m[2,6] = inv[2,2]
    m[3,6] = inv[2,3]
    m[12,8] = -1.0
    m[13,9] = -1.0
    r[17] = m
    
    # O18:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,0] = -ci[0,1]*1
    m[4,1] = -ci[1,1]*1
    m[0,5] = -1.0*inv[0,0]
    m[1,5] = -1.0*inv[0,1]
    m[7,10] = 1
    m[11,14] = 1
    r[18] = m
    
    # O19:  <Bp|  \hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,2] = -ci[2,2]*1
    m[6,3] = -ci[3,2]*1
    m[8,12] = 1
    m[9,13] = 1
    m[2,15] = -1.0*inv[3,2]
    m[3,15] = -1.0*inv[3,3]
    r[19] = m
    
    # O20:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,0] = -ci[0,1]*1
    m[15,1] = -ci[1,1]*1
    m[0,6] = -1.0*inv[0,0]
    m[1,6] = -1.0*inv[0,1]
    m[11,8] = -1.0
    m[13,10] = -1.0
    r[20] = m
    
    # O21:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,2] = ci[2,2]
    m[5,3] = ci[3,2]
    m[2,4] = inv[3,2]
    m[3,4] = inv[3,3]
    m[9,7] = 1
    m[14,12] = 1
    r[21] = m
    
    # O22:  <Bp|  \hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,0] = ci[0,0]
    m[6,1] = ci[1,0]
    m[8,11] = 1
    m[10,13] = 1
    m[0,15] = inv[1,0]
    m[1,15] = inv[1,1]
    r[22] = m
    
    # O23:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}  |Bq> = 
    
    # O24:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[0,0]
    m[1,0] = ci[0,0]*inv[0,1]
    m[0,1] = ci[1,0]*inv[0,0]
    m[1,1] = ci[1,0]*inv[0,1]
    m[2,2] = ci[2,3]*inv[3,2]
    m[3,2] = ci[2,3]*inv[3,3]
    m[2,3] = ci[3,3]*inv[3,2]
    m[3,3] = ci[3,3]*inv[3,3]
    m[5,5] = 1
    m[9,9] = 1
    m[11,11] = 1.0
    m[13,13] = 1.0
    m[14,14] = 1
    m[15,15] = 1.0
    r[24] = m
    
    # O25:  <Bp|  \hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O26:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,2] = ci[2,2]
    m[15,3] = ci[3,2]
    m[2,6] = inv[3,2]
    m[3,6] = inv[3,3]
    m[11,7] = 1
    m[14,10] = 1
    r[26] = m
    
    # O27:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,0] = ci[0,1]
    m[5,1] = ci[1,1]
    m[0,4] = inv[0,0]
    m[1,4] = inv[0,1]
    m[9,8] = 1
    m[13,12] = 1.0
    r[27] = m
    
    # O28:  <Bp|  \hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,2] = -ci[2,3]*1
    m[6,3] = -ci[3,3]*1
    m[7,11] = -1.0
    m[10,14] = -1.0
    m[2,15] = -1.0*inv[2,2]
    m[3,15] = -1.0*inv[2,3]
    r[28] = m
    
    # O29:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,4] = -1.0
    m[5,6] = 1
    m[13,7] = 1
    m[14,8] = -1.0
    r[29] = m
    
    # O30:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,1]*inv[3,2]
    m[3,0] = -ci[0,1]*inv[3,3]
    m[2,1] = -ci[1,1]*inv[3,2]
    m[3,1] = -ci[1,1]*inv[3,3]
    m[0,2] = ci[2,2]*inv[0,0]
    m[1,2] = ci[2,2]*inv[0,1]
    m[0,3] = ci[3,2]*inv[0,0]
    m[1,3] = ci[3,2]*inv[0,1]
    m[9,10] = 1
    m[11,12] = -1.0
    r[30] = m
    
    # O31:  <Bp|  \hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,5] = -1.0
    m[7,13] = -1.0
    m[8,14] = 1
    m[4,15] = 1
    r[31] = m
    
    # O32:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,5] = 1
    m[4,6] = -1.0
    m[11,9] = 1
    m[12,10] = -1.0
    r[32] = m
    
    # O33:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,0]*inv[3,2]
    m[3,0] = -ci[0,0]*inv[3,3]
    m[2,1] = -ci[1,0]*inv[3,2]
    m[3,1] = -ci[1,0]*inv[3,3]
    m[0,2] = ci[2,2]*inv[1,0]
    m[1,2] = ci[2,2]*inv[1,1]
    m[0,3] = ci[3,2]*inv[1,0]
    m[1,3] = ci[3,2]*inv[1,1]
    m[8,7] = 1
    m[14,13] = -1.0
    r[33] = m
    
    # O34:  <Bp|  \hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,4] = 1
    m[9,11] = -1.0
    m[10,12] = 1
    m[5,15] = -1.0
    r[34] = m
    
    # O35:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,2] = -ci[2,2]*1
    m[15,3] = -ci[3,2]*1
    m[2,6] = -1.0*inv[3,2]
    m[3,6] = -1.0*inv[3,3]
    m[11,7] = -1.0
    m[14,10] = -1.0
    r[35] = m
    
    # O36:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,0] = ci[0,0]
    m[4,1] = ci[1,0]
    m[0,5] = inv[1,0]
    m[1,5] = inv[1,1]
    m[8,9] = 1
    m[12,13] = 1.0
    r[36] = m
    
    # O37:  <Bp|  \hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,2] = ci[2,3]
    m[6,3] = ci[3,3]
    m[7,11] = 1
    m[10,14] = 1
    m[2,15] = inv[2,2]
    m[3,15] = inv[2,3]
    r[37] = m
    
    # O38:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}  |Bq> = 
    
    # O39:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[1,0]
    m[1,0] = ci[0,1]*inv[1,1]
    m[0,1] = ci[1,1]*inv[1,0]
    m[1,1] = ci[1,1]*inv[1,1]
    m[2,2] = ci[2,3]*inv[3,2]
    m[3,2] = ci[2,3]*inv[3,3]
    m[2,3] = ci[3,3]*inv[3,2]
    m[3,3] = ci[3,3]*inv[3,3]
    m[4,4] = 1.0
    m[8,8] = 1
    m[11,11] = 1
    m[12,12] = 1.0
    m[14,14] = 1.0
    m[15,15] = 1
    r[39] = m
    
    # O40:  <Bp|  \hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O41:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,0] = ci[0,0]
    m[15,1] = ci[1,0]
    m[0,6] = inv[1,0]
    m[1,6] = inv[1,1]
    m[12,7] = 1
    m[14,9] = 1
    r[41] = m
    
    # O42:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,2] = ci[2,2]
    m[4,3] = ci[3,2]
    m[2,5] = inv[3,2]
    m[3,5] = inv[3,3]
    m[8,10] = 1
    m[11,13] = 1
    r[42] = m
    
    # O43:  <Bp|  \hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,0] = -ci[0,1]*1
    m[6,1] = -ci[1,1]*1
    m[7,12] = -1.0
    m[9,14] = -1.0
    m[0,15] = -1.0*inv[0,0]
    m[1,15] = -1.0*inv[0,1]
    r[43] = m
    
    # O44:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,2] = -ci[2,3]*1
    m[15,3] = -ci[3,3]*1
    m[2,6] = -1.0*inv[2,2]
    m[3,6] = -1.0*inv[2,3]
    m[12,8] = 1
    m[13,9] = 1
    r[44] = m
    
    # O45:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,0] = -ci[0,0]*1
    m[5,1] = -ci[1,0]*1
    m[0,4] = -1.0*inv[1,0]
    m[1,4] = -1.0*inv[1,1]
    m[10,7] = 1
    m[14,11] = 1
    r[45] = m
    
    # O46:  <Bp|  \hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,2] = ci[2,2]
    m[6,3] = ci[3,2]
    m[8,12] = -1.0
    m[9,13] = -1.0
    m[2,15] = inv[3,2]
    m[3,15] = inv[3,3]
    r[46] = m
    
    # O47:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,4] = 1
    m[5,6] = -1.0
    m[13,7] = -1.0
    m[14,8] = 1
    r[47] = m
    
    # O48:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,0]*inv[2,2]
    m[3,0] = ci[0,0]*inv[2,3]
    m[2,1] = ci[1,0]*inv[2,2]
    m[3,1] = ci[1,0]*inv[2,3]
    m[0,2] = -ci[2,3]*inv[1,0]
    m[1,2] = -ci[2,3]*inv[1,1]
    m[0,3] = -ci[3,3]*inv[1,0]
    m[1,3] = -ci[3,3]*inv[1,1]
    m[10,9] = 1
    m[12,11] = -1.0
    r[48] = m
    
    # O49:  <Bp|  \hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,5] = 1
    m[7,13] = 1
    m[8,14] = -1.0
    m[4,15] = -1.0
    r[49] = m
    
    # O50:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[15,0] = -ci[0,0]*1
    m[15,1] = -ci[1,0]*1
    m[0,6] = -1.0*inv[1,0]
    m[1,6] = -1.0*inv[1,1]
    m[12,7] = -1.0
    m[14,9] = -1.0
    r[50] = m
    
    # O51:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,2] = ci[2,3]
    m[5,3] = ci[3,3]
    m[2,4] = inv[2,2]
    m[3,4] = inv[2,3]
    m[10,8] = 1
    m[13,11] = 1
    r[51] = m
    
    # O52:  <Bp|  \hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[6,0] = ci[0,1]
    m[6,1] = ci[1,1]
    m[7,12] = 1
    m[9,14] = 1
    m[0,15] = inv[0,0]
    m[1,15] = inv[0,1]
    r[52] = m
    
    # O53:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}  |Bq> = 
    
    # O54:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[1,0]
    m[1,0] = ci[0,1]*inv[1,1]
    m[0,1] = ci[1,1]*inv[1,0]
    m[1,1] = ci[1,1]*inv[1,1]
    m[2,2] = ci[2,2]*inv[2,2]
    m[3,2] = ci[2,2]*inv[2,3]
    m[2,3] = ci[3,2]*inv[2,2]
    m[3,3] = ci[3,2]*inv[2,3]
    m[5,5] = 1.0
    m[10,10] = 1
    m[12,12] = 1
    m[13,13] = 1
    m[14,14] = 1
    m[15,15] = 1.0
    r[54] = m
    
    # O55:  <Bp|  \hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O56:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O57:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O58:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O59:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,0] = -ci[0,0]*1
    m[7,1] = -ci[1,0]*1
    m[4,11] = -1.0
    m[2,13] = -1.0*inv[2,2]
    m[3,13] = -1.0*inv[2,3]
    m[12,15] = -1.0
    r[59] = m
    
    # O60:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O61:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,4] = -1.0
    m[0,11] = inv[0,0]
    m[1,11] = inv[0,1]
    m[2,12] = -1.0*inv[2,2]
    m[3,12] = -1.0*inv[2,3]
    m[13,15] = 1
    r[61] = m
    
    # O62:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O63:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,2] = -ci[2,2]*1
    m[7,3] = -ci[3,2]*1
    m[4,12] = 1
    m[0,13] = inv[0,0]
    m[1,13] = inv[0,1]
    m[11,15] = -1.0
    r[63] = m
    
    # O64:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,2] = ci[2,2]
    m[13,3] = ci[3,2]
    m[11,4] = 1
    m[0,7] = inv[0,0]
    m[1,7] = inv[0,1]
    m[15,12] = 1
    r[64] = m
    
    # O65:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,0] = ci[0,0]
    m[7,1] = ci[1,0]
    m[4,11] = 1
    m[2,13] = inv[2,2]
    m[3,13] = inv[2,3]
    m[12,15] = 1
    r[65] = m
    
    # O66:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,2] = ci[2,3]
    m[11,3] = ci[3,3]
    m[13,5] = 1
    m[0,9] = inv[0,0]
    m[1,9] = inv[0,1]
    m[15,14] = 1
    r[66] = m
    
    # O67:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O68:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = ci[0,1]
    m[13,1] = ci[1,1]
    m[0,8] = inv[0,0]
    m[1,8] = inv[0,1]
    r[68] = m
    
    # O69:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,2] = -ci[2,3]*1
    m[7,3] = -ci[3,3]*1
    m[2,14] = -1.0*inv[2,2]
    m[3,14] = -1.0*inv[2,3]
    r[69] = m
    
    # O70:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = -ci[0,1]*1
    m[11,1] = -ci[1,1]*1
    m[0,10] = inv[0,0]
    m[1,10] = inv[0,1]
    r[70] = m
    
    # O71:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,5] = -1.0
    m[4,14] = 1
    r[71] = m
    
    # O72:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = -ci[0,0]*1
    m[11,1] = -ci[1,0]*1
    m[12,2] = ci[2,2]
    m[12,3] = ci[3,2]
    m[4,7] = 1
    m[15,13] = -1.0
    r[72] = m
    
    # O73:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,4] = 1
    m[0,11] = -1.0*inv[0,0]
    m[1,11] = -1.0*inv[0,1]
    m[2,12] = inv[2,2]
    m[3,12] = inv[2,3]
    m[13,15] = -1.0
    r[73] = m
    
    # O74:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,5] = 1
    m[4,9] = 1
    r[74] = m
    
    # O75:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,2] = ci[2,3]
    m[7,3] = ci[3,3]
    m[2,14] = inv[2,2]
    m[3,14] = inv[2,3]
    r[75] = m
    
    # O76:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = ci[0,1]
    m[12,1] = ci[1,1]
    m[11,2] = ci[2,3]
    m[11,3] = ci[3,3]
    m[4,8] = 1
    m[15,14] = 1.0
    r[76] = m
    
    # O77:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O78:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,5] = 1.0
    m[4,10] = 1
    r[78] = m
    
    # O79:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,0] = -ci[0,1]*1
    m[7,1] = -ci[1,1]*1
    m[0,14] = -1.0*inv[0,0]
    m[1,14] = -1.0*inv[0,1]
    r[79] = m
    
    # O80:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = -ci[0,0]*1
    m[13,1] = -ci[1,0]*1
    m[12,4] = -1.0
    m[2,7] = inv[2,2]
    m[3,7] = inv[2,3]
    m[15,11] = 1
    r[80] = m
    
    # O81:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,2] = ci[2,2]
    m[7,3] = ci[3,2]
    m[4,12] = -1.0
    m[0,13] = -1.0*inv[0,0]
    m[1,13] = -1.0*inv[0,1]
    m[11,15] = 1
    r[81] = m
    
    # O82:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,2] = -ci[2,3]*1
    m[12,3] = -ci[3,3]*1
    m[2,9] = inv[2,2]
    m[3,9] = inv[2,3]
    r[82] = m
    
    # O83:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,5] = 1
    m[4,14] = -1.0
    r[83] = m
    
    # O84:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,2] = ci[2,3]
    m[13,3] = ci[3,3]
    m[2,8] = inv[2,2]
    m[3,8] = inv[2,3]
    r[84] = m
    
    # O85:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[7,0] = ci[0,1]
    m[7,1] = ci[1,1]
    m[0,14] = inv[0,0]
    m[1,14] = inv[0,1]
    r[85] = m
    
    # O86:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = ci[0,1]
    m[12,1] = ci[1,1]
    m[13,5] = 1.0
    m[2,10] = inv[2,2]
    m[3,10] = inv[2,3]
    m[15,14] = 1
    r[86] = m
    
    # O87:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O88:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,2] = -ci[2,2]*1
    m[13,3] = -ci[3,2]*1
    m[11,4] = -1.0
    m[0,7] = -1.0*inv[0,0]
    m[1,7] = -1.0*inv[0,1]
    m[15,12] = -1.0
    r[88] = m
    
    # O89:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O90:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,2] = -ci[2,3]*1
    m[11,3] = -ci[3,3]*1
    m[13,5] = -1.0
    m[0,9] = -1.0*inv[0,0]
    m[1,9] = -1.0*inv[0,1]
    m[15,14] = -1.0
    r[90] = m
    
    # O91:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,0] = -ci[0,0]*1
    m[9,1] = -ci[1,0]*1
    m[2,11] = -1.0*inv[3,2]
    m[3,11] = -1.0*inv[3,3]
    m[5,13] = -1.0
    m[14,15] = -1.0
    r[91] = m
    
    # O92:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = -ci[0,1]*1
    m[13,1] = -ci[1,1]*1
    m[0,8] = -1.0*inv[0,0]
    m[1,8] = -1.0*inv[0,1]
    r[92] = m
    
    # O93:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,4] = -1.0
    m[5,12] = -1.0
    r[93] = m
    
    # O94:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = ci[0,1]
    m[11,1] = ci[1,1]
    m[0,10] = -1.0*inv[0,0]
    m[1,10] = -1.0*inv[0,1]
    r[94] = m
    
    # O95:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,2] = -ci[2,2]*1
    m[9,3] = -ci[3,2]*1
    m[2,12] = inv[3,2]
    m[3,12] = inv[3,3]
    r[95] = m
    
    # O96:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O97:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,0] = ci[0,0]
    m[9,1] = ci[1,0]
    m[2,11] = inv[3,2]
    m[3,11] = inv[3,3]
    m[5,13] = 1
    m[14,15] = 1
    r[97] = m
    
    # O98:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O99:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O100:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O101:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,2] = -ci[2,3]*1
    m[9,3] = -ci[3,3]*1
    m[0,11] = inv[0,0]
    m[1,11] = inv[0,1]
    m[5,14] = -1.0
    m[13,15] = 1.0
    r[101] = m
    
    # O102:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O103:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,5] = -1.0
    m[0,13] = inv[0,0]
    m[1,13] = inv[0,1]
    m[2,14] = inv[3,2]
    m[3,14] = inv[3,3]
    m[11,15] = -1.0
    r[103] = m
    
    # O104:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,2] = ci[2,2]
    m[14,3] = ci[3,2]
    m[2,7] = inv[3,2]
    m[3,7] = inv[3,3]
    r[104] = m
    
    # O105:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,4] = 1
    m[5,12] = 1
    r[105] = m
    
    # O106:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = -ci[0,0]*1
    m[11,1] = -ci[1,0]*1
    m[14,5] = 1
    m[2,9] = inv[3,2]
    m[3,9] = inv[3,3]
    m[15,13] = -1.0
    r[106] = m
    
    # O107:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,2] = ci[2,3]
    m[9,3] = ci[3,3]
    m[0,11] = -1.0*inv[0,0]
    m[1,11] = -1.0*inv[0,1]
    m[5,14] = 1
    m[13,15] = -1.0
    r[107] = m
    
    # O108:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = ci[0,1]
    m[14,1] = ci[1,1]
    m[11,4] = -1.0
    m[2,8] = inv[3,2]
    m[3,8] = inv[3,3]
    m[15,12] = -1.0
    r[108] = m
    
    # O109:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O110:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,2] = -ci[2,2]*1
    m[11,3] = -ci[3,2]*1
    m[2,10] = inv[3,2]
    m[3,10] = inv[3,3]
    r[110] = m
    
    # O111:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,0] = -ci[0,1]*1
    m[9,1] = -ci[1,1]*1
    m[0,12] = inv[0,0]
    m[1,12] = inv[0,1]
    r[111] = m
    
    # O112:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,4] = -1.0
    m[5,7] = 1
    r[112] = m
    
    # O113:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,2] = ci[2,2]
    m[9,3] = ci[3,2]
    m[2,12] = -1.0*inv[3,2]
    m[3,12] = -1.0*inv[3,3]
    r[113] = m
    
    # O114:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = -ci[0,0]*1
    m[13,1] = -ci[1,0]*1
    m[14,2] = -ci[2,3]*1
    m[14,3] = -ci[3,3]*1
    m[5,9] = 1
    m[15,11] = 1.0
    r[114] = m
    
    # O115:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,5] = 1
    m[0,13] = -1.0*inv[0,0]
    m[1,13] = -1.0*inv[0,1]
    m[2,14] = -1.0*inv[3,2]
    m[3,14] = -1.0*inv[3,3]
    m[11,15] = 1.0
    r[115] = m
    
    # O116:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,4] = -1.0
    m[5,8] = 1
    r[116] = m
    
    # O117:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[9,0] = ci[0,1]
    m[9,1] = ci[1,1]
    m[0,12] = -1.0*inv[0,0]
    m[1,12] = -1.0*inv[0,1]
    r[117] = m
    
    # O118:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = ci[0,1]
    m[14,1] = ci[1,1]
    m[13,2] = -ci[2,2]*1
    m[13,3] = -ci[3,2]*1
    m[5,10] = 1
    m[15,12] = -1.0
    r[118] = m
    
    # O119:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O120:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = ci[0,0]
    m[11,1] = ci[1,0]
    m[12,2] = -ci[2,2]*1
    m[12,3] = -ci[3,2]*1
    m[4,7] = -1.0
    m[15,13] = 1
    r[120] = m
    
    # O121:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O122:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,5] = -1.0
    m[4,9] = -1.0
    r[122] = m
    
    # O123:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,0] = -ci[0,0]*1
    m[8,1] = -ci[1,0]*1
    m[0,13] = -1.0*inv[1,0]
    m[1,13] = -1.0*inv[1,1]
    r[123] = m
    
    # O124:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = -ci[0,1]*1
    m[12,1] = -ci[1,1]*1
    m[11,2] = -ci[2,3]*1
    m[11,3] = -ci[3,3]*1
    m[4,8] = -1.0
    m[15,14] = -1.0
    r[124] = m
    
    # O125:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,4] = -1.0
    m[2,11] = -1.0*inv[3,2]
    m[3,11] = -1.0*inv[3,3]
    m[0,12] = -1.0*inv[1,0]
    m[1,12] = -1.0*inv[1,1]
    m[14,15] = -1.0
    r[125] = m
    
    # O126:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,5] = -1.0
    m[4,10] = -1.0
    r[126] = m
    
    # O127:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,2] = -ci[2,2]*1
    m[8,3] = -ci[3,2]*1
    m[2,13] = -1.0*inv[3,2]
    m[3,13] = -1.0*inv[3,3]
    r[127] = m
    
    # O128:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,2] = -ci[2,2]*1
    m[14,3] = -ci[3,2]*1
    m[2,7] = -1.0*inv[3,2]
    m[3,7] = -1.0*inv[3,3]
    r[128] = m
    
    # O129:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,0] = ci[0,0]
    m[8,1] = ci[1,0]
    m[0,13] = inv[1,0]
    m[1,13] = inv[1,1]
    r[129] = m
    
    # O130:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,0] = ci[0,0]
    m[11,1] = ci[1,0]
    m[14,5] = -1.0
    m[2,9] = -1.0*inv[3,2]
    m[3,9] = -1.0*inv[3,3]
    m[15,13] = 1.0
    r[130] = m
    
    # O131:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O132:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = -ci[0,1]*1
    m[14,1] = -ci[1,1]*1
    m[11,4] = 1.0
    m[2,8] = -1.0*inv[3,2]
    m[3,8] = -1.0*inv[3,3]
    m[15,12] = 1.0
    r[132] = m
    
    # O133:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,2] = -ci[2,3]*1
    m[8,3] = -ci[3,3]*1
    m[4,11] = 1.0
    m[0,14] = -1.0*inv[1,0]
    m[1,14] = -1.0*inv[1,1]
    m[12,15] = 1.0
    r[133] = m
    
    # O134:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[11,2] = ci[2,2]
    m[11,3] = ci[3,2]
    m[2,10] = -1.0*inv[3,2]
    m[3,10] = -1.0*inv[3,3]
    r[134] = m
    
    # O135:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,5] = -1.0
    m[4,13] = 1.0
    r[135] = m
    
    # O136:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O137:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,4] = 1
    m[2,11] = inv[3,2]
    m[3,11] = inv[3,3]
    m[0,12] = inv[1,0]
    m[1,12] = inv[1,1]
    m[14,15] = 1.0
    r[137] = m
    
    # O138:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O139:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,2] = ci[2,3]
    m[8,3] = ci[3,3]
    m[4,11] = -1.0
    m[0,14] = inv[1,0]
    m[1,14] = inv[1,1]
    m[12,15] = -1.0
    r[139] = m
    
    # O140:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O141:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O142:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O143:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,0] = -ci[0,1]*1
    m[8,1] = -ci[1,1]*1
    m[4,12] = 1.0
    m[2,14] = inv[3,2]
    m[3,14] = inv[3,3]
    m[11,15] = -1.0
    r[143] = m
    
    # O144:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = ci[0,0]
    m[14,1] = ci[1,0]
    m[0,7] = inv[1,0]
    m[1,7] = inv[1,1]
    r[144] = m
    
    # O145:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,2] = ci[2,2]
    m[8,3] = ci[3,2]
    m[2,13] = inv[3,2]
    m[3,13] = inv[3,3]
    r[145] = m
    
    # O146:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = -ci[0,0]*1
    m[12,1] = -ci[1,0]*1
    m[0,9] = inv[1,0]
    m[1,9] = inv[1,1]
    r[146] = m
    
    # O147:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,5] = 1
    m[4,13] = -1.0
    r[147] = m
    
    # O148:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,2] = -ci[2,3]*1
    m[14,3] = -ci[3,3]*1
    m[12,4] = -1.0
    m[0,8] = inv[1,0]
    m[1,8] = inv[1,1]
    m[15,11] = 1
    r[148] = m
    
    # O149:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[8,0] = ci[0,1]
    m[8,1] = ci[1,1]
    m[4,12] = -1.0
    m[2,14] = -1.0*inv[3,2]
    m[3,14] = -1.0*inv[3,3]
    m[11,15] = 1
    r[149] = m
    
    # O150:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,2] = -ci[2,2]*1
    m[12,3] = -ci[3,2]*1
    m[14,5] = -1.0
    m[0,10] = inv[1,0]
    m[1,10] = inv[1,1]
    m[15,13] = 1
    r[150] = m
    
    # O151:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O152:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = ci[0,0]
    m[13,1] = ci[1,0]
    m[12,4] = 1
    m[2,7] = -1.0*inv[2,2]
    m[3,7] = -1.0*inv[2,3]
    m[15,11] = -1.0
    r[152] = m
    
    # O153:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O154:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,2] = ci[2,3]
    m[12,3] = ci[3,3]
    m[2,9] = -1.0*inv[2,2]
    m[3,9] = -1.0*inv[2,3]
    r[154] = m
    
    # O155:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,0] = -ci[0,0]*1
    m[10,1] = -ci[1,0]*1
    m[0,11] = inv[1,0]
    m[1,11] = inv[1,1]
    r[155] = m
    
    # O156:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,2] = -ci[2,3]*1
    m[13,3] = -ci[3,3]*1
    m[2,8] = -1.0*inv[2,2]
    m[3,8] = -1.0*inv[2,3]
    r[156] = m
    
    # O157:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,4] = -1.0
    m[5,11] = -1.0
    r[157] = m
    
    # O158:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = -ci[0,1]*1
    m[12,1] = -ci[1,1]*1
    m[13,5] = -1.0
    m[2,10] = -1.0*inv[2,2]
    m[3,10] = -1.0*inv[2,3]
    m[15,14] = -1.0
    r[158] = m
    
    # O159:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,2] = -ci[2,2]*1
    m[10,3] = -ci[3,2]*1
    m[0,12] = -1.0*inv[1,0]
    m[1,12] = -1.0*inv[1,1]
    m[5,13] = -1.0
    m[14,15] = -1.0
    r[159] = m
    
    # O160:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,4] = 1
    m[5,7] = -1.0
    r[160] = m
    
    # O161:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,0] = ci[0,0]
    m[10,1] = ci[1,0]
    m[0,11] = -1.0*inv[1,0]
    m[1,11] = -1.0*inv[1,1]
    r[161] = m
    
    # O162:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,0] = ci[0,0]
    m[13,1] = ci[1,0]
    m[14,2] = ci[2,3]
    m[14,3] = ci[3,3]
    m[5,9] = -1.0
    m[15,11] = -1.0
    r[162] = m
    
    # O163:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O164:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[13,4] = 1.0
    m[5,8] = -1.0
    r[164] = m
    
    # O165:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,2] = -ci[2,3]*1
    m[10,3] = -ci[3,3]*1
    m[2,11] = inv[2,2]
    m[3,11] = inv[2,3]
    r[165] = m
    
    # O166:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = -ci[0,1]*1
    m[14,1] = -ci[1,1]*1
    m[13,2] = ci[2,2]
    m[13,3] = ci[3,2]
    m[5,10] = -1.0
    m[15,12] = 1
    r[166] = m
    
    # O167:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,5] = -1.0
    m[2,13] = inv[2,2]
    m[3,13] = inv[2,3]
    m[0,14] = -1.0*inv[1,0]
    m[1,14] = -1.0*inv[1,1]
    m[12,15] = 1
    r[167] = m
    
    # O168:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,0] = -ci[0,0]*1
    m[14,1] = -ci[1,0]*1
    m[0,7] = -1.0*inv[1,0]
    m[1,7] = -1.0*inv[1,1]
    r[168] = m
    
    # O169:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,4] = 1
    m[5,11] = 1.0
    r[169] = m
    
    # O170:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,0] = ci[0,0]
    m[12,1] = ci[1,0]
    m[0,9] = -1.0*inv[1,0]
    m[1,9] = -1.0*inv[1,1]
    r[170] = m
    
    # O171:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,2] = ci[2,3]
    m[10,3] = ci[3,3]
    m[2,11] = -1.0*inv[2,2]
    m[3,11] = -1.0*inv[2,3]
    r[171] = m
    
    # O172:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[14,2] = ci[2,3]
    m[14,3] = ci[3,3]
    m[12,4] = 1.0
    m[0,8] = -1.0*inv[1,0]
    m[1,8] = -1.0*inv[1,1]
    m[15,11] = -1.0
    r[172] = m
    
    # O173:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O174:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[12,2] = ci[2,2]
    m[12,3] = ci[3,2]
    m[14,5] = 1.0
    m[0,10] = -1.0*inv[1,0]
    m[1,10] = -1.0*inv[1,1]
    m[15,13] = -1.0
    r[174] = m
    
    # O175:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,0] = -ci[0,1]*1
    m[10,1] = -ci[1,1]*1
    m[2,12] = inv[2,2]
    m[3,12] = inv[2,3]
    m[5,14] = 1.0
    m[13,15] = -1.0
    r[175] = m
    
    # O176:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O177:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,2] = ci[2,2]
    m[10,3] = ci[3,2]
    m[0,12] = inv[1,0]
    m[1,12] = inv[1,1]
    m[5,13] = 1.0
    m[14,15] = 1
    r[177] = m
    
    # O178:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O179:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,5] = 1
    m[2,13] = -1.0*inv[2,2]
    m[3,13] = -1.0*inv[2,3]
    m[0,14] = inv[1,0]
    m[1,14] = inv[1,1]
    m[12,15] = -1.0
    r[179] = m
    
    # O180:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O181:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[10,0] = ci[0,1]
    m[10,1] = ci[1,1]
    m[2,12] = -1.0*inv[2,2]
    m[3,12] = -1.0*inv[2,3]
    m[5,14] = -1.0
    m[13,15] = 1
    r[181] = m
    
    # O182:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O183:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O184:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O185:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O186:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O187:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O188:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O189:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O190:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O191:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O192:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,0]*inv[0,0]
    m[1,0] = -ci[0,0]*inv[0,1]
    m[0,1] = -ci[1,0]*inv[0,0]
    m[1,1] = -ci[1,0]*inv[0,1]
    m[11,11] = -1.0
    m[13,13] = -1.0
    m[15,15] = -1.0
    r[192] = m
    
    # O193:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,2]*inv[0,0]
    m[1,2] = -ci[2,2]*inv[0,1]
    m[0,3] = -ci[3,2]*inv[0,0]
    m[1,3] = -ci[3,2]*inv[0,1]
    m[11,12] = 1
    r[193] = m
    
    # O194:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[0,0]
    m[1,0] = ci[0,0]*inv[0,1]
    m[0,1] = ci[1,0]*inv[0,0]
    m[1,1] = ci[1,0]*inv[0,1]
    m[11,11] = 1
    m[13,13] = 1
    m[15,15] = 1
    r[194] = m
    
    # O195:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,3]*inv[0,0]
    m[1,2] = -ci[2,3]*inv[0,1]
    m[0,3] = -ci[3,3]*inv[0,0]
    m[1,3] = -ci[3,3]*inv[0,1]
    m[13,14] = -1.0
    r[195] = m
    
    # O196:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,3]*inv[0,0]
    m[1,2] = ci[2,3]*inv[0,1]
    m[0,3] = ci[3,3]*inv[0,0]
    m[1,3] = ci[3,3]*inv[0,1]
    m[13,14] = 1
    r[196] = m
    
    # O197:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,1]*inv[0,0]
    m[1,0] = -ci[0,1]*inv[0,1]
    m[0,1] = -ci[1,1]*inv[0,0]
    m[1,1] = -ci[1,1]*inv[0,1]
    r[197] = m
    
    # O198:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,2]*inv[0,0]
    m[1,2] = ci[2,2]*inv[0,1]
    m[0,3] = ci[3,2]*inv[0,0]
    m[1,3] = ci[3,2]*inv[0,1]
    m[11,12] = -1.0
    r[198] = m
    
    # O199:  <Bp|  \hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[0,0]
    m[1,0] = ci[0,1]*inv[0,1]
    m[0,1] = ci[1,1]*inv[0,0]
    m[1,1] = ci[1,1]*inv[0,1]
    r[199] = m
    
    # O200:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O201:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,4] = -1.0
    m[11,11] = -1.0
    m[12,12] = -1.0
    m[15,15] = -1.0
    r[201] = m
    
    # O202:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O203:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,5] = -1.0
    r[203] = m
    
    # O204:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,4] = 1
    m[11,11] = 1.0
    m[12,12] = 1
    m[15,15] = 1.0
    r[204] = m
    
    # O205:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O206:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,5] = 1
    r[206] = m
    
    # O207:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O208:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,0]*inv[2,2]
    m[3,0] = -ci[0,0]*inv[2,3]
    m[2,1] = -ci[1,0]*inv[2,2]
    m[3,1] = -ci[1,0]*inv[2,3]
    m[12,11] = 1.0
    r[208] = m
    
    # O209:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,2]*inv[2,2]
    m[3,2] = -ci[2,2]*inv[2,3]
    m[2,3] = -ci[3,2]*inv[2,2]
    m[3,3] = -ci[3,2]*inv[2,3]
    m[12,12] = -1.0
    m[13,13] = -1.0
    m[15,15] = -1.0
    r[209] = m
    
    # O210:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,0]*inv[2,2]
    m[3,0] = ci[0,0]*inv[2,3]
    m[2,1] = ci[1,0]*inv[2,2]
    m[3,1] = ci[1,0]*inv[2,3]
    m[12,11] = -1.0
    r[210] = m
    
    # O211:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,3]*inv[2,2]
    m[3,2] = -ci[2,3]*inv[2,3]
    m[2,3] = -ci[3,3]*inv[2,2]
    m[3,3] = -ci[3,3]*inv[2,3]
    r[211] = m
    
    # O212:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,3]*inv[2,2]
    m[3,2] = ci[2,3]*inv[2,3]
    m[2,3] = ci[3,3]*inv[2,2]
    m[3,3] = ci[3,3]*inv[2,3]
    r[212] = m
    
    # O213:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,1]*inv[2,2]
    m[3,0] = -ci[0,1]*inv[2,3]
    m[2,1] = -ci[1,1]*inv[2,2]
    m[3,1] = -ci[1,1]*inv[2,3]
    m[13,14] = 1.0
    r[213] = m
    
    # O214:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,2]*inv[2,2]
    m[3,2] = ci[2,2]*inv[2,3]
    m[2,3] = ci[3,2]*inv[2,2]
    m[3,3] = ci[3,2]*inv[2,3]
    m[12,12] = 1.0
    m[13,13] = 1.0
    m[15,15] = 1
    r[214] = m
    
    # O215:  <Bp|  \hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,1]*inv[2,2]
    m[3,0] = ci[0,1]*inv[2,3]
    m[2,1] = ci[1,1]*inv[2,2]
    m[3,1] = ci[1,1]*inv[2,3]
    m[13,14] = -1.0
    r[215] = m
    
    # O216:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[0,0]
    m[1,0] = ci[0,0]*inv[0,1]
    m[0,1] = ci[1,0]*inv[0,0]
    m[1,1] = ci[1,0]*inv[0,1]
    m[11,11] = 1.0
    m[13,13] = 1.0
    m[15,15] = 1.0
    r[216] = m
    
    # O217:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,2]*inv[0,0]
    m[1,2] = ci[2,2]*inv[0,1]
    m[0,3] = ci[3,2]*inv[0,0]
    m[1,3] = ci[3,2]*inv[0,1]
    m[11,12] = -1.0
    r[217] = m
    
    # O218:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,0]*inv[0,0]
    m[1,0] = -ci[0,0]*inv[0,1]
    m[0,1] = -ci[1,0]*inv[0,0]
    m[1,1] = -ci[1,0]*inv[0,1]
    m[11,11] = -1.0
    m[13,13] = -1.0
    m[15,15] = -1.0
    r[218] = m
    
    # O219:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,3]*inv[0,0]
    m[1,2] = ci[2,3]*inv[0,1]
    m[0,3] = ci[3,3]*inv[0,0]
    m[1,3] = ci[3,3]*inv[0,1]
    m[13,14] = 1.0
    r[219] = m
    
    # O220:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,3]*inv[0,0]
    m[1,2] = -ci[2,3]*inv[0,1]
    m[0,3] = -ci[3,3]*inv[0,0]
    m[1,3] = -ci[3,3]*inv[0,1]
    m[13,14] = -1.0
    r[220] = m
    
    # O221:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[0,0]
    m[1,0] = ci[0,1]*inv[0,1]
    m[0,1] = ci[1,1]*inv[0,0]
    m[1,1] = ci[1,1]*inv[0,1]
    r[221] = m
    
    # O222:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,2]*inv[0,0]
    m[1,2] = -ci[2,2]*inv[0,1]
    m[0,3] = -ci[3,2]*inv[0,0]
    m[1,3] = -ci[3,2]*inv[0,1]
    m[11,12] = 1.0
    r[222] = m
    
    # O223:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,1]*inv[0,0]
    m[1,0] = -ci[0,1]*inv[0,1]
    m[0,1] = -ci[1,1]*inv[0,0]
    m[1,1] = -ci[1,1]*inv[0,1]
    r[223] = m
    
    # O224:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O225:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O226:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O227:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O228:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O229:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O230:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O231:  <Bp|  \hat{0}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O232:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,0]*inv[3,2]
    m[3,0] = -ci[0,0]*inv[3,3]
    m[2,1] = -ci[1,0]*inv[3,2]
    m[3,1] = -ci[1,0]*inv[3,3]
    m[14,13] = -1.0
    r[232] = m
    
    # O233:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,2]*inv[3,2]
    m[3,2] = -ci[2,2]*inv[3,3]
    m[2,3] = -ci[3,2]*inv[3,2]
    m[3,3] = -ci[3,2]*inv[3,3]
    r[233] = m
    
    # O234:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,0]*inv[3,2]
    m[3,0] = ci[0,0]*inv[3,3]
    m[2,1] = ci[1,0]*inv[3,2]
    m[3,1] = ci[1,0]*inv[3,3]
    m[14,13] = 1
    r[234] = m
    
    # O235:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,3]*inv[3,2]
    m[3,2] = -ci[2,3]*inv[3,3]
    m[2,3] = -ci[3,3]*inv[3,2]
    m[3,3] = -ci[3,3]*inv[3,3]
    m[11,11] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[235] = m
    
    # O236:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,3]*inv[3,2]
    m[3,2] = ci[2,3]*inv[3,3]
    m[2,3] = ci[3,3]*inv[3,2]
    m[3,3] = ci[3,3]*inv[3,3]
    m[11,11] = 1
    m[14,14] = 1
    m[15,15] = 1
    r[236] = m
    
    # O237:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,1]*inv[3,2]
    m[3,0] = -ci[0,1]*inv[3,3]
    m[2,1] = -ci[1,1]*inv[3,2]
    m[3,1] = -ci[1,1]*inv[3,3]
    m[11,12] = -1.0
    r[237] = m
    
    # O238:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,2]*inv[3,2]
    m[3,2] = ci[2,2]*inv[3,3]
    m[2,3] = ci[3,2]*inv[3,2]
    m[3,3] = ci[3,2]*inv[3,3]
    r[238] = m
    
    # O239:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,1]*inv[3,2]
    m[3,0] = ci[0,1]*inv[3,3]
    m[2,1] = ci[1,1]*inv[3,2]
    m[3,1] = ci[1,1]*inv[3,3]
    m[11,12] = 1
    r[239] = m
    
    # O240:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O241:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,4] = -1.0
    r[241] = m
    
    # O242:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O243:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,5] = -1.0
    m[13,13] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[243] = m
    
    # O244:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,4] = 1
    r[244] = m
    
    # O245:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O246:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,5] = 1
    m[13,13] = 1
    m[14,14] = 1.0
    m[15,15] = 1.0
    r[246] = m
    
    # O247:  <Bp|  \hat{0}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O248:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O249:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,4] = 1.0
    m[11,11] = 1
    m[12,12] = 1.0
    m[15,15] = 1
    r[249] = m
    
    # O250:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O251:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,5] = 1.0
    r[251] = m
    
    # O252:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,4] = -1.0
    m[11,11] = -1.0
    m[12,12] = -1.0
    m[15,15] = -1.0
    r[252] = m
    
    # O253:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O254:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[4,5] = -1.0
    r[254] = m
    
    # O255:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O256:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,0]*inv[3,2]
    m[3,0] = ci[0,0]*inv[3,3]
    m[2,1] = ci[1,0]*inv[3,2]
    m[3,1] = ci[1,0]*inv[3,3]
    m[14,13] = 1.0
    r[256] = m
    
    # O257:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,2]*inv[3,2]
    m[3,2] = ci[2,2]*inv[3,3]
    m[2,3] = ci[3,2]*inv[3,2]
    m[3,3] = ci[3,2]*inv[3,3]
    r[257] = m
    
    # O258:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,0]*inv[3,2]
    m[3,0] = -ci[0,0]*inv[3,3]
    m[2,1] = -ci[1,0]*inv[3,2]
    m[3,1] = -ci[1,0]*inv[3,3]
    m[14,13] = -1.0
    r[258] = m
    
    # O259:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,3]*inv[3,2]
    m[3,2] = ci[2,3]*inv[3,3]
    m[2,3] = ci[3,3]*inv[3,2]
    m[3,3] = ci[3,3]*inv[3,3]
    m[11,11] = 1.0
    m[14,14] = 1.0
    m[15,15] = 1.0
    r[259] = m
    
    # O260:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,3]*inv[3,2]
    m[3,2] = -ci[2,3]*inv[3,3]
    m[2,3] = -ci[3,3]*inv[3,2]
    m[3,3] = -ci[3,3]*inv[3,3]
    m[11,11] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[260] = m
    
    # O261:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,1]*inv[3,2]
    m[3,0] = ci[0,1]*inv[3,3]
    m[2,1] = ci[1,1]*inv[3,2]
    m[3,1] = ci[1,1]*inv[3,3]
    m[11,12] = 1.0
    r[261] = m
    
    # O262:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,2]*inv[3,2]
    m[3,2] = -ci[2,2]*inv[3,3]
    m[2,3] = -ci[3,2]*inv[3,2]
    m[3,3] = -ci[3,2]*inv[3,3]
    r[262] = m
    
    # O263:  <Bp|  \hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,1]*inv[3,2]
    m[3,0] = -ci[0,1]*inv[3,3]
    m[2,1] = -ci[1,1]*inv[3,2]
    m[3,1] = -ci[1,1]*inv[3,3]
    m[11,12] = -1.0
    r[263] = m
    
    # O264:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O265:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O266:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O267:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O268:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O269:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O270:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O271:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O272:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,0]*inv[1,0]
    m[1,0] = -ci[0,0]*inv[1,1]
    m[0,1] = -ci[1,0]*inv[1,0]
    m[1,1] = -ci[1,0]*inv[1,1]
    r[272] = m
    
    # O273:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,2]*inv[1,0]
    m[1,2] = -ci[2,2]*inv[1,1]
    m[0,3] = -ci[3,2]*inv[1,0]
    m[1,3] = -ci[3,2]*inv[1,1]
    m[14,13] = 1
    r[273] = m
    
    # O274:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[1,0]
    m[1,0] = ci[0,0]*inv[1,1]
    m[0,1] = ci[1,0]*inv[1,0]
    m[1,1] = ci[1,0]*inv[1,1]
    r[274] = m
    
    # O275:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,3]*inv[1,0]
    m[1,2] = -ci[2,3]*inv[1,1]
    m[0,3] = -ci[3,3]*inv[1,0]
    m[1,3] = -ci[3,3]*inv[1,1]
    m[12,11] = -1.0
    r[275] = m
    
    # O276:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,3]*inv[1,0]
    m[1,2] = ci[2,3]*inv[1,1]
    m[0,3] = ci[3,3]*inv[1,0]
    m[1,3] = ci[3,3]*inv[1,1]
    m[12,11] = 1
    r[276] = m
    
    # O277:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,1]*inv[1,0]
    m[1,0] = -ci[0,1]*inv[1,1]
    m[0,1] = -ci[1,1]*inv[1,0]
    m[1,1] = -ci[1,1]*inv[1,1]
    m[12,12] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[277] = m
    
    # O278:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,2]*inv[1,0]
    m[1,2] = ci[2,2]*inv[1,1]
    m[0,3] = ci[3,2]*inv[1,0]
    m[1,3] = ci[3,2]*inv[1,1]
    m[14,13] = -1.0
    r[278] = m
    
    # O279:  <Bp|  \hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[1,0]
    m[1,0] = ci[0,1]*inv[1,1]
    m[0,1] = ci[1,1]*inv[1,0]
    m[1,1] = ci[1,1]*inv[1,1]
    m[12,12] = 1
    m[14,14] = 1
    m[15,15] = 1
    r[279] = m
    
    # O280:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,0]*inv[2,2]
    m[3,0] = ci[0,0]*inv[2,3]
    m[2,1] = ci[1,0]*inv[2,2]
    m[3,1] = ci[1,0]*inv[2,3]
    m[12,11] = -1.0
    r[280] = m
    
    # O281:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,2]*inv[2,2]
    m[3,2] = ci[2,2]*inv[2,3]
    m[2,3] = ci[3,2]*inv[2,2]
    m[3,3] = ci[3,2]*inv[2,3]
    m[12,12] = 1
    m[13,13] = 1
    m[15,15] = 1.0
    r[281] = m
    
    # O282:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,0]*inv[2,2]
    m[3,0] = -ci[0,0]*inv[2,3]
    m[2,1] = -ci[1,0]*inv[2,2]
    m[3,1] = -ci[1,0]*inv[2,3]
    m[12,11] = 1
    r[282] = m
    
    # O283:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = ci[2,3]*inv[2,2]
    m[3,2] = ci[2,3]*inv[2,3]
    m[2,3] = ci[3,3]*inv[2,2]
    m[3,3] = ci[3,3]*inv[2,3]
    r[283] = m
    
    # O284:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,3]*inv[2,2]
    m[3,2] = -ci[2,3]*inv[2,3]
    m[2,3] = -ci[3,3]*inv[2,2]
    m[3,3] = -ci[3,3]*inv[2,3]
    r[284] = m
    
    # O285:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = ci[0,1]*inv[2,2]
    m[3,0] = ci[0,1]*inv[2,3]
    m[2,1] = ci[1,1]*inv[2,2]
    m[3,1] = ci[1,1]*inv[2,3]
    m[13,14] = -1.0
    r[285] = m
    
    # O286:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,2] = -ci[2,2]*inv[2,2]
    m[3,2] = -ci[2,2]*inv[2,3]
    m[2,3] = -ci[3,2]*inv[2,2]
    m[3,3] = -ci[3,2]*inv[2,3]
    m[12,12] = -1.0
    m[13,13] = -1.0
    m[15,15] = -1.0
    r[286] = m
    
    # O287:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[2,0] = -ci[0,1]*inv[2,2]
    m[3,0] = -ci[0,1]*inv[2,3]
    m[2,1] = -ci[1,1]*inv[2,2]
    m[3,1] = -ci[1,1]*inv[2,3]
    m[13,14] = 1
    r[287] = m
    
    # O288:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O289:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,4] = 1.0
    r[289] = m
    
    # O290:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O291:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,5] = 1.0
    m[13,13] = 1.0
    m[14,14] = 1
    m[15,15] = 1
    r[291] = m
    
    # O292:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,4] = -1.0
    r[292] = m
    
    # O293:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O294:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[5,5] = -1.0
    m[13,13] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[294] = m
    
    # O295:  <Bp|  \hat{1}_{\beta}^{+}\hat{0}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O296:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,0]*inv[1,0]
    m[1,0] = ci[0,0]*inv[1,1]
    m[0,1] = ci[1,0]*inv[1,0]
    m[1,1] = ci[1,0]*inv[1,1]
    r[296] = m
    
    # O297:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,2]*inv[1,0]
    m[1,2] = ci[2,2]*inv[1,1]
    m[0,3] = ci[3,2]*inv[1,0]
    m[1,3] = ci[3,2]*inv[1,1]
    m[14,13] = -1.0
    r[297] = m
    
    # O298:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,0]*inv[1,0]
    m[1,0] = -ci[0,0]*inv[1,1]
    m[0,1] = -ci[1,0]*inv[1,0]
    m[1,1] = -ci[1,0]*inv[1,1]
    r[298] = m
    
    # O299:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = ci[2,3]*inv[1,0]
    m[1,2] = ci[2,3]*inv[1,1]
    m[0,3] = ci[3,3]*inv[1,0]
    m[1,3] = ci[3,3]*inv[1,1]
    m[12,11] = 1.0
    r[299] = m
    
    # O300:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,3]*inv[1,0]
    m[1,2] = -ci[2,3]*inv[1,1]
    m[0,3] = -ci[3,3]*inv[1,0]
    m[1,3] = -ci[3,3]*inv[1,1]
    m[12,11] = -1.0
    r[300] = m
    
    # O301:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = ci[0,1]*inv[1,0]
    m[1,0] = ci[0,1]*inv[1,1]
    m[0,1] = ci[1,1]*inv[1,0]
    m[1,1] = ci[1,1]*inv[1,1]
    m[12,12] = 1.0
    m[14,14] = 1.0
    m[15,15] = 1.0
    r[301] = m
    
    # O302:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,2] = -ci[2,2]*inv[1,0]
    m[1,2] = -ci[2,2]*inv[1,1]
    m[0,3] = -ci[3,2]*inv[1,0]
    m[1,3] = -ci[3,2]*inv[1,1]
    m[14,13] = 1.0
    r[302] = m
    
    # O303:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    m = numpy.zeros(shape=(nbs,nbs)) 
    m[0,0] = -ci[0,1]*inv[1,0]
    m[1,0] = -ci[0,1]*inv[1,1]
    m[0,1] = -ci[1,1]*inv[1,0]
    m[1,1] = -ci[1,1]*inv[1,1]
    m[12,12] = -1.0
    m[14,14] = -1.0
    m[15,15] = -1.0
    r[303] = m
    
    # O304:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O305:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O306:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O307:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{0}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    # O308:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{0}_{\alpha}^{-}  |Bq> = 
    
    # O309:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\alpha}^{-}\hat{1}_{\alpha}^{-}  |Bq> = 
    
    # O310:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{0}_{\beta}^{-}  |Bq> = 
    
    # O311:  <Bp|  \hat{1}_{\beta}^{+}\hat{1}_{\beta}^{+}\hat{1}_{\beta}^{-}\hat{1}_{\beta}^{-}  |Bq> = 
    
    
    return r
    
