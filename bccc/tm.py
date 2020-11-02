#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
# 
# Author: Qingchun Wang @ NJU 
# E-mail: qingchun720@foxmail.com 
# 


from bccc.pub import sign


def tm(np, nt=4, nc=0):
    tmp = {}
    
    for A in range(nc, np):
        _t = {}
        _t[tuple(sorted([(A,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
        tmp[tuple(sorted([(A,1),]))] = _t
        
        _t = {}
        _t[tuple(sorted([(A,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
        tmp[tuple(sorted([(A,2),]))] = _t
        
        _t = {}
        _t[tuple(sorted([(A,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
        tmp[tuple(sorted([(A,3),]))] = _t
        
    for A in range(nc, np):
        for B in range(nc, np):
            if B in {A}: continue
            _t = {}
            _t[tuple(sorted([(A,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,1),])))
            _t[tuple(sorted([(B,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,1),])))
            _t[tuple(sorted([(A,1),(B,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,1),(B,1),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,2),])))
            _t[tuple(sorted([(B,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,1),])))
            _t[tuple(sorted([(A,1),(B,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,1),(B,2),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,1),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,3),])))
            _t[tuple(sorted([(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,1),])))
            _t[tuple(sorted([(A,1),(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,1),(B,3),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,2),])))
            _t[tuple(sorted([(B,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,2),])))
            _t[tuple(sorted([(A,2),(B,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,2),(B,2),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,2),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,3),])))
            _t[tuple(sorted([(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,2),])))
            _t[tuple(sorted([(A,2),(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,2),(B,3),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(B,3),])))
            _t[tuple(sorted([(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([(A,3),])))
            _t[tuple(sorted([(A,3),(B,3),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,3),(B,3),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,4),(B,5),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,4),(B,5),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,6),(B,15),]))] = (sign(sorted([])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,6),(B,15),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,7),(B,13),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,7),(B,13),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,7),(B,14),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,7),(B,14),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,8),(B,13),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,8),(B,13),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,8),(B,14),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,8),(B,14),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,9),(B,11),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,9),(B,11),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,9),(B,12),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,9),(B,12),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,10),(B,11),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,10),(B,11),]))] = _t
            
            _t = {}
            _t[tuple(sorted([(A,10),(B,12),]))] = (sign(sorted([A,B])+sorted([])), tuple(sorted([])))
            tmp[tuple(sorted([(A,10),(B,12),]))] = _t
            
    return tmp
    
