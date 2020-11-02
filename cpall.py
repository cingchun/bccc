#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   cpp.py
#            Des:   cpp
#           Mail:   qingchun720@foxmail.com
#   Created Time:   15:36 六月-28/2019 
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'
    
    
import os, sys
from multiprocessing import Pool


node = None; pathl = sys.argv[1:]
def scp(file):
    cwd = os.getcwd()
    os.system(F'scp {path}/{file} {node}:{cwd}/{path}')


for path in pathl:
    if os.path.isdir(path):
        fl = os.listdir(path)
        pool = Pool()
        pool.map(scp, fl)
        pool.close()
        pool.join()
    elif os.path.splitext(path) == '.py':
        os.system(F'python -m py_compile {path}')
    else:
        print(F'Warning: {path} isnot path or .py file !')
        pass


    

    
    