#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   compall.py
#            Des:   compall
#           Mail:   qingchun720@foxmail.com
#   Created Time:   16:02 六月-28/2019 
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'



import os, sys
from multiprocessing import Pool
Path = os.path
cmd = os.system


def comp(pathl, cwd='./'):
    # for path in pathl:
    #     if os.path.isdir(path):
    #         fl = os.listdir(path)
    #         Pool().map(comp, fl)
    #     elif os.path.splitext(path) == '.py':
    #         os.system(F'python -m py_compile {path}')
    #     else:
    #         # raise RuntimeError(F'Error: {path} isnot path or .py file !')
    #         print(F'Warning: {path} isnot path or .py file !')
    #         pass
    pl,fl = [],[]
    plapp = pl.append; flapp=fl.append
    for p_f in pathl:
        p_f = F'{cwd}{p_f}'
        if Path.isdir(p_f):
            if p_f[-1] != '/': p_f += '/'
            plapp(p_f)
        elif Path.splitext(p_f)[1]=='.py': flapp(F'python -m py_compile {p_f}')

    Pool().map(cmd, fl)
    for p in pl: comp(os.listdir(p), p)




if __name__ == '__main__':
    comp(sys.argv[1:])






