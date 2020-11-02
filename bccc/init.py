#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   init.py
#            Des:   init
#           Mail:   qingchun720@foxmail.com
#   Created Time:   Wed 16:00:27 Jul/17 2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


from itertools import combinations,permutations
from bccc.pub import sign


tsd = { }
def init(np, nfrozen=0, ntuply=2):
    t = {}

    for n in range(1,ntuply+1):
        for bl in combinations(range(nfrozen,np), r=n):
            for sl in tsd[n]:
                key = tuple(zip(bl,sl)); value = 1.0
                t[key] = value
                for _key in permutations(key):
                    t[_key] = sign([b for b,s in _key if 6<s<15])*value

    return t






