#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   derivform.py
#            Des:   project ket
#           Mail:   qingchun720@foxmail.com
#   Created Time:   10:49 二月-25/2019
#
'''
Note: implement mode
        (1) Symbol + list
        (2) Matrix show Expression (row +,  col *)
        (3) Matrix show Expression (row +,  col *, x0 is coeff)
        (4) Binary Structure show Expression (left +, right *)
        (5) Tree show Expression  (operator node, symbol leaf)
        (6) Tree + Structure + Expression
'''

# Variable(or function name) mean:
#   _... :  temp, backup varible

#   s... :  str varible
#   n... :  number
#    ...l:  list, array, set, tuple
#    ...d:  dict
#   e... :  element in ...l
#   k... :  key in ...d
#   a... :  array in Numpy
#   m... :  matrix in Numpy

#   o... :  orbital ...
#   b... :  block ...
#   w... :  wave ...
#   f... :  functor ...
#   c... :  coefficient ...
#   h... :  hamiltonian functor...
#   t... :  T functor ...

'''
TODO:
    1. operator include two mode: symbol and functional 
    2. bfl and 440 repf are expressed as dict
    3. mh (2D dict)
    4. reduce memory

    __hash__ => opt repfl, mh
    *symb => Expr.__init__()
    ws -> wave in detat()
    __str__ -> latex; mean -> __str__(__repr__)
    different __repr__ with __str__

    copy(), deepcopy(), __copy__(), __deepcopy__()
    __bool__() -> Symb0; __bool__() -> Expr0
    str '0' in Block.l -> number 0 and similar 
'''


__version__ = '3.1.26'
__author__ = 'Qingchun Wang'




import os, sys, shutil, subprocess
sys.setrecursionlimit(1000000)
sys.setswitchinterval(10000)
from copy import copy,deepcopy
from datetime import datetime
import logging

from functools import reduce
from operator import mul,add,sub,truediv
from itertools import permutations, combinations, product, chain
from collections import Iterable

import math, random
from fractions import Fraction as Frac
# import numpy, scipy
# from sympy import Symbol

import pickle
from re import compile, findall
numstr = compile(r'^[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?$') # match integer, fraction, exponent
# numstr = compile(r'[-+]?(\b[0-9]+(\.[0-9]+)?|\.[0-9]+)([eE][-+]?[0-9]+\b)?') # eg. w123.345w
# findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", 'A1.45aa, b5., B6.45, F-2.e2ddf')

from multiprocessing import Pool, Manager
# import mpi4py

from bccc.pub import index,sign,prodl,sum_list




# operator
#     include function and operator in writing form
#         function:        ?( , )
#         character:       ( )?( )
#     include unary, binary and ternary from the point of view of the number of operations
#         unary operator belong to function
#         binary operator partly belong to functional, and part belong to symbol
#         ternary operator belong to function
''' unary operator '''
unarycharl = {'+1', '-1',
              '0+', '0-'}
signoperl = {'0-', '0+'}
unaryfuncl = {'abs', 'fabs',
              'ceil', 'floor', 'round',
              'exp', 'sqrt', 'lg', 'ln',
              'sin', 'cos', 'asin', 'acos', 'tan', 'cot', 'atan', 'acot',
              'int', 'float', 'str'}
unaryoperl = unarycharl|unaryfuncl
#             0+, 0-:    sign operator in mathematics
#             +1, -1:    ++, -- in computer science
# ceil, floor, round:    Upward, downward integrate, and rounding of value
#    int, float, str:    int(), float(), str() in code
def isunary(oper):
    '''
    is unary operator ?

    :param oper:

    :return:
    '''

    return oper in unaryoperl
''' binary operator '''
binarycharl = {'=',
               '||', '&&', '!',
               '||=', '&&=', '!=',
               '==', '<>', '<', '>', '>=', '<=',
               '|', '^', '&',
               '|=', '^=', '&=',
               '<<', '>>',
               '<<=', '>>=',
               '+', '-',
               '+=', '-=',
               '*', '/', '//', '%', '++',
               '*=', '/=', '//=', '%=', '++=',
               '~',
               '~='}
swapableoperl = {'+', '*'}
binaryfuncl = {',', ':', '-?', '**', 'pow', '%%', 'log'}
binaryoperl = binarycharl|binaryfuncl
#       !=:    not and equal to, like +=
#       <>:    not eqaul
#       //:    floor divid, with no remainder
#        %:    mod, the remainder
#       ++:    times e.g. x++3 = x+x+x = 3x
#  **, pow:    power function
#  %%, log:    log function
#        ,:    each one in sigma or pi cycle, e.g. i, j
#        ::    ben -> end in sigma or pi cycle,  e.g. 0->n
#       -?:    function for block swap
def isbinary(oper):
    '''
    is binary operator ?

    :param oper:

    :return:
    '''

    return oper in binaryoperl
''' ternary operator '''
def sigma(expr=None, ben=0, end=None):
    '''
    sigma summation

    run function, how to compute
    :param expr:
    :param ben:
    :param end:

    :return:
    '''

    raise NotImplementedError('TODO: sigma() for summation')
def pi(expr=None, ben=0, end=None):
    '''
    pi(II) product

    run function, how to compute
    :param expr:
    :param ben:
    :param end:

    :return:
    '''

    raise NotImplementedError('TODO: pi() for II product')
def proj(fun, bra, ket):
    '''
    proj() for <i| H |j>

    run function, how to compute
    :param bra:
    :param ket:
    :param fun:

    :return:
    '''

    raise NotImplementedError('TODO: proj() for <i| H |j>')
def condit(case, out1, out2):
    '''
    condit() for if-else

    run function, how to compute
    :param case:
    :param out1:
    :param out2:

    :return:
    '''

    return out1 if case else out2
ternaryoperl = {'?:', '+++', '***', 'proj'}
#   ?::    if-then-else, [bool, True-Expr, False-Expr]
#  +++:    sum sigma, [ben, end, expr]
#  ***:    prod pi, [ben, end, expr]
# proj:    put project, [bra, ket, functor]
def isternary(oper):
    '''
    is ternary operator ?

    :param oper:

    :return:
    '''

    return oper in ternaryoperl

''' priority of operator '''
priorityl={    1: {'lambda'},

               2: {'||', '&&', '!'},
               3: {'==', '<>', '<', '>', '>=', '<='},
               4: {'|', '^', '&'},
               5: {'<<', '>>'},
               6: {'+', '-', '+1', '-1'},
               7: {'*', '/', '//', '%', '++'},
               8: {'0+', '0-'},
               9: {'~'},

               10: {',', ':', '-?', '**', 'pow', '%%', 'log',

                'abs', 'fabs', 'ceil', 'floor', 'round', 'exp', 'sqrt', 'lg', 'ln', 'sin',
                'cos', 'asin', 'acos', 'tan', 'cot', 'atan', 'acot', 'int', 'float', 'str',

                '?:', '+++', '***', 'proj'}
           }
def prioritize(oper):
    '''
    prioritize of an operator

    :param oper:

    :return:
    '''

    for i in priorityl:
        if oper in priorityl[i]: return i
    else:
        raise ValueError(F'unrecognized operator {oper} in prioritize()')

''' functional map of operator '''
opermap={
    '+': add,
    '-': sub,
    '*': mul,
    '/': truediv,
  '+++': sigma,
  '***': pi}
def nary(oper):
    '''
    n-nary of operator ?

    :param oper:

    :return:
    '''

    if oper in unaryoperl: return 1
    elif oper in binaryoperl: return 2
    elif oper in ternaryoperl: return 3
    else:
        raise ValueError(F'unrecognized operator {oper} in nary()')




# number
#     int and float
#     include number-Symb(Num), and number-Expr
#     only number-Symb(Num), and number-Expr can be be converted into number
intorfloatl = {int, float}
def isnum(arg):
    '''
    is number (int or float)?

    :param arg:

    :return:
    '''

    return arg.__class__ in intorfloatl
def is0(arg):
    '''
    is number 0 ?

    :param arg:

    :return:
    '''

    return arg.__class__ in intorfloatl and abs(arg)<1.e-12
def is1(arg):
    '''
    is number 1 ?

    :param arg:

    :return:
    '''

    return arg.__class__ in intorfloatl and abs(arg-1.0)<1.e-12
def isnumstr(arg):
    '''
    is str-number ?

    :param arg:

    :return: Mathch object
    '''

    if isinstance(arg, str): return numstr.match(arg)
    else: return False
def str2num(arg):
    '''
    str-number -> number

          Note:     only for number str (eg. '231323.93493') -> number
                    is not Frac (eg. 1/4) -> number
    :param arg:

    :return:
    '''

    return float(arg)
def isstr0(arg):
    '''
    is str-number 0 ?

    :param arg:

    :return:
    '''

    return isinstance(arg, str) and numstr.match(arg) and abs(float(arg))<1.e-12
def isstr1(arg):
    '''
    is str-number 1 ?

    :param arg:

    :return:
    '''

    return isinstance(arg, str) and numstr.match(arg) and abs(float(arg)-1)<1.e-12

def num2str(arg):
    '''
    number (int, float) -> latex

        Note :  only for number -> Frac
                is not str(arg)
    :param arg:

    :return:
    '''

    arg = Frac(arg); numerator, denominator = arg.numerator, arg.denominator
    if denominator is 1 and numerator is 1: return '1'
    elif denominator is 1: return F'{numerator}'
    else: return F'\\frac{{{numerator}}}{{{denominator}}}'
def num2formt(arg):
    '''
    number (int, float) -> str

    :param arg:

    :return:
    '''

    return arg




# Two base class:
#     Symb and Expr
#     For convenience, '' stands for None in all of the following types
'''
在之前版本中，实现了Symb, Expr对 int/float, str 等类型的兼容
但是，在实际中，多数都是 Symb, Expr 的计算，而非 int/float, str 的计算
对 int/float, str 兼容，增加了对 type 的判断，会使行计算做很多无效计算，得不偿失
因此，在这个版本，删除对 int/float, str 兼容
'''
class Symb(object):
    '''
    Class: Mathematical Symbol

           in which every term can be number(int and float), str, Symb or Expr

    NOTE:  c is a very important parameter when initial
           whose default type is str, of course, can be Symb or Expr
           when c = '', 0, None or [],(),{} it is regarded as Symb 0
           None or '':  unkown

    :param c:       script@chief, central, indicates the attributes of Symb
                          S: str                  'xxx'; and support +-* for str class
                          N: number               0, 1.2; and support +-*/ for int and float class
                            Note: The two basic data type(Str, Number), only Number can be add to Symb and Expr

                          O: Orbital              \phi
                          B: Block                \Lambda
                          W: Wave                 \Phi

                          C: Coefficient          C
                          D: Delta                \delta

                          F: Functor              \hat
                          I: Integral             I

                          T: T operator
                          H: H operator
                            Note: T and H can be regarded as subclass of F
    :param l:       script@low,
                    Orbital,            p,q,r,s:   any orbital
                                        i,j,k,l:   occupied orbital
                                        a,b,c,d:   virtual or unoccuped orbital
                    block or geminal,   P,Q,R,S:   any block
                                        I,J,K,L:   occupied block
                                        A,B,C,D:   virtual or unoccuped block
                    wavefunction,          \Phi:
    :param f:       script@foot: state(x, y, z) for Orbital, Block, Wave, Functor and Integral
                    For Orbital,
                                               0:  space
                                               1:  alpha  \alpha
                                              -1:  beta   \beta
                    For Block,
                                               0:  ground state
                                       1, 2, ...:  excited state
                    For Wave,
                                               0:  ground state
                                       1, 2, ...:  excited state
                    For Functor,
                                               0:  ^
                                               1:  +
                                              -1:  -
                    For Integral,
                                               0:  core-hamiltonian integral
                                               1:  fock-hamiltonian integral
                                               2:  coulomb integral
                                               3:  exchange integral
                                               4:  (coclomb - exchange) integral
    '''

    def __init__(self, c=None, l=None, f=None): self.c, self.l, self.f = c, l, f

    # 在这个版本中，删除对0,1以及 特殊运算 的处理
    # 比如：0+s=s+0=s-0=s*1=1*s=s/1=s; (0+a)+(1*b)=(1*a)+(1*b)=(1*a)-(0+b)=a+b
    # 支持 0,1 特殊处理，计算变得杂，实际符号运算中，这并不常见
    # 但是，特别是 1*s=s 保留, 0-e=-e
    def __add__(self, other):
        '''
        +: self + other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            tmp = Expr('+', self, other)
            other.root = tmp
            return tmp
        elif isinstance(other, Symb): return Expr('+', self, other)
        else:
            if abs(other)<1.0e-12: return self
            else: return Expr('+', self, other)
    def __radd__(self, other):
        '''
        +: other + self

        Note: only when type(other)=Num
        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12: return self
            else: return Expr('+', other, self)
        except TypeError: # Symb or Expr
            tmp = Expr('+', other, self)
            if isinstance(other, Expr): other.root = tmp
            return tmp
    def __sub__(self, other):
        '''
        -: self - other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            tmp = Expr('-', self, other)
            other.root = tmp
            return tmp
        elif isinstance(other, Symb): return Expr('-', self, other)
        else:
            if abs(other)<1.0e-12: return self
            else: return Expr('-', self, other)
    def __rsub__(self, other):
        '''
        -: other - self

        Note: only when type(other)=Num and other=0
        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12: return Expr('0-', self)
            else: return Expr('-', other, self)
        except TypeError:
            tmp = Expr('-', other, self)
            if isinstance(other, Expr): other.root = tmp
            return tmp
    def __mul__(self, other):
        '''
        *: self * other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            tmp = Expr('*', self, other)
            other.root = tmp
            return tmp
        elif isinstance(other, Symb): return Expr('*', self, other)
        else:
            if abs(other-1)<1.0e-12: return self
            elif abs(other)<1.0e-12: return 0.0
            else: return Expr('*', self, other)
    def __rmul__(self, other):
        '''
        *: other * self

        Note: only when type(other)=Num and other=1
        :param other:

        :return:
        '''

        try:
            if abs(other-1)<1.0e-12: return self
            elif abs(other)<1.0e-12: return 0.0
            else: return Expr('*', other, self)
        except TypeError:
            tmp = Expr('*', other, self)
            if isinstance(other, Expr): other.root = tmp
            return tmp
    def __truediv__(self, other):
        '''
        /: self / other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            tmp = Expr('/', self, other)
            other.root = tmp
            return tmp
        elif isinstance(other, Symb): return Expr('/', self, other)
        else:
            if abs(other-1)<1.0e-12: return self
            elif abs(other)<1.0e-12: raise ZeroDivisionError('Error:  Symb / 0')
            else: return Expr('/', self, other)
    def __rtruediv__(self, other):
        '''
        /: other / self

        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12: return 0.0
            else: return Expr('/', other, self)
        except TypeError:
            tmp = Expr('/', other, self)
            if isinstance(other, Expr): other.root = tmp
            return tmp
class Expr(object):
    '''
    Class: Mathematical Expression,

            which may a term(even only one Symb), many term
           Tree structure (unary, binary or ternary tree)

    :param mean:

    :param nlayer:      int
    :param nleaf:       int

    :param root:        upper layer Expr
    :param oper:        operator
    :param leafl:       lower layer, list

    In the way, ce can be treated as operator *
                ep can be treated as operator **
    '''

    def __init__(self, oper, *symbl):
        '''
        init: *symb, oper

        :param symb:
        :param oper:
        '''

        self.root, self.oper, self.leafl = None, oper, symbl
        # for leaf in symbl:
        #     if isinstance(leaf, Expr): leaf.root = self

    def __str__(self):
        root, oper = self.root, self.oper; leafl = tuple(F'{e!s}' for e in self.leafl)
        if oper in binaryoperl:
            if oper in binarycharl:
                expr = F'{leafl[0]!s}{oper}{leafl[1]!s}'
                if root:
                    rootoper = root.oper; poper, prootoper = prioritize(oper), prioritize(rootoper)
                    if poper>prootoper or poper is prootoper and rootoper in swapableoperl: return expr
                    else: return F'({expr})'
                else: return expr
            else: # binaryfuncl
                if oper=='-?': return F"(-1)^{{P({','.join(self.leafl[0])})}}"
                else:
                    raise NotImplementedError(F'TODO: Expr.__str__() for operator {oper}')
        elif oper in unaryoperl:
            if oper in unarycharl:
                if oper=='0-': return F'-{leafl[0]}'
                elif oper=='0+': return leafl[0]
                else: return F'({leafl[0]}{oper})'
            else: return F'{oper}({leafl[0]})'
        elif oper in ternaryoperl:  # debugging ...
            if oper=='+++':
                return F"\\displaystyle\\sum_{{{leafl[1]}}}^{{{leafl[2]}}} {leafl[0]}"
            elif oper=='***':
                return F"\\prod_{{{leafl[1]}}}^{{{leafl[2]}}} {leafl[0]}"
            else:
                raise NotImplementedError(F'TODO: Expr.__str__() for operator {oper}')
        else:
            raise ValueError(F'unrecognized operator {oper} in Expr.__str__()')
    __repr__ = __str__
    def __format__(self, format_spec):
        root, oper, leafl = self.root, self.oper, tuple(F'{leaf}' for leaf in self.leafl)
        if oper in binaryoperl:
            if oper in binarycharl:
                expr = F'{leafl[0]}{oper}{leafl[1]}'
                if root:
                    rootoper = root.oper; poper, prootoper = prioritize(oper), prioritize(rootoper)
                    if poper>prootoper or poper is prootoper and rootoper in swapableoperl: return expr
                    else: return F'({expr})'
                else: return expr
            else: # binaryfuncl
                if oper=='-?':
                    # leafl0 = self.leafl[0]
                    # if isinstance(leafl0[0],str): return F"sign(lnew=[{','.join(leafl0)}])"
                    # else: return F'''sign(lnew={'+'.join(F"sorted([{','.join(bl)}])" for bl in leafl0)})'''
                    return F"sign(lnew=[{','.join(self.leafl[0])}])"
                else:
                    raise NotImplementedError(F'TODO: Expr.__format__() for operator {oper}')
        elif oper in unaryoperl:
            if oper in unarycharl:
                if oper=='0-': return F'-{leafl[0]}'
                elif oper=='0+': return leafl[0]
                else: return F'({leafl[0]}{oper})'
            else: return F'math.{oper}({leafl[0]})'
        else:
            raise ValueError(F'unrecognized operator {oper} in Expr.__format__()')

    # 在这个版本中，删除对0,1以及 特殊运算 的处理
    # 比如：0+s=s+0=s-0=s*1=1*s=s/1=s; (0+a)+(1*b)=(1*a)+(1*b)=(1*a)-(0+b)=a+b
    # 支持 0,1 特殊处理，计算变得杂，实际符号运算中，这并不常见
    # 但是，特别是 1*s=s 保留, 0-e=-e
    # 有效的运算包括：s,e +-*/ N,s,e 和 N +-*/ s,e
    def __eq__(self, other):
        '''

        :param other:
        :return:
        '''

        return F'{self}' == F'{other}'

    def __add__(self, other):
        '''
        + plus: self + other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            if other.oper!='0-': tmp = Expr('+', self, other)
            else: tmp = Expr('-', self, other.leafl[0])
            try: self.root = tmp.leafl[1].root = tmp
            except AttributeError: self.root = tmp
            return tmp
        elif isinstance(other, Symb):
            tmp = Expr('+', self, other)
            self.root = tmp
            return tmp
        else:
            if abs(other)<1.0e-12: return self
            else:
                tmp = Expr('+', self, other)
                self.root = tmp
                return tmp
    def __radd__(self, other):
        '''
        + plus: other + self

        # Note: only when type(other)=Num
        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12: return self
            else:
                tmp = Expr('+', other, self)
                self.root = tmp
                return tmp
        except TypeError:
            tmp = Expr('+', other, self)
            try: self.root = other.root = tmp
            except AttributeError: self.root = tmp
            return tmp
    def __sub__(self, other):
        '''
        - decreate: self - other

        :param other:

        :return:
        '''

        if isinstance(other, Expr):
            if other.oper!='0-': tmp = Expr('-', self, other)
            else: tmp = Expr('+', self, other.leafl[0])
            try: self.root = tmp.leafl[1].root = tmp
            except AttributeError: self.root = tmp
            return tmp
        elif isinstance(other, Symb):
            tmp = Expr('-', self, other)
            self.root = tmp
            return tmp
        else:
            if abs(other)<1.0e-12: return self
            else:
                tmp = Expr('-', self, other)
                self.root = tmp
                return tmp
    def __rsub__(self, other):
        '''
        - decreate: other - self

        Note: only when type(other)=Num and other=0
        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12:
                if self.oper!='0-':
                    tmp = Expr('0-', self)
                    self.root = tmp
                    return tmp
                else: # 0 - (0 - E) = E
                    leafl0 = self.leafl[0]
                    if isinstance(leafl0, Expr): leafl0.root = None
                    return leafl0
            else:
                tmp = Expr('-', other, self)
                self.root = tmp
                return tmp
        except TypeError:
            tmp = Expr('-', other, self)
            try: self.root = other.root = tmp
            except AttributeError: self.root = tmp
            return tmp
    def __mul__(self, other):
        '''
        *: self * other

        :param other:

        :return:
        '''

        tmp = Expr('*', self, other)
        try: self.root = other.root = tmp
        except AttributeError: self.root = tmp
        return tmp
    def __rmul__(self, other):
        '''
        *: other * self

        Note: only when type(other)=Num and other=1
        :param other:

        :return:
        '''

        try:
            if abs(other-1)<1.0e-12: return self
            elif abs(other)<1.0e-12: return 0.0
            else:
                tmp = Expr('*', other, self)
                self.root = tmp
                return tmp
        except TypeError:
            tmp = Expr('*', other, self)
            try: self.root = other.root = tmp
            except AttributeError: self.root = tmp
            return tmp
    def __truediv__(self, other):
        '''
        /: self / other

        :param other:

        :return:
        '''

        tmp = Expr('/', self, other)
        try: self.root = other.root = tmp
        except AttributeError: self.root = tmp
        return tmp
    def __rtruediv__(self, other):
        '''
        /: other / self

        :param other:

        :return:
        '''

        try:
            if abs(other)<1.0e-12: return 0.0
            else:
                tmp = Expr('/', other, self)
                self.root = tmp
                return tmp
        except TypeError:
            tmp = Expr('/', other, self)
            try: self.root = other.root = tmp
            except AttributeError: self.root = tmp
            return tmp

class Coeff(Symb):
    '''
    Coeff class

    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='C', l=symbl, f=mode) # default [] has bug

    def __str__(self): return F"{{{self.f}}}_{{{','.join(F'{e!s}' for e in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        symbl,mode = self.l,self.f; nsymb = len(symbl)
        if mode is 't': return F"t[({''.join(F'({b.l},{b.f}),' for b in symbl)})]"
        elif mode is 'r':
            try: return F'r[{symbl[0].l}][{symbl[1]}][{symbl[2]},{symbl[3]}]'
            except AttributeError: return F'r[{symbl[0]}][{symbl[1]},{symbl[2]}]'
        else: return F"{self.f}[{symbl[0]},{symbl[1]}]"

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"C{self.f}{self.l}" == F"C{other.f}{other.l}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"C{self.f}{self.l}")
class Array(Symb):
    '''
    Array class

    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='A', l=symbl, f=mode)

    def __str__(self): return F"{{{self.f}}}_{{{','.join(F'{e!s}' for e in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Array.__format__()')

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"A{self.f}{self.l}" == F"A{other.f}{other.l}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"A{self.f}{self.l}")
class Integral(Symb):
    '''
    Integral class

    Note:  Integral indiscriminate state (orbital ans block)
    :param symbl: Orbital set
    :param mode:   '':  unkown and default
                    0:  core-hamiltonian integral
                    1:  fock-hamiltonian integral
                    2:  coulomb integral
                    3:  exchange integral
                    4:  (coclomb - exchange) integral
    '''

    def __init__(self,  symbl=(), mode=2, upl=()):
        Symb.__init__(self, c='I', l=symbl, f=mode)
        self.upl = upl

    def __str__(self):
        f = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        sl = {i: F'{symb.l}^{{{upl[i]}}}' for i,symb in enumerate(self.l)}
        if f is 2: return F'\\langle{{{sl[0]}{sl[1]}}}\\vert \\hat{{g}} \\vert{{{sl[2]}{sl[3]}}}\\rangle' # 公式输出时会自动换行 ？
        elif f is 0: return F'\\langle{{{sl[0]}}}\\vert \\hat{{h}} \\vert{{{sl[1]}}}\\rangle'
        elif f is 3: return F'\\langle{{{sl[0]}{sl[1]}}}\\vert \\hat{{g}} \\vert{{{sl[3]}{sl[2]}}}\\rangle'
        elif f is 4:
            return '\\langle{{{0}{1}}}\\vert \\hat{{g}} \\vert{{{2}{3}}}\\rangle'\
                   '-\\langle{{{0}{1}}}\\vert \\hat{{g}} \\vert{{{3}{2}}}\\rangle'.format(*sl.values())
        elif f is 1: return F'\\langle{{{sl[0]}}}\\vert \\hat{{f}} \\vert{{{sl[1]}}}\\rangle'
        else:
            raise ValueError(F'unrecognized {f} for Symb.f in Symb.__str__()')
    __repr__ = __str__
    def __format__(self, format_spec):
        f = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        sl = {i: F'p[{upl[i]},{symb.l}]' for i, symb in enumerate(self.l)}
        if f is 2: return F'g[ijkl({sl[0]},{sl[2]},{sl[1]},{sl[3]})]'
        elif f is 0: return F'h[{sl[0]},{sl[1]}]'
        elif f is 3: return F'g[ijkl({sl[0]},{sl[3]},{sl[1]},{sl[2]})]'
        elif f is 4: return F'(g[ijkl({sl[0]},{sl[2]},{sl[1]},{sl[3]})]-g[ijkl({sl[0]},{sl[3]},{sl[1]},{sl[2]})])'
        elif f is 1: return F'f[{sl[0]},{sl[1]}]'
        else:
            raise ValueError(F'unrecognized Symb.f = {f} for Integral in Symb.__format__()')

    def orbsort(self):
        symbl, upl = self.l, self.upl; norb = len(symbl)
        od = {i:(upl[i].l, symbl[i].l) for i in range(norb)} # sort orb: I0，I1, A0, A1
        pr = list(range(0, norb, 2)); qs = list(range(1, norb, 2))
        pr = sorted(pr, key=lambda x:od[x]); qs = sorted(qs, key=lambda x:od[x])
        if [od[i] for i in qs]<[od[i] for i in pr]: pqrs = list(chain(*zip(qs, pr)))
        else: pqrs = list(chain(*zip(pr, qs)))
        self.upl=tuple(upl[i] for i in pqrs); self.l = tuple(symbl[i] for i in pqrs)
class Functor(Symb):
    '''
    Functor class

    :param symbl: Orbital set
    :param pdl: script@hat or head, act on Orbital, Block
                 0,^: Operator
                 1,+: plus, creation operator
                -1,-: delete, annihilation operator
    '''

    def __init__(self, symbl=(), pdl='', upl=()):
        Symb.__init__(self, c='F', l=symbl, f=pdl)
        self.upl = upl

    def __str__(self):
        pdl = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        if upl: return ''.join(F'\\hat{{{s.l}}}_{{{s.f}}}^{{({upl[i]}){pdl[i]}}}' for i,s in enumerate(self.l))
        else: return ''.join(F'\\hat{{{s.l}}}_{{{s.f}}}^{{{pdl[i]}}}' for i,s in enumerate(self.l))
    __repr__ = __str__

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'F{self.l}{self.upl}{self.f}' == F'F{other.l}{other.upl}{other.f}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F'F{self.l}{self.upl}{self.f}')

    @property
    def nsymb(self): return len(self.l)
    @property
    def ne(self):
        if self.l:
            if isinstance(self.l[0], Orbital):
                return sum(1 if e is '+' else -1 for e in self.f)
            elif isinstance(self.l[0], Block):
                ml = self.f
                return sum(b.ne if ml[i] is '+' else -b.ne for i,b in enumerate(self.l))
            else:
                raise TypeError(F'unrecognited type {self.l[0].__class__} for Functor.l in Functor.ne()')
        else: return 0
    @property
    def sz(self):
        if self.l:
            if isinstance(self.l[0], Orbital):
                return sum(1 if o.f == '\\alpha' else -1 for o in self.l)
            elif isinstance(self.l[0], Block):
                return sum(b.sz for b in self.l)
            else:
                raise TypeError(F'unrecognited type {self.l[0].__class__} for Functor.l in Functor.sz()')
        else: return 0
    def blockize(self):
        '''
        blockize Functor

        :param case: one of in, state, ordinal and '+-'
                        up:   block
                     state:   orbital spin
                   ordinal:   orbital ordinal
                        +-:
        :return:
        '''

        Al = tuple(sorted(set(self.upl))) # 只可能在不同block, 不可有相同block(不同state)
        nb = len(Al)

        if nb is 1:
            self.upl = ()
            return 1, (self,), Al
        else:
            ill = tuple([] for A in Al); symbll = deepcopy(ill); pdll = ['' for A in Al]
            symbl = self.l; upl = self.upl; pdl = self.f
            for isymb in range(self.nsymb):
                A = index(upl[isymb], Al)
                ill[A].append(isymb)
                symbll[A].append(symbl[isymb])
                pdll[A] += pdl[isymb]

            # functorl
            functorl = tuple(Functor(symbl=tuple(symbll[A]), pdl=pdll[A]) for A in range(nb))
            # sign of repermutation of functorl  [f1f2][f3] AB...
            s_reperm = sign(list(chain.from_iterable(ill)))
            # sign of movation to block for functorl [f1f2]A [f3]B ...
            nswap = 0
            for ib in range(nb-1): nswap += Al[ib].ne*sum(functorl[ibf].nsymb for ibf in range(ib+1, nb))
            if nswap%2: s_reperm *= -1

            return s_reperm, functorl, Al




# At the level, Orbital < Block < Wave
# Orbital, block, wave is Symb, and also possess an implemented Expr in principle
#
# But in fact, the bottom Orbital is only just Symb, has no Expr
#         and  the top floor Wave isn't a Symb, but Expr
''' Orbital level '''
class Orbital(Symb):
    '''
    Orbital class

    :param ordinal: str
                              '':  unkown
                      p, q, r, s:  arbitrary block
                      i, j, k, l:  occopied block
                      a, b, c, d:  unoccopied block
    :param spin: str or int
                              '':  unkown
                               0:  spin state
                            1, a:  alpha space state
                           -1, b:  beta space state
    :param block: Symb
    '''

    def __init__(self, ordinal='p', spin=''):
        Symb.__init__(self, c='O', l=ordinal, f=spin)
        # # 在 =,>,<,hash,cmp 比较中，每个都在调用 __format__(), 所以这里能够进一步优化，比如定义一个format属性
        # # 但是注意对象是一个可变对象，改变某个属性，这个format属性就必须改变
        # self.str = F'{{{self.l}}}_{{{self.f}}}'
        # self.format = F'{self.l}{self.f}'

    def __str__(self): return F'{{{self.l}}}_{{{self.f}}}'
    __repr__ = __str__
    def __format__(self, format_spec): return F'{self.l}{self.f}'

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'O{self.l}{self.f}' == F'O{other.l}{other.f}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'O{self.l}{self.f}' > F'O{other.l}{other.f}'
nbo_ = 2 # 2 space orbitals in block
bo0 = Orbital(ordinal='0') # 0 occupied
bo1 = Orbital(ordinal='1') # 1 unoccpied or virtual
bo_l = (bo0, bo1)
nbo = 4 # 4 spin orbitals, nbe is 2 in block
bo0a = Orbital(ordinal='0', spin='\\alpha')
bo0b = Orbital(ordinal='0', spin='\\beta')
bo1a = Orbital(ordinal='1', spin='\\alpha')
bo1b = Orbital(ordinal='1', spin='\\beta')
bol = (bo0a, bo0b, bo1a, bo1b)  # block orbital list
''' block level '''
class Block(Symb):
    '''
    Block class

    :param ordinal:str
                              '':  unkown
                      P, Q, R, S:  any Block
                      I, J, K, L:  occopied Block
                      A, B, C, D:  unoccopied Block
    :param state: str or int
                              '':  unkown state
                               0:  ground state
                          1, ...:  state index
    :param wave: Phi(0) or Phi(a)
                 B(0) may be in Phi(0), and may be Phi(a) when excitation on a Block
    '''

    def __init__(self, ordinal='P', state=0): Symb.__init__(self, c='B', l=ordinal, f=state)

    def __str__(self): return F'{{{self.l}}}_{{{self.f}}}'
    __repr__ = __str__
    def __format__(self, format_spec): return F'({self.l},{self.f})'

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' == F'B{other.l}{other.f}' # Python Bug: is 有时是 False, 有时是 True
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' > F'B{other.l}{other.f}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F'B{self.l}{self.f}')

    @property
    def ne(self): return bfd[self.f].ne
    @property
    def sz(self): return bfd[self.f].sz
nhb = 4 # any block labeled as PQRS
bhl = 'PQRS'; b_l = set(bhl)
ntb = 4
bt0l = 'IJKLMNUV'; b0l = set(bt0l) # ground block labeled as IJKLMNUV
bt1l = 'ABCD'; b1l = set(bt1l) # excited block labeled as ABCD




''' Symb: 16 block functor '''
nbf = 16
bfd = {
        0: Functor(symbl=(bo0a, bo0b), pdl='++'),
        1: Functor(symbl=(bo1a, bo1b), pdl='++'),
        2: Functor(symbl=(bo0a, bo1b), pdl='++'),
        3: Functor(symbl=(bo0b, bo1a), pdl='++'),
        4: Functor(symbl=(bo0a, bo1a), pdl='++'),
        5: Functor(symbl=(bo0b, bo1b), pdl='++'),
        6: Functor(symbl=()),
        7: Functor(symbl=(bo0a,), pdl='+'),
        8: Functor(symbl=(bo1a,), pdl='+'),
        9: Functor(symbl=(bo0b,), pdl='+'),
        10: Functor(symbl=(bo1b,), pdl='+'),
        11: Functor(symbl=(bo0a, bo0b, bo1a), pdl='+++'),
        12: Functor(symbl=(bo0a, bo1a, bo1b), pdl='+++'),
        13: Functor(symbl=(bo0a, bo0b, bo1b), pdl='+++'),
        14: Functor(symbl=(bo0b, bo1a, bo1b), pdl='+++'),
        15: Functor(symbl=(bo0a, bo0b, bo1a, bo1b), pdl='++++')
      }
# Symb: 16 states in block
nbs = nbf
bsd = {s: Block(state=s) for s in range(nbs)}
''' Expr: 16 Functor in block '''
class Subbfeq(object):
    '''
    subitem block functor Expr

    A term in block functor (unit structure)
    can be a Expr or structure tpye,

    whose combination is block functor


    subbfe = C x B
    bfe = sum( subbfe )
    '''

    def __init__(self, bf=0, bs=0):
        self.ce = Coeff(symbl=(bf, bs), mode='inv') if bs<4 else 1
        self.block = bsd[bs]
    def __str__(self): return F'{{{self.ce!s}{self.block!s}}}'
    __repr__ = __str__
bfeqd = {
            bfd[0]: (Subbfeq(bf=0, bs=0), Subbfeq(bf=0, bs=1)),
            bfd[1]: (Subbfeq(bf=1, bs=0), Subbfeq(bf=1, bs=1)),
            bfd[2]: (Subbfeq(bf=2, bs=2), Subbfeq(bf=2, bs=3)),
            bfd[3]: (Subbfeq(bf=3, bs=2), Subbfeq(bf=3, bs=3)),
            bfd[4]: (Subbfeq(bf=4, bs=4),),
            bfd[5]: (Subbfeq(bf=5, bs=5),),
            bfd[6]: (Subbfeq(bf=6, bs=6),),
            bfd[7]: (Subbfeq(bf=7, bs=7),),
            bfd[8]: (Subbfeq(bf=8, bs=8),),
            bfd[9]: (Subbfeq(bf=9, bs=9),),
            bfd[10]: (Subbfeq(bf=10, bs=10),),
            bfd[11]: (Subbfeq(bf=11, bs=11),),
            bfd[12]: (Subbfeq(bf=12, bs=12),),
            bfd[13]: (Subbfeq(bf=13, bs=13),),
            bfd[14]: (Subbfeq(bf=14, bs=14),),
            bfd[15]: (Subbfeq(bf=15, bs=15),)
        }
''' Expr: 16 states in block '''
class Subbseq(object):
    '''
    subitem block state Expr

    A term in block configure (unit structure)
    can be a Expr or structure tpye,

    whose combination is block state


    subbse = C x F
    bse = sum( subbse )
    '''

    def __init__(self, bs=0, bf=0):
        self.ce = Coeff(symbl=(bs, bf), mode='ci') if bs<4 else 1
        self.functor = bfd[bf]
    def __str__(self): return F'{{{self.ce!s}{self.functor}}}'
    __repr__ = __str__

    def sortorb(self):
        if sign(self.functor.l) is -1: self.ce = 0.0-self.ce
        self.functor.l = tuple(sorted(self.functor.l))
bseqd = {
            0: (Subbseq(bs=0, bf=0), Subbseq(bs=0, bf=1)),
            1: (Subbseq(bs=1, bf=0), Subbseq(bs=1, bf=1)),
            2: (Subbseq(bs=2, bf=2), Subbseq(bs=2, bf=3)),
            3: (Subbseq(bs=3, bf=2), Subbseq(bs=3, bf=3)),
            4: (Subbseq(bs=4, bf=4),),
            5: (Subbseq(bs=5, bf=5),),
            6: (Subbseq(bs=6, bf=6),),
            7: (Subbseq(bs=7, bf=7),),
            8: (Subbseq(bs=8, bf=8),),
            9: (Subbseq(bs=9, bf=9),),
            10: (Subbseq(bs=10, bf=10),),
            11: (Subbseq(bs=11, bf=11),),
            12: (Subbseq(bs=12, bf=12),),
            13: (Subbseq(bs=13, bf=13),),
            14: (Subbseq(bs=14, bf=14),),
            15: (Subbseq(bs=15, bf=15),)
        }
def base2theory():
    '''
    base of blocak state -> tex file

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'base'
    ftitle = 'Foundation of GVB-BCCC formula derivation'
    with open(F'doc{ossep}{fpre}.tex', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write(F'% {ftitle}\n')
        write('\n\n')

        write( '\\documentclass[14pt]{article}\n')
        write( '\\usepackage[a4paper,hmargin=1.25in,vmargin=1in]{geometry}\n')
        write( '\\usepackage{setspace, amsmath, amsfonts}\n')
        write( '\n\n')

        write( '\\begin{document}\n')
        write(F'    \\title{{{ftitle}}}\n')
        write( '    \\author{Qingchun Wang}\n')
        write( '    \\maketitle\n')
        write( '    \\centerline{Nanjing Univ}\n')
        write( '    \\centerline{qingchun720@foxmail.com}\n')
        write( '    \n')

        write( '    \\linespread{2.0}\selectfont\n')
        write( '    \\thispagestyle{empty}\n')
        write( '    \n')

        write( '    \\newpage\n')
        write( '    \\setlength{\parindent}{0pt}\n')
        write( '    \\setcounter{page}{1}\n')
        write( '    \n')

        # block state
        write(F"    A block has 16 states: ${','.join(F'{b!s}' for b in bsd.values())} $ \\\\ \n")
        # block functor
        write( '    Each state corresponds to a functor in essence: \\\\ \n')
        for f in range(nbf):
            write(F'    {f}: ${bfd[f]} $ \\\\ \n')
        write( '    \n')

        # block functor of linear combination of state
        write( '    However, a functor is often expressed as a linear combination of a number of States: \\\\ \n')
        for f in range(nbf): write(F"    {f}: ${bfd[f]} = {'+'.join(F'{cxb}' for cxb in bfeqd[bfd[f]])} $ \\\\ \n")
        write( '    \n')

        # block state of linear combination of functor
        write( '    Therefore, the composition of each state can be written: \\\\ \n')
        for s in range(nbs):
            write(F"    ${bsd[s]!s} = {'+'.join(F'{cxf}' for cxf in bseqd[s])} \\quad \\vert{{vac}}\\rangle $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




''' repressive functor '''
rfeqd = {}
rd = {}
def rf2block(functor, block):
    '''
    repressive functor project into block

    :param functor: Functor class
    :param block: Block class

    :return:
    '''

    forbl=functor.l; nsymb_1,pdl =len(forbl)-1,functor.f
    _cxfl = deepcopy(bseqd[block.f])

    cxfl = []; append = cxfl.append  # C x f
    for cxf in _cxfl:
        f = cxf.functor; sortorb = cxf.sortorb
        fl = list(f.l); remove = fl.remove; insert = fl.insert
        s_annihil = 1 # s_annihil
        for iforb in range(nsymb_1, -1, -1):
            orb = forbl[iforb]
            exist = orb in fl
            if pdl[iforb] is '-':
                if exist:
                    iborb = index(orb, fl)
                    remove(orb)
                    f.f = f.f[:-1]
                    if iborb%2: s_annihil *= -1
                else: break
            else:  # '+'
                if exist: break
                else:
                    insert(0, orb)
                    f.f += '+'
        else:
            f.l = fl
            if s_annihil is -1: cxf.ce = 0.0-cxf.ce
            sortorb()  # s_create by sort
            append(cxf)

    cxbl = {} # C x b
    for cxf in cxfl:
        _ce = cxf.ce
        for cxb in deepcopy(bfeqd[cxf.functor]):
            ce = _ce*cxb.ce
            s = cxb.block.f
            try: cxbl[s] = cxbl[s]+ce # merge C1 Pp, C2 Pp
            except KeyError: cxbl[s] = ce

    return cxbl
def detrfeq():
    '''
    determine repressive functor equal value

    :return:
    '''

    nrf = 0
    # repressive functor status
    nrfsl = 4; rfsll = (('+', '-'), ('++', '+-', '--'), ('++-', '+--'), ('++--',))
    # nrfs=8; rfsl=['+', '-', '++', '+-', '--', '++-', '+--', '++--']
    ab13 = {1, 3}; ab03 = {0, 3}

    for nrfs in range(nrfsl):
        nrfs1 = nrfs+1
        for ol in product(bol, repeat=nrfs1):
            if nrfs1 is 4 and tuple(o.f for o in ol).count('\\alpha') in ab13: continue
            for rfsl in rfsll[nrfs]:
                rf = Functor(symbl=ol, pdl=rfsl)

                # m = [[0.0 for i in range(nbs)] for j in range(nbs)]
                cxbll = {} # {j:{} for j in range(nbs)}
                for q in range(nbs):
                    cxbll[q] = rf2block(rf, bsd[q])
                    for p in cxbll[q]:
                        r = Coeff(symbl=(nrf, p, q), mode='r')
                        rd[r] = cxbll[q][p]
                        cxbll[q][p] = r
                rfeqd[rf] = cxbll

                nrf += 1
def rfeq2theory():
    '''
    repressive functor equal value -> tex file

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'rep'
    ftitle = 'Repressive matrix of functor in GVB-BCCC formula derivation'
    with open(F'doc{ossep}{fpre}.tex', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write(F'% {ftitle}\n')
        write('\n\n')

        write('\\documentclass[14pt]{article}\n')
        write('\\usepackage[a4paper,hmargin=1.25in,vmargin=1in]{geometry}\n')
        write('\\usepackage{setspace, amsmath, amsfonts}\n')
        write('\n\n')

        write( '\\begin{document}\n')
        write(F'    \\title{{{ftitle}}}\n')
        write( '    \\author{Qingchun Wang}\n')
        write( '    \\maketitle\n')
        write( '    \\centerline{Nanjing Univ}\n')
        write( '    \\centerline{qingchun720@foxmail.com}\n')
        write( '    \n')

        write('    \\linespread{2.0}\selectfont\n')
        write('    \\thispagestyle{empty}\n')
        write('    \n')

        write( '    \\newpage\n')
        write( '    \\setlength{\parindent}{0pt}\n')
        write( '    \\setcounter{page}{1}\n')
        write( '    \n')

        write(F'    There are {len(rfeqd)} functors constructed: \\\\ \n')
        for i,rf in enumerate(rfeqd.keys()):
            write(F'    $\\hat{{O}}_{{{i}}}:  \\langle{{P_p}}\\vert {rf} \\vert{{P_q}}\\rangle => $ \\\\ \n')
            for q in rfeqd[rf]:
                write(F"    $ {rf} \\vert{{P_{{{q}}}}}\\rangle = {'+'.join(F'{ce!s}P_{{{p}}}' for p,ce in rfeqd[rf][q].items())} $ \\\\ \n")
                for p,ce in rfeqd[rf][q].items():
                    write(F'    ${ce!s}\\ =\\ {rd[ce]!s} $ \\\\ \n')
            write('    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def rfeq2code(lang='Python'):
    '''
    repressive functor equal value -> code

    :param lang: str, default is Python

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'rep'
    with open(F'bccc{ossep}{fpre}.py', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b %d %Y\n"))
        write('#\n')
        write('\n\n')

        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')

        write('import numpy\n')
        write('\n\n')

        write(F'nbs = {nbs}\n')
        write(F'nrf = {len(rfeqd)}\n')
        write('\n\n')

        write('def rep(ci):\n')
        write('    r = {}\n')
        write('    inv = numpy.linalg.inv(ci)\n')
        write('    \n')

        for i,rf in enumerate(rfeqd.keys()):
            write(F'    # O{i}:  <Bp|  {rf}  |Bq> = \n')
            if any(rfeqd[rf].values()):
                write('    m = numpy.zeros(shape=(nbs,nbs)) \n')
                for q in rfeqd[rf]:
                    for p,ce in rfeqd[rf][q].items():
                        write(F'    m[{p},{q}] = {rd[ce]}\n')
                write(F'    r[{i}] = m\n')
            write('    \n')
        write('    \n')

        write('    return r\n')
        write('    \n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




''' Tn expansion '''
def tcomb(nt):
    '''
    taylor expansion for T

    :param ntuply:

    :return:

    '''

    tnl = []; append = tnl.append
    for nb in range(1,nt+1):
        for tn in product(range(1,ntb+1), repeat=nb):
            if sum(tn)==nt: append(tuple(sorted(tn, reverse=True)))

    return tuple(sorted(set(tnl), reverse=True))
tnd = {nt: tcomb(nt) for nt in range(1,ntb+nhb+1)}
def tamp(nb):
    '''
    t amplify

    :param nb:

    :return:
    '''

    sl = ()
    for sn in product(range(1, nbs), repeat=nb):
        ne = sum(bfd[s].ne for s in sn)
        sz = sum(bfd[s].sz for s in sn)
        if ne is nb*2 and sz is 0:
            sl += (sn,)

    return sl
tsd = {nb: tamp(nb) for nb in range(1,ntb+1)}
tcd = {(): [[Frac(1,1)]]}
def taylor():
    for nb in range(1,ntb+nhb+1):
        for tl in tnd[nb]:
            for sl in product(*[tsd[t] for t in tl]):
                sl = tuple(chain(*sl))
                key = tuple(Block(ordinal=bt0l[i], state=sl[i]) for i in range(nb))
                if key in tcd:
                    keyl = list(tcd.keys())
                    key = keyl[index(key, keyl)]
                cl = []; append = cl.append
                nl = tuple(tl.count(t) for t in set(tl))
                ctt = reduce(mul, [math.factorial(n) for n in nl])
                ctn = reduce(mul, [math.factorial(n) for n in tl])
                append(Frac(1,ctt*ctn))
                t0 = 0
                for t in tl:
                    append(Coeff(symbl=key[t0:t0+t], mode='t'))
                    t0+=t
                try: tcd[key].append(cl)
                except KeyError: tcd[key] = [cl]
def t2thoery():
    '''
    raylor expansion -> tex file

    :return:
    '''

    raise NotImplementedError('TODO: t2theory()')




''' matrix elem for hamiltion '''
hfeqd = {} # 这种简单方式比 Manager.dict 更快
# hfeqd = Manager().dict()  # for parallel
tud = {}
# tud = Manager().dict()
ad = {}
# ad = Manager().dict()
def tu2theory():
    '''
    array of T element u -> latex

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    for bra in tud:
        for ket,hl in tud[bra].items():
            fpre = F"at{''.join(F'{b.l}{b.f}' for b in bra)}_{''.join(F'{b.l}{b.f}' for b in ket)}"
            ftitle = 'Amplitude equation in GVB-BCCC formula derivation'
            with open(F'doc{ossep}{fpre}.tex', 'w') as fout:
                write = fout.write
                # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

                write(F'% {ftitle}\n')
                write( '\n\n')

                write( '\\documentclass[14pt]{article}\n')
                write( '\\usepackage[a4paper,hmargin=1.25in,vmargin=1in]{geometry}\n')
                write( '\\usepackage{setspace, amsmath, amsfonts}\n')
                write( '\n\n')

                write( '\\begin{document}\n')
                write(F'    \\title{{{ftitle}}}\n')
                write( '    \\author{Qingchun Wang}\n')
                write( '    \\maketitle\n')
                write( '    \\centerline{Nanjing Univ}\n')
                write( '    \\centerline{qingchun720@foxmail.com}\n')
                write( '    \n')

                write( '    \\linespread{2.0}\selectfont\n')
                write( '    \\thispagestyle{empty}\n')
                write( '    \n')

                write( '    \\newpage\n')
                write( '    \\setlength{\parindent}{0pt}\n')
                write( '    \\setcounter{page}{1}\n')
                write( '    \n')

                # # bra
                # write( '    Intermediate array: \\\\ \n')
                # wl = []
                # for ixcl in tud[ket]:
                #     for a in ixcl[1]:
                #         if a in wl: continue
                #         else:
                #             wl.append(a)
                #             sumbl = ''.join(F'\\displaystyle\\sum_{{{e}}}' for e in a.f[:index('r', a.f)] if e in b0l)
                #             if isinstance(iad[a][1], Coeff):
                #                 a_ = iad[a][1]
                #                 sa = F"X{index(a_, list(iad.keys()))}({','.join(F'{b.l}' for b in a_.l)})" if a_.l else F'X{index(a_, list(iad.keys()))}'
                #             else: sa = F'{iad[a][1]!s}'
                #             if a.l:
                #                 write(F"    $X{index(a, list(iad.keys()))}({','.join(F'{b.l}' for b in a.l)}) = {sumbl}{iad[a][0]!s}*{sa}$ \\\\ \n")
                #             else:
                #                 write(F'    $X{index(a, list(iad.keys()))} = {sumbl}{iad[a][0]!s}{sa}$ \\\\ \n')
                #             # write(F'    ${a!s} = {sumbl}{iad[a]!s}$ \\\\ \n')
                # write( '    \n')

                # write(F"    $\\langle{{{bra}}} \\vert \\hat{{H}} \\vert ({'+'.join(''.join(F'T_{t}' for t in tl) for tl in tnd[len(ket)])})\\rangle = $ \\\\ \n")
                # bl = [b.l for b in ket if b.l in b0l]
                for ixcl in tud[bra][ket]:
                    sumbl = ''.join(F'\\displaystyle\\sum_{{{e}}}' for e in sorted(set([b.l for b in ixcl[1].upl if b.l in b0l])))
                    # write(F"    $+{ixcl[0]}{sumbl}{reduce(mul, ixcl[1:])!s} $ \\\\ \n")
                    write(F"    $+\\frac{{{ixcl[0].numerator}}}{{{ixcl[0].denominator}}}{sumbl}{reduce(mul, ixcl[1:])!s} $ \\\\ \n")
                write( '    \n')

                write('\\end{document}\n')
                write('\n')

                # sys.stdout = stdout; sys.stderr = stderr
                # fout.close()
def tu2code():
    '''
    array of T element u -> code

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep
    tab = '    '

    for bra in tud:
        _fpre = F"at{''.join(F'{b.l}{b.f}' for b in bra)}"
        s_bra = F"sign([{','.join(b.l for b in bra if 6<b.f<15)}])"
        k_bra = F"tuple(sorted([{','.join(F'({b.l},{b.f})' for b in bra)}]))"
        bbl = [b.l for b in bra]

        for ket,hl in tud[bra].items():
            fpre = F"{_fpre}_{''.join(F'{b.l}{b.f}' for b in ket)}"
            s_ket = F"sign([{','.join(b.l for b in ket if 6<b.f<15)}])"
            k_ket = F"tuple(sorted([{','.join(F'({b.l},{b.f})' for b in ket)}]))"
            kbl = [b.l for b in ket if b.l not in bbl]

            hld = {}
            for h in hl:
                k = tuple(sorted(set(b.l for b in h[1].upl if b.l in b0l and b.l not in kbl)))
                try: hld[k].append(h)
                except KeyError: hld[k] = [h]

            with open(F'bccc{ossep}at{ossep}{fpre}.py', 'w') as fout:
                # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
                write = fout.write

                write( '#!/usr/bin/env python \n')
                write( '# -*- coding: utf-8 -*- \n')
                write( '#\n')
                write( '#      Copyright:   Qingchun Wang @ NJU\n')
                write(F'#      File Name:   {fpre}.py\n')
                write(F'#            Des:   {fpre}\n')
                write( '#           Mail:   qingchun720@foxmail.com\n')
                write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b %d %Y\n"))
                write('#\n')
                write('\n\n\n\n')

                write("__version__ = '1.0'\n")
                write("__author__ = 'Qingchun Wang'\n")
                write('\n\n')

                write('import numpy\n')
                write('from bccc.pub import sign,ijkl\n')
                write('\n\n\n\n')

                write(F'def {fpre}(h, g, r, t, p):\n') # iterate amplitude
                write( '    hdd = {}\n')
                write( '    n = len(p)\n')
                write( '    \n')

                nbb = 0
                for b in bbl:
                    write(tab*nbb+F'    for {b} in range(0, n):\n')
                    if nbb: write(tab*nbb+F"        if {b} in [{','.join(bbl[:nbb])}]: continue\n")
                    nbb += 1
                write(tab*nbb+F'    u = {k_bra}\n')
                write(tab*nbb+ '    hd = {}\n')
                write(tab*nbb+ '    \n')

                bk = bbl+kbl
                nkb = nbb
                for b in kbl:
                    write(tab*nkb+F'    for {b} in range(0, n):\n')
                    if nkb: write(tab*nkb+F"        if {b} in [{','.join(bk[:nkb])}]: continue\n")
                    nkb += 1
                write(tab*nkb+F'    v = {k_ket}\n')
                write(tab*nkb+ '    huv = 0.0 \n')
                write(tab*nkb+ '    \n')

                for Il,hl in hld.items():
                    if Il:
                        bkI = bk+list(Il)
                        # nI = nkb
                        # forl = []; append = forl.append
                        # for b in Il:
                        #     if nI: append(F"for {b} in range(n) if {b} not in [{','.join(bkI[:nI])}]")
                        #     else: append(F"for {b} in range(n)")
                        #     nI += 1
                        for h in hl: write(tab*nkb+F"    huv += sum({float(h[0])*reduce(mul, h[1:])} {' '.join(forl)})\n")
                    else:
                        for h in hl: write(tab*nkb+F'    huv += {float(h[0])}*{reduce(mul, h[1:])}\n')
                write(tab*nkb+ '    \n')

                write(tab*nkb+F"    hd[v] = hd.get(v, 0.0)+huv*{s_bra}*{s_ket}\n")
                write(tab*nkb+ '    \n')

                write(tab*nbb+ "    try: \n")
                write(tab*nbb+ "        for v,huv in hd.items(): hdd[u][v] = hdd[u].get(v, 0.0)+huv\n")
                write(tab*nbb+ "    except KeyError:\n")
                write(tab*nbb+ "        hdd[u] = hd\n")
                write(tab*nbb+ '    \n')

                write(F'    return hdd\n')
                write( '    \n\n\n\n')

                # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
                # write("if __name__ == '__main__':\n")
                # write('    \n\n')
                # write("    print('End successfully')\n")
                # write('    \n\n\n\n')

                # # sys.stdout = stdout; sys.stderr = stderr
                # fout.close()
def t2theory():
    '''
    array of T element u -> latex

    main function
    :return:
    '''

    raise NotImplementedError('TODO: t2theory()')
def t2code():
    '''
    array of T element u -> code

    main function
    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'at'
    with open(F'bccc{ossep}at{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %hdd %H:%M:%S %b/%d %Y\n"))
        write('#\n')
        write('\n\n')

        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')

        write('from datetime import datetime\n')
        write('import logging\n')
        write('from copy import deepcopy\n')
        write('from multiprocessing import Pool, Manager\n')
        write('\n\n')

        '''
        write('def at(t0, h, pair, ntuply=4):\n')
        write('    t = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        write('    from bccc.at.at_ import at_\n')
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"at_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[1])
        for braf in brafl: write(F'    from bccc.at.{braf} import {braf}\n')
        write('    brafl = (\n             ')
        for i,braf in enumerate(brafl):
            write(F'{braf}, ')
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for ntuply in range(2,ntb+1):
            write(F'    # {ntuply}-block excited types\n')
            write(F'    if ntuply>{ntuply-1}: \n')
            brafl = tuple(F"at_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[ntuply])
            for braf in brafl: write(F'        from bccc.at.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F'{braf}, ')
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write("    print('    it            dE               dt               time/s')\n")
        write('    E0 = h[()][()]+at_(t0, h, pair)\n')
        write('    for it in range(200):\n')
        write('        t1 = datetime.now()\n')
        write('        \n')

        write('        pool = Pool()\n')
        write('        vl = tuple(pool.apply_async(braf, (t0, h, pair, E0)) for braf in brafl)\n')
        write('        pool.close()\n')
        write('        pool.join()\n')
        write('        for v in vl: t.update(v.get())\n')
        write('        # for braf in brafl: t.update(braf(t0, h, pair, E0))')
        write('        \n')

        write('        t2 = datetime.now()\n')
        write('        \n')

        write('        Ecorr = at_(t, h, pair); E = h[()][()]+Ecorr\n')
        write('        dE = E-E0; dt = sum(abs(t[u]-t0[u]) for u in t)\n')
        write("        print(F'    {it+1:3d}      {dE:13.10f}    {dt:13.10f}     {t2-t1}')\n")
        write('        if abs(dt)<1.0e-6: \n')
        write("            print(F'Successfully converged: Ecorr = {Ecorr}')\n")
        write('            break\n')
        write('        else: \n')
        write('            t0 = deepcopy(t)\n')
        write('            E0 = E\n')
        write('        \n')
        write('    else:\n')
        write("        print('Error: Convergence failed ')\n")
        write('    \n')

        write('    return Ecorr,t\n')
        write('    \n\n\n\n')
        '''

        modulel=[]; append=modulel.append
        for bra in tud:
            s = F"at{''.join(F'{b.l}{b.f}' for b in bra)}"
            for ket in tud[bra]: append(F"{s}_{''.join(F'{b.l}{b.f}' for b in ket)}")
        for module in modulel:
            write(F'from bccc.at.{module} import {module}\n')
        write( '\n\n')
        write( 'modulel = [')
        n = 10
        for i in range(len(modulel)//n):
            if i: write(F"           {','.join(modulel[i*10:(i+1)*10])},\n")
            else: write(F"{','.join(modulel[i*10:(i+1)*10])},\n")
        else:
            try:  write(F"           {','.join(modulel[(i+1)*10:])}\n")
            except(UnboundLocalError, NameError): write(F"           {','.join(modulel)}\n")
        write( '          ]')
        write( '\n\n')

        # write( 'tsd = {\n')
        # for nb,ts in tsd.items():
        #     write(F'        {nb}: (')
        #     for i in range(len(ts)//n):
        #         if i: write(F"            {','.join(F'{s}' for s in ts[i*10:(i+1)*10])},\n")
        #         else: write(F"{','.join(F'{s}' for s in ts[i*10:(i+1)*10])},\n")
        #     else:
        #         try: write(F"            {','.join(F'{s}' for s in ts[(i+1)*10:])}),\n")
        #         except(NameError, UnboundLocalError): write(F"{''.join(F'{s},' for s in ts)}),\n")
        # write( '       }\n')
        # write( '\n\n')

        write( 'def at(h, g, r, t, p):\n')
        write( '    hdd = {}\n\n')

        write( '    for f in modulel:\n')
        write( '        for u,hd in f(h, g, r, t, p).items():\n')
        write( '            try:\n')
        write( '                for v,huv in hd.items(): hdd[u][v] = hdd[u].get(v, 0.0)+huv\n')
        write( '            except KeyError: hdd[u] = hd\n')
        write( '    \n\n')

        write( '    return hdd')
        write( '\n\n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def h2t(ket):
    '''
    hamiltonnian project into ket

    :param ket: wave state, list

    :return:
    '''

    hld = {}; ud = {} # {u: [c, <>, r..]}
    nb = len(ket); ne = sum(b.ne for b in ket)
    bl = tuple(b.l for b in ket)
    spinl = ('\\alpha', '\\beta')

    # for h
    if nb<=2 and ne>=1:  # nb<2 for h; p+p- ne>=1
        for prodb in product(ket, repeat=2):
            if len(set(prodb)) is not nb: continue
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    inte = Integral(symbl=(bop, boq), mode=0)
                    functor = Functor(symbl=(bop, boq), pdl='+-', upl=prodb)
                    s_blockize,functorl,upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                    _rlld={}; _ud={} # ce {u: ce} 没有与 inte 相乘
                    for a in spinl:
                        bop.f = boq.f = a
                        bfeql = [None]*nb
                        for i,f in enumerate(functorl):
                            bfeq = rfeqd[f][upl[i].f]  # {s: ce, }
                            if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                            else: break
                        else:
                            for prodcxb in product(*bfeql): # 各个 block 的展开，得到系数、波函数态
                                bsl,rl = zip(*prodcxb)
                                key = deepcopy(upl)
                                for i,b in zip(bsl,key): b.f = i
                                if key in hld: key = ud[key] # 取 key of hld, 如果key已经有了
                                elif key in _rlld: key = _ud[key] # 取 key of _rlld, 如果key已经有了
                                _rl = []; append = _rl.append
                                for b,r in zip(key, rl):
                                    _r = deepcopy(r)
                                    _r.l = (b,)+_r.l
                                    append(_r)
                                try: _rlld[key].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                except KeyError: _ud[key],_rlld[key] = key,[_rl]
                    for key,rll in _rlld.items():  # 将 _rlld 更新到 hld 中
                        _inte = deepcopy(inte)
                        _inte.upl = tuple(key[index(b, bl)] for b in _bl)
                        _inte.orbsort()
                        si = [Frac(s_blockize,1), _inte]
                        try:
                            for rl in rll: hld[key].append(si+rl)
                        except KeyError:
                            for i,rl in enumerate(rll): rll[i] = si+rl
                            ud[key],hld[key] = key,rll

    # for g
    if ne>=2: # p+p+s-r- ne>=2
        for prodb in product(ket, repeat=4):
            if len(set(prodb)) is not nb: continue
            print(F'prodb = {prodb}')
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    for bor in deepcopy(bo_l):
                        for bos in deepcopy(bo_l):
                            inte = Integral(symbl=(bop, boq, bor, bos), mode=2)
                            functor = Functor((bop, boq, bos, bor), pdl='++--', upl=(prodb[0], prodb[1], prodb[3], prodb[2]))
                            s_blockize, functorl, upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                            _rlld = {}; _ud={}  # ce 没有与 inte 相乘
                            for a in spinl:
                                bop.f = bor.f = boq.f = bos.f = a
                                bfeql = [None]*nb
                                for i,f in enumerate(functorl):
                                    bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                    if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                    else: break
                                else:
                                    for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                        bsl,rl = zip(*prodcxb)
                                        key = deepcopy(upl)
                                        for i,b in zip(bsl, key): b.f = i
                                        if key in hld: key = ud[key] # 取 key of hld, 如果key已经有了
                                        elif key in _rlld: key = _ud[key] # 取 key of _rlld, 如果key已经有了
                                        _rl = []; append = _rl.append
                                        for b,r in zip(key, rl):
                                            _r = deepcopy(r)
                                            _r.l = (b,)+_r.l
                                            append(_r)
                                        try: _rlld[key].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError:  _ud[key],_rlld[key] = key,[_rl]
                                        # print(F'len(_rl) = {len(_rl)}')
                            for key,rll in _rlld.items():  # 将 _rlld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(key[index(b, bl)] for b in _bl)
                                _inte.orbsort()
                                si = [Frac(s_blockize,2), _inte]
                                try:
                                    for rl in rll: hld[key].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(rll): rll[i] = si+rl
                                    ud[key],hld[key] = key,rll

                            _rlld = {}; _ud = {}  # ce 没有与 inte 相乘
                            a,b = spinl
                            bop.f = bor.f = a; boq.f = bos.f = b
                            bfeql = [None]*nb
                            for i,f in enumerate(functorl):
                                bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                print(F'bfeq = {bfeq}')
                                if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                else: break
                            else:
                                for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                    bsl,rl = zip(*prodcxb)
                                    key = deepcopy(upl)
                                    for i,b in zip(bsl, key): b.f = i
                                    if key in hld: key = ud[key] # 取 key of hld, 如果key已经有了
                                    elif key in _rlld: key = _ud[key] # 取 key of _rlld, 如果key已经有了
                                    _rl = []; append = _rl.append
                                    for b,r in zip(key, rl):
                                        _r = deepcopy(r)
                                        _r.l = (b,)+_r.l
                                        append(_r)
                                    try: _rlld[key].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                    except KeyError: _ud[key],_rlld[key] = key,[_rl]
                            for key,rll in _rlld.items():  # 将 _rlld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(key[index(b, bl)] for b in _bl)
                                _inte.orbsort()
                                si = [Frac(s_blockize,1), _inte]
                                try:
                                    for rl in rll: hld[key].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(rll): rll[i] = si+rl
                                    ud[key],hld[key] = key,rll

    return hld
def dettu(ket, hfeqd, tud):
    '''
    determine amplitude

    :param ket:
    :param hfeqd

    :return:
    '''

    print('ket = ', ket)

    nkb=len(ket); nbe0=nkb-ntb if nkb>ntb else 0
    chbd={nbc: Frac(1, math.factorial(nbc)) for nbc in range(1, nhb+1)}  # 1/chb H
    for nbc in range(1, nhb+1):  # n blocks correlation
        nbe1 = nbc+1 if nkb>nbc else nkb+1
        for nbe in range(nbe0, nbe1): # nbe blocks from ket
        # for nbe in range(nbc+1): # nbe blocks from ket
            for bchose in combinations(ket, r=nbe): # 从激发态取 nbe
                # print(F'bchose = {bchose}')
                brest = tuple(sorted(set(ket)-set(bchose))) # set -> list is random and disordered
                i0=nkb; i1=i0+nbc-nbe
                _ket = bchose+tuple(Block(ordinal=i) for i in bt0l[i0:i1])
                print(F'_ket = {_ket}') 

                ksl = tuple(b.f for b in _ket); kbl = tuple(b.l for b in _ket)
                if ksl in hfeqd:
                    _hld = deepcopy(hfeqd[ksl])
                    try:
                        bl = tuple(b.l for b in list(_hld.keys())[0])
                        if bl!=kbl:
                            for bl in _hld:
                                for b,i in zip(bl, kbl): b.l = i
                    except IndexError: continue # _hld = {} when nbc=nbe, b.f=6
                else:
                    _hld = h2t(_ket)  # Aa I0 J0 Kk
                    hfeqd[ksl] = _hld

                for _key,_hl in _hld.items(): # 将 _hld 更新(多块加和)到 hld 中
                    _keybrest = _key+deepcopy(brest)
                    _bra = tuple(b for b in _keybrest if b.f); n=len(_bra)
                    if n<=ntb:
                        print(F'{ket} -> {_keybrest}')
                        print(F'len(_hl) = {len(_hl)}\n')
                        keybrest = tuple(sorted(_keybrest))
                        bra = tuple(sorted(_bra))
                        bl,tl = deepcopy([ket,tcd[ket]])
                        # print(F'bl,keybrest = {bl},{keybrest}')
                        i0 = i1 = 0
                        for i,b in enumerate(keybrest):
                            if i<nkb:
                                if b.f is 0:
                                    b.l = bl[i].l = bt0l[i0]
                                    i0 += 1
                                else:
                                    b.l = bl[i].l = bt1l[i1]
                                    i1 += 1
                            else:
                                if b.f is 0:
                                    b.l = bt0l[i0]
                                    i0 += 1
                                else:
                                    b.l = bt1l[i1]
                                    i1 += 1
                        # print(F'bl,keybrest = {bl},{keybrest}')
                        # s_bef = Expr('-?', [b1.l for b1,b2 in zip(_keybrest[:nbe]+_keybrest[nbc:], bchose+brest) if 6<b2.f<15])
                        # print(F's_bef.leafl = {s_bef.leafl}')
                        # s_aft = Expr('-?', [b.l for b in _bra if 6<b.f<15])
                        # print(F's_aft.leafl = {s_aft.leafl}\n')
                        hll = []; append = hll.append
                        for i,j in product(range(len(_hl)), range(len(tl))):
                            # hl = [_hl[i][0]*chbd[nbc]*tl[j][0]]
                            hl = [_hl[i][0]*chbd[nbc]]
                            # hl += _hl[i][1:]+tl[j][1:]
                            # hl += _hl[i][1:]+[s_bef,s_aft] # check Hamiltonian matrix
                            hl += _hl[i][1:] # check Hamiltonian matrix
                            append(hl)
                        if bra in tud:
                            try: tud[bra][bl] += hll
                            except KeyError: tud[bra][bl] = hll
                        else: tud[bra] = {bl: hll}
                        
    # # at
    # td = {}  # {ket: t}
    # ketl = set(); append = ketl.add
    # for _ket in hld.keys(): # Aa..Jj.. -> Aa..Ij..
    #     bbl = tuple(b for b in _ket if b.l in b1l and b.f>0); kbl = deepcopy(tuple(b for b in _ket if b.l in b0l and b.f>0))
    #     for b,l in zip(kbl, 'IJKL'): b.l = l
    #     append(bbl+kbl)
    # ketl.remove(ket)
    # for ket in ketl:
    #     bl, bsl = tuple(b.l for b in ket), tuple(b.f for b in ket); nb = len(bsl)
    #
    #     tket_ = []; append = tket_.append # 有效的分组
    #     for tl in tnd[nb]:
    #         for gl in gld[tl]:
    #             if all(tuple(bsl[i] for i in g) in tsd[len(g)] for g in gl): append(gl)
    #
    #     txtl = []; append=txtl.append
    #     for gl in tket_:
    #         if ket in hld:
    #             for ixcl in hld[ket]:
    #                 # ixcl[0] *= sign([i for i in chain(*gl) if 6<bsl[i]<15])
    #                 # append(ixcl+[Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl])
    #
    #                 # build intermediate array
    #                 # print(F'ixcl = {ixcl}')
    #                 rxtl = [ixcl[0]*sign([i for i in chain(*gl) if 6<bsl[i]<15])]
    #
    #                 bll = [[Block(ordinal=bl[i], state=bsl[i]) for i in g] for g in gl]
    #                 _bll = [[b for b in bl if b.l in b0l] for bl in bll]
    #                 tl = [Coeff(symbl=tuple(bl), mode='t') for bl in bll]
    #                 rll = [[r for b in bl for r in ixcl[2:] if r.l[0]==b] for bl in _bll]
    #                 # print(F'bll = {bll}')
    #
    #                 nvac = _bll.count([]) # del []
    #                 for vac in range(nvac):
    #                     ivac = index([], _bll)
    #                     bll.remove(bll[ivac]); _bll.remove([])
    #                     rxtl.append(tl[ivac]); tl.remove(tl[ivac])
    #                     if rll[ivac]: rxtl.append(rll[ivac])
    #                     rll.remove(rll[ivac])
    #
    #                 al = []
    #                 if bll:
    #                     rl=rll[0]; inteupl,inteol=ixcl[1].upl,ixcl[1].l
    #                     irinte=F"t{''.join(F'{b.l}{b.f}' for b in tl[0].l)}r{''.join(F'{r.l[1]}{r.l[0].l}{r.l[0].f}' for r in rl)}i{''.join(F'{inteupl[o].l}{inteol[o].l}' for o in range(len(inteol)))}"
    #                     v = tl[0]*reduce(mul, rll[0]), ixcl[1]
    #                     k = Coeff(symbl=tuple(chain(*_bll[1:])), mode=irinte)
    #                     # if k in iad and v!=iad[k]:
    #                     #     print(F'Exist k=_k & v!=_v')
    #                     #     print(F' k= {k!s}')
    #                     #     print(F' v= {v!s}')
    #                     #     print(F'_v= {iad[k]!s}')
    #                     # for _k,_v in iad.items():
    #                     #     if v==_v and _k!=k:
    #                     #         print(F'Exist v=_v & k!=_k')
    #                     #         print(F' v= {v!s}')
    #                     #         print(F'_v= {_v!s}')
    #                     #         print(F' k= {k!s}')
    #                     #         print(F'_k= {_k!s}')
    #                     al.append(k)
    #                     iad[k] = v
    #                     for i in range(1,len(bll)):
    #                         v = tl[i]*reduce(mul, rll[i]), k
    #                         k = Coeff(symbl=tuple(chain(*_bll[i+1:])), mode=F"t{''.join(F'{b.l}{b.f}' for b in tl[i].l)}r{''.join(F'{r.l[1]}{r.l[0].l}{r.l[0].f}' for r in rll[i])}c{al[i-1].f}")
    #                         # if k in iad and v!=iad[k]:
    #                         #     print(F'Exist k=_k & v!=_v')
    #                         #     print(F' k= {k!s}')
    #                         #     print(F' v= {v!s}')
    #                         #     print(F'_v= {iad[k]!s}')
    #                         # for _k,_v in iad.items():
    #                         #     if v==_v and _k!=k:
    #                         #         print(F'Exist v=_v & k!=_k')
    #                         #         print(F' v= {v!s}')
    #                         #         print(F'_v= {_v!s}')
    #                         #         print(F' k= {k!s}')
    #                         #         print(F'_k= {_k!s}')
    #                         al.append(k)
    #                         iad[k] = v
    #                 # print(F'al = {al}')
    #                 append([rxtl, al])
    #
    #     if txtl: td[ket] = txtl

    # tu2code(ket, hld)
    # mhuv2theory(ket, hld)
    # mhuv2code(ket, hld)




''' Parallel '''
def parallel():
    '''
    parallel

    :param func: dettu() or detmt()

    :return:
    '''

    mp = Pool()
    for u in tcd.keys(): mp.apply_async(dettu, (u, hfeqd, tud))
    mp.close()
    mp.join()





if __name__ == '__main__':
    # base2theory()

    detrfeq()
    # rfeq2theory()
    # rfeq2code()

    taylor()
    print(F'len(tcd) = {len(tcd)}')

    t1 = datetime.now()
    # if sys.platform=='linux':
    #     parallel()  # for parallel
    #     # Pool().map(dettu, ((u, hfeqd, tud) for u in tcd.keys())) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    # else:
    #     # dettu((), hfeqd, tud)
    # for ket in tcd.keys(): dettu(ket, hfeqd, tud)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # tu2theory()
    tu2code()
    # t2theory()
    t2code()


    print('End successfully')

