#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   derivform.py
#            Des:   39 -> diis
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
    2. reduce memory

    __copy__(), __deepcopy__()
    __bool__() -> Symb0; __bool__() -> Expr0
    str '0' in Block.l -> number 0 and similar 
'''


__version__ = '3.1.45'
__author__ = 'Qingchun Wang'




import os,sys,shutil,subprocess
sys.setrecursionlimit(1000000)
sys.setswitchinterval(10000)
from copy import copy,deepcopy
from datetime import datetime
import logging

from functools import reduce
from operator import mul,add,sub,truediv
from itertools import permutations,combinations,product,chain
from collections import Iterable

import math,random
from fractions import Fraction as Frac
# import numpy, scipy
# from sympy import Symbol

import pickle
from re import split,compile,findall
numstr = compile(r'^[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?$') # match integer, fraction, exponent
# numstr = compile(r'[-+]?(\b[0-9]+(\.[0-9]+)?|\.[0-9]+)([eE][-+]?[0-9]+\b)?') # eg. w123.345w
# findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", 'A1.45aa, b5., B6.45, F-2.e2ddf')

from multiprocessing import Pool,Manager,freeze_support
# import mpi4py

# from bccc.addons import
from bccc.pub import iis,iin,sign,group,prodl,sum_list




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
binaryfuncl = {',', ':', '-?', '**', 'pow', '%%', 'log'} # -? 特殊处理
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
def sigma(variable=None, ben=0, expr=None):
    '''
    sigma summation

    run function, how to compute
    :param variable:
    :param ben:
    :param expr:

    :return:
    '''

    tmp = Expr('+++', variable, ben, expr)
    if isinstance(expr, Expr): expr.root = tmp
    
    return tmp
def pi(variable=None, ben=0, expr=None):
    '''
    pi(II) product

    run function, how to compute
    :param expr:
    :param ben:
    :param end:

    :return:
    '''

    tmp = Expr('***', ben, variable, expr)
    if isinstance(expr, Expr): expr.root = tmp

    return tmp
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
# 在之前版本中，实现了Symb, Expr对 int/float, str 等类型的兼容
# 但是，在实际中，多数都是 Symb, Expr 的计算，而非 int/float, str 的计算
# 对 int/float, str 兼容，增加了对 type 的判断，会使行计算做很多无效计算，得不偿失
# 因此，在这个版本，删除对 int/float, str 兼容
class Symb(object):
    '''
    Class: Mathematical Symbol

           in which every term can be number(int and float), str, Symb or Expr
             在 =,>,<,hash,cmp 比较中，每个都在调用 __format__(), 所以这里能够进一步优化，比如定义一个format属性
             但是注意对象是一个可变对象，改变某个属性，这个format属性就必须改变
             self.str = F'{{{self.l}}}_{{{self.f}}}'
             self.format = F'{self.l}{self.f}'

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
                          A: Array                D

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

    def __str__(self):
        '''
        Symb -> str, latex output

        :return:
        '''

        return F'{self.c}_{{{self.l!s}}}^{{{self.f!s}}}'
    __repr__ = __str__
    def __format__(self, format_spec):
        '''
        Symb -> format, code output

        :param format_spec:

        :return:
        '''

        if not format_spec: return F'{self.c}[{self.l}][{self.f}]'
        else:
            # return format_spec.replace('%s', self.c).replace('%s', self.l).replace('%s', self.f)
            formatd = {'-f': '{0.c}[{0.l}]',
                       '+f': '{0.c}[{0.l}][{0.f]'}
            return formatd[format_spec].format(self)

    # def __copy__(self):
    #     '''
    #     shallow copy, shared IP
    #
    #     NOTE: =, copy and deepcopy are different
    #           1)        =:  all value are shared ID
    #                  copy:  list value are shared ID
    #              deepcopy:  all value are not shared
    #           2) overload of = operator don't be supported in Python due to it similar to labeling, shard ID
    #
    #     :return:
    #     '''
    #
    #     return copy(self)
    # def __deepcopy__(self, memodict={}):
    #     '''
    #     deepcopy, not shared ID
    #
    #     :return:
    #     '''
    #
    #     # return deepcopy(self)
    #     tmp = Symb(c=self.c)
    #
    #     from collections import Iterable
    #     if self.l.__class__ in [list, str, tuple]: tmp.l = self.l
    #     elif isinstance(self.l, Iterable):
    #         '''
    #         Bug1: self.l may not be [], let's assume that
    #         Bug2: elem in self.f also may be iterable
    #
    #             Note: fateful bug in Python language
    #                   deepcopy are not real copy
    #                   >>> l1= [2, 3]
    #                   >>> l2 = [l1, 'di', l1]
    #                   >>> l2_ = deepcopy(l2)
    #                   >>> l2_[0][0] = 9
    #                   >>> l2_
    #                       [[9, 3], 'di', [9, 3]]
    #         '''
    #         tmp.l = [deepcopy(elem) for elem in self.l]
    #     else: tmp.l = deepcopy(self.l)
    #     if self.f.__class__ in [list, str, tuple]: tmp.f = self.f
    #     elif isinstance(self.f, Iterable):
    #         # Bug: self.f may not be [], let's assume that
    #         # Bug: elem also may be iterable
    #         tmp.f = [deepcopy(elem) for elem in self.f]
    #     else: tmp.f = deepcopy(self.f)
    #
    #     return tmp

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

           Tree structure (unary, binary or ternary tree)
           which may a term(even only one Symb), many term

            在这个版本中，删除对0,1以及 特殊运算 的处理
            比如：0+s=s+0=s-0=s*1=1*s=s/1=s; (0+a)+(1*b)=(1*a)+(1*b)=(1*a)-(0+b)=a+b
            支持 0,1 特殊处理，计算变得杂，实际符号运算中，这并不常见
            但是，特别是 1*s=s 保留, 0-e=-e
            有效的运算包括：s,e +-*/ N,s,e 和 N +-*/ s,e

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

        self.root, self.oper, self.leafl = None, oper, symbl # 若symbl中有Expr, 请修改symb.root = self

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
                return F"\\displaystyle\\sum_{{{leafl[0]}={leafl[1]}}} {leafl[2]}"
            elif oper=='***':
                return F"\\prod_{{{leafl[0]}={leafl[1]}}} {leafl[2]}"
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
                    # if isinstance(leafl0[0],str): return F"sign([{','.join(leafl0)}])"
                    # else: return F'''sign({'+'.join(F"sorted([{','.join(bl)}])" for bl in leafl0)})'''
                    return F"sign([{','.join(self.leafl[0])}])"
                else:
                    raise NotImplementedError(F'TODO: Expr.__format__() for operator {oper}')
        elif oper in unaryoperl:
            if oper in unarycharl:
                if oper=='0-': return F'-{leafl[0]}'
                elif oper=='0+': return leafl[0]
                else: return F'({leafl[0]}{oper})'
            else: return F'math.{oper}({leafl[0]})'
        elif oper in ternaryoperl: # debugging ...
            if oper=='+++': return F'sum({leafl[2]} for {leafl[0][1]} in range({leafl[1]},np))'
            elif oper=='***': return F'prodl({leafl[2]} for {leafl[0][1]} in range({leafl[1]},np))'
            else:
                raise NotImplementedError(F'TODO: Expr.__format__() for operator {oper}')
        else:
            raise ValueError(F'unrecognized operator {oper} in Expr.__format__()')

    def __eq__(self, other):
        '''

        :param other:
        :return:
        '''

        return F'E{self}' == F'E{other}'
    def __ne__(self, other):

        return F'E{self}' != F'E{other}'
    # def __copy__(self):
    #     '''
    #     shallow copy, shared IP
    #
    #     NOTE: =, copy and deepcopy are different
    #           1)        =:  all value are shared ID
    #                  copy:  list value are shared ID
    #              deepcopy:  all value are not shared
    #           2) overload of = operator don't be supported in Python due to it similar to labeling, shard ID
    #
    #     :return:
    #     '''
    #
    #     return deepcopy(self)
    # def __deepcopy__(self, memodict={}):
    #     '''
    #     deepcopy, not shared IP
    #
    #     :return:
    #     '''
    #
    #     # return deepcopy(self)
    #     tmp = Expr(oper=self.oper)
    #
    #     from collections import Iterable
    #     for i,leaf in enumerate(self.leafl):
    #         if leaf.__class__ in [list, str, tuple]: tmp.leafl[i] = leaf
    #         elif isinstance(leaf, Iterable):
    #             '''
    #             Bug1: leaf may not be [], let's assume that
    #             Bug2: elem in leaf also may be iterable
    #
    #                 Note: fateful bug in Python language
    #                       deepcopy are not real copy
    #                       >>> l1= [2, 3]
    #                       >>> l2 = [l1, 'di', l1]
    #                       >>> l2_ = deepcopy(l2)
    #                       >>> l2_[0][0] = 9
    #                       >>> l2_
    #                           [[9, 3], 'di', [9, 3]]
    #             '''
    #             tmp.leafl[i] = [deepcopy(elem) for elem in leaf]
    #         else: tmp.leafl[i] = deepcopy(leaf)
    #
    #     return tmp
    def __hash__(self):
        return hash(F'E{self}')

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

    def isroot(self):
        '''
        is root ?

        :return:
        '''

        return not self.root
    def istrunk(self):
        '''
        is trunk ?

        :return:
        '''

        if self.root:
            for leaf in self.leafl:
                if isinstance(leaf, Expr): return True
            else: return False
        else: return False
    def isleaf(self):
        '''
        is leaf ?

        :return:
        '''

        # flag = False
        #
        # for leaf in self.leafl:
        #     if not isinstance(leaf, Expr):
        #         flag = True
        #         break
        #
        # return flag
        return all(not isinstance(leaf, Expr) for leaf in self.leafl)
    def nary(self):
        return nary(self.oper)
    def expand(self):
        raise NotImplementedError('TODO: Expr.exppand()')
    def simplify(self):
        raise NotImplementedError('TODO: Expr.simplify()')

class Coeff(Symb):
    '''
    Coeff class

        such as: Configuration interaction ci, Amplitude t
    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='C', l=symbl, f=mode) # default [] has bug

    def __str__(self): return F"{{{self.f}}}_{{{','.join(F'{e!s}' for e in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        symbl,mode = self.l,self.f; nsymb = len(symbl)
        if mode is 't':
            return F"t[({''.join(F'({b.l},{b.f}),' for b in symbl)})]"
            # return F"t[ud[phiu({{{','.join(F'{b.l},{b.f}' for b in symbl)}}})]]"  # C++ code
        elif mode is 'r':
            try:
                return F'r[{symbl[0].l}][{symbl[1]}][{symbl[2]},{symbl[3]}]'
                # return F'r[rdm({symbl[0].l},{symbl[1]},{symbl[2]},{symbl[3]})]'  # C++ code
            except AttributeError: return F'r[{symbl[0]}][{symbl[1]},{symbl[2]}]'
        else: return F"{self.f}[{symbl[0]},{symbl[1]}]"

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"C{self}" == F"C{other}"
    def __ne__(self, other):
        return F"C{self}" != F"C{other}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"C{self}")
class Array(Symb):
    '''
    Array class

        such as: Representive matrix r, Intermediate array a
    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='A', l=symbl, f=mode)

    def __str__(self): return F"{{{self.f}}}_{{{','.join(F'{b!s}' for b in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        n = self.f[1:]
        return F"a[{n}]{''.join(F'[{b.l}]' for b in self.l)}"

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'A{self}' == F'A{other}'
    def __ne__(self, other):
        return F'A{self}' != F'A{other}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"A{self}")
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
        ol = {i: F'{symb.l}^{{{upl[i]}}}' for i,symb in enumerate(self.l)}
        if f is 2: return F'\\langle{{{ol[0]}{ol[1]}}}\\vert \\hat{{g}} \\vert{{{ol[2]}{ol[3]}}}\\rangle' # 公式输出时会自动换行 ？
        elif f is 0: return F'\\langle{{{ol[0]}}}\\vert \\hat{{h}} \\vert{{{ol[1]}}}\\rangle'
        elif f is 3: return F'\\langle{{{ol[0]}{ol[1]}}}\\vert \\hat{{g}} \\vert{{{ol[3]}{ol[2]}}}\\rangle'
        elif f is 4:
            return '\\langle{{{0}{1}}}\\vert \\hat{{g}} \\vert{{{2}{3}}}\\rangle'\
                   '-\\langle{{{0}{1}}}\\vert \\hat{{g}} \\vert{{{3}{2}}}\\rangle'.format(*ol.values())
        elif f is 1: return F'\\langle{{{ol[0]}}}\\vert \\hat{{f}} \\vert{{{ol[1]}}}\\rangle'
        else:
            raise ValueError(F'unrecognized {f} for Symb.f in Symb.__str__()')
    __repr__ = __str__
    def __format__(self, format_spec):
        f = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        ol = {i: F'p[{upl[i]},{symb.l}]' for i, symb in enumerate(self.l)}
        # ol = {i: F'po({upl[i]},{symb.l})' for i, symb in enumerate(self.l)}  # C++ code
        if f is 2:
            return F'g[dei({ol[0]},{ol[2]},{ol[1]},{ol[3]})]'
            # return F'g[dei({ol[0]},{ol[2]},{ol[1]},{ol[3]})]'  # C++ code
        elif f is 0:
            return F'h[{ol[0]},{ol[1]}]'
            # return F'h[sei({ol[0]},{ol[1]})]'  # C++ code
        elif f is 3: return F'g[dei({ol[0]},{ol[3]},{ol[1]},{ol[2]})]'
        elif f is 4: return F'(g[dei({ol[0]},{ol[2]},{ol[1]},{ol[3]})]-g[dei({ol[0]},{ol[3]},{ol[1]},{ol[2]})])'
        elif f is 1: return F'f[{ol[0]},{ol[1]}]'
        else:
            raise ValueError(F'unrecognized Symb.f = {f} for Integral in Symb.__format__()')

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"I{''.join(b.l for b in self.upl)}{''.join(o.l for o in self.l)}{self.f}" \
               == F"I{''.join(b.l for b in other.upl)}{''.join(o.l for o in other.l)}{other.f}"
    def __ne__(self, other):
        return F"I{''.join(b.l for b in self.upl)}{''.join(o.l for o in self.l)}{self.f}" \
               != F"I{''.join(b.l for b in other.upl)}{''.join(o.l for o in other.l)}{other.f}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"I{''.join(b.l for b in self.upl)}{''.join(o.l for o in self.l)}{self.f}")

    @property
    def orbl(self): return self.l
    @property
    def norb(self): return len(self.l)

    def orbl_dedu(self):
        '''
        symbl of de-duplication

        :return:
        '''

        return tuple(sorted(set(self.l)))
    def orbsort(self):
        symbl, upl = self.l, self.upl; norb = len(symbl)
        od = {i:(upl[i].l, symbl[i].l) for i in range(norb)} # sort orb: I0，I1, A0, A1
        pr = list(range(0, norb, 2)); qs = list(range(1, norb, 2))
        pr = sorted(pr, key=lambda x:od[x]); qs = sorted(qs, key=lambda x:od[x])
        if [od[i] for i in qs]<[od[i] for i in pr]: pqrs = list(chain.from_iterable(zip(qs, pr)))
        else: pqrs = list(chain.from_iterable(zip(pr, qs)))
        self.upl=tuple(upl[i] for i in pqrs); self.l = tuple(symbl[i] for i in pqrs)
    def blockl_dedu(self):
        '''
        upl of de-duplication

        :return:
        '''

        return tuple(sorted(set(self.upl)))
    def blockize(self):
        raise NotImplementedError('TODO: Integral.blockize()')
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
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Functor.__format__()')

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'F{self.l}{self.upl}{self.f}' == F'F{other.l}{other.upl}{other.f}'
    def __ne__(self, other):
        return F'F{self.l}{self.upl}{self.f}' != F'F{other.l}{other.upl}{other.f}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F'F{self.l}{self.upl}{self.f}')

    @property
    def symbl(self): return self.l
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

    def symbl_dedu(self):
        '''
        symbl of de-duplication

        :return:
        '''

        return tuple(sorted(set(self.l)))
    def blockl_dedu(self):
        '''
        upl of de-duplication

            NOTE: 只可能在不同block, 不可有相同block(不同state)
        :return:
        '''

        return tuple(sorted(set(self.upl)))
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
                A = iis(upl[isymb], Al)
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


    def __str__(self): return F'{{{self.l}}}_{{{self.f}}}'
    __repr__ = __str__
    def __format__(self, format_spec): return F'{self.l}{self.f}'

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'O{self}' == F'O{other}'
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F'O{self}' != F'O{other}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'O{self}' > F'O{other}'
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F'O{self}' >= F'O{other}'
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:
        :return:
        '''

        return F'O{self}' < F'O{other}'
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F'O{self}' <= F'O{other}'
    def __hash__(self):
        '''
        hash: used in set

        Note: hash函数只可定义给不可变对象，而该对象是一个可变对象，定义hash实际上是不合规矩的，这是只是为了方便集合操作
        :return:
        '''

        return hash(F'O{self}')
# space orbital
nho_ = 4
hop = Orbital(ordinal='p')
hoq = Orbital(ordinal='q')
hor = Orbital(ordinal='r')
hos = Orbital(ordinal='s')
ho_l = (hop, hoq, hor, hos)
# spin orbital
nho = 8; spinl = ('\\alpha', '\\beta')
hopa = Orbital(ordinal='p', spin='\\alpha')
hopb = Orbital(ordinal='p', spin='\\beta')
hoqa = Orbital(ordinal='q', spin='\\alpha')
hoqb = Orbital(ordinal='q', spin='\\beta')
hora = Orbital(ordinal='r', spin='\\alpha')
horb = Orbital(ordinal='r', spin='\\beta')
hosa = Orbital(ordinal='s', spin='\\alpha')
hosb = Orbital(ordinal='s', spin='\\beta')
hol = (hopa, hopb, hoqa, hoqb, hora, horb, hosa, hosb)
# in one Block
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

        return F'B{self}' == F'B{other}' # Python Bug: is 有时是 False, 有时是 True
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F'B{self}' != F'B{other}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'B{self}' > F'B{other}'
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F'B{self}' >= F'B{other}'
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:
        :return:
        '''

        return F'B{self}' < F'B{other}'
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F'B{self}' <= F'B{other}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F'B{self}')

    @property
    def ne(self): return bfd[self.f].ne
    @property
    def sz(self): return bfd[self.f].sz
# ground (general, any) block labeled as IJKL
ngb = 4; bgl = 'IJKL'; b0l = set(bgl)
gbI = Block(ordinal='I')
gbJ = Block(ordinal='J')
gbK = Block(ordinal='K')
gbL = Block(ordinal='L')
gbl = (gbI, gbJ, gbK, gbL)
# excited block labeled as ABCD
neb = 2; bel = 'ABCD'; b1l = set(bel)
ebA = Block(ordinal='A', state=1)
ebB = Block(ordinal='B', state=1)
ebC = Block(ordinal='C', state=1)
ebD = Block(ordinal='D', state=1)
ebl = (ebA, ebB, ebC, ebD)
''' wave level '''
class Wave(Symb):
    '''
    wave(wavefunction) Symb

    :param ordinal: str,
                                '':  unkown
                        X, Y, Z, W:  wave ordinal, Only one in general
    :param state: str or int or block[]
                                '':  unkown
                                 0:  ground state
                           block[]:  excited state
    :param nblock:
    '''

    def __init__(self, ordinal='\\Phi', state=()): Symb.__init__(self, c='W', l=ordinal, f=state)

    def __str__(self): return F"\\Phi_{{{','.join(F'{b!s}' for b in self.f)}}}"
    __repr__ = __str__
    def __format__(self, format_spec): return F"({''.join(F'({b.l},{b.f}),' for b in self.f)})"

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' == F'phi{other.f}'
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' != F'phi{other.f}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' > F'phi{other.f}'
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' >= F'phi{other.f}'
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' < F'phi{other.f}'
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F'phi{self.f}' <= F'phi{other.f}'
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F'phi{self.f}')

    @property
    def ne1(self): return sum(b.ne for b in self.f)
    @property
    def sz(self): return sum(b.sz for b in self.f)
    @property
    def nb1(self): return len(self.f)
    def sortblock(self): self.f = tuple(sorted(self.f))

    def conjugate(self):
        '''
        conjugated form of wavefunction, bra <-> ket

        :return:
        '''

        raise NotImplementedError('TODO: Wave.conjugate()')
    def groundize(self):
        '''
        excited state -> ground state

        :return:
        '''

        raise NotImplementedError('TODO: Wave.groundize()')
    def excite(self):
        '''
        ground state -> excited state

        :return:
        '''

        raise NotImplementedError('TODO: Wave.excite()')
    def deexcite(self):
        '''
        one excited state -> other state

        :return:
        '''

        raise NotImplementedError('TODO: Wave.deexcite()')
wPhi = Wave(state=())




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
''' Symb: 16 states in block '''
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
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Subbfeq.__format__()')
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

    def __str__(self): return F'{{{self.ce!s}{self.functor!s}}}'
    __repr__ = __str__
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Subbseq.__format__()')

    def indexorb(self, orb): return iis(orb, self.functor.l)
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

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 'base','Foundation of GVB-BCCC formula derivation'
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
        write(F"    A block has {nbs} states: ${','.join(F'{b!s}' for b in bsd.values())}. $ \\\\ \n")
        # block functor
        write( '    Each state corresponds to a functor in essence: \\\\ \n')
        for s,f in bfd.items():
            write(F'    {s}: ${f!s} $ \\\\ \n')
        write( '    \n')

        # block functor of linear combination of state
        write( '    However, a functor is often expressed as a linear combination of a number of states: \\\\ \n')
        for i,(f,bfeq) in enumerate(bfeqd.items()): write(F"    {i}: ${f!s} = {'+'.join(F'{e!s}' for e in bfeq)} $ \\\\ \n")
        write( '    \n')

        # block state of linear combination of functor
        write( '    Therefore, the composition of each state can be written: \\\\ \n')
        for s,bseq in bseqd.items():
            write(F"    ${bsd[s]!s} = {'+'.join(F'{e!s}' for e in bseq)} \\quad \\vert{{vac}}\\rangle $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
qd = {i: Frac(1, math.factorial(i)) for i in range(ngb+neb+1)}
nel = [0,2] # for no T1, list(range(neb+1))
ngl = list(range(neb+ngb+1))
sld = {n: set(sl for sl in product(range(1,nbs), repeat=n)
              if sum(bfd[s].ne for s in sl) is n*2 and sum(bfd[s].sz for s in sl) is 0) for n in nel}




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
                    iborb = iis(orb, fl)
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

    # repressive functor status
    nrfsl = 4; rfsll = (('+', '-'), ('++', '+-', '--'), ('++-', '+--'), ('++--',))
    # nrfs=8; rfsl=['+', '-', '++', '+-', '--', '++-', '+--', '++--']
    nrf=0; ab13={1, 3}
    
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

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 'rep','Repressive matrix of functor in GVB-BCCC formula derivation'
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
        for i,(f,r) in enumerate(rfeqd.items()):
            write(F'    $\\hat{{O}}_{{{i}}}:  \\langle{{P_p}}\\vert {f!s} \\vert{{P_q}}\\rangle => $ \\\\ \n')
            for q,cd in r.items():
                write(F"    $ {f!s} \\vert{{P_{{{q}}}}}\\rangle = {'+'.join(F'{c!s}P_{{{p}}}' for p,c in cd.items())} $ \\\\ \n")
                for c in cd.values():
                    write(F'    ${c!s}\\ =\\ {rd[c]!s} $ \\\\ \n')
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

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre = 'rep'
    with open(F'bccc{ossep}{fpre}.py', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        write('import numpy\n')
        write('\n\n')

        write(F'nbs = {nbs}\n')
        write(F'nrf = {len(rfeqd)}\n')
        write('\n\n')

        write('def rep(ci):\n')
        write('    r = {}\n')
        write('    inv = numpy.linalg.inv(ci)\n')
        write('    \n')

        for i,(f,rfeq) in enumerate(rfeqd.items()):
            write(F'    # O{i}:  <Bp|  {f!s}  |Bq> = \n')
            if any(rfeq.values()):
                write('    m = numpy.zeros(shape=(nbs,nbs)) \n')
                for q,cxbld in rfeq.items():
                    for p,c in cxbld.items():
                        write(F'    m[{p},{q}] = {rd[c]}\n')
                write(F'    r[{i}] = m\n')
            write('    \n')
        write('    \n')

        write('    return r\n')
        write('    \n')

        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




''' wave state for n-block excite '''
ul = () # wave state
uld = {} # wave state (n-block excite)
def wexcite():
    '''
    wave excite

    Note:
    :return:
    '''

    global ul,uld

    for n,sl in sld.items():
        sl = set(tuple(sorted(s)) for s in sl)  # less nu
        ul=[]; append=ul.append
        for s in sorted(sl):
            u = deepcopy(ebl[:n])
            for b,i in zip(u,s): b.f = i
            append(u)

        uld[n] = tuple(ul)

    ul = tuple(chain.from_iterable(uld.values()))
def w2theory():
    '''
    wave (ket part) -> tex file

    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 'wave','Reference wave function in GVB-BCCC formula derivation'
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

        write('    \\newpage\n')
        write('    \\setlength{\parindent}{0pt}\n')
        write('    \\setcounter{page}{1}\n')
        write('    \n')

        write( '    To derive the GVB-BCCC formula, all possible excited states need to be constructed firstly. \\\\ \n')
        write(F'    For N blocks system, the number of wavefunction state is \\\\ \n')
        write(F"    ${'+'.join(F'{len(l)}C^N_{n}' for n,l in uld.items())}  $ \\\\ \n")
        write( '    including 1 ground state, and other states are excited one. \\\\ \n')
        write(F'    In particular, For N-block $(N\\le4 )$ system, the number of wavefunction state can be simply recorded as $ (C^N_{{2N}})^2 $. \\\\ \n')
        write( '    \n')

        write( '    Ground state wavefunction: \\\\ \n')
        write(F'    $\\Phi_{{{ul[0]!s}}}$ \\\\ \n')
        write( '    \n')

        write('    Excited state wavefunction: \\\\ \n')
        for nb in range(1, neb+1):
            write(F'    For {nb}-block excitation, there are $ {len(uld[nb])}C^N_{nb} $ excited types: \\\\ \n')
            for s in uld[nb]:
                write(F'    $\\Phi_{{{s!s}}}$ \\\\ \n')
            write('    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# H (Hamiltonian) functor
hd = {} # H types
hnd = {} # H (one- or two-electron)
def deth():
    raise NotImplementedError('TODO: deth()')
def h2theory():
    raise NotImplementedError('TODO: h2theory()')

# Hamiltion matrix
hfeqd = {} # 这种简单方式比 Manager.dict 更快
# hfeqd = Manager().dict()  # for parallel
def mhuv2theory(bra, hld):
    '''
    Matirx of hamiltionian element (i,j) -> latex

    Note:    subfunction
    :param bra:
    :param ket:

    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}",'Hamiltonian matrix in GVB-BCCC formula derivation'
    with open(F'doc{ossep}{fpre}.tex', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write(F'% {ftitle}\n')
        write( '\n\n')

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

        write('    \\newpage\n')
        write('    \\setlength{\parindent}{0pt}\n')
        write('    \\setcounter{page}{1}\n')
        write('    \n')

        write(F'    $\\hat{{H}}\\vert{{\\Phi_{{ {bra} }} }}\\rangle = $ \\\\ \n')
        sketl = tuple(F'c_{{ {bra}, {ket} }}\\Phi_{{{ket}}}' for ket in hld)
        write(F"    $+{'+'.join(sketl)} $ \\\\ \n")
        write( '    \n')

        write( '    Specifically, \\\\ \n')
        for ket in hld:
            if bra==tuple(b for b in ket if b.f>0):
                write(F'    $c_{{{bra}, {ket}}} = \\langle\\Phi_{{{bra}}}\\vert \\hat{{H}} \\vert\\Phi_{{{ket}}}\\rangle = $ \\\\ \n')
                write(F'    ${hld[ket]!s} $ \\\\ \n')
                write( '    \n')
        write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def mhuv2code(bra, hld):
    '''
    Matirx of hamiltionian element (i,j) ->  code

    Note:     subfuction
    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre = F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}"
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        write('import numpy\n')
        write('from bccc.pub import sign,dei\n')
        write('\n\n')

        tab = '    '
        # bbl,kbra = [b.l for b in bra],F"({','.join(F'{b}' for b in bra)})"
        bbl,kbra = [b.l for b in bra],F"tuple(sorted([{','.join(F'{b}' for b in bra)}]))"; sbbl = [b.l for b in bra if 6<b.f<15]
        
        fl=[]; flapp=fl.append
        for v,hl in hld.items():
            f = F"_{''.join(F'{b.l}{b.f}' for b in v)}({''.join(F'{b}, ' for b in bbl)}r, h, g, p, nc)"
            # f = F"_{''.join(F'{b.l}{b.f}' for b in v)}({''.join(F'{b}, ' for b in bbl)}r, mf)"
            flapp(f)
            
            kbl,kket = [b.l for b in v if b.l in b0l],F"tuple(sorted([{','.join(F'{b}' for b in v)}]))"; skbl = [b.l for b in v if 6<b.f<15]
            
            write(F"def {f}:\n")
            # write( '    h,g,p,nc=mf.h,mf.g,mf.p,mf.nc\n')
            write( '    np = len(p)\n')
            write( '    hd = {}\n')
            if sbbl: write(F"    sbra = sign([{','.join(sbbl)}])\n")
            write( '    \n')

            nkt,tabnkt = 0,''
            for b in kbl:
                # if nkt>0: write(tabnkt+F'    for {b} in range({kbl[nkt-1]}+1, np):\n')
                # else: write(F'    for {b} in range(nc, np):\n')
                # write(tabnkt+F"        if {b} in {{{','.join(bbl+kbl[:nkt])}}}: continue\n")
                write(tabnkt+F'    for {b} in range(nc, np):\n')
                if bbl or nkt: write(tabnkt+F"        if {b} in {{{','.join(bbl+kbl[:nkt])}}}: continue\n")
                nkt += 1
                tabnkt += tab
            write(tabnkt+F'    v = {kket}\n')
            write(tabnkt+ '    hv = 0.0\n')
            if skbl: write(tabnkt+F"    sket = sign([{','.join(skbl)}])\n")
            write(tabnkt+ '    \n')

            for h in hl:
                Il = sorted(set([b.l for b in h[1].upl if b.l in b0l and b.l not in kbl]))
                if Il:
                    # r1l,r0l = [h[0]],[h[1]]
                    # for r in h[2:]:
                    #     if r.l[0].l in Il: r0l.append(r)
                    #     else: r1l.append(r)
                    forl=[]; forlapp=forl.append
                    nI = 0
                    for b in Il:
                        # if nI: forlapp(F"for {b} in range({Il[nI-1]},np) if {b} not in {{{','.join(bbl+kbl+Il[:nI])}}}")
                        # else: forlapp(F"for {b} in range(np) if {b} not in {{{','.join(bbl+kbl)}}}")
                        forlapp(F"for {b} in range(np) if {b} not in {{{','.join(bbl+kbl+Il[:nI])}}}")
                        nI += 1
                    write(tabnkt+F"    hv += sum({prodl(h)} {' '.join(forl)})\n")
                    # write(tabnkt+F"    hv += {prodl(r1l)}*sum({prodl(r0l)} {' '.join(forl)})\n") # 并没有节省时间在python, 因为 sum 作个优化
                else: write(tabnkt+F'    hv += {prodl(h)}\n')
            write(tabnkt+ '    \n')

            # write(tabnkt+ '    hd[v] = hd.get(v,0.0)+hv\n')
            if skbl: write(tabnkt+ '    hd[v] = hd.get(v,0.0)+sket*hv\n')
            else: write(tabnkt+ '    hd[v] = hd.get(v,0.0)+hv\n')
            write(tabnkt+ '    \n')
            
            if sbbl:
                write( '    if sbra is -1:\n')
                write( '        for v,h in hd.items(): hd[v] = -h\n')
                write( '    \n')
            
            write(F'    return hd\n')
            write( '    \n')
        
        write(F'def {fpre}(r, h, g, p, nc=0):\n')
        write( '    np=len(p)\n')
        # write(F'def {fpre}(r, mf):\n')
        # write( '    nc,np=mf.nc,len(mf.p)\n')
        write( '    hdd = {}\n')
        write( '    \n')

        nbt,tabnbt = 0,''
        for b in bbl:
            # if nbt>0: write(tabnbt+F'    for {b} in range({bbl[nbt-1]}+1, np):\n')
            # else: write(F'    for {b} in range(nc, np):\n')
            write(tabnbt+F'    for {b} in range(nc, np):\n')
            if nbt: write(tabnbt+F"        if {b} in {{{','.join(bbl[:nbt])}}}: continue\n")
            nbt += 1
            tabnbt += tab
        write(tabnbt+F'    u = {kbra}\n')
        write(tabnbt+ '    \n')

        write(tabnbt+ '    hd = {}\n')
        for f in fl: write(tabnbt+F'    for v,hv in {f}.items(): hd[v] = hd.get(v,0.0)+hv\n')
        write(tabnbt+ '    \n')
        
        write(tabnbt+ '    hdd[u] = hd\n')
        write(tabnbt+ '    \n')

        write('    return hdd\n')
        write('    \n\n')
def h2w(w):
    '''
    hamiltonnian project into wave

    :param w: wave state, list

    :return:
    '''

    hld = {}; vd = {} # {v: [[c, <>, r..]]}
    nb = len(w); ne = sum(b.ne for b in w)
    bl = tuple(b.l for b in w)

    # for h
    if nb<=2 and ne>=1:  # nb<2 for h; p+p- ne>=1
        for prodb in product(w, repeat=2):
            if len(set(prodb)) is not nb: continue
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    inte = Integral(symbl=(bop, boq), mode=0)
                    functor = Functor(symbl=(bop, boq), pdl='+-', upl=prodb)
                    s_blockize,functorl,upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                    _hld={}; _vd={} # ce {v: r} 没有与 inte 相乘
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
                                v = deepcopy(upl)
                                for i,b in zip(bsl,v): b.f = i
                                if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                _rl=[]; append=_rl.append
                                for b,r in zip(v, rl):
                                    _r = deepcopy(r)
                                    _r.l = (b,)+_r.l
                                    append(_r)
                                try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                except KeyError: _vd[v],_hld[v] = v,[_rl]
                    for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                        _inte = deepcopy(inte)
                        _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                        # _inte.orbsort() # call inte.sort() in detmh
                        si = [Frac(s_blockize,1), _inte]
                        try:
                            for rl in _hl: hld[v].append(si+rl)
                        except KeyError:
                            for i,rl in enumerate(_hl): _hl[i] = si+rl
                            vd[v],hld[v] = v,_hl

    # for g
    if ne>=2: # p+p+s-r- ne>=2
        for prodb in product(w, repeat=4):
            if len(set(prodb)) is not nb: continue
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    for bor in deepcopy(bo_l):
                        for bos in deepcopy(bo_l):
                            inte = Integral(symbl=(bop, boq, bor, bos), mode=2)
                            functor = Functor((bop, boq, bos, bor), pdl='++--', upl=(prodb[0], prodb[1], prodb[3], prodb[2]))
                            s_blockize, functorl, upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                            _hld = {}; _vd={}  # ce 没有与 inte 相乘
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
                                        v = deepcopy(upl)
                                        for i,b in zip(bsl, v): b.f = i
                                        if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                        elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                        _rl=[]; append=_rl.append
                                        for b,r in zip(v, rl):
                                            _r = deepcopy(r)
                                            _r.l = (b,)+_r.l
                                            append(_r)
                                        try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError: _vd[v],_hld[v] = v,[_rl]
                            for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                                # _inte.orbsort() # call inte.sort() in detmh
                                si = [Frac(s_blockize, 2), _inte]
                                try:
                                    for rl in _hl: hld[v].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(_hl): _hl[i] = si+rl
                                    vd[v],hld[v] = v,_hl

                            _hld = {}; _vd = {}  # ce 没有与 inte 相乘
                            a,b = spinl
                            bop.f = bor.f = a; boq.f = bos.f = b
                            bfeql = [None]*nb
                            for i,f in enumerate(functorl):
                                bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                else: break
                            else:
                                for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                    bsl,rl = zip(*prodcxb)
                                    v = deepcopy(upl)
                                    for i,b in zip(bsl, v): b.f = i
                                    if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                    elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                    _rl = []; append=_rl.append
                                    for b,r in zip(v, rl):
                                        _r = deepcopy(r)
                                        _r.l = (b,)+_r.l
                                        append(_r)
                                    try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                    except KeyError: _vd[v], _hld[v] = v, [_rl]
                            for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                                # _inte.orbsort() # call inte.sort() in detmh
                                si = [Frac(s_blockize, 1), _inte]
                                try:
                                    for rl in _hl: hld[v].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(_hl): _hl[i] = si+rl
                                    vd[v],hld[v] = v,_hl

    # for v,hl in hld.items():
    #     hld[v] = sum(prodl(h) for h in hl)
    
    return hld
def detmh(bra, hfeqd):
    '''
    determine matrix of hamiltionian (project into wave)

    Note:       Hamiltonian matrix is complex conjugate. This point is often used in the process of deriving formulas
                <W| H  =>  H |W>
                when it stored and output, <bra| H |ket>

    :return:
    '''

    print('bra = ', bra)
    nb1 = len(bra)
    # cbra = qd[nb1]  # cbra (for AB) 本应被乘, 但为减少 nu/iteration, 顺序化 bs (但若 bs 相同, 需在 code 中 特殊处理), 所以不乘

    # mH
    hld = {} # {v: [[c, <>, r..]...]}
    for nbc in range(1, ngb+1):  # n blocks correlation
        for nbe in range(nbc+1):  # nbe blocks from bra
            nb0=nbc-nbe; cket=qd[nb0]  # cbra (for IJ..) 本应被乘, 顺序化 bs 后不乘
            for bchose in combinations(bra, r=nbe):  # 从激发态取 nbe
                brest = tuple(sorted(set(bra)-set(bchose)))  # set -> list is random and disordered
                s_bef = sign([b.l for b in bchose+brest if b.ne%2])  # H 作用前, 激发块提到前面, code 中会设置 A<B 所以无需表达式
                # s_bef = Expr('-?', [b.l for b in bchose+brest if b.ne%2])  # H 作用前, 激发块提到前面, code 中 AB 没有大小限制, 需表达式 (默认A<B, 若不是添负号)
                _u = bchose+gbl[:nb0]
                # if len(_u+brest)>neb: continue

                sl = tuple(b.f for b in _u); bl = tuple(b.l for b in _u)
                if sl in hfeqd:
                    _hld = deepcopy(hfeqd[sl])
                    try:
                        _bl = tuple(b.l for b in list(_hld.keys())[0])
                        if _bl != bl:
                            for v in _hld:
                                for b,i in zip(v, bl): b.l = i
                    except IndexError: continue  # _hld = {} when nbc=nbe, b.f=6
                else:
                    _hld = h2u(_u)  # Aa I0 J0 Kk
                    hfeqd[sl] = _hld
                    _hld = deepcopy(_hld)

                for v,hl in _hld.items():  # 将 _hld 更新(多块加和)到 hld 中
                    nI1 = len([b for b in v[nbe:] if b.f])
                    if nb0>1 and nI1 is not nb0:  # I0Ja -> J0Ia
                        i1 = 0; i0 = nI1
                        for b in v[nbe:]:
                            if b.f:
                                b.l = bgl[i1]
                                i1 += 1
                            else:
                                b.l = bgl[i0]
                                i0 += 1
                    
                    _u = v+deepcopy(brest); u = tuple(sorted(b for b in _u if b.f))
                    # bsl = [b.f for b in u[len(bra):] if b.f]  # bs for IJ...
                    # if bsl != sorted(bsl): continue  # 如果 bs 不顺序化, 丢弃
                    
                    ''' only for LCC '''
                    if len(u)>neb: continue
                    
                    s_lat = sign([b.l for b in _u if b.ne%2])  # H 作用后, 激发块回到原来位置，code 中会设置 A<B,I<J 所以无需表达式
                    # s_lat = Expr('-?', [b.l for b in _u if b.ne%2])  # H 作用后, 激发块回到原来位置，code 中 AB..,IJ.. 没有大小限制, 需表达式 (默认A<B..I<J.., 若不是添负号)
                    for h in hl:
                        # h[0] = s_bef*s_lat*h[0]  # 不乘1/n, 假设 A<B..., I<J.. (或者 顺序化 bs for AB, IJ..)
                        # h[0] = s_bef*s_lat*cbra*cket*h[0]   # 乘1/n (分别对 AB,IJ..) 没有大小限制
                        h[0] = s_bef*s_lat*cket*h[0]  # 为减少 nu/iteration, 不乘 cbra
                        # h[0] = s_bef*s_lat*h[0]  # 假定 A<B..., I<J. 真实的大小 体现在 code
                        # h[1].orbsort()
                    hld[u] = hld.get(u, [])+hl
            
    # mhuv2theory(bra, hld)
    mhuv2code(bra, hld)
def mh2theory():
    '''
    matrix of hamiltonian (project into wave) -> tex file

    Note:    main function output Hamiltionian and H|W> = ?
    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 'hm','Hamiltonian matrix in GVB-BCCC formula derivation'
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

        write( '    \\linespread{2.0}\selectfont\n')
        write( '    \\thispagestyle{empty}\n')
        write( '    \n')

        write( '    \\newpage\n')
        write( '    \\setlength{\parindent}{0pt}\n')
        write( '    \\setcounter{page}{1}\n')
        write( '    \n')

        write( '    To derive the GVB-BCCC formula, all possible types of excited states need to be constructed firstly. \\\\ \n')
        write(F'    In total, there are 1 ground state and {len(ul)-1} excited types. \\\\ \n')
        write( '    \n')

        write( '    When $\\hat{H}$ operator project onto a wavefunction state, it will be converted to a linear combination of other states. that is \\\\ \n')
        write( '    $\\hat{H}\\vert{\\Phi_{w}}\\rangle = \\displaystyle\\sum_{u,v}{c_{u,v}\\Phi_{v}}$ \\\\ \n')
        write( '    where $c_{u,v} $ is Hamiltonian matrix, and listed separately in other files.\\\\ \n')
        write( '    \n')

        write( '    Specifically, state $ u, v = $ \\\\ \n')
        # write(F"    ${','.join(F'{u}' for i,u in ul)} $ \\\\ \n")
        for u in ul: write(F'    ${u} $ \\\\ \n')
        write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def mh2code(lang='Python'):
    '''
    matrix of Hamiltonian (project into wave) -> code

    main function calculate Hamiltionian matrix
    :param lang: str, one of Python, C++ and C. default is Python

    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre = 'hm'
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        write('from multiprocessing import Pool, Manager\n')
        write('from importlib import import_module\n')
        write('\n\n')

        write('def cal(braf, r, h, g, p, nc=0):\n')
        write("    module = F'bccc.hm.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    hd = fun(r, h, g, p, nc)\n')
        write('    return hd\n')
        write('\n\n')

        write('def hm(r, h, g, p, nt=4, nc=0):\n')
        # write('def hm(r, mf):\n')
        write('    hdd = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        write('    from bccc.hm.hm_ import hm_\n')
        write('    brafl = (hm_,)\n')
        for nt in nel[1:]:
            write(F'    # {nt}-block excited types\n')
            brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[nt])
            for braf in brafl: write(F'    from bccc.hm.{braf} import {braf}\n')
            fll = [brafl[10*i:10*(i+1)] for i in range(0, math.ceil(len(brafl)/10))]
            write( '    brafl += (\n')
            for fl in fll: write(F"            {','.join(fl)},\n")
            write( '            )\n')
        write( '    \n')

        # write('    for braf in brafl: hdd.update(calculate(braf, r, h, g, p))\n')
        # write('    for braf in brafl:\n')
        # write('        for u,hu in calculate(braf, r, h, g, p).items():\n')
        # write('            try:\n')
        # write('               for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        # write('            except KeyError: hdd[u] = hu\n')
        write('    pool = Pool()\n')
        write('    pl = tuple(pool.apply_async(braf, (r, h, g, p, nc)) for braf in brafl)\n')
        # write('    pl = tuple(pool.apply_async(braf, (r, mf)) for braf in brafl)\n')
        write('    pool.close()\n')
        write('    pool.join()\n')
        write('    for p in pl: hdd.update(p.get())\n')
        # write('    for p in pl:\n')
        # write('        for u,hu in p.get().items():\n')
        # write('            try:\n')
        # write('                for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        # write('            except KeyError: hdd[u] = hu\n')
        
        write('    \n')
        write('    return hdd\n')
        write('    \n\n')

        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()

# T (Taylor) funcotr
tl = () # T types
tld = {} # T (n-block excite)
def dett():
    '''
    T operator

    :return:
    '''

    global tl, tld

    for n,sl in sld.items():
        ul = []; append=ul.append
        for s in sorted(sl):
            u = deepcopy(gbl[:n])
            for b,i in zip(u,s): b.f = i
            append(u)
        tld[n] = tuple(ul)

    tl = tuple(chain.from_iterable(tld.values()))
def t2theory():
    '''
    T operator -> tex file

    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 't',' T operator '
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

        write( '    Similar to single reference coupled cluster (SRCC) theory, the T operator of GVB-LBCCC theory are block-based \\\\ \n')
        write(F"    $T = {'+'.join(F'T_{i+1}' for i in range(neb))} $ \\\\ \n")
        write( '    where: \\\\ \n')
        for i in tld:
            bl = tuple(b.l for b in tld[i][0]); bl_ = tuple(b.lower() for b in bl)
            sigmal = tuple("\\displaystyle\\sum_{{{}}}\\displaystyle\\sum_{}".format('>'.join(bl[j::-1]), bl_[j]) for j in range(len(bl)))
            operl = tuple(F'{bl[j]}_{bl_[j]}^+{bl[j]}_0^-' for j in range(len(bl)))
            write(F"    $T_{i} = {''.join(sigmal)}{{{''.join(operl)}}} $ \\\\ \n")
            for bl in tld[i]:
                write(F"    ${','.join(F'{b}={bl[ib].f}' for ib,b in enumerate(bl_))} $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()

# LCC
mt = {}
def detmt():
    '''
    determine matrix of T

        typical way:
                    Ii+:  <Phi_p|  Ii+     |Phi_q>
                 Ii+Jj+:  <Phi_p|  Ii+Jj+  |Phi_q>
          novel way:
                    Oo-:  Oo- |Phi_p>
        For T=T1+T2:
                    Ii+:  |Phi_p> = p+| >
                 Ii+Jj+:  |Phi_p> = p+| >
                          |Phi_p> = Ii+|Jj> I,j=1,2,3

        The best scheme:  set way
                          in < p | t | q >, p,t and q is set, and q = p-t
    :return:
    '''

    for nb in nel[1:]:
        for u in uld[nb]:
            at,l = {},list(range(nb))
            for n in nel[1:]:
                for comb in combinations(l, n):
                    k = tuple(u[i] for i in comb)
                    if tuple(b.f for b in k) in sld[n]:
                        v = tuple(b for b in u if b not in k)
                        # at[k] = (sign([b.l for b in k+v if 6<b.f<15]), v)
                        at[k] = v

            mt[tuple(u)] = at
def mt2theory():
    '''
    matrix of T (project into wave) -> tex file

    :return:
    '''

    def mtij2theory(bra, ket):
        raise NotImplementedError('TODO: mtij2theory() in mt2theory()')

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = 'tm','T matrix in GVB-BCCC formula derivation'
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

        for bra,ketl in mt.items():
            write(F'    $\\Phi_{{{bra!s}}} = $ \\\\ \n')
            for t,(s,ket) in ketl.items():
                if s is 1: write(F"    $+ {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{ket}}} $ \\\\ \n")
                else: write(F"    $- {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{ket}}} $ \\\\ \n")
            write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def mt2code(lang='Python'):
    '''
    matrix of T (project into wave) -> code

    :param lang: str, default is Python

    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre = 'tm'
    with open(F'bccc{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')
        
        write('from bccc.pub import sign\n')
        write('\n\n')
        
        bld = {}
        for bra in mt:
            bl = tuple(b.l for b in bra)
            try: bld[bl].append(bra)
            except KeyError: bld[bl] = [bra]
        
        write('def tm(np, nt=4, nc=0):\n')
        write('    tmp = {}\n')
        write('    \n')

        tab = '    '
        for bl in bld:
            nb,tabnb = 0,''
            for b in bl:
                # if nb>0:  write(tabnb+F'    for {b} in range({bl[nb-1]}+1, np):\n')
                # else: write(F'    for {b} in range(nc, np):\n')
                write(tabnb+F'    for {b} in range(nc, np):\n')
                if nb: write(tabnb+F"        if {b} in {{{','.join(bl[:nb])}}}: continue\n")
                nb += 1
                tabnb += tab
            for bra in bld[bl]:
                # kbra = ''.join(F'({b.l},{b.f}),' for b in bra)
                kbra = F"tuple(sorted([{''.join(F'({b.l},{b.f}),' for b in bra)}]))"
                write(tabnb+ '    _t = {}\n')
                for t,ket in mt[bra].items():
                    # kt = ''.join(F'({b.l},{b.f}),' for b in t)
                    kt = F"tuple(sorted([{''.join(F'({b.l},{b.f}),' for b in t)}]))"
                    # kket = ''.join(F'({b.l},{b.f}),' for b in ket)
                    kket = F"tuple(sorted([{''.join(F'({b.l},{b.f}),' for b in ket)}]))"
                    s = F"sign(sorted([{','.join(b.l for b in t if 6<b.f<15)}])+sorted([{','.join(b.l for b in ket if 6<b.f<15)}]))"
                    write(tabnb+F'    _t[{kt}] = ({s}, {kket})\n')
                write(tabnb+F'    tmp[{kbra}] = _t\n')
                write(tabnb+ '    \n')
        
        write( '    return tmp\n')
        write( '    \n')

        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
# CC
nld = {n: tuple(sorted(set([tuple(sorted(tn, reverse=True)) for nb in range(1,n+1)
            for tn in product(nel[1:], repeat=nb) if sum(tn) is n]), reverse=True)) for n in ngl}
gld = {nl: (prodl([qd[ni]*qd[i]**ni for i,ni in zip(sorted(set(nl)), [nl.count(i) for i in sorted(set(nl))])]),
            group(set(range(n)), nl)) for n in ngl for nl in nld[n]}
na = 0
ad = {}; ral = [[] for i in range(ngb+1)]  # rank of array
# na,lock = Manager().Value('l', 0),Manager().Lock()
# ad = Manager().dict(); ral = Manager().list([Manager().list() for i in range(ngb+1)])
def atu2theory(bra, ald):
    '''
    array of T element u -> latex

    Note:     subfuction
    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    fpre,ftitle = F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}",'Hamiltonian matrix in GVB-BCCC formula derivation'
    with open(F'doc{ossep}{fpre}.tex', 'w') as fout:
        write = fout.write
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout

        write(F'% {ftitle}\n')
        write( '\n\n')

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
        write( '    \\linespread{2.0}\selectfont\n')
        write( '    \\thispagestyle{empty}\n')
        write( '    \n\n')
        
        write( '    \\newpage\n')
        write( '    \\setlength{\parindent}{0pt}\n')
        write( '    \\setcounter{page}{1}\n')
        write( '    \n')

        rfeql = list(rfeqd.items())
        for ket,hl in ald.items():
            write(F'    $\\langle{{{bra}}}\\vert \\hat{{H}} \\vert{{{ket}}}\\rangle = $ \\\\ \n')
            for h in hl:
                s = qd[len(bra)]*h[0]
                bl = ''.join(F"\\displaystyle\\sum_{I}" for I in sorted(set(b.l for b in h[1].upl if b.l in b0l)))
                cl = ''.join(F'\\langle{{{c.l[0].l}}}_{{{c.l[2]}}}\\vert {rfeql[c.l[1]][0]!s} \\vert{{{c.l[0].l}}}_{{{c.l[3]}}}\\rangle' if c.f is 'r' else F'{c!s}' for c in h[2:])
                write(F"    $+\\frac{{{s.numerator}}}{{{s.denominator}}}{bl}{h[1]!s}{cl} $ \\\\ \n")
            write( '    \\\\ \n')
            write( '    \n')
        write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def atu2code(bra, ald):
    '''
    array of T element u -> Python code

    Note:     subfuction
    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    bbl,bsl = tuple(b.l for b in bra),tuple(b.f for b in bra); sbbl = [bbl[i] for i,s in enumerate(bsl) if 6<s<15]
    sbra,kbra = ''.join(F'{b}{s}' for b,s in zip(bbl, bsl)),F"({''.join(F'({b},{s}),' for b,s in zip(bbl,bsl))})"
    tab = '    '

    fpre = F'ta_{sbra}'
    with open(F'bccc{ossep}ta{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        # write('import numpy\n')
        write('from bccc.pub import sign,dei\n')
        write('\n\n')

        # h0_
        write(F'def h0_{sbra}(r, h, g, p, nc=0):\n')
        write( '    np = len(p)\n')
        write( '    hd = {}\n')
        write( '    \n')

        ntb,tabntb = 0,''
        for b in bbl:
            # if ntb: write(tabntb+F'    for {b} in range({bbl[ntb-1]}+1,np):\n')
            # else: write(F'    for {b} in range(nc,np):\n')
            write(tabntb+F'    for {b} in range(nc,np):\n')
            if ntb: write(tabntb+F"        if {b} in {{{','.join(bbl[:ntb])}}}: continue\n")
            ntb += 1
            tabntb += tab
        # write(tabntb+F'    u = {kbra}\n')  # A<B
        write(tabntb+F'    u = tuple(sorted({kbra}))\n')  # AB.. 大小不定, 排序固定态, 特别在减少 nu (s 不重复), 这是另一个态, 在 s 重复时, 覆盖旧的, 只保留一份
        if sbbl: write(tabntb+F"    s = sign([{','.join(sbbl)}])\n")  # 排序 AB.. 的符号
        write(tabntb+ '    hu = 0.0\n')
        write(tabntb+ '    \n')

        for h in ald[()]:
            Il = tuple(sorted(set(b.l for b in h[1].upl if b.l in b0l)))
            if Il:
                forl=[]; forlapp=forl.append
                for i,b in enumerate(Il): forlapp(F"for {b} in range(np) if {b} not in [{','.join(bbl+Il[:i])}]")
                write(tabntb+F"    hu += sum({prodl(h)} {' '.join(forl)})\n")
            else: write(tabntb+F'    hu += {prodl(h)}\n')
        write(tabntb+ '    \n')

        if not sbbl: write(tabntb+ '    hd[u] = hu\n')
        else: write(tabntb+ '    hd[u] = s*hu\n')  # 乘上排序符号
        write(tabntb+ '    \n')

        write( '    return hd\n')
        write( '    \n')
        ald.pop(())

        # hu_
        if bra: # bra != ()
            write(F'def hu_{sbra}(r, h, g, p, nc=0):\n')
            write('    np = len(p)\n')
            write('    hd = {}\n')
            write('    \n')
            
            ntb,tabntb = 0,''
            for b in bbl:
                # if ntb: write(tabntb+F'    for {b} in range({bbl[ntb-1]}+1,np):\n')
                # else: write(F'    for {b} in range(nc,np):\n')
                write(tabntb+F'    for {b} in range(nc,np):\n')
                if ntb: write(tabntb+F"        if {b} in {{{','.join(bbl[:ntb])}}}: continue\n")
                ntb += 1
                tabntb += tab
            # write(tabntb+F'    u = {kbra}\n')  # A<B
            write(tabntb+F'    u = tuple(sorted({kbra}))\n')  # AB.. 大小不定, 排序固定态; 注意, 本应有一个符号, 但这里为 hu, H 作用前后 u 相同, 两次 ++/-- 必为 +
            write(tabntb+ '    hu = 0.0\n')
            write(tabntb+ '    \n')

            for h in ald[bra]:
                Il = tuple(sorted(set(b.l for b in h[1].upl if b.l in b0l)))
                if Il:
                    forl=[]; forlapp=forl.append
                    for i,b in enumerate(Il): forlapp(F"for {b} in range(np) if {b} not in [{','.join(bbl+Il[:i])}]")
                    write(tabntb+F"    hu += sum({prodl(h)} {' '.join(forl)})\n")
                else: write(tabntb+F'    hu += {prodl(h)}\n')
            write(tabntb+ '    \n')
            
            tl = []; append = tl.append
            for nl in nld[ntb]:
                for gl in gld[nl][1]:
                    if all(tuple(bsl[i] for i in g) in sld[len(g)] for g in gl):
                        # append((sign([bbl[i] for i in chain.from_iterable(gl) if 6<bsl[i]<15]), ','.join(F"({''.join(F'{bra[i]},' for i in g)})" for g in gl)))
                        append((sign([bbl[i] for i in chain.from_iterable(gl) if 6<bsl[i]<15]), ','.join(F"tuple(sorted([{''.join(F'{bra[i]},' for i in g)}]))" for g in gl)))
            write(tabntb+F"    hd[u] = hu,[{','.join(F'({s},[{t}])' for s,t in tl)}]\n")
            write(tabntb+ '    \n')
            
            write( '    return hd\n')
            write( '    \n')
            ald.pop(bra)
            
        # ta_
        vld = {}  # '', I, IJ, IJK, IJKL
        for v in ald:
            bl = tuple(b.l for b in v if b.l in b0l)
            vld[bl] = vld.get(bl, [])+[v]

        write(F'def {fpre}(t, r, h, g, p, nc=0):\n')  # Amplitude equation
        write( '    td = {}\n')
        write( '    np = len(p)\n')
        write( '    \n')

        ntb,tabntb = 0,''
        for b in bbl:
            # if ntb: write(tabntb+F'    for {b} in range({bbl[ntb-1]}+1,np):\n')
            # else: write(F'    for {b} in range(nc,np):\n')
            write(tabntb+F'    for {b} in range(nc,np):\n')
            if ntb: write(tabntb+F"        if {b} in {{{','.join(bbl[:ntb])}}}: continue\n")
            ntb += 1
            tabntb += tab
        # write(tabntb+F'    u = {kbra}\n')  # A<B
        write(tabntb+F'    u = tuple(sorted({kbra}))\n')  # AB.. 大小不定, 排序固定态, 特别在减少 nu (s 不重复), 这是另一个态, 在 s 重复时, 覆盖旧的, 只保留一份
        if sbbl: write(tabntb+F"    s = sign([{','.join(sbbl)}])\n")  # 排序 AB.. 的符号
        write(tabntb+ '    tu = 0.0\n')
        write(tabntb+ '    tul = []\n')
        write(tabntb+ '    \n')
        
        # na = 0
        for kbl,vl in vld.items():
            ntk,tabntk,bl = ntb,tabntb,bbl+kbl
            for b in kbl:
                write(tabntk+F'    for {b} in range(nc,np):\n')
                write(tabntk+F"        if {b} in {{{','.join(bl[:ntk])}}}: continue\n")
                ntk += 1
                tabntk += tab
            # forl = [F"for {b} in range(nc,np) if {b} not in {{{','.join(bl[:ntb+i])}}}" for i,b in enumerate(kbl)]; ntk=len(bl)

            for v in vl:
                write(tabntk+F'    # {v}\n')
                # write(tabntb+F'    # {v}\n')
                for h in ald[v]:
                    Il = tuple(sorted(set(b.l for b in h[1].upl if b.l not in bl))); bl_=bl+Il
                    if Il:
                        # forl_=[]; forlapp=forl_.append
                        # for i,b in enumerate(Il): forlapp(F"for {b} in range(np) if {b} not in [{','.join(bl+Il[:i])}]")
                        forl_ = [F"for {b} in range(np) if {b} not in {{{','.join(bl_[:ntk+i])}}}" for i,b in enumerate(Il)]
                        write(tabntk+F"    tu += sum({prodl(h)} {' '.join(forl_)})\n")
                        # tu = F"sum({prodl(h)} {' '.join(forl_)})"
                    else:
                        write(tabntk+F'    tu += {prodl(h)}\n')
                        # tu = F'{prodl(h)}'
                    # if forl: write(tabntb+F"    tul.append(sum({tu} {' '.join(forl)}))  # {na}\n")
                    # else: write(tabntb+F"    tul.append({tu})  # {na}\n")
                    # na += 1
            write(tabntk+ '    \n')
            # write(tabntb+ '    \n')
        # write(tabntb+ "    #for i,v in enumerate(tul): print(F'tul[{i}] = {v}')\n")
        # write(tabntb+ '    tu += sum(tul)\n')
        # write(tabntb+ '    \n')

        # write(tabntb+F"    td[u] = td.get(u, 0.0)+{qd[ntb]}*tu*sign([{','.join(b for b,s in zip(bbl, bsl) if 6<s<15)}])\n")
        # write(tabntb+F"    td[u] = td.get(u, 0.0)+tu*sign([{','.join(b for b,s in zip(bbl, bsl) if 6<s<15)}])\n")
        # write(tabntb+F"    td[u] = td.get(u, 0.0)+tu\n")
        if not sbbl: write(tabntb+F"    td[u] = tu\n")  # 原本有个系数，但我们认为 A1B2 != B2A1!, 所以不考虑
        else: write(tabntb+F"    td[u] = s*tu\n")  # 乘上排序符号
        write(tabntb+ '\n')

        write( "    return td\n")
        write( '    \n\n')

        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
# def atu2code(bra, ald):
#     '''
#     array of T element u -> C++ code
#
#     Note:     subfuction
#     :return:
#     '''
#
#     ossep,oslinesep = os.sep,os.linesep # cross platform
#     bbl,bsl = tuple(b.l for b in bra),tuple(b.f for b in bra); sbbl = [bbl[i] for i,s in enumerate(bsl) if 6<s<15]
#     sbra,kbra = ''.join(F'{b}{s}' for b,s in zip(bbl, bsl)), F"{{{','.join(F'{b},{s}' for b,s in zip(bbl, bsl))}}}"
#     tab = '    '
#
#     fpre=F"ta_{sbra}"
#     with open(F'bccc{ossep}ta{ossep}{fpre}.cpp', 'w') as fout:
#         # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
#         write = fout.write
#
#         write( '// Author: Qingchun Wang @ NJU \n')
#         write( '// E-mail: qingchun720@foxmail.com \n')
#         write('\n\n')
#
#         write( '#include %s\n' %('''"../pub.h"'''))
#         write('\n\n')
#
#         # h0_  # Python code
#         '''write(F'def h0_{sbra}(r, h, g, p, nc=0):\n')
#         write( '    np = len(p)\n')
#         write( '    hd = {}\n')
#         write( '    \n')
#
#         ntb,tabntb = 0,''
#         for b in bbl:
#             # if ntb: write(tabntb+F'    for {b} in range({bbl[ntb-1]}+1,np):\n')
#             # else: write(F'    for {b} in range(nc,np):\n')
#             write(tabntb+F'    for {b} in range(nc,np):\n')
#             if ntb: write(tabntb+F"        if {b} in {{{','.join(bbl[:ntb])}}}: continue\n")
#             ntb += 1
#             tabntb += tab
#         # write(tabntb+F'    u = {kbra}\n')  # A<B
#         write(tabntb+F'    u = tuple(sorted({kbra}))\n')  # AB.. 大小不定, 排序固定态, 特别在减少 nu (s 不重复), 这是另一个态, 在 s 重复时, 覆盖旧的, 只保留一份
#         if sbbl: write(tabntb+F"    s = sign([{','.join(sbbl)}])\n")  # 排序 AB.. 的符号
#         write(tabntb+ '    hu = 0.0\n')
#         write(tabntb+ '    \n')
#
#         for h in ald[()]:
#             Il = tuple(sorted(set(b.l for b in h[1].upl if b.l in b0l)))
#             if Il:
#                 forl=[]; forlapp=forl.append
#                 for i,b in enumerate(Il): forlapp(F"for {b} in range(np) if {b} not in [{','.join(bbl+Il[:i])}]")
#                 write(tabntb+F"    hu += sum({prodl(h)} {' '.join(forl)})\n")
#             else: write(tabntb+F'    hu += {prodl(h)}\n')
#         write(tabntb+ '    \n')
#
#         if not sbbl: write(tabntb+ '    hd[u] = hu\n')
#         else: write(tabntb+ '    hd[u] = s*hu\n')  # 乘上排序符号
#         write(tabntb+ '    \n')
#
#         write( '    return hd\n')
#         write( '    \n')'''
#         ald.pop(())
#
#         # hu_  # Python code
#         if bra: # bra != ()
#             '''write(F'def hu_{sbra}(r, h, g, p, nc=0):\n')
#             write('    np = len(p)\n')
#             write('    hd = {}\n')
#             write('    \n')
#
#             ntb,tabntb = 0,''
#             for b in bbl:
#                 # if ntb: write(tabntb+F'    for {b} in range({bbl[ntb-1]}+1,np):\n')
#                 # else: write(F'    for {b} in range(nc,np):\n')
#                 write(tabntb+F'    for {b} in range(nc,np):\n')
#                 if ntb: write(tabntb+F"        if {b} in {{{','.join(bbl[:ntb])}}}: continue\n")
#                 ntb += 1
#                 tabntb += tab
#             # write(tabntb+F'    u = {kbra}\n')  # A<B
#             write(tabntb+F'    u = tuple(sorted({kbra}))\n')  # AB.. 大小不定, 排序固定态; 注意, 本应有一个符号, 但这里为 hu, H 作用前后 u 相同, 两次 ++/-- 必为 +
#             write(tabntb+ '    hu = 0.0\n')
#             write(tabntb+ '    \n')
#
#             for h in ald[bra]:
#                 Il = tuple(sorted(set(b.l for b in h[1].upl if b.l in b0l)))
#                 if Il:
#                     forl=[]; forlapp=forl.append
#                     for i,b in enumerate(Il): forlapp(F"for {b} in range(np) if {b} not in [{','.join(bbl+Il[:i])}]")
#                     write(tabntb+F"    hu += sum({prodl(h)} {' '.join(forl)})\n")
#                 else: write(tabntb+F'    hu += {prodl(h)}\n')
#             write(tabntb+ '    \n')
#
#             tl = []; append = tl.append
#             for nl in nld[ntb]:
#                 for gl in gld[nl][1]:
#                     if all(tuple(bsl[i] for i in g) in sld[len(g)] for g in gl):
#                         # append((sign([bbl[i] for i in chain.from_iterable(gl) if 6<bsl[i]<15]), ','.join(F"({''.join(F'{bra[i]},' for i in g)})" for g in gl)))
#                         append((sign([bbl[i] for i in chain.from_iterable(gl) if 6<bsl[i]<15]), ','.join(F"tuple(sorted([{''.join(F'{bra[i]},' for i in g)}]))" for g in gl)))
#             write(tabntb+F"    hd[u] = hu,[{','.join(F'({s},[{t}])' for s,t in tl)}]\n")
#             write(tabntb+ '    \n')
#
#             write( '    return hd\n')
#             write( '    \n')'''
#             ald.pop(bra)
#
#         # ta_
#         vld = {}  # '', I, IJ, IJK, IJKL
#         for v in ald:
#             bl = tuple(b.l for b in v if b.l in b0l)
#             vld[bl] = vld.get(bl, [])+[v]
#
#         write(F'void {fpre}(real* t, real* r, real* h, real* g, uint* p, uint nc=0)\n')  # Amplitude equation
#         write( '{   \n')
#         write( '    real v = 0.0;\n')
#         write( '    \n')
#
#         ntb,tabntb = 0,''
#         for b in bbl:
#             # if ntb: write(tabntb+F'    for (uint {b}={bbl[ntb-1]}+1; {b}<np; ++{b})\n')
#             # else: write(F'    for (uint {b}=nc; {b}<np; ++{b})\n')
#             write(tabntb+F'    for (uint {b}=nc; {b}<np; ++{b}) {{\n')
#             if ntb: write(tabntb+F"        if ({'||'.join(F'{b}=={ib}' for ib in bbl[:ntb])}) continue;\n")
#             ntb += 1
#             tabntb += tab
#         # write(tabntb+F'    uint u = ud[{kbra}];\n')  # A<B
#         write(tabntb+F"    uint u = ud[{F'phiu({kbra})' if ntb<2 else F'phiu(usort({kbra}))'}];\n")  # AB.. 大小不定, 排序固定态, 特别在减少 nu (s 不重复), 这是另一个态, 在 s 重复时, 覆盖旧的, 只保留一份
#         if sbbl: write(tabntb+F"    int s = sign({{{','.join(sbbl)}}});\n")  # 排序 AB.. 的符号
#         write(tabntb+ '    real tu = 0.0;\n')
#         AB = F"l{''.join(bbl)}"
#         write(tabntb+F"    us {AB} = {{{','.join(bbl)}}};\n")
#         write(tabntb+ '    \n')
#
#         for kbl,vl in vld.items():
#             ntk,tabntk,bl = ntb,tabntb,bbl+kbl
#             sul = [AB]; ABI = AB
#             for b in kbl:
#                 write(tabntk+F'    for (uint {b}=nc; {b}<np; ++{b}) {{\n')
#                 ABi = sul[ntk-ntb]
#                 write(tabntk+F"        if ({ABi}.find({b}) != {ABi}.end()) continue;\n")
#                 ABI = F"l{''.join(bl[:ntk])}{b}"; sul.append(ABI)
#                 write(tabntk+F'        us {ABI} = {ABi}; {ABI}.insert({b});\n')
#                 ntk += 1
#                 tabntk += tab
#
#             for v in vl:
#                 # write(tabntk+F'    // {v}\n')
#                 for h in ald[v]:
#                     h0,h1= float(h[0]),h[1:]
#                     Il = tuple(sorted(set(b.l for b in h[1].upl if b.l not in bl)))  # len(Il) 只可能 <=1, 因为 为 2 时, 必 v = bra 在 hu 中
#                     if Il:
#                         tabntI = tabntk
#                         for i,I in enumerate(Il):
#                             write(tabntI+F'    for (uint {I}=0; {I}<np; ++{I})\n')
#                             if i is 0: write(tabntI+F'        if ({ABI}.find({I})=={ABI}.end())\n')
#                             else: write(tabntI+F"        if ({ABI}.find({I})=={ABI}.end()) && {' '.join(F'({I}!={b})' for b in Il[:i])}\n")
#                             tabntI += tab
#                         write(tabntI+F'        v += {prodl(h1)};\n')
#                         if is1(h0): write(tabntk+ '    tu += v; v = 0.0;\n')
#                         elif is1(h0+2): write(tabntk+ '    tu += -v; v = 0.0;\n')
#                         else: write(tabntk+F'    tu += {h0}*v; v = 0.0;\n')
#                     else:
#                         if is1(h0): write(tabntk+F'    tu += {prodl(h1)};\n')
#                         elif is1(h0+2): write(tabntk+F'    tu += -{prodl(h1)};\n')
#                         else: write(tabntk+F'    tu += {h0}*{prodl(h1)};\n')
#             write(tabntk+F"    {'}'*(ntk-ntb)}\n")
#             write(tabntk+ '    \n')
#
#         if not sbbl: write(tabntb+F"    td[u] = tu;\n")  # 原本有个系数，但我们认为 A1B2 != B2A1!, 所以不考虑
#         else: write(tabntb+F"    td[u] = s*tu;\n")  # 乘上排序符号
#         write(tabntb+F"{'}'*ntb}\n")
#         write(tabntb+ '\n')
#
#         # write( "    return td;\n")
#         write( '}\n')
#         write( '    \n\n')
#
#         # write('\n\n\n\n')
#         # write("if __name__ == '__main__':\n")
#         # write('    \n\n')
#         # write("    print('End successfully')\n")
#         # write('    \n\n')
#
#         # # sys.stdout = stdout; sys.stderr = stderr
#         # fout.close()
def h2u(u):
    '''
    T project into wave

    :param u: wave state, list

    :return:
    '''

    hld = {}; vd = {} # {v: [[c, <>, r..]]}
    nb = len(u); ne = sum(b.ne for b in u)
    bl = tuple(b.l for b in u)

    # for h
    if nb<=2 and ne>=1:  # nb<2 for h; p+p- ne>=1
        for prodb in product(u, repeat=2):
            if len(set(prodb)) is not nb: continue
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    inte = Integral(symbl=(bop, boq), mode=0)
                    functor = Functor(symbl=(bop, boq), pdl='+-', upl=prodb)
                    s_blockize,functorl,upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                    _hld={}; _vd={} # ce {v: r} 没有与 inte 相乘
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
                                v = deepcopy(upl)
                                for i,b in zip(bsl,v): b.f = i
                                if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                _rl=[]; append=_rl.append
                                for b,r in zip(v, rl):
                                    _r = deepcopy(r)
                                    _r.l = (b,)+_r.l
                                    append(_r)
                                try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                except KeyError: _vd[v],_hld[v] = v,[_rl]
                    for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                        _inte = deepcopy(inte)
                        _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                        # _inte.orbsort()
                        si = [Frac(s_blockize,1), _inte]
                        try:
                            for rl in _hl: hld[v].append(si+rl)
                        except KeyError:
                            for i,rl in enumerate(_hl): _hl[i] = si+rl
                            vd[v],hld[v] = v,_hl

    # for g
    if ne>=2: # p+p+s-r- ne>=2
        for prodb in product(u, repeat=4):
            if len(set(prodb)) is not nb: continue
            _bl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    for bor in deepcopy(bo_l):
                        for bos in deepcopy(bo_l):
                            inte = Integral(symbl=(bop, boq, bor, bos), mode=2)
                            functor = Functor((bop, boq, bos, bor), pdl='++--', upl=(prodb[0], prodb[1], prodb[3], prodb[2]))
                            s_blockize, functorl, upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                            _hld = {}; _vd={}  # ce 没有与 inte 相乘
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
                                        v = deepcopy(upl)
                                        for i,b in zip(bsl, v): b.f = i
                                        if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                        elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                        _rl=[]; append=_rl.append
                                        for b,r in zip(v, rl):
                                            _r = deepcopy(r)
                                            _r.l = (b,)+_r.l
                                            append(_r)
                                        try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError: _vd[v],_hld[v] = v,[_rl]
                            for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                                # _inte.orbsort()
                                si = [Frac(s_blockize, 2), _inte]
                                try:
                                    for rl in _hl: hld[v].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(_hl): _hl[i] = si+rl
                                    vd[v],hld[v] = v,_hl

                            _hld = {}; _vd = {}  # ce 没有与 inte 相乘
                            a,b = spinl
                            bop.f = bor.f = a; boq.f = bos.f = b
                            bfeql = [None]*nb
                            for i,f in enumerate(functorl):
                                bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                else: break
                            else:
                                for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                    bsl,rl = zip(*prodcxb)
                                    v = deepcopy(upl)
                                    for i,b in zip(bsl, v): b.f = i
                                    if v in hld: v = vd[v] # 取 v of hld, 如果key已经有了
                                    elif v in _hld: v = _vd[v] # 取 v of _hld, 如果key已经有了
                                    _rl = []; append=_rl.append
                                    for b,r in zip(v, rl):
                                        _r = deepcopy(r)
                                        _r.l = (b,)+_r.l
                                        append(_r)
                                    try: _hld[v].append(_rl)  # spin 内合并 或 alpha与beta 的合并
                                    except KeyError: _vd[v], _hld[v] = v, [_rl]
                            for v,_hl in _hld.items():  # 将 _hld 更新到 hld 中
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(v[iis(b, bl)] for b in _bl)
                                # _inte.orbsort()
                                si = [Frac(s_blockize, 1), _inte]
                                try:
                                    for rl in _hl: hld[v].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(_hl): _hl[i] = si+rl
                                    vd[v],hld[v] = v,_hl

    # for v,hl in hld.items():
    #     hld[v] = sum(prodl(h) for h in hl)
    
    return hld
def detat(bra, hfeqd):
    '''
    determine matrix of hamiltionian (project into wave)

    Note:       Hamiltonian matrix is complex conjugate. This point is often used in the process of deriving formulas
                <W| H  =>  H |W>
                when it stored and output, <bra| H |ket>

    :return:
    '''

    print('bra = ', bra)
    # cbra = qd[len(bra)]  # cbra (for AB) 本应被乘, 但为方便处理, 令 AaBb!=BbAa, 所以不乘
    global na
    
    # mH
    hld = {}  # {v: [[c, <>, r..]...]}
    for nbc in range(1, ngb+1):  # n blocks correlation
        for nbe in range(nbc+1):  # nbe blocks from bra
            nb0=nbc-nbe; cket=qd[nb0]
            for bchose in combinations(bra, r=nbe):  # 从激发态取 nbe
                brest = tuple(sorted(set(bra)-set(bchose)))  # set -> list is random and disordered
                s_bef = sign([b.l for b in bchose+brest if b.ne%2])  # H 作用前, 激发块提到前面 (与 dethm 不同, 符号会体现在 t 中, 所以这里无须体现在表达式上)
                _u = bchose+gbl[:nb0]
                # if len(_u+brest)>neb: continue
                
                sl = tuple(b.f for b in _u); bl = tuple(b.l for b in _u)
                if sl in hfeqd:
                    _hld = deepcopy(hfeqd[sl])
                    try:
                        _bl = tuple(b.l for b in list(_hld.keys())[0])
                        if _bl != bl:
                            for v in _hld:
                                for b,i in zip(v, bl): b.l = i
                    except IndexError: continue   # _hld = {} when nbc=nbe, b.f=6
                else:
                    _hld = h2u(_u)  # Aa I0 J0 Kk
                    hfeqd[sl] = _hld
                    _hld = deepcopy(_hld)
                    
                for v,hl in _hld.items():  # 将 _hld 更新(多块加和)到 hld 中
                    nI1 = len([b for b in v[nbe:] if b.f])
                    if nb0>1 and nI1 is not nb0:  # I0Ja -> J0Ia
                        i1 = 0; i0 = nI1
                        for b in v[nbe:]:
                            if b.f:
                                b.l = bgl[i1]
                                i1 += 1
                            else:
                                b.l = bgl[i0]
                                i0 += 1

                    _u = v+deepcopy(brest); u = tuple(sorted(b for b in _u if b.f))
                    if not len(u) in ngl: continue  # for no T1
                    s_lat = sign([b.l for b in _u if b.ne%2])  # H 作用后, 激发块回到原来位置 (与 dethm 不同, 符号会体现在 t 中, 所以这里无须体现在表达式上)
                    for h in hl:
                        # h[0] = s_bef*s_lat*h[0]
                        h[0] = s_bef*s_lat*cket*h[0]
                        # h[1].orbsort()
                    hld[u] = hld.get(u, [])+hl
                    
    # at
    ald = {(): hld.get((),[]), bra: hld[bra]}  # {v: [[c, <>, r.., t..]...]}
    hld.pop((), []); hld.pop(bra,[])
    for v,hl in hld.items():
        # print(F'bra->v = {bra}->{v}')
        hd = {}  # 有效的 hl (合并 h)
        for c,inte,*rl in hl:
            inte.orbsort()  # reorder h
            rd={r.l[0].l: r for r in rl}; rl=[rd[k] for k in sorted(rd)]  # reorder rl
            kh = F'{prodl([*rl, inte])}'
            if kh in hd: hd[kh][0] += c
            else: hd[kh] = [c, inte, *rl]
        hl = list(h for h in hd.values() if not is0(h[0]))
        bsl = tuple(b.f for b in v)
        nb,gll = len(bsl),[]  # 有效的 gl (排除不合理 gl)
        for nl in nld[nb]:
            # c,gll = gld[nl]  # 已保证只算一次，所以无须乘c
            for gl in gld[nl][1]:
                if all(tuple(bsl[i] for i in g) in sld[len(g)] for g in gl):
                    gll.append(gl)
            
        for gl,h in product(gll, hl):
            c,inte,*rl = h
            v0=v; rbl=[r.l[0] for r in rl if r.l[2]]  # old v, r/inte 未必依赖于 v
            v = tuple(sorted([b for b in v0 if b not in rbl]+rbl))  # new v, rl/inte 依赖于 v, v 中 b 改变必引起 r/inte 改变
            gl,_ = zip(*sorted(zip(gl, [sum(4 if v[i].l in b0l else -1 for i in g) for g in gl]), key=lambda x: x[1]))  # reorder gl
            ul,ng = [tuple(v[i] for i in g) for g in gl],len(gl)
            s,tl = sign([i for i in chain.from_iterable(gl) if 6<bsl[i]<15]),[Coeff(symbl=ul[i], mode='t') for i in range(ng)]
            ald[v] = ald.get(v, [])+[[c*s, inte, *rl, *tl]]
            # a = c*s*prodl(rl+tl)*inte  # non-array mode
            # print(F'a(!a) = {a}')
            
            # bdl=[{b.l: [b] for b in ul[i]} for i in range(ng)]; tdl=[{''.join(bdl[i].keys()):[tl[i]]} for i in range(ng)]
            # bll=[''.join(bd.keys()) for bd in bdl]; bl=''.join(bll)
            # iI = min(iin(b, bll) for b in b0l if b in bl) if b0l&set(bl) else ng
            # rtl=[(bdl[i], {}, tdl[i]) for i in range(iI,ng)]  # [(bd, rd, td)...]; 在 bd 在 g 内保证唯一；g 外可能有重复 AB AI BIJ
            # bd0,rd0,td0 = {b.l: [] for b in bra},{},{}
            # for bi,bl in chain.from_iterable(bd.items() for bd in bdl[:iI]): bd0[bi] = bl
            # for bi,tl in chain.from_iterable(td.items() for td in tdl[:iI]): td0[bi] = tl
            # rI0l = []
            # for r in rl:  # rl
            #     b=r.l[0]; bi=b.l
            #     if bi in b0l:
            #         if b.f: rtl[iin(iis(b,v),gl)-iI][1][bi] = [r]
            #         else: rI0l.append((bi,b,r))
            #     else: rd0[bi] = [r]
            # else:
            #     for bi,b,r in rI0l: rtl.append(({bi: [b]}, {bi: [r]}, {}))
            # rtl.insert(0, (bd0, rd0, td0))
            # # print(F"bd,rd,td = {'  '.join(F'{list(chain(chain(*bd.values()),chain(*rd.values()),chain(*td.values())))}' for bd,rd,td in rtl)}")
            # # print(F"ul = {'-'.join(''.join(F'{b}' for b in chain.from_iterable(bd.values())) for bd,*_ in rtl)}")
            # # print(F'a(sort) = {c},{s},{prodl(chain.from_iterable(chain(chain(*td.values()),chain(*rd.values())) for _,td,rd in rtl))},{inte}')
            #
            # pgl = [(s, rtl, inte)]  # pgl(prod g list) = [(sign, [(bd, rd, td)...], inte/array)...]; IJKL 必会出现在 inte/array 中
            # bll,ng = [''.join(tr[0].keys()) for tr in rtl],len(rtl)
            # # print(F'bll = {bll}')
            # cbll=[''.join(b for b in bl if b in b0l) for bl in bll]; bl=''.join(list(bd0.keys())+cbll)[::-1]  # 收缩指标
            # # print(F'cbll,bl = {cbll},   {bl}')
            # iI = 1 if cbll[1:] else ng
            # # print(F'iI,ng = {iI},{ng}')
            # for ing in range(ng-1,iI-1,-1):  # i=0: prod(thl[i:}); i=ng: sum(L0)
            #     cbl = cbll[ing]; ncb=len(cbl)
            #     # print(F'cbl = {cbl}')
            #
            #     # K -> J, I, B, A
            #     _pglc,_pgl = [],[]  # _pglc 需 build array 的 _pg, 无需 build 直接存于 _pgl; _pgl 本代 ing 结果,
            #     _pgl1,_pgl0 = pgl,[]  # incb,incb-1 次结果; incb=ncb-1, 开始时 _pgl1 = pgl(pgl 上一代 ing+1 结果), 结束时 _pgl1,_pgl0 = _pgl0,[]
            #     for incb in range(ncb-1, 0, -1):  # 0<incb, rtl[-1] 中必有 t(若无t, 只能K0L0, 但这会被分解 [K0][L0])
            #         s,cbi = -1**(incb+1),cbl[incb]
            #         # print(F' cbi={cbi}')
            #         for _pg in _pgl1:
            #             _pglc.append(_pg)
            #             for bi in bl[iis(cbi, bl)+1:]:
            #                 _bd,_rd,_td = _pg[1][-1]
            #                 if not any(cbi in _kt and bi in _kt for _kt in _td.keys()):  # cbi,bi 在同一个 t 内(有可能一个 g 内含多个 t), 略去
            #                     _s,rtl,inte = deepcopy(_pg)
            #
            #                     _bd,_rd,_td = rtl[-1]
            #                     for b in _bd[cbi]: b.l = bi  # 修改 cbi -> bi
            #
            #                     # rd,td 原始源于不同 b, cbi->bi 不易造成重复, 但 rd,td 的移动使 bd 易重复, 为保证唯一, 修改 cbi 后, 由 rd,td 重新构造 bd, 丢弃原来的 bd
            #                     # 有 t, 必不在同一个 g, 此时只有 _rd[cbi] 被移动, _td[cbi] 被保留为 _td[bi]
            #                     bd,rd,_=rtl[iin(bi, bll)]
            #                     il=[id(b) for b in bd[bi]]
            #                     for r in _rd[cbi]:  # bd 添加 bi from _rd[cbi]
            #                         b = r.l[0]
            #                         if id(b) not in il: bd[bi].append(b)
            #                     rd[bi] = rd.get(bi, [])+_rd.pop(cbi)  # 移动 _rd[cbi]
            #                     for _kt in list(_td.keys()):
            #                         if cbi in _kt:
            #                             kt=_kt.replace(cbi, bi)
            #                             for tb in chain.from_iterable(t.l for t in _td[_kt]):  # _bd 添加 bi from _td[cbi]
            #                                 if tb.l is bi: _bd[bi] = _bd.get(bi, [])+[tb]
            #                             _td[kt] = _td.get(kt, [])+_td.pop(_kt)  # 移动 _td[cbi]
            #                     _bd.pop(cbi)  # del _bd[cbi]
            #
            #                     _pgl0.append((s*_s, rtl, inte))
            #         _pgl1=_pgl0; _pgl0=[]
            #     else:  # incb = 0
            #         s,cbi = -1,cbl[0]
            #         # print(F' cbi={cbi}')
            #         for _pg in _pgl1:
            #             _pglc.append(_pg)
            #             for bi in bl[iis(cbi, bl)+1:]:
            #                 _bd,_rd,_td = _pg[1][-1]
            #                 if not any(cbi in _kt and bi in _kt for _kt in _td.keys()):  # cbi,bi 在同一个 t 内(有可能一个 g 内含多个 t)，略去
            #                     _s,rtl,inte = deepcopy(_pg)
            #
            #                     _bd,_rd,_td = rtl[-1]
            #                     for b in _bd[cbi]: b.l = bi  # 修改 cbi -> bi
            #
            #                     # rd,td 原始源于不同 b, cbi->bi 不易造成重复, 但 bd 易重复, 为保证唯一, 修改 cbi 后, 由 rd,td 重新构造 bd, 丢弃原来的 bd
            #                     # cbi = cbl[0] 必不在同一个 g; 移动 r,t (r->trl[bi], t-> trl[max(kt)]), 而不 del rtl[-1], 后面取 rtl[:-1]
            #                     bd,rd,_ = rtl[iin(bi, bll)]
            #                     il = [id(b) for b in bd[bi]]
            #                     for r in _rd[cbi]:  # bd 添加 bi from _rd[cbi]
            #                         b = r.l[0]
            #                         if id(b) not in il: bd[bi].append(b)
            #                     rd[bi] = rd.get(bi, [])+_rd[cbi]  # 添加 _rd[cbi]
            #                     for _kt in list(_td.keys()):
            #                         kt = _kt.replace(cbi, bi)
            #                         bdt,_,tdt = rtl[max(iin(b, bll) for b in kt)]
            #                         idd = {ib:[id(b) for b in bdt.get(ib,[])] for ib in kt}
            #                         for tb in chain.from_iterable(t.l for t in _td[_kt]):  # bdt 添加 bi from _td[cbi]
            #                             itb = tb.l
            #                             if id(tb) not in idd[itb]: bdt[itb] = bdt.get(itb, [])+[tb]
            #                         tdt[kt] = tdt.get(kt,[])+_td[_kt]  # 添加 _td[cbi]
            #
            #                     _pgl.append((s*_s, rtl[:-1], inte))
            #     # print(' _pglc = ')
            #     # for _pg in _pglc: print(F'     {_pg[0]} ...{_pg[1][-1]}  {_pg[2]}')
            #     # print(' _pgl = ')
            #     # for _pg in _pgl: print(F'   {_pg[0]} ...{_pg[1][-1]}  {_pg[2]}')
            #
            #     # Array
            #     al = []
            #     for s,rtl,_inte in _pglc:
            #         _bd,_rd,_td = rtl[-1]
            #         for b in (_inte.upl if isinstance(_inte, Integral) else _inte.l):
            #             ib = b.l
            #             if id(b) not in [id(b) for b in _bd.get(ib, [])]:
            #                 _bd[ib] = _bd.get(ib, [])+[b]
            #         _cbil=''.join(sorted(b for b in cbl if b in _bd)); _abil=''.join(sorted(b for b in _bd.keys() if b not in _cbil))
            #         _cbd={b:_bd[b] for b in _cbil}; _abd={b:_bd.get(b, []) for b in _abil}
            #         # print(F" _abil,_cbil = {_abil},{_cbil}")
            #         # print(F' _abd,_cbd,_rd,_td = {_abd}  {_cbd}  {_rd}  {_td}')
            #
            #         # 欲保证 array 不重复计算, 首先保证 block,orbit(in inte) 顺序, 其次保证 rl,tl(包括block in tu, 这还会涉及到 +/-) 顺序,
            #         abd,cbd,rd,td,inte = deepcopy([_abd, _cbd, _rd, _td, _inte])
            #         i1 = i0 = 0; abd_,cbd_,rd_,td_ = {},{},{},{}
            #         for _ib in _abil+_cbil:  # reorder block (BKJ->AIJ)
            #             if _ib in b1l:
            #                 ib_ = bel[i1]
            #                 i1 += 1
            #             else:
            #                 ib_ = bgl[i0]
            #                 i0 += 1
            #             if _ib in abd:
            #                 for b in abd[_ib]: b.l = ib_
            #                 abd_[ib_] = abd.pop(_ib)
            #             else:
            #                 for b in cbd[_ib]: b.l = ib_
            #                 cbd_[ib_] = cbd.pop(_ib)
            #             if _ib in rd: rd_[ib_] = rd.pop(_ib)
            #         ktl = list(kt for kt in td if _ib in kt)
            #         for kt in ktl:
            #             tl = td.pop(kt)
            #             td_[''.join(b.l for b in tl[0].l)] = tl
            #         if isinstance(inte, Integral): inte.orbsort()  # reorder orbital in inte (block order in array finished above)
            #         abd,cbd,rd,td = abd_,cbd_,rd_,td_
            #         # print(F" abil,cbil =  {''.join(abd.keys())},{''.join(cbd.keys())}")
            #         # print(F' abd,cbd,rd,td = {abd}  {cbd}  {rd}  {td}')
            #
            #         cbil=sorted(cbd.keys(), reverse=True); abl,_abl=tuple(bl[0] for _,bl in abd.items()),tuple(bl[0] for _,bl in _abd.items())
            #         sl = []  # sign list for reorder block in t
            #         va = inte
            #         for cbi in cbil:  # sigma求和只对冗余指标, 中间数组下标并不需要
            #             cb = cbd[cbi][0]
            #             td_,ktl = {},list(_kt for _kt in td if cbi in _kt)  # 合并 tl (AI, IA)
            #             for kt in ktl:
            #                 kt_ = ''.join(sorted(kt))
            #                 td_[kt_] = td_.get(kt_, [])+td.pop(kt)
            #             kt_l = sorted(td_)  # reorder tl
            #             for kt_ in kt_l:
            #                 tl = td_[kt_]
            #                 for t in tl:
            #                     idd = {id(b):i for i,b in enumerate(t.l)}; t.l = sorted(t.l)  # reorder block in t
            #                     sl.append(sign([idd[id(b)] for b in t.l if b.ne%2]))
            #                 va = prodl([*tl,va])
            #             _rl = [(r, F'{r}') for r in rd.pop(cbi, [])]
            #             if _rl: rl,_ = zip(*sorted(_rl, key=lambda x: x[1]))  # reorder rl; r(的 key) 可能重复, 无法用字典实现对 rl 排序(无法以字典同时排序两个 list 时, 这是一种被推荐的方式)
            #             else: rl = []
            #             va = sigma(cb, 'nc' if cb.f else 0, prodl([*rl,va]))
            #         if sl: s *= prodl(sl)
            #
            #         ka = F'{va}'
            #         # print(F' ka = {ka}')
            #         kv=ad.get(ka, None)
            #         if kv is None:
            #             # lock.acquire(); na.value += 1; lock.release()
            #             na += 1
            #             kv = [Array(symbl=abl, mode=F'a{na}'), va]
            #             ad[ka] = kv
            #             # print(F' non-exsit {kv[0]!s} = {kv[1]}')
            #         # else:
            #         #     print(F' exist {kv[0]!s} = {kv[1]}')
            #         #     print(F' na = {na}')
            #         a = deepcopy(kv[0])
            #         a.l = _abl
            #
            #         al.append((s, rtl[:-1], a))
            #
            #     pgl = al+_pgl
            # else: # end ing
            #     for pg in pgl:
            #         s,rtl,a = pg
            #         a = prodl([c*s, *chain.from_iterable(chain(chain(*rd.values()),chain(*td.values())) for _,rd,td in rtl[:iI]), a])
            #         ald[v] = ald.get(v, [])+[a]
            #         # print(F'a = {a}\n')
            
    # atu2theory(bra, ald)
    atu2code(bra, ald)
def at2theory():
    '''
    array of T element u -> latex

    main function
    :return:
    '''

    raise NotImplementedError('TODO: at2theory()')
def at2code():
    '''
    array of T element u -> code

    main function
    :return:
    '''

    ossep,oslinesep = os.sep,os.linesep # cross platform
    kld={n: tuple(tuple((b.l,b.f) for b in u) for u in ul) for n,ul in uld.items()}; sld={n: tuple(''.join(F'{b}{s}' for b,s in u) for u in ul) for n,ul in kld.items()}
    f0ld = {n: tuple(F'h0_{s}' for s in sl) for n,sl in sld}; ftld = {n: tuple(F'ta_{s}' for s in sl) for n,sl in sld}
    fuld = {n: tuple(F'hu_{s}' for s in sl) for n,sl in sld}; fuld.pop(0)

    fpre = 'hm'
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        write( 'from multiprocessing import Pool, Manager\n')
        write( 'from importlib import import_module\n')
        write( 'import sys, gc\n')
        write( '\n\n')

        write( 'def cal(braf, r, h, g, p, nc=0):\n')
        write( "    module = F'bccc.hm.{braf}'\n")
        write( '    lib = import_module(module)\n')
        write( '    fun = getattr(lib, braf)\n')
        write( '    hd = fun(r, h, g, p, nc)\n')
        write( '    del sys.modules[module]\n')
        write( '    gc.collect()\n')
        write( '    return hd\n')
        write( '\n\n')

        write( 'def ta(r, h, g, p, nt=4):\n')
        write( '    hdd = {}\n')
        write( '    \n')

        def outfl(fl):
            write( '    # bra wavefunction list\n')
            write( '    # ground state\n')
            # write(F'    from bccc.ta.ta_ import ta_\n')
            # write(F'    from bccc.ta.ta_ import h0_\n')
            write( "    f0l = ('h0_',)\n")
            write( "    ftl = ('ta_',)\n")
            write( '    \n')

            write( '    # 1-block excited types\n')
            f0l = f0ld[1]; f0ll = [f0l[10*i:10*(i+1)] for i in range(0, math.ceil(len(f0l)/10))]
            # for f0 in f0l: write(F'    from bccc.ta.{f0} import {f0}\n')
            write( '    f0l += (\n')
            # for f0l in f0ll: write(F"            {','.join(f0l)},\n")
            for f0l in f0ll: write(F'''            {''.join(F"'{f0}'," for f0 in f0l)}\n''')
            write( '            )\n')
            write( '    \n')
            ful = fuld[1]; full = [ful[10*i:10*(i+1)] for i in range(0, math.ceil(len(ful)/10))]
            # for fu in ful: write(F'    from bccc.ta.{fu} import {fu}\n')
            write( '    ful += (\n')
            # for ful in full: write(F"            {','.join(ful)},\n")
            for ful in full: write(F'''            {''.join(F"'{fu}'," for fu in ful)}\n''')
            write( '            )\n')
            write( '    \n')
            ftl = ftld[1]; ftll = [ftl[10*i:10*(i+1)] for i in range(0, math.ceil(len(ftl)/10))]
            # for ft in ftl: write(F'    from bccc.ta.{ft} import {ft}\n')
            write( '    ftl += (\n')
            # for ftl in full: write(F"            {','.join(ftl)},\n")
            for ftl in full: write(F'''            {''.join(F"'{ft}'," for ft in ftl)}\n''')
            write( '            )\n')
            write( '    \n')

            for nt in range(2,neb+1):
                write(F'    # {nt}-block excited types\n')
                write(F'    if nt>{nt-1}: \n')
                f0l = f0ld[nt]; f0ll = [f0l[10*i:10*(i+1)] for i in range(0, math.ceil(len(f0l)/10))]
                # for f0 in f0l: write(F'    from bccc.ta.{f0} import {f0}\n')
                write( '        f0l += (\n')
                # for f0l in f0ll: write(F"                {','.join(f0l)},\n")
                for f0l in f0ll: write(F'''                {''.join(F"'{f0}'," for f0 in f0l)}\n''')
                write( '                )\n')
                write( '        \n')
                ful = fuld[nt]; full = [ful[10*i:10*(i+1)] for i in range(0, math.ceil(len(ful)/10))]
                # for fu in ful: write(F'    from bccc.ta.{fu} import {fu}\n')
                write( '        ful += (\n')
                # for ful in full: write(F"                {','.join(ful)},\n")
                for ful in full: write(F'''                {''.join(F"'{fu}'," for fu in ful)}\n''')
                write( '                )\n')
                write( '        \n')
                ftl = ftld[nt]; ftll = [ftl[10*i:10*(i+1)] for i in range(0, math.ceil(len(ftl)/10))]
                # for ft in ftl: write(F'    from bccc.ta.{ft} import {ft}\n')
                write( '        ftl += (\n')
                # for ftl in full: write(F"                {','.join(ftl)},\n")
                for ftl in full: write(F'''                {''.join(F"'{ft}'," for ft in ftl)}\n''')
                write( '                )\n')
                write( '        \n')

        write( '    pool = Pool()\n')
        write( '    pl = tuple(pool.apply_async(calculate, (braf, r, h, g, p, nc)) for braf in ftl)\n')
        write( '    pool.close()\n')
        write( '    pool.join()\n')
        write( '    # for p in pl: hdd.update(p.get())\n')
        write( '    for p in pl:\n')
        write( '        for u,hu in p.get().items():\n')
        write( '            try:\n')
        write( '                for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        write( '            except KeyError: hdd[u] = hu\n')
        write( '    # # for braf in ftl: hdd.update(cal(braf, r, h, g, p))\n')
        write( '    # for braf in ftl:\n')
        write( '    #     for u,hu in calculate(braf, r, h, g, p).items():\n')
        write( '    #         try:\n')
        write( '    #            for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        write( '    #         except KeyError: hdd[u] = hu\n')
        write( '    \n')
        write( '    return hdd\n')
        write( '    \n\n')

    fpre = 'ta'
    with open(F'bccc{ossep}ta{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')

        write('# from multiprocessing import Pool, Manager\n')
        write('from copy import deepcopy\n')
        write('from datetime import datetime\n')
        write('from bccc.pub import sign\n')
        write('from itertools import permutations\n')
        write('\n\n')

        write('def cal(braf, r, h, g, p):\n')
        write("    module = F'bccc.hm.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    hd = fun(r, h, g, p)\n')
        write('    del sys.modules[module]\n')
        write('    gc.collect()\n')
        write('    return hd\n')
        write('\n\n')

        # kld={n: tuple(tuple((b.l,b.f) for b in u) for u in ul) for n,ul in uld.items()}; sld={n: tuple(''.join(F'{b}{s}' for b,s in u) for u in ul) for n,ul in kld.items()}

        write('def ta(t0, h, p, nt=4, nc=0):\n')
        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        write('    from bccc.ta.ta_ import ta_\n')
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[1])
        for braf in brafl: write(F'    from bccc.ta.{braf} import {braf}\n')
        write('    brafl = (\n             ')
        for i,braf in enumerate(brafl):
            write(F'{braf}, ')
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for nt in range(2,neb+1):
            write(F'    # {nt}-block excited types\n')
            write(F'    if nt>{nt-1}: \n')
            brafl = tuple(F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[nt])
            for braf in brafl: write(F'        from bccc.ta.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F'{braf}, ')
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write('    Eref = h[()][()]\n')
        write('    \n')

        write('    t0 = {v: sign([b for b,s in v if 6<s<15])*tu for u,tu in deepcopy(t0).items() for v in permutations(u)}\n')
        write('    Ecorr=ta_(t0, h, p, nc); E0 = Eref+Ecorr\n')
        write("    print(F'Ecorr,E0 = {Ecorr},{E0}')\n")
        write('    \n')

        write("    print('    it            dE               dt               time/s')\n")
        write('    for it in range(200):\n')
        write('        t = {}\n')
        write('        \n')

        write('        t1 = datetime.now()\n')
        write('        # pool = Pool()\n')
        write('        # vl = tuple(pool.apply_async(braf, (t0, h, p, E0, nc)) for braf in brafl)\n')
        write('        # pool.close()\n')
        write('        # pool.join()\n')
        write('        # for v in vl: t.update(v.get())\n')
        write('        for braf in brafl: t.update(braf(t0, h, p, E0, nc))\n')
        write('        t2 = datetime.now()\n')
        write('        \n')

        write('        dt = sum(abs(t[u]-t0[u]) for u in t)\n')
        write('        \n')

        write('        t0 = {v: sign([b for b,s in v if 6<s<15])*tu for u,tu in t.items() for v in permutations(u)}\n')
        write('        Ecorr=ta_(t0, h, p, nc); E = Eref+Ecorr\n')
        write('        \n')

        write('        dE = E-E0\n')
        write('        \n')

        write("        print(F'    {it+1:3d}      {dE:13.10f}    {dt:13.10f}     {t2-t1}', flush=True)\n")
        write('        if abs(dt)<1.0e-7: \n')
        write("            print(F'Successfully converged: Ecorr = {Ecorr}')\n")
        write('            break\n')
        write('        else: \n')
        write('            E0 = E\n')
        write('        \n')
        write('    else: \n')
        write("        print('Error: Convergence failed ')\n")
        write('    \n')

        write('    return Ecorr,t\n')
        write('    \n\n')


        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()

# initial guess
def init2theory():
    raise NotImplementedError('TODO: init2theory()')
def init2code():
    raise NotImplementedError('TODO: init2code()')
# iterative equation
def iter2theory():
    raise NotImplementedError('TODO: iter2theory()')
def iter2code():
    ossep,oslinesep = os.sep,os.linesep # cross platform
    sld = {n: tuple(''.join(F'{b.l}{b.f}' for b in u) for u in ul) for n,ul in uld.items()}
    f0ld = {n: tuple(F'h0_{s}' for s in sl) for n,sl in sld.items()}; ftld = {n: tuple(F'ta_{s}' for s in sl) for n,sl in sld.items()}
    fuld = {n: tuple(F'hu_{s}' for s in sl) for n,sl in sld.items()}; fuld.pop(0)

    # fpre = 'ta'
    # with open(F'bccc{ossep}ta{ossep}{fpre}.h', 'w') as fout:
    #     # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
    #     write = fout.write
    #
    #     write( '// Author: Qingchun Wang @ NJU\n')
    #     write( '// E-mail: qingchun720@foxmail.com \n')
    #     write( '\n\n')
    #     write( '#ifndef TA_H\n')
    #     write( '#define TA_H\n')
    #     write( '\n\n')
    #     write( '#include %s\n' %('''"../pub.h"'''))
    #     write( '\n\n')
    #     write( '%s\n' %'''extern "C" {''')
    #     write( '\n')
    #     write( 'real* ta(real* t, uint* pul, uint* ul, uint nu, real* r, uint np, real* h, uint no, real* g, uint* lT0, uint n2, uint* p, uint nc);\n')
    #     write( '\n')
    #     write( '}\n')
    #     write( '\n\n')
    #     for ft in chain.from_iterable(ftld.values()):
    #         write(F'void {ft}(real* t, real* r, real* h, real* g, uint* p, uint nc=0);\n')
    #     write( '\n\n')
    #
    #     write( '\n\n\n\n')
    #     write( '#endif\n')
    #     write( '\n\n')
    # with open(F'bccc{ossep}ta{ossep}{fpre}.cpp', 'w') as fout:
    #     # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
    #     write = fout.write
    #
    #     write( '// Author: Qingchun Wang @ NJU\n')
    #     write( '// E-mail: qingchun720@foxmail.com \n')
    #     write( '\n\n')
    #     write( '#include %s\n' %('''"ta.h"'''))
    #     write( '#include %s\n' %('''"omp.h"'''))
    #     write( '\n\n')
    #     write( 'real* ta(real* t, uint* pul, uint* ul, uint nu, real* r, uint np, real* h, uint no, real* g, uint* lT0, uint n2, uint* p, uint nc)\n')
    #     write( '{\n')
    #     write( '     //test();\n')
    #     write( '     ::no = no; ::lT0 = lT0; ::n2 = n2;\n')
    #     write( '     ::np = np; ::nc = nc;\n')
    #     write( '     ::pul=pul; ud = l2d(ul, nu); ::nu = nu;\n')
    #     write( '     td = new real[nu]();\n')
    #     write( '     \n')
    #     write( '     #pragma omp parallel sections\n')
    #     write( '     {\n')
    #     for ft in chain.from_iterable(ftld.values()):
    #         write( '          #pragma omp section\n')
    #         write(F'          {ft}(t, r, h, g, p, nc);\n')
    #     write( '     }\n')
    #     write( '     \n')
    #     write( '     return td;\n')
    #     write( '}\n')
    #     write( '\n\n')
    #
    #     write( '\n\n\n\n')
    #     write( 'int main(int narg, char* argl[])\n')
    #     write( '{\n')
    #     write( '     /* test input */\n')
    #     write( '     %s\n' %'''printf(" command line: %s\\n", join(argl, narg).c_str());''')
    #     write( '     %s\n' %'''printf(" \\n");''')
    #     write( '     \n')
    #     write( '     test();\n')
    #     write( '     \n\n')
    #     write( '}\n')
    #     write( '\n\n')

    fpre = 'iter'
    with open(F'bccc{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write
        
        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '# \n')
        write( '# Author: Qingchun Wang @ NJU \n')
        write( '# E-mail: qingchun720@foxmail.com \n')
        write( '# \n')
        write( '\n\n')
        
        write('from copy import deepcopy\n')
        write('from datetime import datetime\n')
        write('from numpy import zeros,ones,array,matrix,ctypeslib,linalg,uint32,float64\n')
        write('from ctypes import c_uint\n')
        write('\n\n')
        
        write('def iter(t, r, h, g, p, nt=4, nc=0, thd=1.0e-7, ncyc=200, xE=0.015):\n')
        write( '    # fun list\n')
        write( '    # ground state\n')
        write(F'    # from bccc.ta.ta_ import ta_ # Ecorr\n')
        write( "    # ftl = (ta_,)\n")
        write( "    ftl = ()\n")
        write(F'    from bccc.ta.ta_ import h0_\n')
        write( "    f0l = (h0_,)\n")
        write( "    ful = ()\n")
        for nt in nel[1:]:
            write(F'    # {nt}-block excited types\n')
            ftl = ftld[nt]; ftll = [ftl[10*i:10*(i+1)] for i in range(0, math.ceil(len(ftl)/10))]
            for ft in ftl: write(F'    from bccc.ta.{ft} import {ft}\n')
            write( '    ftl += (\n')
            for fl in ftll: write(F"            {','.join(fl)},\n")
            write( '            )\n')
            f0l = f0ld[nt]; f0ll = [f0l[10*i:10*(i+1)] for i in range(0, math.ceil(len(f0l)/10))]
            for f0,ft in zip(f0l,ftl): write(F'    from bccc.ta.{ft} import {f0}\n')
            write( '    f0l += (\n')
            for fl in f0ll: write(F"            {','.join(fl)},\n")
            write( '            )\n')
            ful = fuld[nt]; full = [ful[10*i:10*(i+1)] for i in range(0, math.ceil(len(ful)/10))]
            for fu,ft in zip(ful,ftl): write(F'    from bccc.ta.{ft} import {fu}\n')
            write( '    ful += (\n')
            for fl in full: write(F"            {','.join(fl)},\n")
            write( '            )\n')
        write( '    \n')
        
        write( '    from multiprocessing import Pool\n')
        write( '    from bccc.pub import sign,prodl,phiu,lowTri0\n')
        write( '    # h0_\n')
        write( '    h0d = {}\n')
        write( '    pool = Pool()\n')
        write( '    l = tuple(pool.apply_async(f0, (r, h, g, p, nc)) for f0 in f0l)\n')
        write( '    pool.close()\n')
        write( '    pool.join()\n')
        write( '    for h0 in l: h0d.update(h0.get())\n')
        write( '    # hu_\n')
        write( '    hud = {}\n')
        write( '    pool = Pool()\n')
        write( '    l = tuple(pool.apply_async(fu, (r, h, g, p, nc)) for fu in ful)\n')
        write( '    pool.close()\n')
        write( '    pool.join()\n')
        write( '    for hu in l: hud.update(hu.get())\n')
        write( '    \n')
        
        write( "    print(F'threshold,max_cycle,dEshift = {thd},{ncyc},{xE}')\n")
        write( "    print('initial from LCC')\n")
        # write( '    t = {u:h0d.get(u, 0.0)/(Eref-hu) for u,(hu,_) in hud.items()}\n')  # like-bcpt2 guess
        write( '    from itertools import permutations,product\n')
        # write( '        Eref,Ecorr=h0d[()],ta_(t0, r, h, g, p, nc)[()]; E = Eref+Ecorr\n')  # correlation energy
        write( '    Eref,Ecorr=h0d[()],sum(h0d.get(u, 0.0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu) for u,(_,tu) in hud.items()); E0=Eref+Ecorr\n')
        write( "    print(F'Ecorr = {Ecorr}')\n")
        write( '    \n')
        
        # write( '    # Create data for C++\n')
        # write( '    no,n2=h.shape[0],g.shape[0]; lT0=array(lowTri0(no),dtype=uint32)\n')
        # write( '    np,nr,nbs=len(p),312,16; rll = zeros(shape=(np, nr, nbs, nbs))\n')
        # write( '    for P in range(np):\n')
        # write( '        for R,rR in r[P].items():\n')
        # write( '            rll[P,R] = rR\n')
        # write( '    pul=ones(shape=(nt,), dtype=uint32); pul[1]=pu=np*nbs\n')
        # write( '    for i in range(2,nt): pul[i] = pul[i-1]*pu\n')
        write( '    udd = {u:{u:1.0} for u in t}  # all actiove t\n')
        write( '    for u,ud in udd.items():\n')
        write( '        for v in list(permutations(u))[1:]:\n')
        write( '            ud[v] = sign([b for b,s in v if 6<s<15])\n')
        write( '    t0 = {v:s*t[u] for u,ud in udd.items() for v,s in ud.items()}  # unordered t\n')
        # write( '    phid={u:phiu(u, pul) for u in t0}; nu=len(phid)\n')
        # write( '    id,il={u:i for i,u in enumerate(phid.keys())},array(tuple(phid.values()), dtype=uint32)\n')
        # write( '    \n')
        #
        # write( '    from os import path as Path\n')
        # write( '    path = Path.dirname(__file__)\n')
        # write( "    libta = ctypeslib.load_library('libta.so', path)\n")
        # write( '    libta.ta.argtypes = [\n')
        # write( "        ctypeslib.ndpointer(dtype=float64, ndim=1, flags='C_CONTIGUOUS'),  # tl\n")
        # write( "        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # pul\n")
        # write( "        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # il\n")
        # write( "        c_uint,                                                            # nu\n")
        # write( "        ctypeslib.ndpointer(dtype=float64, ndim=4, flags='C_CONTIGUOUS'),  # r\n")
        # write( "        c_uint,                                                            # np\n")
        # write( "        ctypeslib.ndpointer(dtype=float64, ndim=2, flags='C_CONTIGUOUS'),  # h\n")
        # write( "        c_uint,                                                            # no\n")
        # write( "        ctypeslib.ndpointer(dtype=float64, ndim=2, flags='C_CONTIGUOUS'),  # g\n")
        # write( "        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # lT0\n")
        # write( "        c_uint,                                                            # n2\n")
        # write( "        ctypeslib.ndpointer(dtype=uint32, ndim=2, flags='C_CONTIGUOUS'),   # p\n")
        # write( "        c_uint,                                                            # nc\n")
        # write( '                   ]\n')
        # # write( '    libta.ta.restype = POINTER(c_double)\n')
        # # write( "    libta.ta.restype = ctypeslib.ndpointer(dtype=float64, ndim=1, flags='C_CONTIGUOUS')\n")  # bug
        # write( '    libta.ta.restype = ctypeslib.ndpointer(dtype=float64, shape=(nu,))\n')
        # write( '    \n')

        write( "    print('    it          Ecorr              dE               dt               time/s')\n")
        # write( '    al(t0, r, h, g, p, nc)\n')
        # write( '    td_ = ta_(t0, r, h, g, p, nc)\n')
        # write( "    print(F'td_ = {td_}')\n")
        # write( "    exit()\n")
        write( '    ndiis=6; idiis = [*list(range(1,ndiis)), 0, ndiis]\n')
        write( '    objl,errl = [],[]\n')
        write( '    A = zeros((ndiis+1,ndiis+1))-1; A[ndiis,ndiis] = 0\n')
        write( '    b=zeros((ndiis+1, )); b[ndiis]=-1\n')
        write( '    for it in range(ncyc):\n')
        write( '        # ta_\n')
        write( '        t1 = datetime.now()\n')
        write( '        pool = Pool()\n')
        write( '        l = tuple(pool.apply_async(ft, (t0, r, h, g, p, nc)) for ft in ftl)\n')
        write( '        pool.close()\n')
        write( '        pool.join()\n')
        write( '        td = {u: t for e in l for u,t in e.get().items()}\n')
        # write( '        tl = array(tuple(t0.values()), dtype=float64)\n')
        # write( '        l = libta.ta(tl, pul, il, nu, rll, np, h, no, g, lT0, n2, p, nc)\n')
        # write( '        td = {u:l[id[u]] for u in hud}\n')
        write( '        t2 = datetime.now()\n')
        write( '        \n')
        
        write( '        dd = {u: (h0d.get(u, 0.0)+td[u]+(hu-E0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu))/(xE+hu-E0) for u,(hu,tu) in hud.items()}\n')  # new iteraion
        write( '        t,dt = {u:t0[u]-du for u,du in dd.items()},sum([abs(du) for du in dd.values()])\n')
        # write( '        t={u: (td[u]+h0d.get(u, 0.0))/(E0-hu)-sum(s*prodl([t[v] for v in vl]) for s,vl in tu) for u,(hu, tu) in hud.items()}\n')  # old iteration
        # write( '        dd = {u: t[u]-t0[u] for u in t}\n')
        write( '        if dt>thd:\n')
        write( '            objl.append(t); errl.append(dd)\n')
        write( '            if it+1>=ndiis:\n')
        write( '                if it+1>ndiis:\n')
        write( '                    objl.remove(objl[0]); errl.remove(errl[0])\n')
        write( '                    A[:,:] = A[:,idiis]; A[:,:] = A[idiis,:]\n')
        write( '                    for i in range(ndiis):\n')
        write( '                        A[ndiis-1,i] = A[i,ndiis-1] = sum(errl[i][u]*du for u,du in dd.items())\n')
        write( '                else:\n')
        write( '                    for i in range(ndiis):\n')
        write( '                        for j in range(i,ndiis):\n')
        write( '                            A[i,j] = A[j,i] = sum(errl[i][u]*errl[j][u] for u in errl[i])\n')
        write( '                x=linalg.solve(A,b)\n')
        write( '                t = {u:sum(x[i]*objl[i][u] for i in range(ndiis)) for u in t.keys()}\n')
        write( '        Ecorr=sum(h0d.get(u, 0.0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu) for u,(_,tu) in hud.items())\n')
        write( '        E=Eref+Ecorr; dE=E-E0\n')
        write( "        print(F'    {it+1:3d}      {Ecorr:.8f}      {dE:13.10f}    {dt:13.10f}     {t2-t1}', flush=True)\n")
        write( '        if dt<thd: \n')
        write( "            print(F'Successfully converged: Ecorr = {Ecorr}')\n")
        write( '            break\n')
        write( '        else: E0,t0 = E,{v:s*t[u] for u,ud in udd.items() for v,s in ud.items()}\n')
        write( '        \n')
        
        write( "    else: print('Error: Convergence failed ')\n")
        write( '    \n')
        
        write( '    return Ecorr,t\n')
        write( '    \n\n')
        
        # write('\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# Parallel
def parallel(func):
    '''
    parallel

    :param func: detmh() or detat()

    :return:
    '''
    
    # for u in ul: func(u, hfeqd)
    pool = Pool()
    for u in ul: pool.apply_async(func, (u, hfeqd))
    pool.close()
    pool.join()






if __name__ == '__main__':
    print(F'neb = {neb}')
    print(F'nel = {nel}')
    print( 'code= Python')
    # base2theory()
    
    detrfeq()
    # rfeq2theory()
    # rfeq2code()
    
    wexcite()
    # w2theory()
    
    ''' initial from GVB-BCLCC '''
    # deth()
    # h2theory()
    detmt()
    # mt2theory()
    mt2code()
    t1 = datetime.now()
    if sys.platform=='linux':
        parallel(detmh)  # for parallel
        # Pool().map(func, ((u, hfeqd) for u in ul)) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    else:
        for u in ul: detmh(u, hfeqd)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # mh2theory()
    mh2code()
    # init2theory()
    # init2code()
    
    ''' iteration for GVB-BCCC '''
    # dett()
    # t2theory()
    t1 = datetime.now()
    if sys.platform=='linux':
        parallel(detat)  # for parallel
        # Pool().map(func, ((u, hfeqd) for u in ul)) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    else:
        for u in ul: detat(u, hfeqd)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # at2theory()
    # at2code()
    # iter2theory()
    iter2code()
    
    print('End successfully')

