#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   derivform.py
#            Des:   I0Ja+IaJ0 for intermediate array
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


__version__ = '3.1.28'
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
from re import compile, findall
numstr = compile(r'^[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?$') # match integer, fraction, exponent
# numstr = compile(r'[-+]?(\b[0-9]+(\.[0-9]+)?|\.[0-9]+)([eE][-+]?[0-9]+\b)?') # eg. w123.345w
# findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", 'A1.45aa, b5., B6.45, F-2.e2ddf')

from multiprocessing import Pool,Manager
# import mpi4py

# from bccc.addons import
from bccc.pub import index,sign,group,prodl,sum_list




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
        elif oper in ternaryoperl: # debugging ...
            if oper=='+++':
                I0l = [b.l for b in self.leafl[0][0]]
                Al=[]; Il=[]
                for b in self.leafl[0][1]:
                    if b.l in set('IJKL'): Il.append(b.l)
                    else: Al.append(b.l)
                AIl = Al+Il; _AIl = Al+list('IJKL'[:len(Il)])

                func = F"sigma_{''.join(I0l)}_{''.join(AIl)}{len(Expr.ternfunl)}"
                if AIl:
                    argv = F"({','.jion(AIl)}, r, h, g, p)"
                    para = F"({','.join(_AIl)}, r, h, g, p)"
                else: argv = '(r, h, g, p)'

                funclinel = []
                funclinel.append(F'def {func+argv}:')
                funclinel.append(F"    # {','.join(list(chain.from_iterable(self.leafl[0])))}")
                funclinel.append( '    npair = len(p)')
                funclinel.append( '    tmp = 0.0')
                tab = '    '; ntab=1
                for ib,b in enumerate(I0l):
                    if ib: funclinel.append(tab*ntab +F'for {b} in range({I0l[ib-1]}+1,{self.leafl[1][1]}):')
                    else: funclinel.append(tab*ntab +F'for {b} in range({self.leafl[1][0]},{self.leafl[1][1]}):')
                    ntab+=1
                    # if I1ll: funclinel.append(tab*ntab + '%s = %s' % (strize(I1ll+I0ll), strize(list(bIJKL[:len(I0ll)])+I1ll)))
                    if AIl: funclinel.append(tab*ntab +F"if {b} in [{','.join(AIl)}]: continue")

                funclinel.append(tab*ntab+'tmp+='+self.leafl[2].code())
                funclinel.append('    return tmp')

                Expr.ternfunl.append('\n'.join(funclinel)+'\n')

                if AIl: return func+para
                else: return func+argv
            elif oper=='***': return '\\prod_{{{}}}^{{{}}}{{{}}}'.format(*self.leafl)
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

    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='A', l=symbl, f=mode)

    def __str__(self): return F"{{{self.f}}}_{{{','.join(b.l for b in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Array.__format__()')

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

        return hash(F"A{self.f}")
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
        if f is 2: return F'g[ijkl({ol[0]},{ol[2]},{ol[1]},{ol[3]})]'
        elif f is 0: return F'h[{ol[0]},{ol[1]}]'
        elif f is 3: return F'g[ijkl({ol[0]},{ol[3]},{ol[1]},{ol[2]})]'
        elif f is 4: return F'(g[ijkl({ol[0]},{ol[2]},{ol[1]},{ol[3]})]-g[ijkl({ol[0]},{ol[3]},{ol[1]},{ol[2]})])'
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
neb = 4; bel = 'ABCD'; b1l = set(bel)
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

    def indexorb(self, orb): return index(orb, self.functor.l)
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
cgd = {i: Frac(1, math.factorial(i)) for i in range(ngb+1)}
sld = {n: set(sl for sl in product(range(1,nbs), repeat=n)
              if sum(bfd[s].ne for s in sl) is n*2 and sum(bfd[s].sz for s in sl) is 0) for n in range(neb+1)}



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

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
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

    global ul, uld

    # # ground state
    # uld[0] = ((),)
    # # excited state
    # for n in range(1, neb+1):
    #     bll=[]; append = bll.append
    #     for sl in product(range(1,nbs), repeat=n):
    #         if sum(bfd[s].ne for s in sl) is n*2 and sum(bfd[s].sz for s in sl) is 0:
    #             bl = deepcopy(ebl[:n])
    #             for b,s in zip(bl,sl): b.f = s
    #             append(bl)
    #     uld[n] = tuple(bll)
    for n,sl in sld.items():
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

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'wave'
    ftitle = 'Reference wave function in GVB-BCCC formula derivation'
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
        write(F"    ${'+'.join(F'{len(ul)}C^N_{n}' for n,ul in uld.items())}  $ \\\\ \n")
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
def mhuv2code(bra, ced):
    '''
    Matirx of hamiltionian element (i,j) ->  code

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}"
    kbar = ''.join([F'({b.l}, {b.f}),' for b in bra])
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
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

        ibld = {}
        for ket in ced:
            lfl = tuple((b.l, b.f) for b in ket)

            sket = ''.join('{}{}'.format(*lf) for lf in lfl)
            ibl = tuple(lf[0] for lf in lfl)
            ketf = F"{sket}({', '.join(ibl)}, r, h, g, p)"

            key = tuple(filter(lambda lf: lf[1], lfl))
            ib1l = tuple(lf[0] for lf in key)

            try: ibld[ibl][ket] = (ib1l, key, ketf)
            except KeyError: ibld[ibl] = {ket: (ib1l, key, ketf)}

            write(F"def {ketf}:\n")
            write(F'    return ({ced[ket]})\n')
            write( '    \n')
        write('\n\n\n\n')


        write(F'def {fpre}(r, h, g, p):\n')
        write( '    npair = len(p)\n')
        write( '    hdd = {}\n')
        write( '    \n')

        tab = '    '

        ibAl = tuple(b.l for b in bra)
        nbbra = 0
        for b in ibAl:
            if nbbra>0:  write(tab*nbbra+F'    for {b} in range({ibAl[nbbra-1]}+1, npair):\n')
            else: write(F'    for {b} in range(npair):\n')
            nbbra += 1
        write(tab*nbbra+F"    bbl = {{{','.join(ibAl)}}}\n")
        write(tab*nbbra+ '    hd = {}\n')
        write(tab*nbbra+ '    \n')

        for ibl in ibld:
            ibIl = tuple(b for b in ibl if b in b0l)

            nbket = 0
            for b in ibIl:
                if nbket>0: write(tab*(nbbra+nbket)+F'    for {b} in range({ibIl[nbket-1]}+1, npair):\n')
                else: write(tab*(nbbra+nbket)+F'    for {b} in range(npair):\n')
                if ibAl: write(tab*(nbbra+nbket)+F"        if {b} in bbl: continue\n")
                nbket += 1

            nb = nbbra+nbket
            for ket in ibld[ibl]:
                _ib1l, _key, ketf = ibld[ibl][ket]
                key = ''.join('({}, {}),'.format(*lf) for lf in _key)
                ib1l = set(_ib1l)

                if ib1l&b1l and ib1l&b0l: write(tab*nb+F'    key = tuple(sorted(({key})))\n')
                else: write(tab*nb+F'    key = ({key})\n')
                # write(tab*nb+F'    try: hd[key] += {ketf}\n')
                # write(tab*nb+F'    except KeyError: hd[key] = {ketf}\n')
                write(tab*nb+F'    hd[key] = hd.get(key, 0.0)+{ketf}\n')
                write(tab*nb+ '    \n')

        write(tab*nbbra+F"    hdd[({kbar})] = hd\n")
        write(tab*nbbra+ '    \n')

        write('    return hdd\n')
        write('    \n\n\n\n')


        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n\n\n')

        # # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def mhuv2theory(bra, ced):
    '''
    Matirx of hamiltionian element (i,j) -> latex

        subfunction
    :param bra:
    :param ket:

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}"
    ftitle = 'Hamiltonian matrix in GVB-BCCC formula derivation'
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
        sketl = tuple(F'c_{{ {bra}, {ket} }}\\Phi_{{{ket}}}' for ket in ced)
        write(F"    $+{'+'.join(sketl)} $ \\\\ \n")
        write( '    \n')

        write( '    Specifically, \\\\ \n')
        for ket in ced:
            if bra==tuple(b for b in ket if b.f>0):
                write(F'    $c_{{{bra}, {ket}}} = \\langle\\Phi_{{{bra}}}\\vert \\hat{{H}} \\vert\\Phi_{{{ket}}}\\rangle = $ \\\\ \n')
                write(F'    ${ced[ket]!s} $ \\\\ \n')
                write( '    \n')
        write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def h2w(w):
    '''
    hamiltonnian project into wave

    :param w: wave state, list

    :return:
    '''

    hd = {}; keyd = {} # {v: h={<>:ce}} {v: v}
    nb = len(w); ne = sum(b.ne for b in w)
    ibl = tuple(b.l for b in w)

    # for h
    if nb<=2 and ne>=1:  # nb<2 for h; p+p- ne>=1
        for prodb in product(w, repeat=2):
            if len(set(prodb)) is not nb: continue
            _ibl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    inte = Integral(symbl=(bop, boq), mode=0)
                    functor = Functor(symbl=(bop, boq), pdl='+-', upl=prodb)
                    s_blockize,functorl,upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                    _hd={}; _keyd={} # ce {u: ce} 没有与 inte 相乘
                    for a in spinl:
                        bop.f = boq.f = a
                        bfeql = [None]*nb
                        for i,f in enumerate(functorl):
                            bfeq = rfeqd[f][upl[i].f]  # {s: ce, }
                            if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                            else: break
                        else:
                            for prodcxb in product(*bfeql): # 各个 block 的展开，得到系数、波函数态
                                bsl,cel = zip(*prodcxb)
                                key = deepcopy(upl)
                                for i,b in zip(bsl,key): b.f = i
                                if key in hd: key = keyd[key] # 取 key of hd, 如果key已经有了
                                elif key in _hd: key = _keyd[key] # 取 key of _hd, 如果key已经有了
                                ce = 1.0
                                for b,_ce in zip(key, cel):
                                    _ce = deepcopy(_ce)
                                    _ce.l = (b,)+_ce.l
                                    ce = ce*_ce
                                try: _hd[key] = _hd[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                except KeyError: _keyd[key],_hd[key] = key,ce

                    for key in _hd:  # 将 _hd 更新到 hd 中
                        ce = _hd[key]
                        _inte = deepcopy(inte)
                        _inte.upl = tuple(key[index(b, ibl)] for b in _ibl)
                        _inte.orbsort()
                        if s_blockize is -1: ce = 0.0-ce
                        if key in hd:
                            try: hd[key][_inte] = hd[key][_inte]+ce
                            except KeyError: hd[key][_inte] = ce
                        else: keyd[key],hd[key] = key,{_inte: ce}

    # for g
    if ne>=2: # p+p+s-r- ne>=2
        for prodb in product(w, repeat=4):
            if len(set(prodb)) is not nb: continue
            _ibl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    for bor in deepcopy(bo_l):
                        for bos in deepcopy(bo_l):
                            inte = Integral(symbl=(bop, boq, bor, bos), mode=2)
                            functor = Functor((bop, boq, bos, bor), pdl='++--', upl=(prodb[0], prodb[1], prodb[3], prodb[2]))
                            s_blockize, functorl, upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                            _hd = {}; _keyd={}  # ce 没有与 inte 相乘
                            for a in spinl:
                                bop.f = bor.f = boq.f = bos.f = a
                                bfeql = [None]*nb
                                for i,f in enumerate(functorl):
                                    bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                    if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                    else: break
                                else:
                                    for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                        bsl,cel = zip(*prodcxb)
                                        key = deepcopy(upl)
                                        for i,b in zip(bsl, key): b.f = i
                                        if key in hd: key = keyd[key] # 取 key of hd, 如果key已经有了
                                        elif key in _hd: key = _keyd[key] # 取 key of _hd, 如果key已经有了
                                        ce = 1.0
                                        for b,_ce in zip(key, cel):
                                            _ce = deepcopy(_ce)
                                            _ce.l = (b,)+_ce.l
                                            ce = ce*_ce
                                        try: _hd[key] = _hd[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError:  _keyd[key],_hd[key] = key,ce
                            for key in _hd: _hd[key] = 0.5*_hd[key]
                            for a,b in (spinl, ):  # a = alpha, b = beta
                                bop.f = bor.f = a; boq.f = bos.f = b
                                bfeql = [None]*nb
                                for i,f in enumerate(functorl):
                                    bfeq = rfeqd[f][upl[i].f] # {s: ce, }
                                    if bfeq: bfeql[i] = tuple((s,ce) for s,ce in bfeq.items())
                                    else: break
                                else:
                                    for prodcxb in product(*bfeql): # 各个 block 的展开，构成系数, 波函数态
                                        bsl,cel = zip(*prodcxb)
                                        key = deepcopy(upl)
                                        for i,b in zip(bsl, key): b.f = i
                                        if key in hd: key = keyd[key] # 取 key of hd, 如果key已经有了
                                        elif key in _hd: key = _keyd[key] # 取 key of _hd, 如果key已经有了
                                        ce = 1.0
                                        for b,_ce in zip(key, cel):
                                            _ce = deepcopy(_ce)
                                            _ce.l = (b,)+_ce.l
                                            ce = ce*_ce
                                        try: _hd[key] = _hd[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError: _keyd[key],_hd[key] = key,ce

                            for key in _hd:  # 将 _hd 更新到 hd 中
                                ce = _hd[key]
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(key[index(b, ibl)] for b in _ibl)
                                _inte.orbsort()
                                if s_blockize is -1: ce = 0.0-ce
                                if key in hd:
                                    try: hd[key][_inte] = hd[key][_inte]+ce
                                    except KeyError: hd[key][_inte] = ce
                                else: keyd[key],hd[key] = key,{_inte: ce}

    for key in hd:
        hd[key] = sum(ce*inte for inte,ce in hd[key].items())

    return hd
def detmh(bra, hfeqd):
    '''
    determine matrix of hamiltionian (project into wave)

    Note:       Hamiltonian matrix is complex conjugate. This point is often used in the process of deriving formulas
                <W| H  =>  H |W>
                when it stored and output, <bra| H |ket>

    :return:
    '''

    print('bra = ', bra)

    # mH
    hd = {} # {ket: ce}
    for nbcorr in range(1, ngb+1):  # n blocks correlation
        for nbexcite in range(nbcorr+1): # nbexcite blocks from bra
            for bchose in combinations(bra, r=nbexcite): # 从激发态取 nbexcite
                brest = tuple(sorted(set(bra)-set(bchose))) # set -> list is random and disordered
                s_bef = sign([b.l for b in bchose+brest if b.ne%2]) # 激发块提到前面
                _ws = bchose+gbl[:nbcorr-nbexcite]

                key = tuple(b.f for b in _ws); iws = tuple(b.l for b in _ws)
                if key in hfeqd:
                    _hd = deepcopy(hfeqd[key])
                    try:
                        bl = tuple(b.l for b in list(_hd.keys())[0])
                        if bl != iws:
                            for bl in _hd:
                                for b, i in zip(bl, iws): b.l = i
                    except IndexError: continue   # _hd = {} when nbcorr=nbexcite, b.f=6
                else:
                    _hd = h2w(_ws)  # Aa I0 J0 Kk
                    hfeqd[key] = _hd

                for key,ce in _hd.items(): # 将 _hd 更新(多块加和)到 hd 中
                    ket = key+brest

                    leaf = [b.l for b in ket if b.ne%2]; ib1l = [b.l for b in ket if b.ne%2 and b.l in b1l]
                    _leaf = set(leaf); _ib1l = sorted(ib1l)
                    if _leaf&b1l and (_leaf&b0l or _ib1l!=ib1l): # 还原激发块
                        s_lat = Expr('-?', leaf)
                        ce = s_lat*ce
                    if s_bef is -1: ce = 0.0-ce

                    ket = tuple(sorted(ket))
                    try: hd[ket] = hd[ket]+ce
                    except KeyError: hd[ket] = ce

    # mhuv2theory(bra, hd)
    mhuv2code(bra, hd)
def mh2theory():
    '''
    matrix of hamiltonian (project into wave) -> tex file

    Note:    main function output Hamiltionian and H|W> = ?
    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'hm'
    ftitle = 'Hamiltonian matrix in GVB-BCCC formula derivation'
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
        write( '    $\\hat{H}\\vert{\\Phi_{u}}\\rangle = \\displaystyle\\sum_{u,v}{c_{u,v}\\Phi_{v}}$ \\\\ \n')
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

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'hm'
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b/%d %Y\n"))
        write('#\n')
        write('\n\n')

        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')

        write('from multiprocessing import Pool, Manager\n')
        write('from importlib import import_module\n')
        write('import sys, gc\n')
        write('\n\n')

        write('def calculate(braf, r, h, g, p):\n')
        write("    module = F'bccc.hm.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    hd = fun(r, h, g, p)\n')
        write('    del sys.modules[module]\n')
        write('    gc.collect()\n')
        write('    return hd\n')
        write('\n\n')

        write('def hm(r, h, g, p, ntuply=4):\n')
        write('    hdd = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        # write(F'    from bccc.hm.hm_ import hm_\n')
        write("    brafl = ('hm_',)\n")
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[1])
        # for braf in brafl: write(F'    from bccc.hm.{braf} import {braf}\n')
        write('    brafl += (\n             ')
        for i,braf in enumerate(brafl):
            write(F"'{braf}', ")
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for ntuply in range(2,neb+1):
            write(F'    # {ntuply}-block excited types\n')
            write(F'    if ntuply>{ntuply-1}: \n')
            brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[ntuply])
            # for braf in brafl: write(F'        from bccc.hm.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F"'{braf}', ")
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write('    pool = Pool()\n')
        write('    vl = tuple(pool.apply_async(calculate, (braf, r, h, g, p)) for braf in brafl)\n')
        write('    pool.close()\n')
        write('    pool.join()\n')
        write('    for v in vl: hdd.update(v.get())\n')
        write('    # for braf in brafl: hdd.update(calculate(braf, r, h, g, p))\n')
        write('    \n')
        write('    return hdd\n')
        write('    \n\n\n\n')


        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
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

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 't'
    ftitle = ' T operator '
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
            sigmal = tuple("\\displaystyle\\sum_{{{}}}\\displaystyle\\sum_{}".format('>'.join(reversed(bl[:j+1])), bl_[j]) for j in range(len(bl)))
            operl = tuple(F'{bl[j]}_{bl_[j]}^+{bl[j]}_0^-' for j in range(len(bl)))
            write(F"    $T_{i} = {''.join(sigmal)}{{{''.join(operl)}}} $ \\\\ \n")
            for bl in tld[i]:
                write(F"    ${','.join(F'{b}={bl[ib].f}' for ib,b in enumerate(bl_))} $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()

# T matrix for LCC
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

    for nbbra in range(1, neb+1):
        for bra in uld[nbbra]:
            tmp = {}
            bbl = set(bra)
            for nbt in range(1, nbbra+1):
                for t in combinations(bra, nbt):
                    if sum(b.ne for b in t) is 2*len(t) and sum(b.sz for b in t) is 0:
                        ket = tuple(sorted(bbl - set(t)))
                        s = sign(tuple(b.l for b in t+ket if 6<b.f<15))
                        tmp[t] = (s, ket)
            mt[tuple(bra)] = tmp
def mt2theory():
    '''
    matrix of T (project into wave) -> tex file

    :return:
    '''

    def mtij2theory(bra, ket):
        raise NotImplementedError('TODO: mtij2theory() in mt2theory()')

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'tm'
    ftitle = 'T matrix in GVB-BCCC formula derivation'
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

        for bra in mt:
            write(F'    $\\Phi_{{{bra!s}}} = $ \\\\ \n')
            for t in mt[bra]:
                if mt[bra][t][0] is 1: write(F"    $+ {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{mt[bra][t][1]}}} $ \\\\ \n")
                else: write(F"    $- {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{mt[bra][t][1]}}} $ \\\\ \n")
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

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    def mtij2code(bar, ket):
        raise NotImplementedError('TODO: mtij2code() in mt2code()')

    fpre = 'tm'
    with open(F'bccc{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b/%d %Y\n"))
        write('#\n')
        write('\n\n')


        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')


        bld = {}
        for bra in mt:
            bl = tuple(b.l for b in bra)
            try: bld[bl].append(bra)
            except KeyError: bld[bl] = [bra]

        write('def tm(npair, ntuply=4):\n')
        write('    tmp = {}\n')
        write('    \n')

        tab = '    '
        for bl in bld:
            nb = 0
            for b in bl:
                if nb>0:  write(tab*nb+F'    for {b} in range({bl[nb-1]}+1, npair):\n')
                else: write(F'    for {b} in range(npair):\n')
                nb += 1

            for bra in bld[bl]:
                kbra = ''.join(F'({b.l},{b.f}),' for b in bra)
                write(tab*nb+ '    _t = {}\n')
                for t in mt[bra]:
                    kt = ''.join(F'({b.l},{b.f}),' for b in t)
                    kket = ''.join(F'({b.l},{b.f}),' for b in mt[bra][t][1])
                    write(tab*nb+F'    _t[({kt})] = ({mt[bra][t][0]}, ({kket}))\n')
                write(tab*nb+F'    tmp[({kbra})] = _t\n')
                write(tab*nb+ '    \n')

            if nb<4: write(F'    if ntuply is {nb}: return tmp\n')
            else: write( '    return tmp\n')
            write( '    \n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
# Amplitude fof CC
nld = {n: tuple(sorted(set([tuple(sorted(tn, reverse=True)) for nb in range(1,n+1)
            for tn in product(range(1,neb+1), repeat=nb) if sum(tn) is n]), reverse=True)) for n in range(neb+ngb+1)}
gld = {nl: group(set(range(n)), nl) for n in range(neb+ngb+1) for nl in nld[n]}
iad = {}
# ia = Manager().dict()
def atu2theory(bra, td):
    '''
    array of T element u -> latex

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = F"at_{''.join(F'{b.l}{b.f}' for b in bra)}"
    ftitle = 'Hamiltonian matrix in GVB-BCCC formula derivation'
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

        # bra
        bsl = tuple(b.f for b in bra); nb = len(bsl)
        tbra_ = []; append = tbra_.append
        for tl in nld[nb][1:]:
            for gl in gld[tl]:
                if all(tuple(bsl[i] for i in g) in tsd[len(g)] for g in gl):
                    append(Coeff(symbl=tuple(tuple(bra[i] for i in g) for g in gl), mode='t'))

        if bra:
            ce = (Coeff(symbl=((), bra))+sum(t*Coeff(symbl=(bra,ket), mode='mh') for ket,t in td.items()))/(Coeff(mode='E')-Coeff(symbl=(bra,bra), mode='mh'))-sum(tbra_)
            write(F"    ${Coeff(symbl=bra, mode='t')!s} = {ce!s}$ \\\\ \n")
        else:
            ce = sum(t*Coeff(symbl=(bra,ket), mode='mh') for ket,t in td.items())
            write(F"    $E_{{corr}} = {ce!s}$ \\\\ \n")
        write( '    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
def atu2code(bra, hld):
    '''
    array of T element u -> code

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname,ossep = os.name,os.sep
    oslinesep = os.linesep

    fpre = F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}"
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
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

        tab = '    '
        bbl=[b.l for b in bra]; kbra=F"tuple(sorted([{''.join(F'{b},' for b in bra)}]))"

        fl=[]; flapp=fl.append
        for v,hl in hld.items():
            f = F"_{''.join(F'{b.l}{b.f}' for b in v)}({''.join(F'{b}, ' for b in bbl)}r, h, g, p)"
            flapp(f)

            kbl=[b.l for b in v if b.l in b0l]; kket=F"tuple(sorted([{''.join(F'{b},' for b in v)}]))"

            write(F"def {f}:\n")
            write( '    np = len(p)\n')
            write( '    hd = {}\n')
            write( '    \n')

            nkt=0
            for b in kbl:
                write(tab*nkt+F'    for {b} in range(np):\n')
                write(tab*nkt+F"        if {b} in {{{','.join(bbl+kbl[:nkt])}}}: continue\n")
                nkt += 1
            tabnkt = tab*nkt
            write(tabnkt+F"    v = tuple(sorted([{','.join(F'{b}' for b in v)}]))\n")
            write(tabnkt+ '    hv = 0.0\n')
            write(tabnkt+ '    \n')

            for h in hl:
                Il = sorted(set([b.l for b in h[1].upl if b.l in b0l and b.l not in kbl]))
                if Il:
                    forl=[]; forlapp=forl.append
                    nI = 0
                    for b in Il:
                        forlapp(F"for {b} in range(np) if {b} not in [{','.join(bbl+kbl+Il[:nI])}]")
                        nI += 1
                    write(tabnkt+F"    hv += sum({reduce(mul, h)} {' '.join(forl)})\n")
                else: write(tabnkt+F'    hv += {reduce(mul, h)}\n')
            write(tabnkt+ '    \n')

            write(tabnkt+ '    hd[v] = hd.get(v,0.0)+hv\n')
            write(tabnkt+ '    \n')

            write(F'    return hd\n')
            write( '    \n')

        write(F'def {fpre}(r, h, g, p):\n')
        write( '    np = len(p)\n')
        write( '    hdd = {}\n')
        write( '    \n')

        nbt = 0
        for b in bbl:
            write(tab*nbt+F'    for {b} in range(np):\n')
            if nbt: write(tab*nbt+F"        if {b} in {{{','.join(bbl[:nbt])}}}: continue\n")
            nbt += 1
        tabnbt = tab*nbt
        write(tabnbt+F'    u = {kbra}\n')
        write(tabnbt+ '    \n')

        for f in fl:
            write(tabnbt+F'    hd = {f}\n')
            write(tabnbt+ '    try:\n')
            write(tabnbt+F'        for v,hv in hd.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
            write(tabnbt+ '    except KeyError: hdd[u] = hd\n')
            write(tabnbt+ '    \n')
        write(tabnbt+ '    \n')

        write('    return hdd\n')
        write('    \n\n\n\n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n\n\n')

        # # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
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
                        _inte.upl = tuple(v[index(b, bl)] for b in _bl)
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
                                _inte.upl = tuple(v[index(b, bl)] for b in _bl)
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
                                _inte.upl = tuple(v[index(b, bl)] for b in _bl)
                                # _inte.orbsort()
                                si = [Frac(s_blockize, 1), _inte]
                                try:
                                    for rl in _hl: hld[v].append(si+rl)
                                except KeyError:
                                    for i,rl in enumerate(_hl): _hl[i] = si+rl
                                    vd[v],hld[v] = v,_hl

    # for v,hl in hld.items():
    #     hld[v] = sum(reduce(mul, h) for h in hl)
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
    nb = len(bra); cbra = cgd[nb]

    # mH
    hld = {} # {ket: ce}
    for nbc in range(1, ngb+1):  # n blocks correlation
        for nbe in range(nbc+1): # nbe blocks from bra
            nb0=nbc-nbe; cket=cbra*cgd[nb0]
            for bchose in combinations(bra, r=nbe): # 从激发态取 nbe
                brest = tuple(sorted(set(bra)-set(bchose))) # set -> list is random and disordered
                s_bef = Expr('-?', [b.l for b in bchose+brest if b.ne%2]) # 激发块提到前面
                _u = bchose+gbl[:nb0]

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

                for v,hl in _hld.items(): # 将 _hld 更新(多块加和)到 hld 中
                    nI1 = len([b for b in v[nbe:] if b.f])
                    if nb0>1 and nI1 is not nb0: # I0Ja -> J0Ia
                        i1 = 0; i0 = nI1
                        for b in v[nbe:]:
                            if b.f:
                                b.l = bgl[i1]
                                i1 += 1
                            else:
                                b.l = bgl[i0]
                                i0 += 1

                    _u = v+deepcopy(brest); u = tuple(sorted(b for b in _u if b.f))
                    s_lat = Expr('-?', [b.l for b in _u if b.ne%2])
                    for h in hl:
                        h[0] = h[0]*cket*s_bef*s_lat
                        h[1].orbsort()
                    hld[u] = hld.get(u, [])+hl

    # # at
    # td = {}  # {ket: t}
    # ketl = set(); append = ketl.add
    # for _ket in hld.keys(): # Aa..Jj.. -> Aa..Ij..
    #     bbl = tuple(b for b in _ket if b.l in b1l and b.f>0); kbl = deepcopy(tuple(b for b in _ket if b.l in b0l and b.f>0))
    #     for b,l in zip(kbl, 'IJKL'): b.l = l
    #     append(bbl+kbl)
    # ketl.remove(bra)
    # for ket in ketl:
    #     bl, bsl = tuple(b.l for b in ket), tuple(b.f for b in ket); nb = len(bsl)
    #
    #     tket_ = []; append = tket_.append # 有效的分组
    #     for tl in nld[nb]:
    #         for gl in gld[tl]:
    #             if all(tuple(bsl[i] for i in g) in tsd[len(g)] for g in gl): append(gl)
    #     txtl = []; append=txtl.append
    #     for gl in tket_:
    #         oddbl = [bl[i] for i in chain(*gl) if 6<bsl[i]<15]; _oddbl = set(oddbl)
    #         if _oddbl&b0l and _oddbl&b1l:
    #             oddbl = []
    #             for g in gl:
    #                 _t = tuple(bl[i] for i in g if 6<bsl[i]<15)
    #                 if _t: oddbl.append(_t)
    #             s = Expr('-?', oddbl)
    #             append(s*reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))
    #         else:
    #             s = sign(oddbl)
    #             if s is -1: append(0.0-reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))
    #             else: append(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))
    #
    #     if txtl:
    #         # td[ket] = sum(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)) for gl in tket_)
    #         # td[ket] = sum(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)) for gl in tket_)
    #         # td[ket] = sum(txtl)
    #         td[ket] = txtl
    #         # td[ket] = tuple(tuple(tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g) for g in gl) for gl in tket_)

    # atu2theory(bra, hld)
    atu2code(bra, hld)
    # mhuv2theory(bra, hld)
    # mhuv2code(bra, hld)
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

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    # fpre = 'at'
    # with open(F'bccc{ossep}at{ossep}{fpre}.py', 'w') as fout:
    #     # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
    #     write = fout.write
    #
    #     write( '#!/usr/bin/env python \n')
    #     write( '# -*- coding: utf-8 -*- \n')
    #     write( '#\n')
    #     write( '#      Copyright:   Qingchun Wang @ NJU\n')
    #     write(F'#      File Name:   {fpre}.py\n')
    #     write(F'#            Des:   {fpre}\n')
    #     write( '#           Mail:   qingchun720@foxmail.com\n')
    #     write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b/%d %Y\n"))
    #     write('#\n')
    #     write('\n\n')
    #
    #     write("__version__ = '1.0'\n")
    #     write("__author__ = 'Qingchun Wang'\n")
    #     write('\n\n')
    #
    #     write('from multiprocessing import Pool, Manager\n')
    #     write('from copy import deepcopy\n')
    #     write('from datetime import datetime\n')
    #     write('from importlib import import_module\n')
    #     write('import sys, gc\n')
    #     write('\n\n')
    #
    #     write('def calculate(braf, t0, h, p, E0):\n')
    #     write("    module = F'bccc.hm.{braf}'\n")
    #     write('    lib = import_module(module)\n')
    #     write('    fun = getattr(lib, braf)\n')
    #     write('    t1 = fun(t0, h, p, E0)\n')
    #     write('    del sys.modules[module]\n')
    #     write('    gc.collect()\n')
    #     write('    return t1\n')
    #     write('\n\n')
    #
    #     write('def at(t0, h, p, ntuply=4):\n')
    #     write('    t = {}\n')
    #     write('    \n')
    #
    #     write('    # bra wavefunction list\n')
    #     write('    # ground state\n')
    #     write('    from bccc.at.at_ import at_\n')
    #     write('    \n')
    #
    #     write('    # 1-block excited types\n')
    #     brafl = tuple(F"at_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[1])
    #     for braf in brafl: write(F'    from bccc.at.{braf} import {braf}\n')
    #     write('    brafl = (\n             ')
    #     for i,braf in enumerate(brafl):
    #         write(F'{braf}, ')
    #         if (i+1)%5 is 0: write('\n             ')
    #     write('             \n')
    #     write('             )\n')
    #     write('    \n')
    #
    #     for ntuply in range(2,neb+1):
    #         write(F'    # {ntuply}-block excited types\n')
    #         write(F'    if ntuply>{ntuply-1}: \n')
    #         brafl = tuple(F"at_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[ntuply])
    #         for braf in brafl: write(F'        from bccc.at.{braf} import {braf}\n')
    #         write( '        brafl += (\n                 ')
    #         for i,braf in enumerate(brafl):
    #             write(F'{braf}, ')
    #             if (i+1)%5 is 0: write('\n                 ')
    #         write('                 \n')
    #         write('                 )\n')
    #         write('        \n')
    #
    #     write("    print('    it            dE               dt               time/s')\n")
    #     write('    E0 = h[()][()]+at_(t0, h, p)\n')
    #     write('    for it in range(200):\n')
    #     write('        t1 = datetime.now()\n')
    #     write('        \n')
    #
    #     write('        pool = Pool()\n')
    #     write('        vl = tuple(pool.apply_async(braf, (t0, h, p, E0)) for braf in brafl)\n')
    #     write('        pool.close()\n')
    #     write('        pool.join()\n')
    #     write('        for v in vl: t.update(v.get())\n')
    #     write('        # for braf in brafl: t.update(braf(t0, h, p, E0))')
    #     write('        \n')
    #
    #     write('        t2 = datetime.now()\n')
    #     write('        \n')
    #
    #     write('        Ecorr = at_(t, h, p); E = h[()][()]+Ecorr\n')
    #     write('        dE = E-E0; dt = sum(abs(t[u]-t0[u]) for u in t)\n')
    #     write("        print(F'    {it+1:3d}      {dE:13.10f}    {dt:13.10f}     {t2-t1}')\n")
    #     write('        if abs(dt)<1.0e-6: \n')
    #     write("            print(F'Successfully converged: Ecorr = {Ecorr}')\n")
    #     write('            break\n')
    #     write('        else: \n')
    #     write('            t0 = deepcopy(t)\n')
    #     write('            E0 = E\n')
    #     write('        \n')
    #     write('    else:\n')
    #     write("        print('Error: Convergence failed ')\n")
    #     write('    \n')
    #
    #     write('    return Ecorr,t\n')
    #     write('    \n\n\n\n')

    fpre = 'hm'
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b/%d %Y\n"))
        write('#\n')
        write('\n\n')

        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')

        write('from multiprocessing import Pool, Manager\n')
        write('from importlib import import_module\n')
        write('import sys, gc\n')
        write('\n\n')

        write('def calculate(braf, r, h, g, p):\n')
        write("    module = F'bccc.hm.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    hd = fun(r, h, g, p)\n')
        write('    del sys.modules[module]\n')
        write('    gc.collect()\n')
        write('    return hd\n')
        write('\n\n')

        write('def hm(r, h, g, p, ntuply=4):\n')
        write('    hdd = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        # write(F'    from bccc.hm.hm_ import hm_\n')
        write("    brafl = ('hm_',)\n")
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[1])
        # for braf in brafl: write(F'    from bccc.hm.{braf} import {braf}\n')
        write('    brafl += (\n             ')
        for i,braf in enumerate(brafl):
            write(F"'{braf}', ")
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for ntuply in range(2,neb+1):
            write(F'    # {ntuply}-block excited types\n')
            write(F'    if ntuply>{ntuply-1}: \n')
            brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[ntuply])
            # for braf in brafl: write(F'        from bccc.hm.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F"'{braf}', ")
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write('    pool = Pool()\n')
        write('    pl = tuple(pool.apply_async(calculate, (braf, r, h, g, p)) for braf in brafl)\n')
        write('    pool.close()\n')
        write('    pool.join()\n')
        write('    # for p in pl: hdd.update(p.get())\n')
        write('    for p in pl:\n')
        write('        for u,hu in p.get().items():\n')
        write('            try:\n')
        write('                for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        write('            except KeyError: hdd[u] = hu\n')
        write('    # # for braf in brafl: hdd.update(calculate(braf, r, h, g, p))\n')
        write('    # for braf in brafl:\n')
        write('    #     for u,hu in calculate(braf, r, h, g, p).items():\n')
        write('    #         try:\n')
        write('    #            for v,hv in hu.items(): hdd[u][v] = hdd[u].get(v, 0.0)+hv\n')
        write('    #         except KeyError: hdd[u] = hu\n')
        write('    \n')
        write('    return hdd\n')
        write('    \n\n\n\n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()






# CC iteration
def iterate2theory():
    raise NotImplementedError('TODO: iterate2theory()')
def iterate2code():
    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = 'iterate'
    with open(F'bccc{ossep}hm{ossep}{fpre}.py', 'w') as fout:
        # stdout = sys.stdout; stderr = sys.stderr; sys.stdout = fout; sys.stderr = fout
        write = fout.write

        write( '#!/usr/bin/env python \n')
        write( '# -*- coding: utf-8 -*- \n')
        write( '#\n')
        write( '#      Copyright:   Qingchun Wang @ NJU\n')
        write(F'#      File Name:   {fpre}.py\n')
        write(F'#            Des:   {fpre}\n')
        write( '#           Mail:   qingchun720@foxmail.com\n')
        write(datetime.now().strftime("#   Created Time:   %a %H:%M:%S %b/%d %Y\n"))
        write('#\n')
        write('\n\n')

        write("__version__ = '1.0'\n")
        write("__author__ = 'Qingchun Wang'\n")
        write('\n\n')

        # write('from multiprocessing import Pool, Manager\n')
        # write('\n\n')

        write('def corr(mh, ntuply=4):\n')
        write('    tmp = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        brafl = ('hm_',)
        write('    from bccc.hm.hm_ import hm_\n')
        write('    \n')

        write('    # 1-block excited types\n')
        brafl += tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[1])
        for braf in brafl: write(F'    from bccc.hm.{braf} import {braf}\n')
        write('    brafl = (\n             ')
        for i,braf in enumerate(brafl):
            write(F'{braf}, ')
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for ntuply in range(2,neb+1):
            write(F'    if ntuply>{ntuply-1}:')
            write(F'        # {ntuply}-block excited types\n')
            brafl += tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in uld[ntuply])
            for braf in brafl: write(F'        from bccc.hm.{braf} import {braf}\n')
            write( '        brafl = (\n                 ')
            for i,braf in enumerate(brafl):
                write(F'{braf}, ')
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        # write('    pool = Pool()\n')
        # write('    vl = {}\n')
        # write('    for i,braf in enumerate(brafl):  vl[i] = pool.apply_async(braf, (r, h, g, p))\n')
        # write('    pool.close()\n')
        # write('    pool.join()\n')
        # write('    for i in vl:  tmp.update(vl[i].get())\n')
        write('    \n')
        write('    return tmp\n')
        write('    \n\n\n\n')


        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# Parallel
def parallel():
    '''
    parallel

    :param func: detmh() or detmt()

    :return:
    '''

    mp = Pool()
    # for u in ul: mp.apply_async(detmh, (u, hfeqd))
    for u in ul: mp.apply_async(detat, (u, hfeqd))
    mp.close()
    mp.join()




if __name__ == '__main__':
    # base2theory()

    detrfeq()
    # rfeq2theory()
    # rfeq2code()

    wexcite()
    # w2theory()

    t1 = datetime.now()
    # if sys.platform=='linux':
    #     parallel()  # for parallel
    #     # Pool().map(detmh, ((u, hfeqd) for u in ul)) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    # else:
    #     for u in ul: detmh(u, hfeqd)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # mh2theory()
    # mh2code()

    # dett()
    # t2theory()
    # detmt()
    # mt2theory()
    # mt2code()
    t1 = datetime.now()
    if sys.platform=='linux':
        parallel()  # for parallel
        # Pool().map(detmh, ((u, hfeqd) for u in ul)) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    else:
        for u in ul: detat(u, hfeqd)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # at2theory()
    # at2code()

    print('End successfully')

