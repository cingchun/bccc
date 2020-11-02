#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   derivform.py
#            Des:   Simplify expression of block state
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
    ws -> wave in detmh()
    __str__ -> latex; mean -> __str__(__repr__)
    different __repr__ with __str__
    
    copy(), deepcopy(), __copy__(), __deepcopy__()
    __bool__() -> Symb0; __bool__() -> Expr0
    str '0' in Block.l -> number 0 and similar 
'''


__version__ = '3.1.23'
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




'''
operator 
    include function and operator in writing form
        function:        ?( , )
        character:       ( )?( )
    include unary, binary and ternary from the point of view of the number of operations
        unary operator belong to function
        binary operator partly belong to functional, and part belong to symbol
        ternary operator belong to function
'''
# unary operator
unarycharl = {'+1', '-1',
              '0+', '0-'}
signoperl = {'0-', '0+'}
unaryfuncl = {'abs', 'fabs',
              'ceil', 'floor', 'round',
              'exp', 'sqrt', 'lg', 'ln',
              'sin', 'cos', 'asin', 'acos', 'tan', 'cot', 'atan', 'acot',
              'int', 'float', 'str'}
unaryoperl = unarycharl|unaryfuncl
'''
            0+, 0-:    sign operator in mathematics
            +1, -1:    ++, -- in computer science
ceil, floor, round:    Upward, downward integrate, and rounding of value
   int, float, str:    int(), float(), str() in code
'''
def isunary(oper):
    '''
    is unary operator ?

    :param oper:

    :return:
    '''

    return oper in unaryoperl
# binary operator
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
'''
                 !=:    not and equal to, like += 
                 <>:    not eqaul
                 //:    floor divid, with no remainder
                  %:    mod, the remainder
                 ++:    times e.g. x++3 = x+x+x = 3x
            **, pow:    power function
            %%, log:    log function
                  ,:    each one in sigma or pi cycle, e.g. i, j
                  ::    ben -> end in sigma or pi cycle,  e.g. 0->n
                 -?:    function for block swap
'''
def isbinary(oper):
    '''
    is binary operator ?

    :param oper:

    :return:
    '''

    return oper in binaryoperl
# ternary operator
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
'''
                 ?::    if-then-else, [bool, True-Expr, False-Expr] 
                +++:    sum sigma, [ben, end, expr]
                ***:    prod pi, [ben, end, expr]
               proj:    put project, [bra, ket, functor]
'''
def isternary(oper):
    '''
    is ternary operator ?

    :param oper:

    :return:
    '''

    return oper in ternaryoperl

# priority of operator
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

# functional map of operator
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




'''
number
    int and float
    include number-Symb(Num), and number-Expr
    only number-Symb(Num), and number-Expr can be be converted into number
'''
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




'''
Two base class:
    Symb and Expr
    For convenience, '' stands for None in all of the following types
'''
# 在之前版本中，实现了Symb, Expr对 int/float, str 等类型的兼容
# 但是，在实际中，多数都是 Symb, Expr 的计算，而非 int/float, str 的计算
# 对 int/float, str 兼容，增加了对 type 的判断，会使行计算做很多无效计算，得不偿失
# 因此，在这个版本，删除对 int/float, str 兼容
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

    # def __new__(cls, *args, **kwargs): pass
    def __init__(self, c=None, l=None, f=None): self.c, self.l, self.f = c, l, f
    # def __del__(self): pass

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

    # def copy(self):
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
    # def __copy__(self):
    #     return self.copy()
    # def deepcopy(self):
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
    # def __deepcopy__(self, memodict={}):
    #     return self.deepcopy()

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

    # def __new__(cls, *args, **kwargs): pass
    def __init__(self, oper, *symbl):
        '''
        init: *symb, oper

        :param symb:
        :param oper:
        '''

        self.root, self.oper, self.leafl = None, oper, symbl
        # for leaf in symbl:
        #     if isinstance(leaf, Expr): leaf.root = self
    # def __del__(self): pass

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
                if oper=='-?': return F"(-1)^{{P({','.join(leafl[0])})}}"
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
                return F"\\displaystyle\\sum_{{{leafl[1]}}}^{{{leafl[2]}}} {leafl[0]} \\\\ \n"
            elif oper=='***':
                return F"\\prod_{{{leafl[1]}}}^{{{leafl[2]}}} {leafl[0]} \\\\ \n"
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
                    leafl0 = self.leafl[0]
                    if isinstance(leafl0[0],str): return F"sign([{','.join(leafl0)}])"
                    else: return F'''sign({'+'.join(F"sorted([{','.join(bl)}])" for bl in leafl0)})'''
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
                    argv = F"({','.jion(AIl)}, r, h_mo, g_mo, p)"
                    para = F"({','.join(_AIl)}, r, h_mo, g_mo, p)"
                else: argv = '(r, h_mo, g_mo, p)'

                funclinel = []
                funclinel.append(F'def {func+argv}:')
                funclinel.append(F"    # {','.join(list(chain.from_iterable(self.leafl[0])))}")
                funclinel.append( '    np = len(p)')
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

    # def copy(self):
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
    # def __copy__(self):
    #     return self.copy()
    # def deepcopy(self):
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
    # def __deepcopy__(self, memodict={}):
    #     return self.deepcopy()

    # 在这个版本中，删除对0,1以及 特殊运算 的处理
    # 比如：0+s=s+0=s-0=s*1=1*s=s/1=s; (0+a)+(1*b)=(1*a)+(1*b)=(1*a)-(0+b)=a+b
    # 支持 0,1 特殊处理，计算变得杂，实际符号运算中，这并不常见
    # 但是，特别是 1*s=s 保留, 0-e=-e
    # 有效的运算包括：s,e +-*/ N,s,e 和 N +-*/ s,e
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
        if mode is 't':
            bl = [b.l for b in symbl]; _bl = set(bl)
            if _bl&iwb0l and _bl&iwb1l: return F"t[tuple(sorted(({','.join(F'({b.l},{b.f})' for b in symbl)})))]"
            else: return F"t[({''.join(F'({b.l},{b.f}),' for b in symbl)})]"
        elif mode is 'h':
            bra,ket = symbl; nbb,nbk = len(bra),len(ket)
            kbra = F"({''.join(F'({b.l},{b.f}),' for b in bra)})"
            bl = [b.l for b in ket]; _bl = set(bl)
            if _bl&iwb0l and _bl&iwb1l: kket = F"tuple(sorted(({','.join(F'({b.l},{b.f})' for b in ket)})))"
            else: kket = F"({''.join(F'({b.l},{b.f}),' for b in ket)})"
            return F"h[{kbra}][{kket}]"
        else:
            if nsymb is 4: return F'r[{symbl[0].l}][{symbl[1]}][{symbl[2]},{symbl[3]}]'
            elif nsymb is 2: return F'{self.f}[{symbl[0]},{symbl[1]}]' # for rep.py
            elif nsymb is 3: return F'r[{symbl[1].l}][{symbl[0]}][{symbl[1].f},{symbl[2].f}]'
            else:
                raise NotImplementedError(F'TODO: Coeff.__format__() for len(Symb.l) = {len(symbl)}')
class Delta(Symb):
    '''
    Delta class

    :param symbl: symbol set
    '''

    def __init__(self, symbl=(), mode=''): Symb.__init__(self, c='D', l=symbl, f=mode)

    def __str__(self): return F"\\delta^{{{self.f}}}_{{{','.join(F'{e!s}' for e in self.l)}}}"
    __repr__ = __str__
    def __format__(self, format_spec):
        raise NotImplementedError('TODO: Delta.__format__()')
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
        norb = len(self.l); f = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        if upl: sl = {i: F'{symb.l}^{{{upl[i]}}}' for i,symb in enumerate(self.l)}
        else: sl = {i: symb.l for i, symb in enumerate(self.l)}
        if f is 2: return F'\\langle{{{sl[0]}{sl[1]}}}\\vert \\hat{{v}} \\vert{{{sl[2]}{sl[3]}}}\\rangle \\\\ \n    ' # 公式输出时会自动换行 ？
        elif f is 0: return F'\\langle{{{sl[0]}}}\\vert \\hat{{h}} \\vert{{{sl[1]}}}\\rangle \\\\ \n    '
        elif f is 3: return F'\\langle{{{sl[0]}{sl[1]}}}\\vert \\hat{{v}} \\vert{{{sl[3]}{sl[2]}}}\\rangle \\\\ \n'
        elif f is 4:
            return '\\langle{{{0}{1}}}\\vert \\hat{{v}} \\vert{{{2}{3}}}\\rangle ' \
                   '- \\langle{{{0}{1}}}\\vert \\hat{{v}} \\vert{{{3}{2}}}\\rangle \\\\ \n'.format(*sl.values())
        elif f is 1: return F'\\langle{{{sl[0]}}}\\vert \\hat{{f}} \\vert{{{sl[1]}}}\\rangle \\\\ \n'
        else:
            raise ValueError(F'unrecognized {f} for Symb.f in Symb.__str__()')
    __repr__ = __str__
    def __format__(self, format_spec):
        norb = len(self.l); f = self.f; upl = {i: up.l for i, up in enumerate(self.upl)}
        sl = {i: F'p[{upl[i]},{symb.l}]' for i, symb in enumerate(self.l)}
        if f is 2: return F'g_mo[ijkl({sl[0]},{sl[2]},{sl[1]},{sl[3]})] \n                '
        elif f is 0: return F'h_mo[{sl[0]},{sl[1]}] \n                '
        elif f is 3: return F'g_mo[ijkl({sl[0]},{sl[3]},{sl[1]},{sl[2]})] \n                '
        elif f is 4: return F'(g_mo[ijkl({sl[0]},{sl[2]},{sl[1]},{sl[3]})] - g_mo[ijkl({sl[0]},{sl[3]},{sl[1]},{sl[2]})]) \n                '
        elif f is 1: return F'f_mo[{sl[0]},{sl[1]}] \n                '
        else:
            raise ValueError(F'unrecognized Symb.f = {f} for Integral in Symb.__format__()')

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"I_{''.join(o.l for o in self.l)}_{''.join(b.l for b in self.upl)}_{self.f}" \
               == F"I_{''.join(o.l for o in other.l)}_{''.join(b.l for b in other.upl)}_{other.f}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"I_{''.join(o.l for o in self.l)}_{''.join(b.l for b in self.upl)}_{self.f}")

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
        symbl, upl = self.l, self.upl
        orbd = {i:(symbl[i].l,upl[i].l) for i in range(self.norb)}
        pr = list(range(0, self.norb, 2)); qs = list(range(1, self.norb, 2))
        pr = sorted(pr, key=lambda x:orbd[x]); qs = sorted(qs, key=lambda x:orbd[x])
        if [orbd[i] for i in qs]<[orbd[i] for i in pr]: pqrs = list(chain(*zip(qs, pr)))
        else: pqrs = list(chain(*zip(pr, qs)))
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

        return F'{self.l}{self.upl}{self.f}' == F'{other.l}{other.upl}{other.f}'
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




'''
At the level, Orbital < Block < Wave
Orbital, block, wave is Symb, and also possess an implemented Expr in principle

But in fact, the bottom Orbital is only just Symb, has no Expr
        and  the top floor Wave isn't a Symb, but Expr 
'''
# Orbital level
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

        return F'{self.l}{self.f}' == F'{other.l}{other.f}'
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F'{self.l}{self.f}' != F'{other.l}{other.f}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'{self.l}{self.f}' > F'{other.l}{other.f}'
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F'{self.l}{self.f}' >= F'{other.l}{other.f}'
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:
        :return:
        '''

        return F'{self.l}{self.f}' < F'{other.l}{other.f}'
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F'{self.l}{self.f}' <= F'{other.l}{other.f}'
    def __hash__(self):
        '''
        hash: used in set

        Note: hash函数只可定义给不可变对象，而该对象是一个可变对象，定义hash实际上是不合规矩的，这是只是为了方便集合操作
        :return:
        '''

        return hash(F'O{self.l}{self.f}')
# space orbital
nho_ = 4
hop = Orbital(ordinal='p')
hoq = Orbital(ordinal='q')
hor = Orbital(ordinal='r')
hos = Orbital(ordinal='s')
ho_l = (hop, hoq, hor, hos)
# spin orbital
nho = 8
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
nbo = 4  # 4 spin orbitals, nbe is 2 in block
bo0a = Orbital(ordinal='0', spin='\\alpha')
bo0b = Orbital(ordinal='0', spin='\\beta')
bo1a = Orbital(ordinal='1', spin='\\alpha')
bo1b = Orbital(ordinal='1', spin='\\beta')
bol = (bo0a, bo0b, bo1a, bo1b)  # block orbital list
# block level
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
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' != F'B{other.l}{other.f}'
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' > F'B{other.l}{other.f}'
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' >= F'B{other.l}{other.f}'
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:
        :return:
        '''

        return F'B{self.l}{self.f}' < F'B{other.l}{other.f}'
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F'B{self.l}{self.f}' <= F'B{other.l}{other.f}'
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
# any block labeled as PQRS
nhb = 4 # for debug
hbP = Block(ordinal='P')
hbR = Block(ordinal='R')
hbQ = Block(ordinal='Q')
hbS = Block(ordinal='S')
hbl = (hbP, hbQ, hbR, hbS)
ntb = 4
tbA = Block(ordinal='A', state=1)
tbB = Block(ordinal='B', state=1)
tbC = Block(ordinal='C', state=1)
tbD = Block(ordinal='D', state=1)
tbl = (tbA, tbB, tbC, tbD)
# ground block labeled as IJKL
iwb0l = set('IJKL')
wb0l = tuple(Block(ordinal=i) for i in 'IJKL')
# excited block labeled as ABCD
iwb1l = set('ABCD')
wb1l = tuple(Block(ordinal=i, state=1) for i in 'ABCD')
# wave level
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

    def __str__(self): return F"\\Phi_{{{self.f}}}"
    __repr__ = __str__
    def __format__(self, format_spec): return F"({','.join(F'({b.l},{b.f})' for b in self.f)})"

    def __eq__(self, other):
        '''
        self == other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" == F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __ne__(self, other):
        '''
        !=: sefl != other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" != F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __gt__(self, other):
        '''
        >: self > other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" > F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __ge__(self, other):
        '''
        >=: self >= other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" >= F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __lt__(self, other):
        '''
        <: self < other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" < F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __le__(self, other):
        '''
        <=: self <= other ?

        :param other:

        :return:
        '''

        return F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}" <= F"phi_{''.join(F'{b.l}{b.f}' for b in other.f)}"
    def __hash__(self):
        '''
        hash: used in set

        :return:
        '''

        return hash(F"phi_{''.join(F'{b.l}{b.f}' for b in self.f)}")

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




# Symb: 16 block functor
nbf = 16
bfd = {
        0: Functor(symbl=(bo0a, bo0b), pdl='++'),
        1: Functor(symbl=(bo1a, bo1b), pdl='++'),
        2: Functor(symbl=(bo0a, bo1b), pdl='++'),
        3: Functor(symbl=(bo0b, bo1a), pdl='++'),
        4: Functor(symbl=(bo0a, bo1a), pdl='++'),
        5: Functor(symbl=(bo0b, bo1b), pdl='++'),
        6: Functor(symbl=()),
        7:  Functor(symbl=(bo0a,), pdl='+'),
        8:  Functor(symbl=(bo1a,), pdl='+'),
        9:  Functor(symbl=(bo0b,), pdl='+'),
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
# Expr: 16 Functor in block
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
            bfd[7]: (Subbfeq(bf=7,  bs=7),),
            bfd[8]: (Subbfeq(bf=8,  bs=8),),
            bfd[9]: (Subbfeq(bf=9,  bs=9),),
            bfd[10]: (Subbfeq(bf=10, bs=10),),
            bfd[11]: (Subbfeq(bf=11, bs=11),),
            bfd[12]: (Subbfeq(bf=12, bs=12),),
            bfd[13]: (Subbfeq(bf=13, bs=13),),
            bfd[14]: (Subbfeq(bf=14, bs=14),),
            bfd[15]: (Subbfeq(bf=15, bs=15),)
        }
# Expr: 16 states in block
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
        write(F"    A block has 16 states: ${','.join(F'{b!s}' for b in bsd.values())} $ \\\\ \n")
        # block functor
        write( '    Each state corresponds to a functor in essence: \\\\ \n')
        for f in range(nbf):
            write(F'    {f}: ${bfd[f]!s} $ \\\\ \n')
        write( '    \n')

        # block functor of linear combination of state
        write( '    However, a functor is often expressed as a linear combination of a number of States: \\\\ \n')
        for f in range(nbf): write(F"    {f}: ${bfd[f]!s} = {'+'.join(F'{bf!s}' for bf in bfeqd[bfd[f]])} $ \\\\ \n")
        write( '    \n')

        # block state of linear combination of functor
        write( '    Therefore, the composition of each state can be written: \\\\ \n')
        for s in range(nbs):
            write(F"    ${bsd[s]!s} = {'+'.join(F'{bs!s}' for bs in bseqd[s])} \\quad \\vert{{vac}}\\rangle $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# repressive functor
rfeqd = {}
rd = {}
def rf2block(functor, block):
    '''
    repressive functor project into block

    :param functor: Functor class
    :param block: Block class

    :return:
    '''

    forbl = functor.l; nsymb_1 = functor.nsymb-1; pdl = functor.f
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

                m = [[0.0 for i in range(nbs)] for j in range(nbs)]
                cxbll = {j:{} for j in range(nbs)}
                for jbs in range(nbs):
                    cxbll[jbs] = rf2block(rf, bsd[jbs])
                    for ibs in cxbll[jbs]:
                        m[ibs][jbs] = cxbll[jbs][ibs]  # Note: [jbs][ibs] order in cxbll
                        cxbll[jbs][ibs] = Coeff(symbl=(nrf, ibs, jbs), mode='r')
                rd[nrf] = m
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
        for irf,rf in enumerate(rfeqd.keys()):
            write(F'    $\\hat{{O}}_{{{irf}}}:  \\langle{{P_p}}\\vert {rf!s} \\vert{{P_q}}\\rangle => $ \\\\ \n')
            for q in rfeqd[rf]:
                write(F"    $ {rf!s} \\vert{{P_{{{q}}}}}\\rangle = {'+'.join(F'{ce!s}P_{{{p}}}' for p,ce in rfeqd[rf][q].items())} $ \\\\ \n")
                for p in rfeqd[rf][q]:
                    write(F'    ${rfeqd[rf][q][p]!s}\\ =\\ {rd[irf][p][q]!s} $ \\\\ \n')
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
        write('    tmp = {}\n')
        write('    inv = numpy.linalg.inv(ci)\n')

        for irf,rf in enumerate(rfeqd.keys()):
            write('    \n')
            write(F'    # O{irf}:  <Bp|  {rf!s}  |Bq> = \n')
            if any(rfeqd[rf].values()):
                write('    m = numpy.zeros(shape=(nbs,nbs)) \n')
                for q in rfeqd[rf]:
                    for p in rfeqd[rf][q]:
                        write(F'    m[{p},{q}] = {rd[irf][p][q]}\n')
                write(F'    tmp[{irf}] = m\n')
        write('    \n')

        write('    return tmp\n')
        write('    \n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# possible wave state
sl = () # wave state
snd = {b: () for b in range(ntb+1)} # wave state (n-block excite)
def taylor(nt):
    tmp = {}

    for n in range(nt+nhb+1):
        tnl = []; append = tnl.append
        for nb in range(1,n+1):
            for tn in product(range(1,nt+1), repeat=nb):
                if sum(tn)==n: append(tuple(sorted(tn, reverse=True)))
        tmp[n] = tuple(sorted(set(tnl), reverse=True))

    return tmp
tld = taylor(ntb)
def group(l, nl):
    '''
    group elements

    :param l: set
    :param n: int, number of element

    :return:
    '''

    n = nl[0]
    if n is 1: return (tuple(sorted((i,) for i in l)),)
    elif n is len(l): return ((tuple(sorted(l)),),)
    else: return tuple(sorted(set(tuple(sorted((tuple(sorted(coi)),)+coj)) for coi in combinations(l, r=n) for coj in group(l-set(coi), nl[1:]))))
gld = {}
for n in range(ntb+nhb+1):
    for nl in tld[n]:
        gld[nl] = group(set(range(n)), nl)
tsd = {b: set() for b in range(ntb+1)}
def wexcite():
    '''
    wave excite

    Note:
    :return:
    '''

    global sl, snd, tsd
    # gound state
    sl += ((),)
    snd[0] += ((),)
    tsd[0].add(())
    # excited state
    for nb in range(1, ntb+1):
        for sn in product(range(1, nbs), repeat=nb):
            bl = deepcopy(wb1l[:nb])
            ne = 0; sz = 0
            for i in range(nb):
                bl[i].f = sn[i]
                ne += bl[i].ne
                sz += bl[i].sz
            if ne is nb*2 and sz is 0:
                tsd[nb].add(sn)
                sl += (bl,)
                snd[nb] += (bl,)
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
        write(F"    ${'+'.join(F'{len(ul)}C^N_{n}' for n,ul in snd.items())}  $ \\\\ \n")
        write( '    including 1 ground state, and other states are excited one. \\\\ \n')
        write(F'    In particular, For N-block $(N<4 )$ system, the number of wavefunction state can be simply recorded as $ (C^N_{{2N}})^2 $. \\\\ \n')
        write( '    \n')

        write( '    Ground state wavefunction: \\\\ \n')
        write(F'    $\\Phi_{{{sl[0]!s}}}$ \\\\ \n')
        write( '    \n')

        write('    Excited state wavefunction: \\\\ \n')
        for nb in range(1, ntb+1):
            write(F'    For {nb}-block excitation, there are $ {len(snd[nb])}C^N_{nb} $ excited types: \\\\ \n')
            for s in snd[nb]:
                write(F'    $\\Phi_{{{s!s}}}$ \\\\ \n')
            write('    \n')

        write('\\end{document}\n')
        write('\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# T (Taylor) funcotr
tl = () # T types
tnd = {b: () for b in range(1, ntb+1)} # T (n-block excite)
def dett():
    '''
    T operator

    :return:
    '''

    global tl,tnd
    for nb in range(1, ntb+1):
        for sl in product(range(1, nbs), repeat=nb):
            I0l = deepcopy(wb0l[:nb])
            I1l = deepcopy(wb0l[:nb])
            ne = 0; sz = 0
            for i in range(nb):
                I1l[i].f = sl[i]
                ne += I1l[i].ne
                sz += I1l[i].sz
            if ne is nb*2 and sz is 0:
                tl += (I1l,)
                tnd[nb] += (I1l,)
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
        write(F"    $T = {'+'.join(F'T_{i+1}' for i in range(ntb))} $ \\\\ \n")
        write( '    where: \\\\ \n')
        for i in tnd:
            bl = tuple(b.l for b in tnd[i][0]); bl_ = tuple(b.lower() for b in bl)
            sigmal = tuple("\\displaystyle\\sum_{{{}}}\\displaystyle\\sum_{}".format('>'.join(reversed(bl[:j+1])), bl_[j]) for j in range(len(bl)))
            operl = tuple(F'{bl[j]}_{bl_[j]}^+{bl[j]}_0^-' for j in range(len(bl)))
            write(F"    $T_{i} = {''.join(sigmal)}{{{''.join(operl)}}} $ \\\\ \n")
            for bl in tnd[i]:
                write(F"    ${','.join(F'{b}={bl[ib].f}' for ib,b in enumerate(bl_))} $ \\\\ \n")
        write( '    \n')

        write( '\\end{document}\n')
        write( '\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()

# matrix elem for T
tfeqd = {}
def t2wave(ws):
    '''
    T project into wave

    :param ws: wave state, list

    :return:
    '''

    raise NotImplementedError('TODO: t2wave() for T project into wave')
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

    for nbbra in range(1, ntb+1):
        for bra in snd[nbbra]:
            tmp = {}
            bbl = set(bra)
            for nbt in range(1, nbbra+1):
                for t in combinations(bra, nbt):
                    if sum(b.ne for b in t) is 2*len(t) and sum(b.sz for b in t) is 0:
                        ket = tuple(sorted(bbl - set(t)))
                        s = sign(tuple(b.l for b in t+ket if b.ne%2))
                        tmp[t] = (s, ket)
            tfeqd[tuple(bra)] = tmp
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

        for bra in tfeqd:
            write(F'    $\\Phi_{{{bra!s}}} = $ \\\\ \n')
            for t in tfeqd[bra]:
                if tfeqd[bra][t][0] is 1: write(F"    $+ {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{tfeqd[bra][t][1]}}} $ \\\\ \n")
                else: write(F"    $- {''.join(F'{b!s}^+' for b in t)} \\Phi_{{{tfeqd[bra][t][1]}}} $ \\\\ \n")
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
        for bra in tfeqd:
            bl = tuple(b.l for b in bra)
            try: bld[bl].append(bra)
            except KeyError: bld[bl] = [bra]

        write('def tm(np, nt=4, nc=0):\n')
        write('    tmp = {}\n')
        write('    \n')

        tab = '    '
        for bl in bld:
            nb = 0
            for b in bl:
                if nb>0:  write(tab*nb+F'    for {b} in range({bl[nb-1]}+1, np):\n')
                else: write(F'    for {b} in range(np):\n')
                nb += 1

            for bra in bld[bl]:
                kbra = ''.join(F'({b.l},{b.f}),' for b in bra)
                write(tab*nb+ '    _t = {}\n')
                for t in tfeqd[bra]:
                    kt = ''.join(F'({b.l},{b.f}),' for b in t)
                    kket = ''.join(F'({b.l},{b.f}),' for b in tfeqd[bra][t][1])
                    write(tab*nb+F'    _t[({kt})] = ({tfeqd[bra][t][0]}, ({kket}))\n')
                write(tab*nb+F'    tmp[({kbra})] = _t\n')
                write(tab*nb+ '    \n')

            if nb<4: write(F'    if nt is {nb}: return tmp\n')
            else: write( '    return tmp\n')
            write( '    \n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()


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

    fpre = F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}"
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
        for tl in tld[nb][1:]:
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
def atu2code(bra, td):
    '''
    array of T element u -> code

    Note:     subfuction

    :return:
    '''

    # cross platform
    osname = os.name
    ossep = os.sep
    oslinesep = os.linesep

    fpre = F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}"
    with open(F'bccc{ossep}ta1{ossep}{fpre}.py', 'w') as fout:
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

        kbd = {}
        for ket in td:
            kb = tuple(b.l for b in ket if b.l in iwb0l)
            try: kbd[kb].append(ket)
            except KeyError: kbd[kb] = [ket]

        if bra:
            write(F'def {fpre}(t, h, p, E):\n') # iterate amplitude
            write('    t1 = {}\n')
        else: write(F'def {fpre}(t, h, p):\n') # calculate Energy
        write( '    np = len(p)\n')
        write( '    \n')

        tab = '    '

        bl,bsl = tuple(b.l for b in bra),tuple(b.f for b in bra)
        nb = len(bsl); kbra = F"({''.join([F'({b.l},{b.f}),' for b in bra])})"
        tbra_ = []; append = tbra_.append
        for tl in tld[nb][1:]:
            for gl in gld[tl]:
                if all(tuple(bsl[i] for i in g) in tsd[len(g)] for g in gl):
                    s = sign([bl[i] for i in chain(*gl) if 6<bsl[i]<15])
                    if s is 1: append(reduce(mul, tuple(Coeff(symbl=tuple(bra[i] for i in g), mode='t') for g in gl)))
                    else: append(0.0-reduce(mul, tuple(Coeff(symbl=tuple(bra[i] for i in g), mode='t') for g in gl)))
                    # append(Coeff(symbl=tuple(tuple(bra[i] for i in g) for g in gl), mode='t'))
        nbbra = 0
        for b in bl:
            if nbbra>0:  write(tab*nbbra+F'    for {b} in range({bl[nbbra-1]}+1, np):\n')
            else: write(F'    for {b} in range(np):\n')
            nbbra += 1
        write(tab*nbbra+ '    _t = 0.0\n')
        write(tab*nbbra+ '    \n')

        for kb in kbd:
            il = tuple(b for b in kb if b in iwb0l)

            nbket = 0
            for b in il:
                if nbket>0: write(tab*(nbbra+nbket)+F'    for {b} in range({il[nbket-1]}+1, np):\n')
                else: write(tab*(nbbra+nbket)+F'    for {b} in range(np):\n')
                if bl: write(tab*(nbbra+nbket)+F"        if {b} in {{{', '.join(bl)}}}: continue\n")
                nbket += 1

            nb = nbbra+nbket
            for ket in kbd[kb]:
                # if ib1l&iwb1l and ib1l&iwb0l: write(tab*nb+F'    key = tuple(sorted(({key})))\n')
                # else: write(tab*nb+F'    key = ({key})\n')
                # write(tab*nb+F"    _t += {td[ket]*Coeff(symbl=(bra, ket), mode='h')} \n")
                txtl = td[ket]; ntxt = len(txtl)
                if ntxt>1:
                    write(tab*nb+'    tmp = (    \n')
                    for txt in txtl: write(tab*nb+F'          +{txt}    \n')
                    write(tab*nb+'           )\n')
                    write(tab*nb+F"    _t += tmp*{Coeff(symbl=(bra, ket), mode='h')} \n")
                else: write(tab*nb+F"    _t += {td[ket][0]*Coeff(symbl=(bra, ket), mode='h')} \n")
            write(tab*nb+ '    \n')

        if bra: write(tab*nbbra+F"    t1[{kbra}] = (_t+h[()].get({kbra},0))/(E-h[{kbra}][{kbra}])-({sum(tbra_)})\n")
        write(tab*nbbra+ '    \n')

        if bra: write(F"    return t1\n")
        else: write(F"    return _t\n")
        write('    \n\n\n\n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n\n\n')

        # # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()
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

    fpre = 'ta'
    with open(F'bccc{ossep}ta1{ossep}{fpre}.py', 'w') as fout:
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
        write('from copy import deepcopy\n')
        write('from datetime import datetime\n')
        write('from importlib import import_module\n')
        write('import sys, gc\n')
        write('\n\n')

        write('def calculate(braf, t0, h, p, E0):\n')
        write("    module = F'bccc.ta1.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    t1 = fun(t0, h, p, E0)\n')
        write('    del sys.modules[module]\n')
        write('    gc.collect()\n')
        write('    return t1\n')
        write('\n\n')

        write('def ta(t0, h, p, nt=4, nc=0):\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        write('    from bccc.ta1.ta_ import ta_\n')
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[1])
        for braf in brafl: write(F'    from bccc.ta1.{braf} import {braf}\n')
        write('    brafl = (\n             ')
        for i,braf in enumerate(brafl):
            write(F'{braf}, ')
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for nt in range(2,ntb+1):
            write(F'    # {nt}-block excited types\n')
            write(F'    if nt>={nt}: \n')
            brafl = tuple(F"ta_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[nt])
            for braf in brafl: write(F'        from bccc.ta1.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F'{braf}, ')
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write("    t,Eref = {},h[()][()]\n")
        write('    Ecorr=ta_(t0, h, p); E0=Eref+Ecorr\n')
        write("    print(F'Ecorr,E0 = {Ecorr},{E0}')\n")
        write('    \n')

        write("    print('    it            dE               dt               time/s')\n")
        write('    for it in range(200):\n')
        write('        t1 = datetime.now()\n')
        write('        \n')

        write('        pool = Pool()\n')
        write('        vl = tuple(pool.apply_async(braf, (t0, h, p, E0)) for braf in brafl)\n')
        write('        pool.close()\n')
        write('        pool.join()\n')
        write('        for v in vl: t.update(v.get())\n')
        # write('        for braf in brafl: t.update(braf(t0, h, p, E0))\n')
        write('        \n')

        write('        t2 = datetime.now()\n')
        write('        \n')

        write('        Ecorr=ta_(t, h, p); E=Eref+Ecorr\n')
        write('        dE=E-E0; dt=sum(abs(t[u]-t0[u]) for u in t)\n')
        write("        print(F'    {it+1:3d}      {dE:13.10f}    {dt:13.10f}     {t2-t1}', flush=True)\n")
        write('        if dt<1.0e-7: \n')
        write("            print(F'Successfully converged: Ecorr = {Ecorr}')\n")
        write('            break\n')
        write('        else: \n')
        write('            t0,E0 = deepcopy(t),E\n')
        write('        \n')
        write('    else:\n')
        write("        print('Error: Convergence failed ')\n")
        write('    \n')

        write('    return Ecorr,t\n')
        write('    \n\n\n\n')

        # write('\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        # write("if __name__ == '__main__':\n")
        # write('    \n\n')
        # write("    print('End successfully')\n")
        # write('    \n\n')

        # sys.stdout = stdout; sys.stderr = stderr
        # fout.close()




# H (Hamiltonian) functor
hd = {} # H types
hnd = {} # H (one- or two-electron)
def deth():
    raise NotImplementedError('TODO: deth()')
def h2theory():
    raise NotImplementedError('TODO: h2theory()')

# matrix elem for hamiltion
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
    with open(F'bccc{ossep}hm1{ossep}{fpre}.py', 'w') as fout:
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
            ketf = F"{sket}({', '.join(ibl)}, r, h_mo, g_mo, p)"

            key = tuple(filter(lambda lf: lf[1], lfl))
            ib1l = tuple(lf[0] for lf in key)

            try: ibld[ibl][ket] = (ib1l, key, ketf)
            except KeyError: ibld[ibl] = {ket: (ib1l, key, ketf)}

            write(F"def {ketf}:\n")
            write(F'    return ({ced[ket]})\n')
            write( '    \n')
        write('\n\n\n\n')


        write(F'def {fpre}(r, h_mo, g_mo, p):\n')
        write( '    np = len(p)\n')
        write( '    tmp = {}\n')
        write( '    \n')

        tab = '    '

        ibAl = tuple(b.l for b in bra)
        nbbra = 0
        for b in ibAl:
            if nbbra>0:  write(tab*nbbra+F'    for {b} in range({ibAl[nbbra-1]}+1, np):\n')
            else: write(F'    for {b} in range(np):\n')
            nbbra += 1
        write(tab*nbbra+F"    bbl = {{{','.join(ibAl)}}}\n")
        write(tab*nbbra+ '    _t = {}\n')
        write(tab*nbbra+ '    \n')

        for ibl in ibld:
            ibIl = tuple(b for b in ibl if b in iwb0l)

            nbket = 0
            for b in ibIl:
                if nbket>0: write(tab*(nbbra+nbket)+F'    for {b} in range({ibIl[nbket-1]}+1, np):\n')
                else: write(tab*(nbbra+nbket)+F'    for {b} in range(np):\n')
                if ibAl: write(tab*(nbbra+nbket)+F"        if {b} in bbl: continue\n")
                nbket += 1

            nb = nbbra+nbket
            for ket in ibld[ibl]:
                _ib1l, _key, ketf = ibld[ibl][ket]
                key = ''.join('({}, {}),'.format(*lf) for lf in _key)
                ib1l = set(_ib1l)

                if ib1l&iwb1l and ib1l&iwb0l: write(tab*nb+F'    key = tuple(sorted(({key})))\n')
                else: write(tab*nb+F'    key = ({key})\n')
                write(tab*nb+F'    try: _t[key] += {ketf}\n')
                write(tab*nb+F'    except KeyError: _t[key] = {ketf}\n')
                write(tab*nb+ '    \n')

        write(tab*nbbra+F"    tmp[({''.join([F'({b.l},{b.f}),' for b in bra])})] = _t\n")
        write(tab*nbbra+ '    \n')

        write('    return tmp\n')
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
def h2wave(ws):
    '''
    hamiltonnian project into wave

    :param ws: wave state, list

    :return:
    '''

    ced = {}; keyd = {} # {u: {<>: ce}} {u: u}
    nb = len(ws); ne = sum(b.ne for b in ws)
    ibl = tuple(b.l for b in ws)
    spinl = ('\\alpha', '\\beta')

    # for h
    if nb<=2 and ne>=1:  # nb<2 for h; p+p- ne>=1
        for prodb in product(ws, repeat=2):
            if len(set(prodb)) is not nb: continue
            _ibl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    inte = Integral(symbl=(bop, boq), mode=0)
                    functor = Functor(symbl=(bop, boq), pdl='+-', upl=prodb)
                    s_blockize,functorl,upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                    _ced={}; _keyd={} # ce {u: ce} 没有与 inte 相乘
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
                                if key in ced: key = keyd[key] # 取 key of ced, 如果key已经有了
                                elif key in _ced: key = _keyd[key] # 取 key of _ced, 如果key已经有了
                                ce = 1.0
                                for b,_ce in zip(key, cel):
                                    _ce = deepcopy(_ce)
                                    _ce.l = (b,)+_ce.l
                                    ce = ce*_ce
                                try: _ced[key] = _ced[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                except KeyError: _keyd[key],_ced[key] = key,ce

                    for key in _ced:  # 将 _ced 更新到 ced 中
                        ce = _ced[key]
                        _inte = deepcopy(inte)
                        _inte.upl = tuple(key[index(b, ibl)] for b in _ibl)
                        _inte.orbsort()
                        if s_blockize is -1: ce = 0.0-ce
                        if key in ced:
                            try: ced[key][_inte] = ced[key][_inte]+ce
                            except KeyError: ced[key][_inte] = ce
                        else: keyd[key],ced[key] = key,{_inte: ce}

    # for v
    if ne>=2: # p+p+s-r- ne>=2
        for prodb in product(ws, repeat=4):
            if len(set(prodb)) is not nb: continue
            _ibl = tuple(b.l for b in prodb)
            for bop in deepcopy(bo_l):
                for boq in deepcopy(bo_l):
                    for bor in deepcopy(bo_l):
                        for bos in deepcopy(bo_l):
                            inte = Integral(symbl=(bop, boq, bor, bos), mode=2)
                            functor = Functor((bop, boq, bos, bor), pdl='++--', upl=(prodb[0], prodb[1], prodb[3], prodb[2]))
                            s_blockize, functorl, upl = functor.blockize()  # 算符分块化 sign，f1b1,f2b2

                            _ced = {}; _keyd={}  # ce 没有与 inte 相乘
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
                                        if key in ced: key = keyd[key] # 取 key of ced, 如果key已经有了
                                        elif key in _ced: key = _keyd[key] # 取 key of _ced, 如果key已经有了
                                        ce = 1.0
                                        for b,_ce in zip(key, cel):
                                            _ce = deepcopy(_ce)
                                            _ce.l = (b,)+_ce.l
                                            ce = ce*_ce
                                        try: _ced[key] = _ced[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError:  _keyd[key],_ced[key] = key,ce
                            for key in _ced: _ced[key] = 0.5*_ced[key]
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
                                        if key in ced: key = keyd[key] # 取 key of ced, 如果key已经有了
                                        elif key in _ced: key = _keyd[key] # 取 key of _ced, 如果key已经有了
                                        ce = 1.0
                                        for b,_ce in zip(key, cel):
                                            _ce = deepcopy(_ce)
                                            _ce.l = (b,)+_ce.l
                                            ce = ce*_ce
                                        try: _ced[key] = _ced[key]+ce  # spin 内合并 或 alpha与beta 的合并
                                        except KeyError: _keyd[key],_ced[key] = key,ce

                            for key in _ced:  # 将 _ced 更新到 ced 中
                                ce = _ced[key]
                                _inte = deepcopy(inte)
                                _inte.upl = tuple(key[index(b, ibl)] for b in _ibl)
                                _inte.orbsort()
                                if s_blockize is -1: ce = 0.0-ce
                                if key in ced:
                                    try: ced[key][_inte] = ced[key][_inte]+ce
                                    except KeyError: ced[key][_inte] = ce
                                else: keyd[key],ced[key] = key,{_inte: ce}

    # ced {u: ce}
    for key in ced:
        ced[key] = sum(ce*inte for inte,ce in ced[key].items())

    return ced
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
    ced = {} # {ket: ce}
    for nbcorr in range(1, nhb+1):  # n blocks correlation
        for nbexcite in range(nbcorr+1): # nbexcite blocks from bra
            # print('nbexcite = ', nbexcite)
            for bchose in combinations(bra, r=nbexcite): # 从激发态取 nbexcite
                # print('bchose = ', bchose)
                brest = tuple(sorted(set(bra)-set(bchose))) # set -> list is random and disordered
                # print('brest = ', brest)
                s_bef = sign([b.l for b in bchose+brest if b.ne%2]) # 激发块提到前面
                # print('s_bef = ', s_bef)
                _ws = bchose+wb0l[:nbcorr-nbexcite]
                # if len(_ws+brest)>ntb: continue # for debug
                # print('                             _ws = ', _ws)

                key = tuple(b.f for b in _ws); iws = tuple(b.l for b in _ws)
                if key in hfeqd:
                    _ced = deepcopy(hfeqd[key])
                    try:
                        bl = tuple(b.l for b in list(_ced.keys())[0])
                        if bl != iws:
                            for bl in _ced:
                                for b, i in zip(bl, iws): b.l = i
                    except IndexError: continue   # _ced = {} when nbcorr=nbexcite, b.f=6
                else:
                    _ced = h2wave(_ws)  # Aa I0 J0 Kk
                    hfeqd[key] = _ced

                for key, ce in _ced.items(): # 将 _ced 更新(多块加和)到 ced 中
                    ket = key+brest

                    leaf = [b.l for b in ket if b.ne%2]; ib1l = [b.l for b in ket if b.ne%2 and b.l in iwb1l]
                    _leaf = set(leaf); _ib1l = sorted(ib1l)
                    if _leaf&iwb1l and (_leaf&iwb0l or _ib1l!=ib1l): # 还原激发块
                        s_lat = Expr('-?', leaf)
                        ce = s_lat*ce
                    if s_bef is -1: ce = 0.0-ce

                    ket = tuple(sorted(ket))
                    try: ced[ket] = ced[ket]+ce
                    except KeyError: ced[ket] = ce

    # at
    td = {}  # {ket: t}
    ketl = set(); append = ketl.add
    for _ket in ced.keys(): # Aa..Jj.. -> Aa..Ij..
        bbl = tuple(b for b in _ket if b.l in iwb1l and b.f>0); kbl = deepcopy(tuple(b for b in _ket if b.l in iwb0l and b.f>0))
        for b,l in zip(kbl, 'IJKL'): b.l = l
        append(bbl+kbl)
    ketl.remove(bra)
    for ket in ketl:
        bl, bsl = tuple(b.l for b in ket), tuple(b.f for b in ket); nb = len(bsl)

        tket_ = []; append = tket_.append # 有效的分组
        for tl in tld[nb]:
            for gl in gld[tl]:
                if all(tuple(bsl[i] for i in g) in tsd[len(g)] for g in gl): append(gl)
        txtl = []; append=txtl.append
        for gl in tket_:
            oddbl = [bl[i] for i in chain(*gl) if 6<bsl[i]<15]; _oddbl = set(oddbl)
            if _oddbl&iwb0l and _oddbl&iwb1l:
                oddbl = []
                for g in gl:
                    _t = tuple(bl[i] for i in g if 6<bsl[i]<15)
                    if _t: oddbl.append(_t)
                s = Expr('-?', oddbl)
                append(s*reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))
            else:
                s = sign(oddbl)
                if s is -1: append(0.0-reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))
                else: append(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)))

        if txtl:
            # td[ket] = sum(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)) for gl in tket_)
            # td[ket] = sum(reduce(mul, tuple(Coeff(symbl=tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g), mode='t') for g in gl)) for gl in tket_)
            # td[ket] = sum(txtl)
            td[ket] = txtl
            # td[ket] = tuple(tuple(tuple(Block(ordinal=bl[i], state=bsl[i]) for i in g) for g in gl) for gl in tket_)

    # atu2theory(bra, td)
    # atu2code(bra, td)
    # mhuv2theory(bra, ced)
    mhuv2code(bra, ced)
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
        write(F'    In total, there are 1 ground state and {len(sl)-1} excited types. \\\\ \n')
        write( '    \n')

        write( '    When $\\hat{H}$ operator project onto a wavefunction state, it will be converted to a linear combination of other states. that is \\\\ \n')
        write( '    $\\hat{H}\\vert{\\Phi_{u}}\\rangle = \\displaystyle\\sum_{u,v}{c_{u,v}\\Phi_{v}}$ \\\\ \n')
        write( '    where $c_{u,v} $ is Hamiltonian matrix, and listed separately in other files.\\\\ \n')
        write( '    \n')

        write( '    Specifically, state $ u, v = $ \\\\ \n')
        # write(F"    ${','.join(F'{u}' for i,u in sl)} $ \\\\ \n")
        for u in sl: write(F'    ${u} $ \\\\ \n')
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
    with open(F'bccc{ossep}hm1{ossep}{fpre}.py', 'w') as fout:
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

        write('def calculate(braf, r, h_mo, g_mo, p):\n')
        write("    module = F'bccc.hm1.{braf}'\n")
        write('    lib = import_module(module)\n')
        write('    fun = getattr(lib, braf)\n')
        write('    tmp = fun(r, h_mo, g_mo, p)\n')
        write('    del sys.modules[module]\n')
        write('    gc.collect()\n')
        write('    return tmp\n')
        write('\n\n')

        write('def hm(r, h_mo, g_mo, p, nt=4, nc=0):\n')
        write('    h = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        write('    from bccc.hm1.hm_ import hm_\n')
        write("    brafl = (hm_,)\n")
        write('    \n')

        write('    # 1-block excited types\n')
        brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[1])
        for braf in brafl: write(F'    from bccc.hm1.{braf} import {braf}\n')
        write('    brafl += (\n             ')
        for i,braf in enumerate(brafl):
            write(F"{braf}, ")
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for nt in range(2,ntb+1):
            write(F'    # {nt}-block excited types\n')
            write(F'    if nt>={nt}: \n')
            brafl = tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[nt])
            for braf in brafl: write(F'        from bccc.hm1.{braf} import {braf}\n')
            write( '        brafl += (\n                 ')
            for i,braf in enumerate(brafl):
                write(F"{braf}, ")
                if (i+1)%5 is 0: write('\n                 ')
            write('                 \n')
            write('                 )\n')
            write('        \n')

        write('    pool = Pool()\n')
        write('    vl = tuple(pool.apply_async(braf, (r, h_mo, g_mo, p)) for braf in brafl)\n')
        write('    pool.close()\n')
        write('    pool.join()\n')
        write('    for v in vl: h.update(v.get())\n')
        write('    \n')
        write('    return h\n')
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
    with open(F'hm{ossep}{fpre}.py', 'w') as fout:
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

        write('def corr(hm, nt=4):\n')
        write('    tmp = {}\n')
        write('    \n')

        write('    # bra wavefunction list\n')
        write('    # ground state\n')
        brafl = ('hm_',)
        write('    from bccc.hm.hm_ import hm_\n')
        write('    \n')

        write('    # 1-block excited types\n')
        brafl += tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[1])
        for braf in brafl: write(F'    from bccc.hm.{braf} import {braf}\n')
        write('    brafl = (\n             ')
        for i,braf in enumerate(brafl):
            write(F'{braf}, ')
            if (i+1)%5 is 0: write('\n             ')
        write('             \n')
        write('             )\n')
        write('    \n')

        for nt in range(2,ntb+1):
            write(F'    if nt>{nt-1}:')
            write(F'        # {nt}-block excited types\n')
            brafl += tuple(F"hm_{''.join(F'{b.l}{b.f}' for b in bra)}" for bra in snd[nt])
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
        # write('    for i,braf in enumerate(brafl):  vl[i] = pool.apply_async(braf, (r, h_mo, g_mo, p))\n')
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
    for u in sl: mp.apply_async(detmh, (u, hfeqd))
    mp.close()
    mp.join()






if __name__ == '__main__':
    # base2theory()

    detrfeq()
    # rfeq2theory()
    # rfeq2code()

    wexcite()
    # w2theory()

    # dett()
    # t2theory()
    # detmt()
    # mt2theory()
    # mt2code()

    t1 = datetime.now()
    if sys.platform=='linux':
        parallel()  # for parallel
        # Pool().map(detmh, ((u, hfeqd) for u in sl)) # map中的函数只能接受一个参数, 即detmh(bh)，所以必须把detmh的变量包装
    else:
        for ws in sl: detmh(ws, hfeqd)
    t2 = datetime.now()
    print(F'derivate formula {t2-t1}')
    # at2theory()
    # at2code()
    # mh2theory()
    # mh2code()



    print('End successfully')



