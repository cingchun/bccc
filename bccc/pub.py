#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   pubfun.py
#            Des:   public function for BCCC package are written in  this module
#           Mail:   qingchun720@foxmail.com
#   Created Time:   16:58 二月-25/2019 
#


__version__ = '3.1.2'
__author__ = 'Qingchun Wang'


from itertools import combinations, combinations_with_replacement, product, permutations, chain
from copy import deepcopy
from functools import reduce
from operator import mul


''' globle variables '''
lT0 = None
def lowTri0(no):
    lT0 = [0]*no
    
    for i in range(1, no):
        lT0[i] = lT0[i-1]+i
        
    return lT0
def dei(i, j, k, l):
    '''
    Index of Lower Triangular

    :param i:
    :param j:
    :param k:
    :param l:

    :return:
    '''

    ij = lT0[j]+i if i<=j else lT0[i]+j
    kl = lT0[l]+k if k<=l else lT0[k]+l

    return ij, kl
def phiu(u, pul):
    '''
    wave state u -> integer
    
    :param u:
    
    :return:
    '''

    return sum(pul[i]*(b*16+s) for i,(b, s) in enumerate(u))


# public funciton for derivate formula
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
    else:
        return tuple(sorted(set(tuple(sorted((tuple(sorted(coi)),)+coj)) for coi in combinations(l, r=n) for coj in group(l-set(coi), nl[1:]))))

def iis(elem, l):
    '''
    index of elem/list is elem/list of list(l)

    :param elem:
    :param l:

    :return:
    '''

    for i, e in enumerate(l):
        if e==elem: return i
    # else: return -1
def iin(elem, ll):
    '''
    index of elem/list in list of list(ll)
    
    :param elem:
    :param ll:
    :return:
    '''
    
    for i,l in enumerate(ll):
        if elem in l: return i
    # else: return None
def sign(lnew, lold=None):
    '''
    sign of swap(lold -> lnew)

    :param lnew:
    :param lold:

    :return:
    '''

    if len(lnew)<2: return 1 # int, for is 1
    else:
        if lold: il = tuple(lold.index(inew) for inew in lnew)
        else: il = lnew
        # il=[]
        # for new in lnew:
        #     for iold, old in enumerate(lold):
        #         if index(iold, il): continue
        #         else:
        #             if new==old:
        #                 il.append(iold)
        #                 break

        n = 1
        for i, e in enumerate(il):
            for j in range(i):
                if il[j] > e: n += 1

        if n%2: return 1
        else: return -1

def suml(l):
    '''
    sum(l)

    :param l:
    :return:
    '''

    return sum(l)
def prodl(l):
    '''
    prod(l)

    :param l:
    :return:
    '''

    return reduce(mul, l)
def sum_list(l):
    '''
    sum(list)

    :param l:
    :return:
    '''

    return list(chain.from_iterable(l))
def prod_list(*l, **kwargs):
    '''
    pord(list): Cartesian products of list

                e.g product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
                    product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    :param l:
    :param repeat:

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in product(*l, **repeat)]
    return [list(i) for i in product(*l, **kwargs)]
# Permutation and combination
def perm(l, m):
    '''
    Permutations of n elements in l, P(n,m)

    :param l: n=len(l)
    :param m: m<n

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in permutations(l, r=m)]
    return [list(i) for i in permutations(l, r=m)]
def perm_replace(l, m):
    '''
    Permutations with replacement of n elements in l, P(n,m) = n**m

    :param l: n=len(l)
    :param m:

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in product(l, repeat=m)]
    return [list(i) for i in product(l, repeat=m)]
def perm_replace_(l, m):
    '''
    Cartesian products or Permutations with replacement of n elements in l, P(n,m)

    :param l: m=len(l)
    :param m:

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in product(l, repeat=m) if len(set(ll))==len(l)]
    return [list(i) for i in product(l, repeat=m) if len(set(i))==len(l)]
def comb(l, m):
    '''
    Combinations of n elements in l, C(n,m)

    :param l: n=len(l)
    :param m: m<n

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in combinations(l, r=m)]
    return [list(i) for i in combinations(l, r=m)]
def comb_replace(l, m):
    '''
    Combinations with replacement of n elements in l, C(n,m)

    :param l: n=len(l)
    :param m:

    :return:
    '''

    # return [[copy.deepcopy(e) for e in ll] for ll in combinations_with_replacement(l, r=m)]
    return [list(i) for i in combinations_with_replacement(l, r=m)]


# calculate matrix of Hamiltion
def backblock(bl, nel):
    '''
    back block

        Note: default [] has bug in Python, Not recommended
              But unmodify [], also acceptable.
    :param bl: block list
    :param nel: number of electorn list

    :return:
    '''

    if len(bl)<2: return 1
    else:
        tmp = 1
        for ib in range(len(bl)):
            ne = 0
            for jb in range(ib+1, len(bl)):
                if bl[ib] > bl[jb]:
                    ne += nel[jb]
            tmp += nel[ib]*ne

        if tmp%2: return 1
        else: return -1
def sortblock(ll):
    '''
    list of tuple (P, p)

    :param ll:

    :return:
    '''

    tmp = []
    for l in ll:
        if isinstance(l, tuple):
            for i,_l in enumerate(tmp):
                if l[0]<_l[0]:
                    tmp.insert(i, l)
                    break
            else: tmp.append(l)
        else: tmp.append(l)

    return tuple(tmp)


# calculate energy
def bcpt2(epair, h0d):
    '''
    BCPT2 energy

    :param epair:
    :param mh:

    :return:
    '''

    # tmp = 0.0
    # for ket in list(mh[()].keys())[1:]:
    #     # # print('ket=', ket)
    #     # if ket[0]==() and ket[1]!=():
    #     #     # print('ket=', ket[1])
    #     #     # print('ket = ', ket)
    #     #     Ediff = 0.0
    #     #     if isinstance(ket[1][0], tuple):
    #     #         # print('ket=', ket[1])
    #     #         for Pp in ket[1]:
    #     #             Ediff += epair[Pp[0]][0]-epair[Pp[0]][Pp[1]]
    #     #     else:
    #     #         Ediff += epair[ket[1][0]][0]-epair[ket[1][0]][ket[1][1]]
    #     #
    #     #     # print('mh[ket]=', mh[ket])
    #     #     tmp += mh[ket]**2/Ediff
    #     # diff = 0.0
    #     # if isinstance(ket[0], tuple):
    #     #     for ws in ket:
    #     #         diff += epair[ws[0]][0]-epair[ws[0]][ws[1]]
    #     # else: diff = epair[ket[0]][0]-epair[ket[0]][ket[1]]
    #     diff = sum(epair[b][0]-epair[b][s] for b,s in ket)
    #
    #     if abs(diff) < 1.0e-10:
    #         print(F'ket = {ket}')
    #         raise ZeroDivisionError('Denominator is 0 in BCPT2 calculation')
    #     tmp += mh[()][ket]**2/diff
    #
    # return tmp
    return sum(h**2/sum(epair[b][0]-epair[b][s] for b,s in u) for u,h in list(h0d.items())[1:])
def bclcc(npair, mh, nt=4, nc=0):
    from bccc.tm import tm
    mt = tm(npair, nt, nc)
    id = {u: i for i,u in enumerate(mt.keys())}; n = len(id)
    import numpy
    mT = numpy.zeros(shape=(n, n)); b = numpy.zeros(n)
    print('mT.shape = ', mT.shape)
    for u,i in id.items():
        b[i] = -mh[()].get(u, 0)
        for v,j in id.items():
            mT[i,j] = mh[u].get(v,0)
            if v in mt[u]:
                s,k = mt[u][v]
                if s is 1: mT[i,j] -= mh[()].get(k,0)
                else: mT[i][j] += mh[()].get(k,0)
    t = numpy.linalg.solve(mT, b)
    # return sum(t[i]*mh[()].get(u,0.0) for u,i in id.items())
    return sum(t[i]*mh[()].get(u,0.0) for u,i in id.items()), {u:t[i] for u,i in id.items()}






if __name__ == '__main__':


    print('End successfully')

