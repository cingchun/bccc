#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
# 
# Author: Qingchun Wang @ NJU 
# E-mail: qingchun720@foxmail.com 
# 


from copy import deepcopy
from datetime import datetime
from numpy import zeros,ones,array,matrix,ctypeslib,linalg,uint32,float64
from ctypes import c_uint


def iter(t, r, h, g, p, nt=4, nc=0, thd=1.0e-7, ncyc=100, xE=0.015):
    # fun list
    # ground state
    # from bccc.ta.ta_ import ta_ # Ecorr
    # ftl = (ta_,)
    ftl = ()
    from bccc.ta.ta_ import h0_
    f0l = (h0_,)
    ful = ()
    # 1-block excited types
    from bccc.ta.ta_A1 import ta_A1
    from bccc.ta.ta_A2 import ta_A2
    from bccc.ta.ta_A3 import ta_A3
    ftl += (
            ta_A1,ta_A2,ta_A3,
            )
    from bccc.ta.ta_A1 import h0_A1
    from bccc.ta.ta_A2 import h0_A2
    from bccc.ta.ta_A3 import h0_A3
    f0l += (
            h0_A1,h0_A2,h0_A3,
            )
    from bccc.ta.ta_A1 import hu_A1
    from bccc.ta.ta_A2 import hu_A2
    from bccc.ta.ta_A3 import hu_A3
    ful += (
            hu_A1,hu_A2,hu_A3,
            )
    # 2-block excited types
    from bccc.ta.ta_A1B1 import ta_A1B1
    from bccc.ta.ta_A1B2 import ta_A1B2
    from bccc.ta.ta_A1B3 import ta_A1B3
    from bccc.ta.ta_A2B2 import ta_A2B2
    from bccc.ta.ta_A2B3 import ta_A2B3
    from bccc.ta.ta_A3B3 import ta_A3B3
    from bccc.ta.ta_A4B5 import ta_A4B5
    from bccc.ta.ta_A6B15 import ta_A6B15
    from bccc.ta.ta_A7B13 import ta_A7B13
    from bccc.ta.ta_A7B14 import ta_A7B14
    from bccc.ta.ta_A8B13 import ta_A8B13
    from bccc.ta.ta_A8B14 import ta_A8B14
    from bccc.ta.ta_A9B11 import ta_A9B11
    from bccc.ta.ta_A9B12 import ta_A9B12
    from bccc.ta.ta_A10B11 import ta_A10B11
    from bccc.ta.ta_A10B12 import ta_A10B12
    ftl += (
            ta_A1B1,ta_A1B2,ta_A1B3,ta_A2B2,ta_A2B3,ta_A3B3,ta_A4B5,ta_A6B15,ta_A7B13,ta_A7B14,
            ta_A8B13,ta_A8B14,ta_A9B11,ta_A9B12,ta_A10B11,ta_A10B12,
            )
    from bccc.ta.ta_A1B1 import h0_A1B1
    from bccc.ta.ta_A1B2 import h0_A1B2
    from bccc.ta.ta_A1B3 import h0_A1B3
    from bccc.ta.ta_A2B2 import h0_A2B2
    from bccc.ta.ta_A2B3 import h0_A2B3
    from bccc.ta.ta_A3B3 import h0_A3B3
    from bccc.ta.ta_A4B5 import h0_A4B5
    from bccc.ta.ta_A6B15 import h0_A6B15
    from bccc.ta.ta_A7B13 import h0_A7B13
    from bccc.ta.ta_A7B14 import h0_A7B14
    from bccc.ta.ta_A8B13 import h0_A8B13
    from bccc.ta.ta_A8B14 import h0_A8B14
    from bccc.ta.ta_A9B11 import h0_A9B11
    from bccc.ta.ta_A9B12 import h0_A9B12
    from bccc.ta.ta_A10B11 import h0_A10B11
    from bccc.ta.ta_A10B12 import h0_A10B12
    f0l += (
            h0_A1B1,h0_A1B2,h0_A1B3,h0_A2B2,h0_A2B3,h0_A3B3,h0_A4B5,h0_A6B15,h0_A7B13,h0_A7B14,
            h0_A8B13,h0_A8B14,h0_A9B11,h0_A9B12,h0_A10B11,h0_A10B12,
            )
    from bccc.ta.ta_A1B1 import hu_A1B1
    from bccc.ta.ta_A1B2 import hu_A1B2
    from bccc.ta.ta_A1B3 import hu_A1B3
    from bccc.ta.ta_A2B2 import hu_A2B2
    from bccc.ta.ta_A2B3 import hu_A2B3
    from bccc.ta.ta_A3B3 import hu_A3B3
    from bccc.ta.ta_A4B5 import hu_A4B5
    from bccc.ta.ta_A6B15 import hu_A6B15
    from bccc.ta.ta_A7B13 import hu_A7B13
    from bccc.ta.ta_A7B14 import hu_A7B14
    from bccc.ta.ta_A8B13 import hu_A8B13
    from bccc.ta.ta_A8B14 import hu_A8B14
    from bccc.ta.ta_A9B11 import hu_A9B11
    from bccc.ta.ta_A9B12 import hu_A9B12
    from bccc.ta.ta_A10B11 import hu_A10B11
    from bccc.ta.ta_A10B12 import hu_A10B12
    ful += (
            hu_A1B1,hu_A1B2,hu_A1B3,hu_A2B2,hu_A2B3,hu_A3B3,hu_A4B5,hu_A6B15,hu_A7B13,hu_A7B14,
            hu_A8B13,hu_A8B14,hu_A9B11,hu_A9B12,hu_A10B11,hu_A10B12,
            )
    
    from multiprocessing import Pool
    from bccc.pub import sign,prodl,phiu,lowTri0
    # h0_
    h0d = {}
    pool = Pool()
    l = tuple(pool.apply_async(f0, (r, h, g, p, nc)) for f0 in f0l)
    pool.close()
    pool.join()
    for h0 in l: h0d.update(h0.get())
    # hu_
    hud = {}
    pool = Pool()
    l = tuple(pool.apply_async(fu, (r, h, g, p, nc)) for fu in ful)
    pool.close()
    pool.join()
    for hu in l: hud.update(hu.get())
    
    print(F'threshold,max_cycle,dEshift = {thd},{ncyc},{xE}')
    print('initial from LCC')
    from itertools import permutations,product
    Eref,Ecorr=h0d[()],sum(h0d.get(u, 0.0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu) for u,(_,tu) in hud.items()); E0=Eref+Ecorr
    print(F'Ecorr = {Ecorr}')
    
    # Create data for C++
    no,n2=h.shape[0],g.shape[0]; lT0=array(lowTri0(no),dtype=uint32)
    np,nr,nbs=len(p),312,16; rll = zeros(shape=(np, nr, nbs, nbs))
    for P in range(np):
        for R,rR in r[P].items():
            rll[P,R] = rR
    pul=ones(shape=(nt,), dtype=uint32); pul[1]=pu=np*nbs
    for i in range(2,nt): pul[i] = pul[i-1]*pu
    udd = {u:{u:1.0} for u in t}  # all actiove t
    for u,ud in udd.items():
        for v in list(permutations(u))[1:]:
            ud[v] = sign([b for b,s in v if 6<s<15])
    t0 = {v:s*t[u] for u,ud in udd.items() for v,s in ud.items()}  # unordered t
    phid={u:phiu(u, pul) for u in t0}; nu=len(phid)
    id,il={u:i for i,u in enumerate(phid.keys())},array(tuple(phid.values()), dtype=uint32)
    
    from os import path as Path
    path = Path.dirname(__file__)
    libta = ctypeslib.load_library('libta.so', path)
    libta.ta.argtypes = [
        ctypeslib.ndpointer(dtype=float64, ndim=1, flags='C_CONTIGUOUS'),  # tl
        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # pul
        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # il
        c_uint,                                                            # nu
        ctypeslib.ndpointer(dtype=float64, ndim=4, flags='C_CONTIGUOUS'),  # r
        c_uint,                                                            # np
        ctypeslib.ndpointer(dtype=float64, ndim=2, flags='C_CONTIGUOUS'),  # h
        c_uint,                                                            # no
        ctypeslib.ndpointer(dtype=float64, ndim=2, flags='C_CONTIGUOUS'),  # g
        ctypeslib.ndpointer(dtype=uint32, ndim=1, flags='C_CONTIGUOUS'),   # lT0
        c_uint,                                                            # n2
        ctypeslib.ndpointer(dtype=uint32, ndim=2, flags='C_CONTIGUOUS'),   # p
        c_uint,                                                            # nc
                   ]
    libta.ta.restype = ctypeslib.ndpointer(dtype=float64, shape=(nu,))
    sl = [(1, 1), (1, 2), (2, 2), (3, 3), (6, 15), (7, 13), (7, 14), (8, 13), (8, 14)]
    tud = [tuple(sorted(zip(bl,sl))) for bl,sl in product(permutations(list(range(np)),r=2), sl)]  # independent t
    hud0 = {u:hu for u,hu in hud.items() if u in tud}
    
    print('    it          Ecorr              dE               dt               time/s')
    ndiis=6; idiis = [*list(range(1,ndiis)), 0, ndiis]
    objl,errl = [],[]
    A = zeros((ndiis+1,ndiis+1))-1; A[ndiis,ndiis] = 0
    b=zeros((ndiis+1, )); b[ndiis]=-1
    for it in range(ncyc):
        # ta_
        t1 = datetime.now()
        tl = array(tuple(t0.values()), dtype=float64)
        l = libta.ta(tl, pul, il, nu, rll, np, h, no, g, lT0, n2, p, nc)
        td = {u:l[id[u]] for u in hud}
        t2 = datetime.now()
        
        dd = {u: (h0d.get(u, 0.0)+td[u]+(hu-E0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu))/(xE+hu-E0) for u,(hu,tu) in hud0.items()}
        t,dt = {u:t0[u]-du for u,du in dd.items()},sum([abs(du) for du in dd.values()])
        if dt>thd:
            objl.append(t); errl.append(dd)
            if it+1>=ndiis:
                if it+1>ndiis:
                    objl.remove(objl[0]); errl.remove(errl[0])
                    A[:,:] = A[:,idiis]; A[:,:] = A[idiis,:]
                    for i in range(ndiis):
                        A[ndiis-1,i] = A[i,ndiis-1] = sum(errl[i][u]*du for u,du in dd.items())
                else:
                    for i in range(ndiis):
                        for j in range(i,ndiis):
                            A[i,j] = A[j,i] = sum(errl[i][u]*errl[j][u] for u in errl[i])
                x=linalg.solve(A,b)
                t = {u:sum(x[i]*objl[i][u] for i in range(ndiis)) for u in hud0}
        for P in range(nc,np):
            for Q in range(P+1,np):
                t[((P,1),(Q,3))]=t[((P,2),(Q,3))]=t[((P,3),(Q,1))]=t[((P,3),(Q,2))]=0.0
                t[((P,4),(Q,5))]=t[((P,5),(Q,4))]=-t[((P,3),(Q,3))]
                t[((P,9),(Q,11))]=-t[((P,7),(Q,13))]; t[((P,11),(Q,9))]=-t[((P,13),(Q,7))]
                t[((P,9),(Q,12))]=-t[((P,7),(Q,14))]; t[((P,12),(Q,9))]=-t[((P,14),(Q,7))]
                t[((P,10),(Q,11))]=-t[((P,8),(Q,13))]; t[((P,11),(Q,10))]=-t[((P,13),(Q,8))]
                t[((P,10),(Q,12))]=-t[((P,8),(Q,14))]; t[((P,12),(Q,10))]=-t[((P,14),(Q,8))]
        Ecorr=sum(h0d.get(u, 0.0)*sum(s*prodl([t[v] for v in vl]) for s,vl in tu) for u,(_,tu) in hud.items())
        E=Eref+Ecorr; dE=E-E0
        print(F'    {it+1:3d}      {Ecorr:.8f}      {dE:13.10f}    {dt:13.10f}     {t2-t1}', flush=True)
        if dt<thd: 
            print(F'Successfully converged: Ecorr = {Ecorr}')
            break
        else: E0,t0 = E,{v:s*t[u] for u,ud in udd.items() for v,s in ud.items()}
        
    else: print('Error: Convergence failed ')
    
    return Ecorr,t
    

