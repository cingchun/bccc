#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   gvbbccc.py
#            Des:   for bccc
#           Mail:   qingchun720@foxmail.com
#   Created Time:   15:18 八月-24/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


import subprocess,datetime
from collections import defaultdict as ddict

from bccc import pub
import gvb

from functools import reduce
import numpy,math


class GVBBCCC(object):
    '''
    GVB-BCCC

    Attributes:


    Saved results：
    '''

    def __init__(self, mf, nt=2):
        # general parameters from mol for mf object
        self.ref = mf
        self.mol = mf.mol
        self.mo_coeff = mf.mo_coeff
        self.max_memory = mf.max_memory
        self.stdout = mf.stdout
        self.verbose = mf.verbose

        # parameters for GVB reference
        # self.nobt,self.nbf = mf.nobt,mf.nbs
        self.nobt,self.nbf = mf.nocc+mf.np,mf.nbs
        self.nc,self.nocc,self.nvir = mf.nco,mf.nocc,mf.nvir
        self.npa,self.np = mf.npa,mf.np
        self.p,self.ci = mf.pair[:self.nocc],mf.ci
        self.h,self.g = mf.h_mo,mf.g_mo

        # parameters for BCCC
        self.nt,self.r = nt,None  # r: matrix representation
        self.thd,self.ncyc,self.xE = 1.0e-5,100,0.015
        self.converged,self.t = False,None  # t: T amplitude
        self.e_nuc,self.e_ele = mf.e_nuc,0.0
        self.e_corr,self.e_tot = 0.0,0.0

    # from memory_profiler import profile
    # @profile(precision=4,stream=open('memory_profiler.dat','w+'))
    def kernel(self):
        mol = self.mol
        output = mol.output; fpre = output[:output.rfind('_')]
        nobt,nocc = self.nobt,self.nocc
        pub.lT0 = pub.lowTri0(nobt)
        nc=self.nc; nt,np=self.nt,self.np
        thd,ncyc,xE = self.thd,self.ncyc,self.xE

        # mo integral and repressive matrix
        print('get mo integral and repressive matrix')
        datfile=F'{fpre}_gvb_init.dat'; fchkfile=F'{datfile}.fchk'
        from pyscf import scf,ao2mo
        hcore = scf.hf.get_hcore(mol)
        eri = mol.intor('int2e', aosym='s8')
        obts = self.mo_coeff[:,:nobt]
        self.h = reduce(numpy.dot, (obts.T, hcore, obts))
        self.g = ao2mo.kernel(eri, obts)
        from bccc.rep import rep
        r = [None]*nocc
        for P in range(nocc):
            mc = numpy.identity(16)
            if nc<=P:
                mc[0,0:2] = self.ci[P]
                mc[1,0:2] = [-mc[0,1], mc[0,0]]
            a = 1/math.sqrt(2.)
            mc[2,2:4] = [a, -a]
            mc[3,2:4] = [a, a]
            r[P] = rep(mc)
        self.r = r
        print( '\n\n\n\n')

        # check GVB energy
        print('check GVB energy')
        from bcpt2.read_pair import read_eng_dat
        E_gvb,e_nuc = read_eng_dat(datfile)
        print(F'e_nuc = {e_nuc}')
        print(F'E(GVB) = {E_gvb:.10f}   (.dat)')
        mf = self.ref
        mf.h_mo,mf.g_mo = self.h,self.g; gvb.rgvb.lT0 = pub.lT0
        mf.e_ele = mf.energy_elec()
        E_gvb = mf.e_nuc+mf.e_ele
        print(F'E(GVB) = {E_gvb:.10f}   (Ref.)')
        from bccc.hm.hm_ import hm_
        h0d = hm_(r, self.h, self.g, self.p, nc)[()]
        # h0d = hm_(self)[()]
        E_gvb = self.e_nuc+h0d[()]
        print(F'E(GVB) = {E_gvb:.10f}')
        print('\n\n\n\n')

        # BCPT2
        print('check GVB-BCPT2 energy')
        from bcpt2 import gvbpt2
        ael_pt2,E_bcpt2 = gvbpt2.gvb_bcpt2(mol, datfile, fchkfile, fixed_core=nc, fixed_virt=self.nvir-self.np)
        ael = [None]*nocc
        for P in range(nocc):
            ael[P] = numpy.zeros(16)
            if nc<=P:
                ael[P][:15] = ael_pt2[P-nc]
                ael[P][:] = ael[P][[4, 5, 6, 7, 8, 9, 15, 0, 1, 2, 3, 10, 11, 12, 13, 14]]
        print(F'E(GVB-BCPT2) = {E_bcpt2:.10f}   (Ref.)')
        from bccc.pub import bcpt2
        corr_bcpt2 = bcpt2(ael, h0d)
        print(F'corr(GVB-BCPT2) = {corr_bcpt2}')
        E_bcpt2 = E_gvb+corr_bcpt2
        print(F'E(GVB-BCPT2) = {E_bcpt2}')
        print('\n\n\n\n')
        
        # check FCI/MCSCF/DMRG
        print('check FCI/MCSCF/DMRG', flush=True)
        hf = scf.RHF(self.mol)
        hf.max_cycle = 1
        hf.kernel()
        hf.mo_coeff = self.mo_coeff
        from pyscf import mcscf
        mc = mcscf.CASCI(hf, 2*np, 2*np, ncore=nc)
        #mc.fcisolver.nroots = nu
        mc.fcisolver.memory = 30  # GB
        mc.fix_spin_(ss=0.0)
        #mc.fcisolver.level_shift = 0.2
        #mc.fcisolver.pspace_size = 1200
        if np>6:  # DMRG
            from pyscf import dmrgscf
            mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=1000)
            mc.fcisolver.block_extra_keyword = ['memory, 30, g']
        #mc.kernel()
        print('\n\n\n\n')
        
        # initilize t
        print('initialize t')
        t = self.t
        if t is None:
            print('t0 from GVB-LBCCC')
            from bccc.hm.hm import hm
            t1 = datetime.datetime.now()
            mh = hm(r, self.h, self.g, self.p, nt, nc)
            # mh = hm(self)
            t2 = datetime.datetime.now()
            print(F'calculate H matrix {t2-t1}')
            from bccc.pub import bclcc
            corr_bclcc,t = bclcc(nocc, mh, nt, nc)
            print(F'corr(GVB-BCLCC{nt}) = {corr_bclcc}')
            E_bclcc = E_gvb+corr_bclcc
            print(F'E(GVB-BCLCC{nt}) = {E_bclcc}')
            import gc
            del mh
            gc.collect()
        else:
            finit=t  # initial t0 file
            print(F't0 from {finit}')
            t = ddict(float)
            with open(finit, 'r') as input:
               linel=input.readlines()
               linel=linel[pub.iis('t(final) =\n', linel)+1:]
               linel=linel[:pub.iis('\n', linel)]
               from re import findall
               for line in linel:
                   sl = findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", line.strip('\r\n'))
                   t[tuple((int(sl[i]), int(sl[i+1])) for i in range(0,len(sl)-1,2))] = float(sl[-1])
        print('\n\n\n\n')

        # GVB-BCCC
        print('calculate GVB-BCCC energy')
        # from bccc.ta1.ta import ta
        # print('From ta1')
        from bccc.iter import iter
        # corr_bccc,t = ta(t, mh, self.p, nt, nc)
        corr_bccc,t = iter(t, r, self.h, self.g, self.p, nt, nc, thd, ncyc, xE)
        print(F'corr(GVB-BCCC{nt}) = {corr_bccc}')
        E_bccc = E_gvb+corr_bccc
        print(F'E(GVB-BCCC{nt}) = {E_bccc}')
        print('t(final) =')
        for u,tu in t.items():
            if abs(tu)>1.0e-6: print(F'u,tu = {u},{tu}')
        #print('t symmetry =')
        #ul = list(t.keys()); ul_,nu_ = [],len(t)
        #for i,u in enumerate(ul):
        #    if u not in ul_:
        #        tu=t[u]
        #        print(F'{u}: {tu}')
        #        ul_.append(u)
        #        tu_=abs(tu)
        #        for v in ul[i+1:]:
        #            if v not in ul_:
        #                tv=t[v]
        #                if abs(tu_-abs(tv))<1.0e-7:
        #                    print(F'{v}: {tv}')
        #                    ul_.append(v)
        #    print()
        print('\n\n\n\n')

