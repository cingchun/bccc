#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   gvbbccc.py
#            Des:   for H matrix check
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
        self.ref = mf
        self.mol = mf.mol
        self.mo_coeff = mf.mo_coeff
        self.max_memory = mf.max_memory
        #self.stdout = mf.stdout
        #self.verbose = mf.verbose

        self.nobt,self.nbf = mf.nobt,mf.nbs
        self.nocc,self.nvir = mf.nocc,mf.nvir
        self.npa,self.np = mf.npa,mf.np
        self.p,self.ci = mf.pair[:self.nocc],mf.ci
        self.h,self.g = mf.h_mo,mf.g_mo

        self.nc = mf.nco
        self.nt = nt
        self.thd,self.ncyc,self.xE = 1.0e-7,100,0.015

        self.converged,self.t = False,None # amplitude
        self.e_nuc,self.e_ele = mf.e_nuc,0.0
        self.e_corr,self.e_tot = 0.0,0.0

    # from memory_profiler import profile
    # @profile(precision=4,stream=open('memory_profiler.dat','w+'))
    def kernel(self):
        mol = self.mol
        output = mol.output; fpre = output[:output.rfind('_')]
        nobt,nocc = self.nobt,self.nocc
        # lowTri = [0]*nobt
        # for i in range(1, nobt): lowTri[i] = lowTri[i-1]+i
        # pub.lowTri=lowTri; gvb.rgvb.i_lowTri=lowTri
        pub.lT0 = pub.lowTri0(nobt)
        nc=self.nc; nt,np=self.nt,self.np
        thd,ncyc,xE = self.thd,self.ncyc,self.xE

        # mo integral and repressive matrix
        print('get mo integral and repressive matrix')
        datfile=F'{fpre}_gvb_init.dat'; fchkfile=F'{datfile}.fchk'
        from pyscf import scf,ao2mo
        hcore = scf.hf.get_hcore(mol)
        eri = mol.intor('int2e', aosym='s8')
        #onp=18; gvbfile=F'{fpre}_gvb{onp}.out'  # reload CI
        #with open(gvbfile, 'r') as input:
        #    linel=input.readlines()
        #    linel=linel[pub.iis('GVB CI:\n', linel)+1:]
        #    linel=linel[:pub.iis('\n', linel)]
        #    from re import findall
        #    for p in range(self.np):
        #        p1=-p-1
        #        sl = findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", linel[p1].strip('\r\n'))
        #        self.ci[p1] = [float(sl[1]), float(sl[3])]
        #        print(F'{self.p[p1, 0] + 1:2d}({self.ci[p1, 0]:6f}) <-> {self.p[p1, 1] + 1:2d}({self.ci[p1, 1]:6f})')
        #from gvb.rgvb import load_mo  # reload mo
        #load_mo(self, F'{fpre}_gvb{onp}.fchk')
        obts = self.mo_coeff
        self.h = reduce(numpy.dot, (obts.T, hcore, obts))
        self.g = ao2mo.kernel(eri, obts)
        from bccc.rep import rep
        # from bcpt2 import gvbpt2
        # ael_pt2, mcl_pt2 = gvbpt2.gvb_bcpt2(mol, datfile, fchkfile, fixed_core=nc)
        # print(F'mcl_pt2.shape = {mcl_pt2.shape}')
        nbs = 16
        rl = [None]*nocc
        # print( 'mcl = ')
        for P in range(nocc):
            # ael[P] = numpy.zeros(nbs)
            # if nc<=P:
            #     ael[P][:15] = ael_pt2[P-nc]
            #     ael[P][:] = ael[P][[4, 5, 6, 7, 8, 9, 15, 0, 1, 2, 3, 10, 11, 12, 13, 14]]
            # # mc = numpy.zeros(shape=(nbs, nbs))
            # # mc[:15, :15] = mcl_pt2[P]; mc[15][15] = 1.
            # # mc[:,:] = mc[[4, 5, 6, 7, 8, 9, 15, 0, 1, 2, 3, 10, 11, 12, 13, 14], :]
            # # mc[:,:] = mc[:, [4, 6, 7, 5, 8, 9, 15, 0, 1, 2, 3, 10, 11, 12, 13, 14]]
            # # mcl[P] = mc
            # # rl[P] = mrfeqd(mc)
            mc1 = numpy.identity(nbs)
            if nc<=P:
                # self.ci[P] = mcl_pt2[P-nc][4][4:6]
                # mc1[0,0:2] = [mcl_pt2[P-nc][4,4], mcl_pt2[P-nc][4,5]]
                mc1[0,0:2] = self.ci[P]
                mc1[1,0:2] = [-mc1[0,1], mc1[0,0]]
            a = 1/math.sqrt(2.)
            mc1[2,2:4] = [a, -a]
            mc1[3,2:4] = [a, a]
            # print(F'P = {P}')
            # print(F'    0: {mc1[0][:4]}')
            # print(F'    1: {mc1[1][:4]}')
            # print(F'    2: {mc1[2][:4]}')
            # print(F'    3: {mc1[3][:4]}')
            # print(F'    7: {mc1[7][7:9]}')
            # print(F'    8: {mc1[8][7:9]}')
            # print(F'    9: {mc1[9][9:11]}')
            # print(F'   10: {mc1[10][9:11]}')
            # print(F'   11: {mc1[11][11:13]}')
            # print(F'   12: {mc1[12][11:13]}')
            # print(F'   13: {mc1[13][13:15]}')
            # print(F'   14: {mc1[14][13:15]}')
            rl[P] = rep(mc1)
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
        h0d = hm_(rl, self.h, self.g, self.p, nc)[()]
        # h0d = hm_(rl, self)[()]
        E_gvb = self.e_nuc+h0d[()]
        print(F'E(GVB) = {E_gvb:.10f}')
        print('\n\n\n\n')

        # BCPT2
        print('check GVB-BCPT2 energy')
        from bcpt2 import gvbpt2
        ael_pt2,E_bcpt2 = gvbpt2.gvb_bcpt2(mol, datfile, fchkfile, fixed_core=nc, fixed_virt=self.nvir-self.np)
        ael = [None]*nocc
        for P in range(nocc):
            ael[P] = numpy.zeros(nbs)
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
        # from bccc.hm1.hm import hm
        # print('From hm1')
        #from bccc.hm.hm import hm
        #t1 = datetime.datetime.now()
        #mh = hm(rl, self.h, self.g, self.p, self.nt, nc)
        #t2 = datetime.datetime.now()
        #print(F'calculate H matrix {t2-t1}')
        #nu=len(mh); id={u: i for i, u in enumerate(mh.keys())}
        #mH = numpy.zeros(shape=(nu, nu))
        #for bra,i in id.items():
        #    # print(F'bra = {bra}')
        #    for ket,j in id.items():
        #        mH[i,j] = mh[bra].get(ket, 0.0)
        #        # if abs(mH[i,j])>1.0e-6: print(F'   ket = {ket}: {mH[i,j]}')
        #eig,v = numpy.linalg.eigh(mH)
        ## print(F'eig = {eig}')
        #E_fci = self.e_nuc+eig[0]
        #print(F'E(FCI) = {E_fci}')
        hf = scf.RHF(self.mol)
        hf.max_cycle = 1
        hf.kernel()
        # gvb.rgvb.load_mo(hf, fchkfile)
        # hf.mo_coeff[:,nc:nocc+np] = self.mo_coeff[:, list(range(nc, 2*nocc-nc, 2))+list(range(2*nocc-1-nc, nc, -2))]
        hf.mo_coeff = self.mo_coeff
        # FCI
        # from pyscf import fci
        # mc = fci.FCI(hf)
        # mc.kernel()
        # CASCI
        from pyscf import mcscf
        mc = mcscf.CASCI(hf, 2*np, 2*np, ncore=nc)
        #mc.fcisolver.nroots = nu
        mc.fcisolver.memory = 30  # GB
        #e_fci,eig_fci,*_ = mc.kernel()
        ## print(F'eig_fci = {eig_fci}')
        #E_fci = e_fci[0]
        #print(F'E(FCI) = {E_fci} (Ref.)')
        # CASSCF
        # mc = mcscf.CASSCF(hf, 2*np, 2*np, ncore=nc)
        mc.fix_spin_(ss=0.0)
        mc.fcisolver.level_shift = 0.2
        mc.fcisolver.pspace_size = 1200
        if np>6: # DMRG
            from pyscf import dmrgscf
            mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=1000)
            mc.fcisolver.block_extra_keyword = ['memory, 30, g']
        mc.kernel()
        #err = [abs(a-b) for a,b in zip(eig, eig_fci)]
        #emax = max(err)
        #print(F'max error = {emax}')
        #emean = sum(err)/nu
        #print(F'mean error = {emean}')
        print('\n\n\n\n')
        
        # initilize t
        print('initialize t')
        # from bccc.init import init
        # t0 = init(nocc)
        from bccc.hm.hm import hm
        t1 = datetime.datetime.now()
        mh = hm(rl, self.h, self.g, self.p, self.nt, nc)
        # mh = hm(rl, self)
        t2 = datetime.datetime.now()
        print(F'calculate H matrix {t2-t1}')
        from bccc.pub import bclcc
        corr_bclcc,t0 = bclcc(nocc, mh, nt, nc)
        print(F'corr(GVB-BCLCC{nt}) = {corr_bclcc}')
        E_bclcc = E_gvb+corr_bclcc
        print(F'E(GVB-BCLCC{nt}) = {E_bclcc}')
        import gc
        del mh
        gc.collect()
        #t0 = ddict(float)
        #finit='C4H10_4_6-31g_bccc2b.out'
        #print(F'finit = {finit}')
        #with open(finit, 'r') as input:
        #   linel=input.readlines()
        #   linel=linel[pub.iis('t(final) =\n', linel)+1:]
        #   linel=linel[:pub.iis('\n', linel)]
        #   from re import findall
        #   for line in linel:
        #       sl = findall(r"[-+]?\d+\.?\d*[eE]?[-+]?\d*", line.strip('\r\n'))
        #       t0[tuple((int(sl[i]), int(sl[i+1])) for i in range(0,len(sl)-1,2))] = float(sl[-1])
        #print('t(initial) =')
        #for u,tu in t0.items():
        #    if abs(tu)>1.0e-6: print(F'u,tu = {u},{tu}')
        print('\n\n\n\n')

        # GVB-BCCC
        print('calculate GVB-BCCC energy')
        # from bccc.ta1.ta import ta
        # print('From ta1')
        from bccc.iter import iter
        # corr_bccc,t = ta(t0, mh, self.p, nt, nc)
        corr_bccc,t = iter(t0, rl, self.h, self.g, self.p, nt, nc, thd, ncyc, xE)
        print(F'corr(GVB-BCCC{nt}) = {corr_bccc}')
        E_bccc = E_gvb+corr_bccc
        print(F'E(GVB-BCCC{nt}) = {E_bccc}')
        print('t(final) =')
        for u,tu in t.items():
            if abs(tu)>1.0e-6: print(F'u,tu = {u},{tu}')
        print('\n\n\n\n')

