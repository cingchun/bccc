#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   __init__.py.py
#            Des:   __init__.py
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:04 8月-11/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


def BCCC(mf, nt=2, df=None):
    '''
    Block Correlation Coupled Cluster

    :param mf: multi-reference wavefunction. one of GVB, CAS, APSG, DMRG, UGVB, UCAS, UAPSG, UDMRG
    :param nt: ntuply of T
    :param nc: number of core orbitals
    :param thd: threshold
    :param ncyc: number of cycle
    :param xE: x-factor of Energy
    :param df: Density fitting, (mo_coeff, mo_occ)

    :return:
    '''

    from bccc import gvbbccc,ugvbbccc,gvbbcccdf,ugvbbcccdf
    from bccc import casbccc,ucasbccc

    from pyscf import gto
    import gvb
    if isinstance(mf,gto.Mole):
        na,nb = mf.nelec
        if na is nb: mf = gvb.RGVB(mf, pop='meta_Lowdin', weight='trans_dipo', auxbasis='sto-6g', sch='sch-ii')
        else: mf = gvb.UGVB(mf, pop='meta_Lowdin', weight='trans_dipo', auxbasis='sto-6g', sch='sch-ii')
        mf.kernel()

    if isinstance(mf, gvb.rgvb.GVB):
        if df: return gvbbcccdf.GVBBCCCDF(mf, nt, df)
        else: return gvbbccc.GVBBCCC(mf, nt)
    elif isinstance(mf, gvb.ugvb.UGVB):
        na,nb = mf.mol.nelec
        if df: return ugvbbcccdf.UGVBBCCCDF(mf, na-nb, nt, df)
        else:
            return ugvbbccc.UGVBBCCC(mf, na-nb, nt)
    else:
        raise NotImplementedError('TODO: other --BCCC-- mode')







