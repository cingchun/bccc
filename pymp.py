#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   pyic.py
#            Des:   python invoke(call) c
#           Mail:   qingchun720@foxmail.com
#   Created Time:   16:48 10月-22/2019
#


__version__ = '1.0'
__author__ = 'Qingchun Wang'


import os
from multiprocessing import Pool,Manager
# from mpi4py import MPI
from ctypes import CDLL,c_double,c_long
from datetime import datetime
# import numpy,scipy,math

lock = Manager().Lock()
d,nd = Manager().dict(),Manager().Value('l', 0)
def f1(i):
    pid = os.getpid()
    print(F' {i}: Hello world from os.getpid() = {pid}')
    lock.acquire()
    nd.value += 1
    lock.release()
    d[i%4] = d.get(i%4, [])+[pid]
    return pid
N = 100
def func(i):
    n = 0
    
    # ithread = os.getpid()
    # print("from OMP thread %d: i = %d" %(ithread, i))
    for j in range(N):
        for k in range(N):
            n += (i+1)*(j+1)*(k+1)
    
    return n
def pymp():
    '''
    python multiprocess(parallel)
    
    Threading多线程，multiprocess多进程，两者用法相似, 都须将任务封装成函数
    multiprocess.Process: 创建一个进程，不能有返回值，但若想返回可借用Queue来保存
    multiprocess.Pool: 创建多个进程，可以有返回值
        multiprocess.Pool.map：同map, 将一个函数作用一系列值上，返回多个结果
        multiprocess.Pool.apply 同步 multiprocess.Pool.apply_async 异步：
             一次创建一个进程，返回一个结果
             所以要借助推导式创建一系列，返回多个
             可以是多个函数并行，如果有多个函数要并行，将函数迭代传入
        
        multiprocess.close 不继续添加
        multiprocess.join 等待所有子进程执行完毕
        
    https://blog.csdn.net/weixin_38611497/article/details/81490960
    
    :return:
    '''

    ''' mulitprocessing '''
    pool = Pool()
    print(F' os.cpu_count() = {os.cpu_count()}')
    # pool.map(f1, range(12))
    pl = tuple(pool.apply_async(f1, (i,)) for i in range(11))
    pool.close()
    pool.join()
    print(F' nd.value = {nd.value}')
    for pid,i in d.items(): print(F' d[{pid}] = {i}')

    ''' mpi4py '''
    # comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()
    # node_name = MPI.Get_processor_name() # get the name of the node
    # print('Hello world from process %d at %s.' %(rank, node_name))
    
    v,nl = 0.0,list(range(N))
    pool = Pool()
    # pl = tuple(pool.apply_async(func, (n,)) for n in nl)
    # pool.close()
    # pool.join()
    # vl = [p.get() for p in pl]
    vl = pool.map(func, nl) #
    # for _ in vl: v += _
    v = sum(vl)

    return v
libcmp = CDLL('./cmp.so')
cmp = libcmp.cmp
cmp.restype = c_long



def cpymp():
    n = pymp()



if __name__ == '__main__':
    print("*** Test python call parallel c ***")
    t1 = datetime.now()
    n = cmp()
    # n = pymp()
    t2 = datetime.now()
    print(F' n = {n}')
    print(F' Elapsed time: {t2-t1} sec. ')
    print("*** cloze testing ***")
    
    