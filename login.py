#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#         Author:   Qingchun Wang @ NJU
#         E-Mail:   qingchun720@foxmail.com
#


import requests
import os,time

addr="https://p.nju.edu.cn/portal/login"
username='dg1624065'
password='072088_w'

header={
    'Accept': '*/*',
    'Accept-Encoding': 'gzip, deflate, br',
    'Accept-Language': 'zh-CN,zh;q=0.9',
    'Connection': 'keep-alive',
    'Content-Length': '107',
    'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8',
    'Cookie': 'login=bQ0pOyR6IXU7PJaQQqRAcBPxGAvxAcroYpuUxcq5od7dlmpltnEal5DQ0gjD6r1n3%252Fhz5Ndv3l%252FxDuNn8jsHuEDCr2BFRDfRYRw2lSpGv8mAsB%252FTG6xFGqlgUw0Xjk65OPxQNFGhLmZ24drwZxp8kv8nzffCTVZo9pEs7xzVqNwVNbU64ooymQU%253D; _ga=GA1.3.908443897.1522253864; login=bQ0pOyR6IXU7PJaQQqRAcBPxGAvxAcroYpuUxcq5od7dlmpltnEal5DQ0gjD6r1n3%252Fhz5Ndv3l%252FxDuNn8jsHuEDCr2BFRDfRYRw2lSpGv8mAsB%252FTG6xFGqlgUw0Xjk65OPxQNFGhLmZ24drwZxp8kv8nzffCTVZo9pEs7xzVqNwVNbU64ooymQU%253D',
    'Host': 'p.nju.edu.cn',
    'Origin': 'https://p.nju.edu.cn',
    'Referer': 'https://p.nju.edu.cn/srun_portal_pc.php?ac_id=1&',
    'User-Agent': 'Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/67.0.3396.99 Safari/537.36',
    'X-Requested-With': 'XMLHttpRequest'
}
data = {
    'action': 'login',
    'username': username,
    'password': password,
    'ac_id': '1',
    'user_ip': '',
    'nas_ip': '',
    'user_mac': '',
    'save_me': '0',
    'ajax': '1'
}


while True:
    if os.system('ping www.baidu.com > nul'):
        response = requests.post(addr, data=data, headers=header)
    time.sleep(1800)
    
    