#! /bin/bash
#
# Author: Qingchun Wang @ NJU
# E-mail: qingchun720@foxmail.com
#

username='1805607'
password='zzdbl1314520'


while true
do
	ping -c 1 www.baidu.com > /dev/null
	if [ $? -ne 0 ];then
		curl -d "username=${username}&password=${password}" http://p.nju.edu.cn/portal_io/login </dev/null 2>/dev/null
		#curl -d "username=${username}&password=${password}&Onclick=" http://p.nju.edu.cn/portal_io/login
  fi
	sleep 30m
done

