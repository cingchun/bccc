:: login on win64
::
:: Author: Qingchun Wang @ NJU
:: E-mail: qingchun720@foxmail.com
::


set username=1805607
set password=zzdbl1314520

cd /d D:\Documents\PyCharmProjects\bccc
bash login.sh
:loop
    curl http://p.nju.edu.cn/portal_io/login -d "username=%username%&password=%password%&Onclick="
    ::for %%i in (1,2) do timeout /t 86400 /nobreak > nul
    timeout /t 21600 /nobreak > nul
goto :loop

::pause

