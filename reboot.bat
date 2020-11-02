:: reboot win64
::
:: Author: Qingchun Wang @ NJU
:: E-mail: qingchun720@foxmail.com
::

for /l %%i in (1,1,7) do timeout /t 86400 /nobreak > nul
shutdown -r -f -t 5

