set MYRSCRIPT=\\ufomcchead\cc\sasha\clusterhlpr.r
set MYRTERM=\\ufomcchead\cc\_programme\R\R-3.1.0\bin\x64\Rterm
set MYRMEM=2048M
%MYRTERM% --max-mem-size=%MYRMEM% --no-restore --no-save < %MYRSCRIPT%