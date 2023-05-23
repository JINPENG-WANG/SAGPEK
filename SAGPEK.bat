@echo off
Rem The arg1 stores TYPE. The arg2 stores tag file name containing anchoring sequences. The arg3 stores the orientation parameter r or nothing. 
set arg1=%1
set arg2=%2
set arg3=%3

Rscript Src/Script01.get.amPeak.R %arg1%

perl Src/Script02.get.genotypes.pl %arg1% %arg2% %arg3%
echo "$SAGPEK: Program completed"
