@echo off
set arg1=%1
set arg2=%2
set arg3=%3

Rscript Src/Script01.get.amPeak.R %arg1%

perl Src/Script02.get.genotypes.pl %arg1% %arg2% %arg3%
echo "$SAGPEK: Program completed"
