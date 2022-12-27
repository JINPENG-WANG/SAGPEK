echo param %1 
echo param %2
set arg1=%1
set arg2=%2


Rscript Src/Script01.get.amPeak.R

perl Src/Script02.get.genotypes.pl %arg1% %arg2%
