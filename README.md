如有使用问题，请联系负责人： 王金鹏，博士，副研究员 wangjinpeng0225@163.com

山东省农业科学院畜牧兽医研究所

研究方向：奶牛基因组选择，重要经济性状关键基因挖掘



# SAGPEK
  Software packages for Sanger sequencing have shown little progress over time, although it is still widely used in human genetic disease diagnosis and livestock animal breeding. Determining the genotypes of tens to hundreds of sites in hundreds to thousands of samples simultaneously in Sanger sequencing data by manual visual conformation with traditional software is still time-consuming and error-prone. Here we present SAGPEK, which automatically identifies genotypes of target loci from hundreds to thousands of ABI-format Sanger sequencing data and directly reports genotypes to the output file. 
  SAGPEK is suitable for human genetic disease screening, drug-resistant mutation identification, livestock animal functional mutation detection, and other similar situations.

# Installation
For the Windows platform, users should install R and Perl interpreters, like the recent releases of R from the comprehensive R archive network (CRAN) website (https://www.r-project.org/) and Perl from the ActiveState website (https://www.activestate.com/products/perl/). After the installation, the directories in which executable programs are located should be added to the PATH environment variable. Before running the package, users should put the ABI-format files generated by the Sanger sequencing technology into the "ABI" file folder of the SAGPEK package. Then open the Command Prompt window, and change the working directory to the SAGPEK package directory with the "cd" command. Run the program with the "SAGPEK.bat TYPE" command, where the TYPE should be replaced with HBV, PAH, cattle_CN, cattle_DUMPS, TEST, or custom according to the actual situations. If the "custom" type is used, users should also provide the name of a file storing the customized anchoring sequences located in the "Custom" directory and add the file name to the command as "SAGPEK.bat custom FILE". The file must consist of two columns. The first column stores the names of anchoring sequences, and the second column stores the sequences. After the completion, the result file consisting of genotypes or amino acid alterations will be generated in the Genotype directory. 

For Unix-like systems, users can install SAGPEK directly, as R and Perl are built-in software along with the system. Users should put the ABI-format files in the "ABI" directory and run the SAGPEK.sh with the "bash SAGPEK.sh -t TYPE" command. The TYPE should be replaced with HBV, PAH, cattle_CN, cattle_DUMPS, TEST, or custom. If the "custom" type is used, users should provide a file name of anchoring sequences with the -g option, for example, "bash SAGPEK.sh -t custom -g FILE". Result files are saved under the Genotype directory.

# Usage Examples
We list a set of common uses of SAGPEK below. We provide ABI format input files in the ABI folder, which can be used to test the SAGPEK. 
For Windows platforms, the SAGPEK can be used like:

$  SAGPEK.bat TEST

$  SAGPEK.bat HBV

$  SAGPEK.bat PAH

$  SAGPEK.bat cattle_CN

$  SAGPEK.bat cattle_DUMPS

$  SAGPEK.bat custom ab.tags.txt

$  SAGPEK.bat custom ab.tags.txt r


For Unix-like systems, the SAGPEK can be used like:

$ bash SAGPEK.bash -h 
       #This will print help page.
 
$ bash SAGPEK.bash -u 
       #This will print the examples about how to use SAGPEK.

$ bash SAGPEK.bash -t TEST

$ bash SAGPEK.bash -t HBV

$ bash SAGPEK.bash -t PAH

$ bash SAGPEK.bash -t cattle_CN

$ bash SAGPEK.bash -t cattle_DUMPS

$ bash SAGPEK.bash -t custom -g ab.tags.txt

$ bash SAGPEK.bash -t custom -g ab.tags.txt -r

For the TEST, HBV, PAH, cattle_CN, and cattle_DUMPS types, SAGPEK has built-in anchoring sequences. For the custom type, anchoring sequences will be obtained from the "ab.tags.txt" file, which SAGPEK provided as an example file under the "Custom" folder.
