# ATPdock: 
A template-based method for ATP-specific protein-ligand docking.

## Pre-requisite:
    - Python3, Java
    - MGLTools software (http://mgltools.scripps.edu/downloads)
    - OpenBabel software (http://openbabel.org/wiki/Category:Installation)
    - Linux system

## Installation:

*Install and configure the softwares of Python3, Java, MGLTools and OpenBabel in your Linux system.
*Download this repository at https://github.com/brightrao/ATPdock.git. Then, uncompress it and run the following command lines on Linux System.

Please make sure that python3 includes the modules of 'os', 'math', 'numpy', 'random', 'subprocess', 'sys', and 'shutil'. If any one modules does not exist, please using 'pip3 install XXXX' command install the python revelant module. Here, "XXXX" is one module name.

~~~
  $ chmod -R 777 ./software
  $ chmod -R 777 ./basefile
  $ java -jar FileUnion.jar ./PLDB/ ./PLDB.tar.gz
  $ rm -rf ./PLDB
  $ tar zxvf PLDB.tar.gz
  $ cp [THE INSTALLATION PATH OF MGLTools]/bin/pythonsh ./basefile
  $ cp [THE INSTALLATION PATH OF MGLTools]/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py ./basefile
  $ cp [THE INSTALLATION PATH OF MGLTools]/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py ./basefile
~~~

NOTE THAT, [THE INSTALLATION PATH OF MGLTools] SHOULD BE THE MGLTools INSTALLATION PATH.

*The file of "Config.properties" should be set as follows:
~~~
PPSALIGN_EXE_ABOSOLTE_PATH=XXX/software/PPSalign
APOC_EXE_ABOSOLTE_PATH=XXX/software/apoc
PPSSEARCH_DATABASE_ABOSOLTE_PATH=XXX/PLDB
~~~

NOTE THAT, THE ABOVE "XXX" SHOULD BE THE ABSOLUTE PATH OF THE DOWNLOADED ATPDOCK PACKAGE, e.g., /home/junh/ATPdock-main


## Prepare the query folder

*The query folder have to contain three files, i.e., "pdb.pdb", "tem.txt", "pdb.site".
~~~
"pdb.pdb" contains the 3D structure information of the query receptor protein in PDB format.

"tem.txt" has two lines: 1) the first line is sequence identity cutoff, ranging from 0.3 to 1 (suggested 0.3); 2) the second line is searched ligand type. If only search for template proteins that bind to ATP and ADP, the second line is ~ATP~~ADP~. If all ligand types are allowed, the second line is NULL. More ligand type information could be found at https://zhanglab.ccmb.med.umich.edu/BioLiP/ligand.html. For example,
        0.3
        ~ATP~~ADP~~AMP~~GTP~~GDP~

"pdb.site" is the information of binding residue type and index. In this file, each line means one pocket. User can define multiple binding pockets of the protein. For example, if one protein has two binding pockets, there are two lines, i.e.,
        V10 G38 C39 G40 R61 P89 G90 D91 G92 K93
        E65 I66 V67 N104 R106
~~~

NOTE THAT, IF YOU DONOT KWON THE ATP-BINDING POCKET INFORMATION OF THE QUERY PROTEIN, WE STRONGLY SUGGESTED THAT USING THE WEB-SERVER OF ATPbind (https://zhanglab.ccmb.med.umich.edu/ATPbind/) TO PREDICTED ITS BINDING POCKETS. BASED ON THE ATPbind-PREDICTED POCKET, PREPARING THE FILE OF "pdb.site".


## Run example
~~~
  $ python3 ATPdock.py yy/example
~~~
NOTE THAT, "yy" should be the absolute path of the example folder.
The outputted "ATPx" folder, "x" is a number. In the folder, "final.pdb" is the docking result of ATP.

## Update History:

- First release     2021-05-31 (see https://github.com/brightrao/ATPdock/)
- Update PPS-search 2021-07-28

## References

[1] Liang Rao, Ning-Xin Jia, Jun Hu, Dong-Jun Yu, and Gui-Jun Zhang. ATPdock: a template-based method for ATP-specific protein-ligand docking. Bioinformatics. sumitted.
