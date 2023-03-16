TAMIPAMI: Software to find protospacer-adjacent motif (PAM) sites
======
Author
------
* Adam Rivers, US Department of Agriculture, Agricultural Research Service
* Lei Ma, US Department of Agriculture, Agricultural Research Service

Citation
------
TBD

Introduction
------
TBD

Installation
------
TBD

Usage
------

positional arguments:

  {library,custom}      Choose whether to use the predetermined library or your own custom spacer/orientation combination

    library             Toggles predetermined library mode, exclusive with custom. Specify the -lib-name with one of ["RTW544",
                        "RTW555", "RTW572", "RTW574"]

    custom              Toggles the custom mode, exclusive with library. Use -spacer SPACER and -orientation ["5prime","3prime"]

options:

  -h, --help            show this help message and exit

  --cont1 CONT1, -c CONT1
                        A forward .fastq, .fq, .fastq.gz or .fq.gz file. .

  --cont2 CONT2, -c2 CONT2
                        A reverse .fastq, .fq, .fastq.gz or .fq.gz file.

  --exp1 EXP1, -e EXP1  A forward .fastq, .fq, .fastq.gz or .fq.gz file. .

  --exp2 EXP2, -e2 EXP2
                        A reverse .fastq, .fq, .fastq.gz or .fq.gz file.

  --log LOG             Filename for log file

  --length [1-10]       The length of the PAM or TAM sequences

  --out OUT, -o OUT     Filename to save full dataframe of sequence frequencies

Prebuilt configurations for spacer/orientation pairs are:

+--------+-----------------------+--------------+
| Name | Spacer | Orientation |
+========+=======================+==============+
| RTW572 | GGAATCCCTTCTGCAGCACCTGG| 3prime |
+--------+-----------------------+--------------+
| RTW554 | GGGCACGGGCAGCTTGCCGG | 3prime |
+--------+-----------------------+--------------+
| RTW555 | GTCGCCCTCGAACTTCACCT | 5prime |
+--------+-----------------------+--------------+
| RTW574 | CTGATGGTCCATGTCTGTTACTC| 5prime |
+--------+-----------------------+--------------+

Example
------
Finds TAMs of length 4 on the 3 prime end of spacer 'GGAATCCCTTCTGCAGCACCTGG'
.. code-block:: bash
    tamipami -c control_R1.fastq -c2 control_R2.fastq -e exp_R1.fastq -e2 exp_R2.fastq --length 4 library -lib-name RTW472 

Finds TAMs of length 5 on the 5 prime end of custom spacer 'GGAATCCCTTCTGCAGCACCTGG'
.. code-block:: bash
    tamipami -c control_R1.fastq -c2 control_R2.fastq -e exp_R1.fastq -e2 exp_R2.fastq custom -spacer GGAATCCCTTCTGCAGCACCTGG -orientation 5prime

