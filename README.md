# ScePlantCells
This repository contains C++ code for modeling plant cells using the subcellular element model.

Below are the compilation instructions for the model and the supporting programs 
for bulk job submission on UCR's HPCC.

====================================================

## main.cpp / makefile / program.out   

====================================================

compile with: make all 
run with: sbatch <SBATCH_SCRIPT.sh> 
(.sh files are automatically generated with
batchmaker.cpp/batchGenerator.out)

IN: Folder names (Formatted in batch script), parameter flags.
Parameters:
-WR <double> 
	This parameter corresponds to alpha_{WUS} in the SI of 
	the main text.
-CKR <double> 
	This parameter corresponds to alpha_{CK} in the SI of 
	the main text.
-div <int>
	This signifies a choice of division plane mechanism.
	1: CAE-E Division mechanism
	2: CED Division mechanism
	3: CAE-M Division mechanism
	4: Layer Specific Division Mechanism
-WUS_change <boolean>
	This flag allows for misexpression of WUSCHEL.
	true: Underexpression of WUS in the expression domain
	(Reduced [WUS]_0 as in the supporting information for the main text).
	false: Wildtype calibrated values for WUS.

====================================================

Below are listed the components of submitting bulk jobs for the SAM SCE model.

The model is made to work with 
the batchmaker.cpp and CSV_interpreter programs, so compile by

make clean
make all
g++ batchmaker.cpp -o batchGenerator.out
g++ CSV_interpreter.cpp -o CSV_interpreter.out

and run it with

./CSV_interpreter.out <csv file name>.csv



====================================================

## CSV_interpreter.(cpp/out)

====================================================

Compile with:  g++ CSV_interpreter.cpp -o CSV_interpreter.out
run with: ./CSV_interpreter.out foo.csv
This file takes a CSV foo with simulation constants detailed on
line 1, variable parameters on line 2, and the values that 
those variables take on line 3 and beyond.  Each row of the 
CSV will correspond to one submitted simulation. 

NOTE:  YOU MUST COMPILE THE BATCH MAKER FIRST OR THIS WILL NOT RUN

IN: K = index of first test name, 
CSV formatted as
-batchmaker_flag_1,flag_value_1, ...,-batchmaker_flag_M,flag_value_M,-test,NAME
>Parameter flags (-P1, ... ,-PM) (M flags taken from main.cpp)
xK1,x12,...,xKP
x(K+1)1,x(k+2)2,...,x(K+1)P
...
x(K+N)1,X(K+N)2,...,x(K+N)P
>x_N1,x_N2,...,x_NM (where N is the number of simulations to be submitted.)

OUTPUT: N simulations submitted with parameter vectors (Xk1,...,XkP)
submitted via slurm as NAME_K, NAME_(K+1), ... , NAME_(K+N). 

List of batchmaker_flags (Optional marked with *)
-test <NAME>
-hours <##> 
-cores <##>
-nodes <##>
-partition <short/batch/intel>
*-mem <#> (# is in GB, default 2)
*-nodes <#> (defaults to 1).

Example usage: (Note that in this example, we started with TEST_NAME_3)

bash-4.2$ ./CSV_interpreter.out TEST.csv

Please enter integer k for row 1's index: 3
Parameter set 3
3 3
Command: ./batchGenerator.out -p batch -hours 24 -cores 12 -test TEST_NAME_3 -par -WR
0.7 -par -CK 1 -par -div 3 -par -TC 3 -flag OOP_off
Submitted batch job 539277
Parameter set 4
2 4
Command: ./batchGenerator.out -p batch -hours 24 -cores 12 -test TEST_NAME_4 -par -WR
0.7 -par -CK 1 -par -div 2 -par -TC 4 -flag OOP_off
Submitted batch job 539278 
>>>This is the console output of simulations TEST_NAME_3 and TEST_NAME_4 being submitted 
>>>via drawing parameters from TEST.csv.



====================================================

## batchmaker.cpp / batchGenerator.out

====================================================

Compile with: g++ batchmaker.cpp -o batchGenerator.out
run with: ./batchGenerator.out -batch_flags -par -parameter_flag -parameter_value -flag <flag>

IN: Flag values:
mandatory:
-p <partition>  //e.g. batch,short, etc.
-hours <int>  //Two digit. Mandatory unless minutes are given.
-cores <int>  //CPUs per task
-test <str> //Name


optional: 
-mem <int> (memory in gigabytes per node. Defaults to 2.)
-mins <int> (Number of minutes.  Can replace hours)
-nodes <int> (Number of nodes.  Defaults to 1.)
-bigdata //Sends output to /bigdata/wchenlab/shared/Plant_SCE_output/, local otherwise

OUT: sbatch script entitled AUTO_BATCH.sh

Example: ./batchGenerator.out -p batch -hours 24
-cores 12 -test TEST_NAME_4 -par -WR 0.7 -par -CK 1 -par -div 2 -par -TC 4 -flag OOP_off

Example output: AUTO_BATACH.sh containing...

 #!/bin/bash -l
 #SBATCH --nodes=1
 #SBATCH --ntasks=1
 #SBATCH --cpus-per-task=12
 #SBATCH --mem-per-cpu=2G
 #SBATCH --time=0-24:00:00
 #SBATCH --output=myTEST_NAME_4.stdout
 #SBATCH --job-name="TEST_NAME_4"
 #SBATCH -p batch 
 export OMP_NUM_THREADS 12
 mkdir Animate_Cyt_TEST_NAME_4
 mkdir Nematic_test_1
 mkdir Locations_test_1
 mkdir Animate_No_Cyt_TEST_NAME_4
 ./program Animate_Cyt_TEST_NAME_4 Locations_test_TEST_NAME_4 Nematic_test_TEST_NAME_4 Animate_No_Cyt_TEST_NAME_4 -WR 0.7 -CK 1 -div 2 -TC 4
