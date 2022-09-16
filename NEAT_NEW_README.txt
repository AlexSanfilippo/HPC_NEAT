#######README#######
26 August, 2022.  Alexander Sanfilippo
updated 16th Sept., 2022
HPC MSc, Trinity College Dublin

#final readme.  Need specifics on setting up miniconda, conda environment, etc.  


FIRST TIME SET-UP
After downloading the files from the github page, it would be best to move them to their
own folder in a directory of your chosing.  Next, a python environment will need to be set
up that can handle several modules including tensorflow, keras, numpy, and gym.  I used
miniconda to set this up, and would reccomend using the same methods described below.
After installing miniconda, create a new environment.
	$conda create --name myenv
	Where  "myenv" is a name of your choosing.  
Next, open the conda environment by entering
	$conda activate myenv
The following python modules (and tensorflow)  are needed and should be 
installed in the following order:
	1. numpy
	2. tensorflow (2.5 on chuck)
	2. keras
	3. gym
	4. pydot*
	5. pydotplus*
	6. python-dev


*pydot and pydotplus are optional to the main program  as they only provide support for 
creating node graphs of the neural networks.

LINKING
In order to use a conda environemnt alongside our C++ program, we must set the LD_LIBRARY_PATH
variable to the location of the library of the  environment we have created.  For me this address was:

	/home/users/mschpc/2021/sanfilia/miniconda3/envs/summer_project/lib/  
and so I needed to call
	$export LD_LIBRARY_PATH=/home/users/mschpc/2021/sanfilia/miniconda3/envs/summer_project/lib/
in order to tell the linker where to find the proper python library.  

COMPILATION
Once the conda environment is created and activate, a couple more commands must be called. 
	$module load cports
	$module load intel
Next, the program can be compiled by simply calling the make file:
	$make

EXECUTION
Once compiled, the program can be called from the command line with the following command:
	$mpirun -np <N> ./testingGrounds
"N" is a variable that stands for the number of processors to run the program on.  Do not include the 
"<>" characters.  



Note:  After initial setup, when re-logging onto the server the falling four commands must be re-entered: 
	$conda activate myenv
	$export LD_LIBRARY_PATH=<Your Path>
	$module load cports
	$module load intel

Cleaning:
	the make file can also clean the executable (and other associated files) by running: 
	$make clean
