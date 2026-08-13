WARNING: Please copy cimipack to your own directory. DO NOT run cimi in share directory

# to compile cimi

# make 		: compile cimi and copy cimi.exe in the run directory

# make sat	: compile cimi's satellite output program and copy sat_flux.exe in directory "run_sat"

# make clean	: remove the run directory and cimi.exe

# make rundir 	: create a run directory

# make debug	: debug mode

# make test	: compile cimi and conduct a test run for 2 minutes


## to use gfortran or ifortran (ifort), open src/Makefile and comment "ifx" and uncomment "gfort ~" or "ifort~"
