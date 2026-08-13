# compile options
ifx=ifx -r8 -132 -save -O3 -mcmodel=medium -Wl,--no-relax 
ifx1=ifx -r8 -132 -save -O3 
#ifort=ifort -r8 -132 -save -O3 -qopenmp
debug_ifx=ifx -r8 -132 -save -O0 -traceback -g -check -debug -mcmodel=medium -Wl,--no-relax
gfort=gfortran -fdefault-real-8 -ffixed-line-length-132 -fno-automatic -fno-align-commons -O4 


# EMIC wave diffusion test for H
EMIC_DIFFUSE_TEST:
	${ifx} -o unit_test_diffuse_ion.exe ModIonDiff.f90 unit_test_diffuse_ion.f90

clean:
	rm -rf *.mod
	rm -rf *.exe
