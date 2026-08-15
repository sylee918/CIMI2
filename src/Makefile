default: CIMI

CIMI:
	make EXE

CIMI_NO_FLC:
	make EXE_NO_FLC

MODULES=\
	ModCurvScatt.f90 \
	ModIonDiff.f90 \
	ModTwinsHden.f90 \
	ModMate.f90 \
	w05sc_fok.f90

OBJECT=\
       	plasmasphere.f90 \
       	T96_01.f TS04c.f \
       	geopack_2008.f \
       	Lstar2.f90 \
	rcm6.f90 \
	gmresm.f90 

ifx=ifx -r8 -132 -save -O3 -xHost -qopenmp -mcmodel=medium -Wl,--no-relax 
ifx1=ifx -r8 -132 -save -O3 -xHost -qopenmp
ifort=ifort -r8 -132 -save -O3 -xHost -qopenmp -mcmodel=medium 
#ifort=ifort -r8 -132 -save -O3 -qopenmp
debug_ifx=ifx -r8 -132 -save -O0 -traceback -g -check -debug -qopenmp -mcmodel=medium -Wl,--no-relax
debug_ifort=ifort -r8 -132 -save -O0 -traceback -g -check -debug -mcmodel=medium
gfort=gfortran -fdefault-real-8 -ffixed-line-length-132 -fno-automatic -fno-align-commons -O4 


FLC_TEST:
	@cp ModCurvScatt.f90 ./src_unit_test/
	@cd src_unit_test;	make EXE_FLC
	@cd src_unit_test;	rm ModCurvScatt.f90
	cp src_unit_test/unit_test_curv_radius.exe ./run_test/
	cp src_unit_test/unit_test_curv_scatt.exe ./run_test/
	cp tools/plot_d_curv_scatt.pro ./run_test/
	@cd run_test;	./unit_test_curv_radius.exe
	@cd run_test;	./unit_test_curv_scatt.exe

EXE:
	${ifx} -o cimi.exe ${MODULES} cimi.f90 ${OBJECT}
	#${ifort} -o cimi.exe ${MODULES} cimi.f90 ${OBJECT}
	#${gfort} -o cimi.exe ${MODULES} cimi.f90 ${OBJECT}

sat_compile:
	${ifx1} -o sat_flux.exe sat_flux.f90

debug:
	${debug_ifx} -o cimi_debug.exe ${MODULES} cimi.f90 ${OBJECT}
	#${debug_ifort} -o cimi_debug.exe ${MODULES} cimi.f90 ${OBJECT}

# EMIC wave diffusion test for H
EMIC_DIFFUSE_TEST:
	${ifx} -o unit_test_diffuse_ion.exe ModIonDiff.f90 unit_test_diffuse_ion.f90

clean:
	rm -rf *.mod
	rm -rf *.exe
