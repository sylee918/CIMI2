default: CIMI

current_dir="$(realpath .)"
share_dir="/data/share/cimi/cimipack_babyblue"
rundir_path=$(current_dir)'/run/'
inputdir_path=$(current_dir)'/inputs/'
runsatdir_path=$(current_dir)'/run_sat/'

CIMI:
	make check
	@if test -d $(rundir_path); then\
		echo "the run directory exists";\
	else \
		make rundir; \
	fi
	@cd src; make
	cp src/cimi.exe $(rundir_path)
	# create .exe that calculate CIMI's output along the satellite position
	make sat

test:
	make check
	@if test -d $(rundir_path); then\
		echo "the run directory exists";\
	else \
		make rundir; \
	fi
	@cd src; make
	cp src/cimi.exe $(rundir_path)
	cd $(rundir_path); rm ./cimi.dat
	cd $(rundir_path); cp ./cimi.dat.test ./cimi.dat
	cd $(rundir_path); ./cimi.exe 	

sat:
	make check
	@if test -d $(runsatdir_path); then\
		echo "the run_sat directory exists";\
		rm $(runsatdir_path)*.fls; ln -s $(rundir_path)*.fls $(runsatdir_path) ;\
		rm $(runsatdir_path)*.pot; ln -s $(rundir_path)*.pot $(runsatdir_path) ;\
		rm $(runsatdir_path)*.rtp; ln -s $(rundir_path)*.rtp $(runsatdir_path) ;\
	else \
		make runsatdir; \
	fi
#	@if test -f $(runsatdir_path)sat_flux.exe ; then\
#		echo "sat_flux.exe already exists";\
#	else \
#		cd src; make sat_compile; \
#	fi
	cd src; make sat_compile
	cp src/sat_flux.exe $(runsatdir_path)
	cp $(runsatdir_path)/* $(rundir_path)
	rm -rf $(runsatdir_path)

rundir:
	make check
	rm -rf $(rundir_path)
	mkdir $(rundir_path)
	cp -r $(inputdir_path)/* $(rundir_path)
	rm $(rundir_path)/D_*.dat 
	cd $(rundir_path); ln -s $(inputdir_path)/D_* ./
	rm $(rundir_path)/cimi.dat
	cd $(rundir_path); ln -s ./cimi.dat.mate ./cimi.dat
	cp tools/*.pro $(rundir_path)
	mkdir $(rundir_path)/src_cimi
	cp src/cimi.f90 $(rundir_path)/src_cimi
	cp src/plasmasphere.f90 $(rundir_path)/src_cimi
	cp src/sat_flux.f90 $(rundir_path)/src_cimi
	#inputs/* tools/*.pro run

runsatdir:
	make check
	rm -rf $(runsatdir_path)
	mkdir $(runsatdir_path)
	cp -r $(inputdir_path)/*.info $(runsatdir_path)
	cp -r $(inputdir_path)/*.pos $(runsatdir_path)
	cp -r $(inputdir_path)/sat_flux.dat $(runsatdir_path)
	#cd $(runsatdir_path); ln -s $(rundir_path)*.fls ./
	#cd $(runsatdir_path); ln -s $(rundir_path)*.pot ./
	#cd $(runsatdir_path); ln -s $(rundir_path)*.rtp ./
	cp tools/plot_sat_flux.pro $(runsatdir_path)
	cp tools/map_sat_flux.pro $(runsatdir_path)


debug:
	make check
	@cd src; make debug
	echo " ... a debug executable was created, cimi_debug.exe" 
	cp src/cimi_debug.exe $(rundir_path)
	echo " ... a debug executable is running" 
	cd $(rundir_path); ./cimi_debug.exe 	

EMIC:
	make check
	make rundir
	cp src/ModIonDiff.f90 srcTest/
	@cd srcTest; make EMIC_DIFFUSE_TEST
	cp srcTest/unit_test_diffuse_ion.exe $(rundir_path)
	rm srcTest/ModIonDiff.f90
	cd $(rundir_path); ./unit_test_diffuse_ion.exe	

# check if the code is comiled in the share directory.
# this check works only in gs673-babyblue2
check:
	@if [[ $(share_dir) == $(current_dir) ]]; then\
		echo " ";\
		echo "** Error: ";\
		echo "** the code cannot be compiled in the share directory.";\
		echo "** copy $(share_dir) to your home directory and compile the code.";\
		echo " ";\
		exit 1; \
	fi


clean:
	rm -rf $(rundir_path)
	rm -rf $(runsatdir_path)
	@cd src;	make clean
	@cd srcTest;	make clean
