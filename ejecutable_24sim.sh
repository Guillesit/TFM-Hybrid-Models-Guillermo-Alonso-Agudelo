#!/bin/bash



#Deter
cp Parameters/pop12 Deterministic/Probando/pop12
cp Parameters/sim_params Deterministic/Probando/sim_params
cp Parameters/eq_params Deterministic/Probando/eq_params
cd Deterministic


	mkdir 1
	cp EXE 1/EXE
	cp Probando/pop12 1/pop12
	cp Probando/sim_params 1/sim_params
	cp Probando/eq_params 1/eq_params
	# cp test.sh $i/test.sh
	cd 1;
	mkdir OutputValuesOxygen
    mkdir OutputValuesPopulation
	./EXE pop12 pop12 sim_params eq_params a &
	cd ..
	cd ..
	wait
	



#Hybrid
cp Parameters/pop12 BasicHybrid/Probando/pop12
cp Parameters/sim_params BasicHybrid/Probando/sim_params
cp Parameters/eq_params BasicHybrid/Probando/eq_params
cd BasicHybrid


	i=1
	while [ $i -lt 11 ]
     do
        mkdir $i
		cp EXE $i/EXE
		cp Probando/pop12 $i/pop12
		cp Probando/sim_params $i/sim_params
		cp Probando/eq_params $i/eq_params
		# cp test.sh $i/test.sh
		cd $i
		mkdir OutputValuesOxygen
    	mkdir OutputValuesPopulation
		./EXE pop12 pop12 sim_params eq_params a &
		cd ..
		i=$(($i+1))

     done
	cd ..
	wait




#FreeBoundaryHybrid
cp Parameters/pop12 FreeBoundaryHybrid/Probando/pop12
cp Parameters/sim_params FreeBoundaryHybrid/Probando/sim_params
cp Parameters/eq_params FreeBoundaryHybrid/Probando/eq_params
cd FreeBoundaryHybrid


	i=1
	while [ $i -lt 11 ]
     do
        mkdir $i
        cp EXE $i/EXE
        cp Probando/pop12 $i/pop12
        cp Probando/sim_params $i/sim_params
        cp Probando/eq_params $i/eq_params
        # cp test.sh $i/test.sh
        cd $i
        mkdir OutputValuesOxygen
        mkdir OutputValuesPopulation
        ./EXE pop12 pop12 sim_params eq_params a &
        cd ..
        i=$(($i+1))

     done
	cd ..
	wait


