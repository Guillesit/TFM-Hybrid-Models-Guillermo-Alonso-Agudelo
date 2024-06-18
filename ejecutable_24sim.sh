#!/bin/bash



#Deter
cp Parameters/pop12 Fisher-Kolgomorov-20240417_NataliaBP/Probando/pop12
cp Parameters/sim_params Fisher-Kolgomorov-20240417_NataliaBP/Probando/sim_params
cp Parameters/eq_params Fisher-Kolgomorov-20240417_NataliaBP/Probando/eq_params
cd Fisher-Kolgomorov-20240417_NataliaBP


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
cp Parameters/pop12 BuenoHybrid/Probando/pop12
cp Parameters/sim_params BuenoHybrid/Probando/sim_params
cp Parameters/eq_params BuenoHybrid/Probando/eq_params
cd BuenoHybrid


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




#Laia
cp Parameters/pop12 Laia/Probando/pop12
cp Parameters/sim_params Laia/Probando/sim_params
cp Parameters/eq_params Laia/Probando/eq_params
cd Laia


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


