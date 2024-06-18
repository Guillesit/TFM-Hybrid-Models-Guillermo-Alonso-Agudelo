#Stoch
cp Parameters/pop12 F-K\ Stoch/Probando/pop12
cp Parameters/sim_params F-K\ Stoch/Probando/sim_params
cp Parameters/eq_params F-K\ Stoch/Probando/eq_params
cd F-K\ Stoch/


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

