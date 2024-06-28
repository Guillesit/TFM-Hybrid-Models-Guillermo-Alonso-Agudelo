i=1
while [ $i -lt 24 ]
     do
        mkdir $i
	cp EXE $i/EXE
	cp Probando/pop12 $i/pop12
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
	sleep 2
     done
