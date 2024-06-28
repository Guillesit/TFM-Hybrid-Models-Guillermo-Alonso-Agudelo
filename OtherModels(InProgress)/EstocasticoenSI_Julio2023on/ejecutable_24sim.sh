i=21
while [ $i -lt 24 ]
     do
        mkdir $i
	cp EXE $i/EXE
	cp Probando/pop12 $i/pop12
	cp Probando/oxy12 $i/oxy12
	cp Probando/pars_pop $i/pars_pop
	cp Probando/pars_sim $i/pars_sim
	# cp test.sh $i/test.sh
	cd $i
	mkdir OutputValuesOxygen
    	mkdir OutputValuesPopulation
	./EXE pop12 oxy12 pars_pop pars_sim a &
	cd ..
	i=$(($i+1))
	sleep 2
     done
