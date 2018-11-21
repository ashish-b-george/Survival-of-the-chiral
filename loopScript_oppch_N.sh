#!/bin/sh
#$ -l h_rt=48:00:00
#$ -N myjob
#$ -j y
home="/usr3/graduate/ashishge/"
path="/projectnb/qedksk/ashishge/data/new_sim/oppch_with_N/gdnoise/try5/"

if [ ! -d $path ]
then
mkdir $path
fi


# Varies N and averages over many runs
declare -i no_runs
no_runs=1 ### number of runs for averaging


p00=0.00 ## constant diffusion/ denisty independent migration
N_arr=(50 100 200 400 800 1600 3200 6400 12800) ## carrying capacity/ inverse noise strength
p1=0.10 ## al+ar
g=0.1 ## growth probability per step
p2=0.065 ##chirality
h=0.1 ## ratio of strains

echo Hello World 
## makes directories with Nvalue as name , subdirectories numbered from 1 to no_runs has data from each simulation
for N in ${N_arr[*]}
do
newdir=$path'/N'$N
if [ ! -d $newdir ]
    then
        mkdir $newdir
fi
cd $newdir
for ((i=1; i<=no_runs;i++))
	do
	rundir=$newdir'/'$i
	if [ ! -d $rundir ]
	then
		mkdir $rundir
	fi
    ## each folder contains the code used to run the simulations and a child bash script that runs the code. (copied below)
	cp -rf $home'/model_simulations.c' $rundir'/line_bash.c'
    cp -rf $home'/qsubScript_oppch_N.sh' $rundir'/qsubS.sh'
	cd $rundir
    ## child bash script takes in parameters used to run the simulation via a string replacement operation
	sed -i "s/Nval/$N/" qsubS.sh
	sed -i "s/p00val/$p00/" qsubS.sh
	sed -i "s/p2val/$p2/" qsubS.sh
	echo $p2
    sed -i "s/p1val/$p1/" qsubS.sh
    sed -i "s/hval/$h/" qsubS.sh
    sed -i "s/gval/$g/" qsubS.sh

    echo $h
    echo $N

    qsub qsubS.sh
	cd $home	
	

done
done

echo Hello World 	






