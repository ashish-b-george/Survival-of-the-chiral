#!/bin/sh
#$ -l h_rt=48:00:00
#$ -N chiral_test
#$ -j y
home="/usr3/graduate/ashishge/off_lattice/"
path="/projectnb/qedksk/ashishge/data/off_lattice/systematic_analysis/oppch_inv_wm/F0_6try4smallerK_msp5Additional/"
if [ ! -d $path ]
then
mkdir $path
fi

cp $home"off_lattice_oppch63_f0.sh" $path"off_lattice_oppch63_f0.sh"  ## just to have the file as a record
cp $home"submit_python_analysis.sh" $path"submit_python_analysis.sh"  ## just to have the file as a record
jobname="oppch63F"
parameter_name="F"
counter=0
##########starting at i=2!!
declare -i no_runs
no_runs=6
echo Hello World
max_population_size=320000
box_length=25
domain_width=500.
spread_in_concentration_sensing=1.0
migration_rate=1.0
#chirality_arr=(0.0 -0.5 -1.0 -2.0 -4.0 -8.0 -16.0 -100.0 0.5 1.0 2.0 4.0 8.0 16.0 100.0)
chirality_arr=(8.0 10.0)
#chirality=8.0
BC_flag=0
per_deme_carrying_capacity=6.
migration_step_size=0.5
sensing_flag=0
growth_rate=0.2
f0_arr=(0.2 0.3 0.4 0.5 0.6 0.7 0.8)
IC_length=5.0
wmflag=1
IC_population_size=10000
Time_steps=3000
Allee_coefficient=0.0 ### negative or 0  means no allee effect
migration_func_flag_arr=(6)

for migration_func_flag in ${migration_func_flag_arr[*]}
do
    mfdir=$path'/Mfunc'$migration_func_flag
    if [ ! -d $mfdir ]
    then
    mkdir $mfdir
    fi

for chirality in ${chirality_arr[*]}
do
    chdir=$mfdir'/ch'$chirality
    if [ ! -d $chdir ]
        then
        mkdir $chdir
        fi
for f0 in ${f0_arr[*]}
do
    Kdir=$chdir'/F'$f0
    if [ ! -d $Kdir ]
    then
    mkdir $Kdir
    fi


    for ((i=2; i<=no_runs;i++))
    do
        rundir=$Kdir'/'$i
        if [ ! -d $rundir ]
        then
        mkdir $rundir
        fi
    cd $rundir
    cp -rf $home'/off_lattice_main.cpp' $rundir'/off_lattice_main.cpp'
    cp -rf $home'/off_lattice_population.hpp' $rundir'/off_lattice_population.hpp'

    g++ -Wall -l /usr/local/include -c off_lattice_main.cpp
    #g++ -std=c++0x -L/usr/local/lib off_lattice_main.o -lgsl -lgslcblas -lm
    g++ -I/usr/local/include -L/usr/local/lib off_lattice_main.o -lgsl -lgslcblas -lm

    chirality_B=$chirality
    minus1=-1.0
    chirality_A=$(bc <<< "scale=2;$chirality*$minus1")
    qsub -N $jobname$counter -l h_rt=96:00:00 -j y -b y a.out -F "-N $max_population_size -W $domain_width -L $box_length -v $spread_in_concentration_sensing -z $migration_step_size -m $migration_rate -g $growth_rate -a $chirality_A -b $chirality_B -K $per_deme_carrying_capacity -B $BC_flag -C $sensing_flag -n $IC_population_size -f $f0 -l $IC_length -w $wmflag -T $Time_steps -A $Allee_coefficient -M $migration_func_flag"
    counter=$((counter+1))
    done
done
done
done

echo "counter is"$counter
cd $path
########## | is used instead of / in sed  so that path name can be passed on. sed is flexible with character one can use! ######
sed -i "s|folder_path|$path|" submit_python_analysis.sh
sed -i "s|parameter_name|$parameter_name|" submit_python_analysis.sh

qsub -N 'plot'$jobname -hold_jid $jobname'*' submit_python_analysis.sh



