#!/bin/sh
#$ -l h_rt=48:00:00
#home="/usr3/graduate/ashishge"
#path="/usr3/graduate/ashishge/programrun"

l=800   #600
b=3000  #9000
r=20
T=120000
N=Nval
g=gval
p00=p00val #p
p0=0.0 #q
p1=p1val #e
p2=p2val #n
h=hval
wmflag=1


#  try1.sh
#  2D C SIMULATIONS
#
#  Created by Ashish Bino George on 4/4/15.
#  Copyright (c) 2015 Ashish Bino George. All rights reserved.
#$ -N myjob
#$ -j y
echo Hello World 

#gcc -std=c99 /usr3/graduate/ashishge/code.c
#gcc -std=c99 -Wall -l /usr/local/include -c line_bash.c


gcc -std=c99 -Wall -l /usr/local/include -c line_bash.c
gcc -L/usr/local/lib line_bash.o -lgsl -lgslcblas -lm
./a.out -l $l -b $b -r $r -N $N -T $T -g $g -p $p00 -q $p0 -e $p1 -n $p2 -w $h -F $wmflag
#gcc usr3/graduate/ashishge/particle-hole/2D-C-SIMULATIONS/squishedball/Squishedball_ptr_ran_line_holeWF3.1.c


