#!/bin/sh -l
#$ -l h_rt=10:00:00
#$ -j y
home="/usr3/graduate/ashishge/off_lattice/"
path="/projectnb/qedksk/ashishge/data/off_lattice/systematic_analysis/oppch_inv_wm/F0_6try4smallerK_msp5Additional/"
parameter="F"
module load python
echo $path


echo python "$home"off_lattice_cluster_plot.py -i "$path"
python "$home"off_lattice_cluster_plot.py -i "$path"

python "$home"off_lattice_parameter_analysis.py -i "$path" -p "$parameter"

#echo python $home off_lattice_cluster_plot.py -i $path
#python $home off_lattice_cluster_plot.py -i $path



#python $home"off_lattice_cluster_plot.py" -i "$path"
#python "$home"off_lattice_cluster_plot.py -i "$path"
