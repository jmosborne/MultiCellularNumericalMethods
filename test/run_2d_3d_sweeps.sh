#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestCellSorting.hpp
#

min_step=6; #6
max_step=12; #18

methods[0]="FE";
methods[1]="RK4";
# methods[2]="RK3";
# methods[3]="MP";
# methods[0]="BE_Tol_5";
# methods[3]="AM2_Tol_5";
# methods[1]="BE_Tol_10";
# methods[7]="AM2_Tol_10";

#simulation_type="2d_node";
#simulation_type="2d_mesh";
simulation_type="2d_vertex";
#simulation_type="3d_node";

compression=1.0; #so no compression
ccd=12.0;
end_time=50;
is_stochastic=1;

min_amt=0;

num_sims=1;
num_sweeps=1;

for (( i=0 ; i<${num_sweeps} ; i++))
do
	start_sim=`expr $i \* $num_sims`;
	end_sim=`expr $start_sim + $num_sims - 1`
	echo $start_sim
	echo $end_sim
	echo $end_time

	for (( k=0 ; k<${#methods[*]} ; k++))
	do
		echo "method " ${methods[$k]};
		# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
		# ">" directs std::cout to the file.
		# "2>&1" directs std::cerr to the same place.
		# "&" on the end lets the script carry on and not wait until this has finished.
		nice -20 /home/chaste/lib/projects/MultiCellularNumericalMethods/test/TestNumerics \
														-simulation_type $simulation_type \
														-step_range_lower $min_step \
														-step_range_upper $max_step \
														-seed_range_lower $start_sim \
														-seed_range_upper $end_sim \
														-method ${methods[$k]} \
														-min_amt $min_amt \
														-ccd $ccd \
														-stochastic $is_stochastic \
														-end_time $end_time \
														-compression $compression \
														$> /home/chaste/Chaste/projects/MultiCellularNumericalMethods/test/output/${simulation_type}_Run_${i}_${methods[$k]}_Output.txt 2>&1 &
	done
done

echo "Jobs submitted"
