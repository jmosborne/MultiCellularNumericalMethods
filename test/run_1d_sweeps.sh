#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestCellSorting.hpp
#

min_step=6;
max_step=12;

methods[0]="FE";
methods[1]="MP";
methods[2]="RK3";
methods[3]="RK4";
#methods[0]="BE_Tol_5";
#methods[1]="AM2_Tol_5";
#methods[2]="BE_Tol_10";
#methods[3]="AM2_Tol_10";



compressions[0]=1.0;
#compressions[1]=0.5;
#compressions[2]=0.9;

simulation_type="1d_compression";

end_time=50;
ccd=12.0; 
is_stochastic=1;

min_amt=0;

for (( j=0 ; j<${#compressions[*]} ; j++))
do
	echo "compression " ${compressions[$j]};

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
														-seed_range_lower 0 \
														-seed_range_upper 0 \
														-method ${methods[$k]} \
														-min_amt $min_amt \
														-ccd $ccd \
														-stochastic $is_stochastic \
														-end_time $end_time \
														-compression ${compressions[$j]} \
														$> /home/chaste/Chaste/projects/MultiCellularNumericalMethods/test/output/${simulation_type}_Run_${i}_${compressions[$j]}_${methods[$k]}_Output.txt 2>&1 &
	done
done

echo "Jobs submitted"
