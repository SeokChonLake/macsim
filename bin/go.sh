#! /bin/bash

TRACE_ROOT=/home/jinku/macsim_traces

ip="160810"

declare -a CPU
CPU=(
		dealII
		calculix
	)



declare -a WRKL
WRKL=(
		16MB
     )

#$k=1

for i in ${CPU[*]}
do
#	$k=$k+1
	if [ "$i" = "calculix" ]
		then	
			cp a2.txt a.txt
	fi


	for j in ${WRKL[*]}
	do
#		TRACE_PATH_C=$TRACE_ROOT/$i/trace.txt
	
#		mkdir $i
#		rm trace_file_list
#		touch trace_file_list
#		echo "1" >> trace_file_list
#		echo $TRACE_PATH_C >> trace_file_list
#		mv trace_file_list $i
#		./macsim --l3_repl==PDP <a.txt  | tee -i PDP.out
	
#		mv *.out $i
#		mv *.csv $i







		TRACE_PATH_C=$TRACE_ROOT/$i/trace.txt
		TRACE_PATH_I=$TRACE_ROOT/$ip/$j/trace.txt
	
		mkdir $i"_"$j
		rm trace_file_list
		touch trace_file_list
		echo "2" >> trace_file_list
		echo $TRACE_PATH_C >> trace_file_list
		echo $TRACE_PATH_I >> trace_file_list
#		mv trace_file_list $i
		
		./macsim --l3_repl==PDP < a.txt | tee -i PDP.out
	
		mv *.out $i"_"$j
		mv *.csv $i"_"$j
	done

done




