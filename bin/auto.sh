#! /bin/bash
#rm macsim
#./build.py -j5
m="PDP_WORKLOAD"

for i in $(seq 1 5)
do
	rm trace_file_list
	cp trace_file_list$i trace_file_list
	rm a.txt
	cp a"$i".txt a.txt

	if [ "$i" = "1" ];
		then
			./macsim --l3_repl==PDP --num_sim_cores=1 --num_sim_large_cores=1 <a.txt | tee -i PDP.out
			cat a"$i".txt
	fi

	
	if [ "$i" = "2" ];
		then
			./macsim --l3_repl==PDP --num_sim_cores=2 --num_sim_large_cores=2 <a.txt | tee -i PDP.out
			cat a"$i".txt
	fi

	
	if [ "$i" = "3" ];
		then
			./macsim --l3_repl==PDP --num_sim_cores=2 --num_sim_large_cores=2 <a.txt | tee -i PDP.out
			cat a"$i".txt
	fi
	
	
	if [ "$i" = "4" ];
		then
			./macsim --l3_repl==PDP --num_sim_cores=3 --num_sim_large_cores=3 <a.txt | tee -i PDP.out
			cat a"$i".txt
	fi


	if [ "$i" = "5" ];
		then
			./macsim --l3_repl==PDP --num_sim_cores=3 --num_sim_large_cores=3 <a.txt | tee -i PDP.out
			cat a"$i".txt
	fi

	mkdir $m$i
	mv *.out* $m$i
	mv *.csv $m$i
done


