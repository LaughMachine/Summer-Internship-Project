#!/bin/bash

for i in 20 50 80; 
    do for j in 2 5 8 10 12; 
        do for k in 92 96; 
            do for l in 'seq 11 20'; 
                sge_run --grid_submit=batch --grid_quiet --grid_mem=30G "./multi_thread_sim_2class_1.py experiment_parameters.csv ${j} ${i} .${k} ${l}";
            done
        done
    done
done