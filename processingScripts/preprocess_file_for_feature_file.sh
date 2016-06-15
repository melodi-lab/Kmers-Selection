#!/bin/bash

inputFile=$1 # specify the file to be processed. 

if [ -f ${inputFile}.ranked ]; then 
    echo "The file ${inputFile} has already been processed and move on"
    exit 0
fi


# process to obtain the list of K-mers
less $inputFile | head -n 1 | awk '{$1=""; print}' > ${inputFile}.K_mers_list



# process to obtain the list of region names
awk '{print $1}' $inputFile | tail -n +2 > ${inputFile}.region_list



# process to obtain the feature graph
tail -n +2 $inputFile  | awk '{for(idx=2;idx<=NF;idx++){if ($idx!=0) {printf idx-1 " " $idx " "}}; printf "\n"}' > ${inputFile}.featureGraph



# run the greedy algorithm to generate the rank of the genome regions

../SubmodularOptimization/greedyAlgorithmOptimization/submod_selection  -graph ${inputFile}.featureGraph -n `less ${inputFile}.region_list | wc -l` -per 100 -out ${inputFile}.ranked -ref ${inputFile}.region_list 





