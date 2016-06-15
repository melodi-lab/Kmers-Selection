# Kmers-Selection


(1) cd to directory "Kmers-Selection/SubmodularOptimization/greedyAlgorithmOptimization", and compile the greedy optimization code by typing the command  "make"
(2) cd to diretory "Kmers-Selection/processingScripts/" and run the following command to generate the ranking of the genome regions:
./preprocess_file_for_feature_file.sh ../dataFile/hg38.width\=1000000.offset\=1000000.k\=3.txt

The generated ranking of the regions will be stored in 
"dataFile/hg38.width\=1000000.offset\=1000000.k\=3.txt.ranked.ref"
