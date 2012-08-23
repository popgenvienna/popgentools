### Allele frequencies, etc.

## this is a collection of scripts which allow you to calculate allele frequencies and allele frequency changes. These scripts were specifically written for time series data. Therefore, they might be too specific for general purposes, but anyway, I put them here and explain how they work!:


###################################################### AF.py #############################################################


## AF.py can be used to calculate the allele frequencies (AF) of several populations based on the allele which is increasing beteween two timepoints. E.g you have a sync file input.sync with three populations at three timepoints 0,15 and 27 and you want to know the AF of the allele rising between 0 and 15 you need to put following commandline: 

python /popgentools/users/martin/Allele_frequencies/AF.py --input input.cmh --allpop 1,2,3 --start 1 --end 2 --replicate 1,2,3 --names 0,15,27 > output.af

## here, the program uses all pops --allpop 1,2,3 to detect the two alternative alleles (or the two most common alternative alleles), picks the allele rising between 1 and 2 and prints the frequencies of this allele for the populations 1,2,3 with the names 0,15 and 27. Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH. If you have replicates for each timepoint, you can incorporate them for the detection of the alleles and to pick the allele between the timepoints. e,g. Three Replicates for 0, 15 and 27:

python /popgentools/users/martin/Allele_frequencies/AF.py --input input.cmh --allpop 1,2,3,4,5,6,7,8,9 --start 1,2,3 --end 4,5,6 --replicate 1,4,7 --names 0,15,27 > output.af


###################################################### AF_threshold.py #############################################################


## AF_threshold.py only prints SNPs for which the major allele is above a certain allele frequency threshold in a given (or multiple) population(s). The Input is a sync or a cmh file.:

python /popgentools/users/martin/Allele_frequencies/AF_threshold.py  --input input.sync --pops 1,2,3 --threshold 0.95 > input_095.cmh

###################################################### AFC.py #############################################################


## AFC.py is similar to AF.py and has the same command line. It can be used to calculate the allele frequencies changes (AFC) between several populations based on the allele which is increasing beteween two timepoints. Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH. E.g you have a sync file input.sync with three populations at three timepoints 0,15 and 27 and you want to know the all pairwise differences in allele frequency of the allele rising between 0 and 15 you need to put following commandline: 

python /popgentools/users/martin/Allele_frequencies/AFC.py --input input.cmh --type CMH --allpop 1,2,3 --start 1 --end 2 --replicate 1,2,3 --names 0,15,27 > output.afc

## here, the program uses all pops --allpop 1,2,3 to detect the two alternative alleles (or the two most common alternative alleles), picks the allele rising between 1 and 2 and substracts the AF of two populations in all combinations (In this case: 15-0, 27-0 and 27-15). If you have replicates for each timepoint, you can incorporate them for the detection of the alleles and to pick the allele between the timepoints. e,g. Three Replicates for 0, 15 and 27:

python /popgentools/users/martin/Allele_frequencies/AFC.py --input input.cmh --type CMH --allpop 1,2,3,4,5,6,7,8,9 --start 1,2,3 --end 4,5,6 --replicate 1,4,7 --names 0,15,27 > output.afc


###################################################### AFC_cov.py #############################################################

## Similarily to the above scripts AFC_cov.py uses all pops --allpop 1,2,3 to detect the two alternative alleles (or the two most common alternative alleles) and picks the allele rising between the start and end. However, in contrast to AFC.py this script calculates the average AFC between the start and end populations (if multiple replicates are available) and appends it at the end of each line. Additionally it appends the cummulative coverage of the start and end populations at the end. 

python /popgentools/users/martin/Allele_frequencies/AFC_cov.py  --input input.sync --all 1,2,3,4,5,6,7,8,9 --start 1,2,3 --end 7,8,9 > input.ac
