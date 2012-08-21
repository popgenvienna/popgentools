####### this pipeline describes the detection of contamination of PoolSeq data with a foreign species. In our case we detected that pools of melanogaster flies contained simulans reads within them.

### note that rpy2 needs to be installed on the computer to run the scripts 

## at first you need to create a dataset consiting of positions being divergent between the species. This means that alleles in one species are not present in the other species (above a minimum count threshold). Cases, where one species is polymorphic for alleles, which are not found in species are also considered as divergence. For the first step, you need one or more CLEAN populations of each species mapped against the same reference and stored in the sync file format. 

## The script detect_divergence.py extracts divergent positions and stores the divergent alleles as well as the synced populations in the output. You have to define the populations of the two different species ( --species1, --species2) as well as the mincount for each of the populations per species (--species1-mincount, --species2-mincount ). Additionally you have to define a global min-coverage and max-coverage. Make sure to put the populations of the species which is causing the contamination to species2. Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH.

## E.g. if you have three populations of melanogaster (1,2,3) and three of simulans (4,5,6) , your command line would look like this: 

python /popgentools/users/martin/simulans_contamination/detect_divergence.py --input input.sync --output sim_mel.div --species1 1,2,3 --species2 4,5,6 --species1-mincount 1,1,1 --species2-mincount 1,1,1 --min-coverage 50 --max-coverage 200

## this dataset can now be used to caluclate the allele frequency of the divergent allele in the presumably contaminated dataset. Apart from flase psotives in the divergence dataset, contamination should be even across the whole genome and only sampling error due to coverage should cause variance in allele frequency. The median of the frequencies of the divergent alleles in the contaminated dataset should represent the approximate level of contamination. These frequencies can be calculated with the script compare_divergence.py

## The script compare_divergence.py needs a sync file with the contaminated population(s) and the divergence dataset as inputs. It then calculate the allele frequency of the divergent allele of the species causing the contamination for each of the populations in defined in the commandline.
## E.g. you have a sync file with four populations and expect three (1,2,4) to be contaminated. Then your commandline should look like this:

python /popgentools/users/martin/simulans_contamination/compare_divergence.py --input contaminated.sync --div output_cov100.div --pops 1,2,4 --min-count 2,2,3 --out output.af

## The results from the above script can be presented as a histograms with R:

echo '''
data=read.table("output.af")
par(mfrow=c(1,3))
for (i in 2:length(data)){
hist(data[,i])
} 
''' > output_hist.r

Rscript output_hist.r

## if you found contamination in your reads, you can use a program called bam2fastq (located inside this folder) to extract pair-end reads from your BAM files from the contaminated populations. see http://www.hudsonalpha.org/gsl/software/bam2fastq.php for more details. The commnandline looks as following: set -s to tell the program that your BAM file does not contain errors -o specifies the output file. Note that the "#" at the end will be replaced by _1 and _2 in the final output. The input BAM file needs to be put at last : 

/popgentools/users/martin/simulans_contamination/bam2fastq-1.1.0/bam2fastq -s -o output# contaminated.bam

## now you can use these fastq files as the input for Ram's python script find_simulans_reads.py , which uses gsnap to map the reads against the simulans and the melanogaster genome. The output of this script consists of four fastq files: Fwd and Rev of melanogaster specific reads and FWd and Rev of simulans specific reads. These datasets can be used to purge the originally contaminated BAM file from the simulans contamination using the script fix_bam.py

## fix_bam.py needs the pysam package to be installed. To do so type:

sudo easy_install pysam

## E.g. you have a contaminated BAM file named contamination.bam and ran Ram's script find_simulans_reads.py  using the reads from this BAM file. Your commnadline to purge the contamination.bam should look like this:

python /popgentools/users/martin/simulans_contamination/fix_bam.py --input contamination.bam --sim simulans_1.fastq --mel melanogaster_1.fastq

## you will get three outputs contamination.bam_mel, contamination.bam_sim and contamination.bam_missed. contamination.bam_mel can now be used to create a pileup and synced files. After that compare_divergence.py can be used again to test how well the purging worked