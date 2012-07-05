README:

Description: 

The aim of this analysis pipeline is to test, whether SNPs in the proximity of candidate loci (showing strong allele frequency differences between compared populations) are showing a similar pattern even though the signal (p-value from FET or CHM) is too low for being included into the candidate SNP list. 
The pipeline identifies the proximate SNPs up- and downstream of candidate SNPs (from FET or CHM), identifies haplotypes from reads containing both the adjacent and the candidate SNPs and allows the comparison of haplotype changes between populations. 

Workflow:

## I would STRONGLY recommend to read the descriptions of the scripts before using them by typing: python script.py -h ##

1) As a first step one needs an *.igv file for all SNPs identified in the observed populations (e.g. the list of SNPs extracted from the output of Roberts FST script). Additionally, one needs another *.igv file containing the candidate SNPs of interest. For the first script, it is important, that the three first columns in *.igv files are comprised by "Chormosome", "Start" and "End". All other entries can be individual. Extract_haplo.py will extract haplotypes comprised of candidate SNP and adjacent SNP if the adjacent SNP is located within the range defined in the script parameters. The longer the range, the less reads will contain both SNPs and thus the coverage for those "long" haplotypes can become quite low. If there is no SNP within the range, the script will only output the candidate position.

usage:

python extract_haplo.py --cand cand.igv --snp snp.igv --range 50 --out haplo

The script will output two files: haplo_u and haplo_d. haplo_u contains haplotypes where the candidate is located upstream to the next SNP, and the other one haplo_d, the other case, where the candidate is downstream. 

2) The second script parses the sam file and identifies reads, where both SNPs in the haplotypes from the input haplo_u or haplo_d are covered. The haplotypic states are recorded and two output files are created which contain counts and frequencies of all identified states/haplotype. This script has to be repeated for all populations, which should be compared.

usage: 

python find_haplo.py --inp haplo --sam Basepop.sam --out Basepop
python find_haplo.py --inp haplo --sam selectedpop1.sam --out selectedpop1
python find_haplo.py --inp haplo --sam selectedpop2.sam --out selectedpop2

This script also produces two outputs for the _u and the _d haplotypes, respectively

3) The last script allows to compare the different populations and test whether haplotype frequencies are significantly different between the populations by applying a Fisher's exact test. You have to repeat the script for the _u and the _d files.

## sync_haplo.py requires a special python module, which allows to use R from within python. You need to install the module by downloading the *.tar file from
http://sourceforge.net/projects/rpy/files/rpy2/ and installing it in Terminal by going to the directory where the compressed file was unpacked, e.g.:

cd ~/Downloads/rpy2
sudo python setup.py install

##

usage:

python --hap Basepop_u,selectedpop1_u,selectedpop2_u --inp haplo_u -o haplotypes_u.out
python --hap Basepop_d,selectedpop1_d,selectedpop2_d --inp haplo_d -o haplotypes_d.out

4) in order to compare the number of significant haplotypes for pairwise comparisons and different contigs, use haplo_summary.py to produce summary tables. (read the description by typing -h to see a comprehensive description) 

usage: 

python haplo_summary.py -i haplotypes_out_u -j haplotypes_out_d > haplotypes.summary

enjoy!

