README:

Description:

If your candidates are base on the FST script of Robert: READ THIS !!!!


This analysis pipeline allows to identify the effects of SNPs in the genome using a program called "SNPeff". The pipeline was written for the analysis of D. melanogaster populations, but can be also used for D. simulans (please ask me in this case!). If SNPeff is not installed on your computer go to: http://snpeff.sourceforge.net/ and download the latest version. No further installation is necessary (if you intend to use the preinstalled genomes, e.g. dmel5.31). To start the pipeline it is necessary to have an *.igv file with all SNPs identified in the observed populations (e.g. the list of SNPs extracted from the output of Roberts FST script). Additionally, one needs another *.igv file containing the candidate SNPs of interest. Furthermore, a synchronized file containing all positions of the genome is needed.

## I would STRONGLY recommend to read the descriptions of the scripts before using them by typing: python script.py -h ##

1) As a first step, you have to extract SNP positions from the sync file. You can use an IGV file filtered for SNP position (e.g. from the output of Robert's FST_sliding.pl script, t) or by extracting all SNPs from a sync file with Roberts snp_frequency-diff.pl script. 
a) Therefore use your full sync file as an input and type: 

perl SNP-frequency-diff.pl --input Full.sync --output-prefix SNPs --min-count 2 --min-coverage 10

b)Use the output with ending "_rc" as the input for igv2sync.py:

usage: 

python igv2sync.py -i SNPs.igv -s Full.sync -t igv > SNPs.sync
or 
python igv2sync.py -i SNPs_rc -s Full.sync -t snp_diff > SNPs.sync

The resulting sync file is a subset of the total sync file and will only contain positions present in the *.igv file

2) Next you will need to extract the reference and alternative alleles from the sync file and produce a proper input for SNPeff. (read the description of sync2SNPs.py by typing -h for further information)

usage: 

python sync2SNPs.py -s SNPs.sync > SNPs.SNPs

3) go to the SNPeff directory and run SNPeff on your input of all SNPs in the analysis. This will produce a huge output providing information of SNP effects of all SNPs in your dataset. Although it is possible to run the analysis also on smaller datasets (e.g. only your candidate SNPs) I would recommend to run it on all, as you will be able to reuse this output for every subset of SNPs afterwards.

usage:

java -Xmx10g -jar ./snpEff.jar dmel5_31 SNPs.SNPs (-if 1 -of 1)* -upDownStreamLen 200 > SNPs.snpeff

* if the version of SNPeff is below 1.7 also put these parameters

4) in the next step you will extract the SNPeff information for the SNPs of your interest and merge the data. (read the description of link_snpeff.py by typing -h for further information). The information about effects (synonymous, nonsynonymous, stop_codon, etc.) will be appended at the end of each line. If you want to compare the effects of your candidate SNPs to all SNPs in the dataset, you have to repeat this step for all SNPs.

usage: 

python link_snpeff.py -s SNPs.snpeff -i candidates.igv > candidates.genes
python link_snpeff.py -s SNPs.snpeff -i SNPs.igv > SNPs.genes

5) compare the effects and distribution of the SNPs in the candidates and the total number of SNPs with snp_summary.py. This will produce summary tables and genelists which can further used for GO Term analyses, etc. (read the description of snp_summary.py by typing -h for further information)


usage: 

python snp_summary.py -i SNPs.genes -j candidates.genes -o out -v yes -g yes -f yes 

depending on the options you have chosen you will get several outputs (out.summary, out.genelist1, out.genelist2 and out.overlap)

enjoy!