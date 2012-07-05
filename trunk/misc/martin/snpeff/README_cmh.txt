If you use the output of the cmh test remove the last column (the p-values), e.g. by printing all columns except the last with auk. The output will be a normal sync file (I would suggest to use *.sync as a file ending, e.g. full_data.sync)
or 

README:

Description:

If your candidates are base on the CMH script of Ram: READ THIS !!!!


This analysis pipeline allows to identify the effects of SNPs in the genome using a program called "SNPeff". The pipeline was written for the analysis of D. melanogaster populations, but can be also used for D. simulans (please ask me in this case!). If SNPeff is not installed on your computer go to: http://snpeff.sourceforge.net/ and download the latest version. No further installation is necessary (if you intend to use the preinstalled genomes, e.g. dmel5.31). To start the pipeline it is necessary to have an *.igv file with all SNPs identified in the observed populations (e.g. the list of SNPs extracted from the output of Roberts FST script). Additionally, one needs another *.igv file containing the candidate SNPs of interest. Furthermore, a synchronized file containing all positions of the genome is needed.

## I would STRONGLY recommend to read the descriptions of the scripts before using them by typing: python script.py -h ##

1) As a first step, you have to extract SNP positions from the CMH output, by combining the "discarded.cmh" and the "considered.cmh" if you used a p-value filter in the cmh-test-v2.pl or by using the full cmh output if you did not set an p-value for the cmh test. 
You will need to extract a sync file by removing the last column (awk '{print $1"\t"$2… [until column before last column]}SNPs.cmh > SNPs.sync). Furthermore, I would suggest to make another igv file by removing the sync columns and keeping the p-values (e.g. awk '{print $1"\t"$2-1"\t"$2"\tfeature\t"$x}' SNPs.cmh > SNPs.igv [where "$x" is the last column ]). 

Repeat the whole procedure for your candidates. (candidates.cmh > candidates.sync and canidates.cmh > candidates.igv)


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