### prepare data for a SNPeff annotation 

## in this folder you will find a copy of SNPeff 2.0.3, with the annotation of D. melanogaster 5.40. This version is tested and working.

# at first you will need to create a VCF file from a sync (e.g. CMH test output), which is needed as an input for SNPeff:

python /popgentools/users/martin/SNPeff/sync2vcf.py --sync input.sync > input.vcf

# then you need to run SNPeff. You will get two outputs: 1) a *.snpeff file with the annotations and 2) a *.htmal file with a nice visual summary. With -ud you can define how many basepairs up- and downstream will be included to a gene as UPSTREAMS and DOWNSTREAMS. with -s you define the output summary file and with -o TXT you define that the annotation output is in TEXT style rather than VCF style. This needs to be set to TXT, otherwise the output cannot be parsed by the downstream software :

cd popgentools/users/martin/SNPeff/snpEff_2_0_3

java -Xmx10g -jar ./snpEff.jar dm5.40 input.vcf -inOffset 1 -outOffset 1 -ud 200 -s input_SUMMARY.html -o TXT > input.snpeff

# with link_snpeff.py this annotation file can be used to append the annotation as the last columns to any text file with at least two columns defining the chromsome and the position. See the help of link_snpeff.py for more details. By convention I use the extension *.genes for the output.

python /popgentools/users/martin/SNPeff/link_snpeff.py --snpeff input.snpeff --input candidates.sync > candidates.genes 