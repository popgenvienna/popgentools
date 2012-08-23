## one of the key results in our experimental evolution paper was the observation that CMH P-Values and Allele frequencies decay rapidly around candidate SNPs. 


## The script used for this is called LD_dist.py. This script averages (or calculates the median) P-Values or Allele frequencies of SNPs adjacent to candidates up to a defined distances in bins of a defined sizes and then averages these bins across all candidates. See help witin the script for more information.

python /popgentools/users/martin/P-value_decay/LD_dist.py --candidates candidates.cmh --snps full.cmh --bin-size 10 --maxdist 1000 --measure median > candidates.ld

