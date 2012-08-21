########################################################### identification of alleles of individuals from indivdual sequencing ###############################


### The following pipeline describes my strategy to extract the sire allele of F1 hybrids between an male of a population of interest an a female of the isogenic strain mel36. The input has to be a sync file with the reference strain mel36 as first population followed by the iF1 indivduals. The script extract_consensus.py is used to create an output, which loosly resembles a pileup file and contains the alleles of all individuals. The procedure is rather complex and I would recommend reading the extensive help of the script before the first usage. 
## E.g. you have a sync file with the reference strain and 11 individuals sequenced and you want to extract the sire alleles, your command line should look like this:

python /popgentools/users/martin/individual_data/extract_consensus.py --input individuals.sync --min-coverage 20 --max-coverage 0.05 --min-count 20 --CI-mode individual --output output_file

## the outputfile output_file.af can for example be used to visualize the distribution of the allele frequencies in the indivduals (which should cluster around 0.5) by plotting histograms in R like this:

echo '''
data=read.table("output_file.af")
par(mfrow=c(2,6))
for (i in 2:length(data)){
hist(data[,i])
} 
''' > output_hist.r

Rscript output_hist.r

## the outputfile: output_file.consensus can now be used for further downstream analyses:


############################################################################### FST and pi ####################################################################

## For example, I wrote a script very specifically for my inversion project, but it might be also useful for some of you. The script FST4inversions.py calculates FST and pi for different combination of individuals. In my case, I knew for each indivdual, which inversion it was carrying. Therefore, I was able to group them according to the inversion. E.g. from the 11 indivduals 4 are carrying inversion In(2L)t. Therefore, I can calculate pi for this group and for the other 7, which are not carrying the inversion. Then, I can calculate FST between these two groups similarily to PoPoolation2. For some SNPs allelic information is not available for all individuals. Therefore, you have to define the minimum number individuals for which allelic information has to be available to perform the calculations (e.g. 3). See the help in the script for further information. In our example the command line would look like this:

python /popgentools/users/martin/individual_data/FST4inversions.py --input output_file.consensus --min-count 3 --in2lt 2,3,6,10 --all 0,1,2,3,4,6,7,8,10 --names 52,53,80,100,117,136,150,168,106 --output output

## I wrote another script to bin the values of the former analysis in non-overlapping windows for a better visual representation. This script can be used either on the FST or the pi output. You will need the header of a SAM file (usually called inh.sam) to provide total length of the chromosomal arms. Again, this script is very specific, but perhaps helpful for somebody....

python /popgentools/users/martin/individual_data/binning.py --input output_fstpi.fst --length inh.sam --data fst --window-size 10000 --output output_10k_fst

## and 

python /popgentools/users/martin/individual_data/binning.py --input output_fstpi.pi --length inh.sam --data pi --window-size 10000 --output output_10k_pi

## now this data can be for example visualized using R, assuming that the FST value for In(3R)Mo is in column 3 in the output_10k_fst_values.txt and the pi values for In(3R)Mo are in column 3 and for all other indivduals in column 4 in output_10k_pi_values.txt:

echo '''
data=read.table("output_10k_fst_values.txt")
data1=read.table("output_10k_pi_values.txt")
pdf("output_10k_3R.pdf",width=15,height=10)
par(cex=2)
plot(data[data$V1=="3R",2],data[data$V1=="3R",6],type="l",main="3R",ylim=c(0,1),xlim=c(0,30000000),xlab="distance",ylab="Pi/FST",lwd=2)
points(data1[data1$V1=="3R",2],data1[data1$V1=="3R",7],type="l",col="red")
points(data1[data1$V1=="3R",2],data1[data1$V1=="3R",9],type="l",col="blue")
rect(15922589,0,28234649,1,col="transparent",border="black",lwd=3,lty=3)
rect(17058538,0,24780914,1,col="transparent",border="black",lwd=3,lty=1)
rect(12401910,0,20580518,1,col="transparent",border="black",lwd=3,lty=4)
dev.off()
''' > output_10k_3R.r

Rscript output_10k_3R.r

## here pi for In(3R)Mo is red, for all other indviduals blue, FST is black and the three inversions on 3R are black boxes with different line types

######################################################################### LD ##########################################################################################

## finally there is a script to calculate pairwise r^2 (as a measurement of LD) for a set of SNPs along a chromsomal arm. This script named LD_heatmap.py produces a tabular output of the distance matrix and PNG file containing a Heatmap of all pairwise comparisons with r^2 highlighted in color. The LDheatmap R Package needs to be installed. It can be found here: http://cran.r-project.org/web/packages/LDheatmap/index.html.see the help within the script for more details. Lets assume we want to test the LD on 2L for 4 indivduals carrying In(2L)t using 500 SNPs randomly picked along the chromosome. Then your commandline should look like this:

python /popgentools/users/martin/individual_data/LD_heatmap.py --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L








