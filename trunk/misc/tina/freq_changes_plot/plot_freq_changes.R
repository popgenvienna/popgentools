args = commandArgs()

options("show.error.messages" = TRUE)

counter = 1
arg =  args[1]

max = length(args)


inFile=""
pop_freq_reference=c()
pop_freq_treatment1=c()
pop_freq_treatment2=c()
CMH_columns=c()
CMH_ref_1=0
CMH_ref_2=0
CMH_1_2=0

cutoff = 0.05
col_base = hsv(1,0,0)

col_treatment1_sig = hsv(1,1,0.7)
col_treatment2_sig = hsv(1,1,0.4)

col_treatment1_nonSig = hsv(1,0,0.9)
col_treatment2_nonSig = hsv(1,0,0.6)

chromosomes=c("X", "2L", "2R", "3L", "3R", "4")

printNonSig = 0
width.png = 72*30;
height.png = 72*20;

xlab=""
ylab=""

ylim=c(-0.35,0.6)
#ylab="distance of light/dark freq from reference"

cex_text=1.7;

highlight=0;
highlightedGenes = data.frame();
index = 1;

while (counter <= max){
	if (arg == "--input"){
		counter=counter+1;
		arg = args[counter];
		inFile = arg; 		
	}
	if ( arg == "--pop-freq-reference"){
		counter=counter+1;
		arg = args[counter];
		pop_freq_reference = as.numeric(strsplit(arg, ",")[[1]])
	}
	if ( arg == "--pop-freq-treatment1"){
		counter=counter+1;
		arg = args[counter];
		pop_freq_treatment1 = as.numeric(strsplit(arg, ",")[[1]])
	}
	if ( arg == "--pop-freq-treatment2"){
		counter=counter+1;
		arg = args[counter];
		pop_freq_treatment2 = as.numeric(strsplit(arg, ",")[[1]])
	}
	
	if (arg == "--CMH-columns"){
		counter=counter+1;
		arg = args[counter];
		CMH_columns = as.numeric(strsplit(arg, ",")[[1]])		
		CMH_ref_1= CMH_columns[1]
		CMH_ref_2= CMH_columns[2]
		CMH_1_2= CMH_columns[3]
	}
	
	if (arg == "--cutoff"){
		counter=counter+1;
		arg = args[counter];
		cutoff = as.numeric(arg); 		
	}
	if (arg == "--col-base"){
		counter=counter+1;
		arg = args[counter];
		col_base = arg; 		
	}
	if (arg == "--col-treatment1"){
		counter=counter+1;
		arg = args[counter];
		col_treatment1_sig = arg; 		
	}
	if (arg == "--col-treatment2"){
		counter=counter+1;
		arg = args[counter];
		col_treatment2_sig = arg; 		
	}
	if (arg == "--chromosomes"){
		counter=counter+1;
		arg = args[counter];
		chromosomes = strsplit(arg, ",")[[1]]; 		
	}
	if (arg == "--print-non-significant"){
		printNonSig=1;
	}
	if (arg == "--width"){
		counter=counter+1;
		arg = args[counter];
		width.png = as.numeric(arg); 		
	}
	if (arg == "--height"){
		counter=counter+1;
		arg = args[counter];
		height.png = as.numeric(arg); 		
	}
	if (arg == "--cex-text"){
		counter=counter+1;
		arg = args[counter];
		cex_text = as.numeric(arg); 		
	}
	
	
	if (arg == "--xlab"){
		counter=counter+1;
		arg = args[counter];
		while ((substr(arg,1,2)!="--") & (counter <=max)){
			xlab=paste(xlab,arg, collapse=" ");
			counter=counter+1;
			arg = args[counter];
		}
		counter=counter-1;
		arg = args[counter];		
	}
	if (arg == "--ylab"){
		counter=counter+1;
		arg = args[counter];
		while ((substr(arg,1,2)!="--") & (counter <=max)){
			ylab=paste(ylab,arg, collapse=" ");
			counter=counter+1;
			arg = args[counter];
		}
		counter=counter-1;
		arg = args[counter];
	}
	if (arg == "--ylim"){
		counter=counter+1;
		arg = args[counter];
		ylim = as.numeric(strsplit(arg, ",")[[1]]);
	}
	
	if (arg == "--highlight"){
		counter=counter+1;
		arg = args[counter];
		highlight = as.numeric(arg);	
		highlightedGenes=data.frame(
									label=vector(mode="character", length = highlight), 
									start = vector(mode="numeric", length = highlight),
									end = vector(mode="numeric", length = highlight),
									chromosome = vector(mode="character", length = highlight),	
									textx = vector(mode="numeric", length=highlight),
									texty = vector(mode="numeric", length=highlight),
									min = vector(mode="numeric", length=highlight),
									max = vector(mode="numeric", length=highlight),
									stringsAsFactors=FALSE
								)
	}
	if (arg == "--region"){
		counter=counter+1;
		arg = args[counter];
		info = strsplit(arg,",")[[1]];
		highlightedGenes$label[index] = paste(highlightedGenes$label[index],as.character(info[1]));
		highlightedGenes$start[index] = as.numeric(info[2]);
		highlightedGenes$end[index] = as.numeric(info[3]);
		highlightedGenes$chromosome[index] = as.character(info[4]);
		#highlightedGenes$textx[index] = as.numeric(info[5]);
		#highlightedGenes$texty[index] = as.numeric(info[6]);
		highlightedGenes$min[index] = -1;
		index=index+1;	
	}
	counter = counter+1;
	arg=args[counter];	
}

means_by_row=function(df){
	means=c()
	
	len = length(df[,1])
	col = length(df[1,])
	for (i in 1:len){
		sum = 0
		for (j in 1:col){
			sum = sum + df[i,j]
		}
		means=c(means, sum/col)
	}
	return(means)
}



data <-read.table(inFile, header=FALSE, sep="\t")

meansTreatment1 = c()
meansTreatment2 = c()

treatment1freq <- data[,pop_freq_treatment1] 
meansTreatment1freq <- means_by_row(treatment1freq)

treatment2freq <- data[,pop_freq_treatment2] 
meansTreatment2freq <- means_by_row(treatment2freq)

referenceFreq <- data[,pop_freq_reference]

meansReferenceFreq=c()
if (is.vector(referenceFreq)){
	meansReferenceFreq=referenceFreq
}else{
	meansReferenceFreq <- means_by_row(referenceFreq)
}

data$increasedFreqInDark = abs(meansTreatment2freq - meansReferenceFreq)
direction = 2*(meansTreatment2freq > meansReferenceFreq) -1
data$increasedFreqInLight = direction * (meansTreatment1freq - meansReferenceFreq)


#sort according to chr/pos
data =  data[order( factor(data[,1], levels=chromosomes), data[,2]), ]

ticks_labels = c()
ticks_diff = c(0)
lastbase = 0

for (i in 1:length(chromosomes)){
	chr = chromosomes[i]
	length = length(subset(data, data[,1]==chr)[,2])
	lastbase = length + lastbase
	tick = lastbase - ceiling(length/2)
	ticks_labels=c(ticks_labels, tick)
	ticks_diff = c(ticks_diff, lastbase)
}

len = length(data[,1])

#default, generic, do separately for bickels fwith highlighting
outFile = paste(sep=".", inFile,"png")
oldpar<-par()
png(outFile, height=height.png, width=width.png)

par(cex.axis=2, mex = 3,cex.lab=2.3)
plot(seq(1:len), seq(from=0, to=0, length.out=len), xlim=c(0, len+1), type = "l",ylim=ylim, xaxt="n", xlab=xlab, ylab=ylab, col=col_base)
axis(1, at = ticks_labels, lab = chromosomes, tick=FALSE)

for (i in 1:length(data[,2])){
	
	nonSig1 = 0
	nonSig2 = 0
	
	if (data[i,CMH_ref_1] > cutoff){nonSig1 = 1}
	if (data[i,CMH_ref_2] > cutoff){nonSig2 = 1}
	
		
	if (!nonSig1 || !nonSig2){
		
		with(data, lines(c(i,i), c(increasedFreqInLight[i],increasedFreqInDark[i]), lty = 3 ))	
		
		if (data[i, CMH_ref_1] <= cutoff){
		#light significant
			with(data, points(i, increasedFreqInLight[i], col=col_treatment1_sig, pch=19))
		}else{
		#light non significant
			with(data, points(i, increasedFreqInLight[i], col=col_treatment1_nonSig, pch=19))	
		}
		if (data[i, CMH_ref_2] <= cutoff){
		#dark significant
			with(data, points(i, increasedFreqInDark[i], col=col_treatment2_sig, pch=19))	
		}else{
		#dark non significant	
			with(data, points(i, increasedFreqInDark[i], col=col_treatment2_nonSig, pch=19))	
		}
	}else{
		#print nonsignificant changes, if significant
		if(printNonSig == 1){
			with(data, lines(c(i,i), c(increasedFreqInLight[i],increasedFreqInDark[i]), lty = 3 ))		
			with(data, points(i, increasedFreqInLight[i], col=col_treatment1_nonSig, pch=19))	
			with(data, points(i, increasedFreqInDark[i], col=col_treatment2_nonSig, pch=19))	
		}
	}
	

	for (j in 1:highlight){
		if (as.character(data[i,1]) == highlightedGenes$chromosome[j]){
			if ((as.numeric(data[i,2]) <= highlightedGenes$end[j])&(as.numeric(data[i,2])>=highlightedGenes$start[j])){
				if (highlightedGenes$min[j] == -1){highlightedGenes$min[j]=i}
				highlightedGenes$max[j]=i
			}
		}
	}	
}

#highlight genes
m = min(abs(ylim[1]), abs(ylim[2]))
for (j in 1:highlight){
	lines(c(highlightedGenes$min[j]-0.5, highlightedGenes$max[j]+0.5, highlightedGenes$max[j]+0.5, highlightedGenes$min[j]-0.5, highlightedGenes$min[j]-0.5), c(-m,-m,m,m,-m))
#	text(highlightedGenes$textx[j], highlightedGenes$texty[j], labels=highlightedGenes$label[j], cex=cex_text)
	text(highlightedGenes$min[j]-0.5, m+0.05, labels=highlightedGenes$label[j], cex=cex_text, pos=4)

}


#
# legend
#

legend(x=40, y=0.7, 
	legend=c(paste(sep=" ", "cutoff p-value", cutoff), 
		"light flies, significant", 
		"light flies, others", 
		"dark flies, significant", 
		"dark flies, others"), 
	pch=c(19,19, 19, 19, 19), 
	col = c(hsv(1,0,1),
		   col_treatment1_sig,
		   col_treatment1_nonSig,
		   col_treatment2_sig,
		   col_treatment2_nonSig), 
	lwd=c(0,0,0,0,0), 
	lty=c(0,1,2,1,0), 
	pt.lwd=c(1,1,1,1),
	cex = cex_text
)

dev.off()	















#plot(seq(1:len), seq(from=0, to=0, length.out=len), xlim=c(0, len+1), type = "l", ylim=c(-0.35,0.6), xaxt="n", xlab="chromosome", ylab="distance of light/dark freq from reference")
##tan: X:9106807-9122238
##cip4: 3L:4322616-4364120

##ebony: 3R:17050772-17067798
##bab1+bab2+trio: 3L:991207-1191114 (all bickels positions included)

#tan_min = 0
#tan_max = 0
#cip4_min = 0
#cip4_max = 0
#ebony_min = 0
#ebony_max = 0
#bab_min = 0
#bab_max = 0


#for (i in 1:length(data$position)){
	
#	nonSigLight = 0
#	nonSigDark = 0

#	if (data$lightCMH[i] > cutoff){nonSigLight = 1}
#	if (data$darkCMH[i] > cutoff){nonSigDark = 1}



#	if (!nonSigLight || !nonSigDark){

#		with(data, lines(c(i,i), c(increasedFreqInLight[i],increasedFreqInDark[i]), lty = 3 ))	
	
#		if (data$lightCMH[i] <= cutoff_p){
#			#light significant
#			with(data, points(i, increasedFreqInLight[i], col=lightColSig, pch=19))
#		}else{
#			#light non significant
#			with(data, points(i, increasedFreqInLight[i], col=lightColNonSig, pch=19))	
#		}
#		if (data$darkCMH[i] <= cutoff_p){
#			#dark significant
#			with(data, points(i, increasedFreqInDark[i], col=darkColSig, pch=19))	
#		}else{
#			#dark non significant	
#			with(data, points(i, increasedFreqInDark[i], col=darkColNonSig, pch=19))	
#		}
#	
##		if (1-data$lightCMH[i] <= cutoff_1_min_p){
##			#light nondistinguishable from early-late
##			with(data, points(i, increasedFreqInLight[i], col=lightColNonDist, pch=19))
##		}
###			#dark nondistinguishable from early-late
##			with(data, points(i, increasedFreqInDark[i], col=darkColNonDist, pch=19))	
##						
##		}
#	
##tan
#		if (data$chromosome[i] == "X"){
#			if ((data$position[i]<=9122238)&&(data$position[i]>=9106807)){
#				if (tan_min == 0){
#					tan_min=i	
#				}
#				tan_max=i		
#			}
#	
#		}	
#
##bab
#		if (data$chromosome[i] == "3L"){
#			if ((data$position[i]<=1191114)&&(data$position[i]>=991207)){
#				if (bab_min == 0){
#					bab_min=i	
#				}
#				bab_max=i	
#			}
#		}	
#
##cip4
#		if (data$chromosome[i] == "3L"){
#			if ((data$position[i]<=4364120)&&(data$position[i]>=4322616)){
#				if (cip4_min == 0){
#					cip4_min=i	
#				}
#				cip4_max=i	
#			}
#		}	
#
##ebony
#		if (data$chromosome[i] == "3R"){
#			if ((data$position[i]<=17067798)&&(data$position[i]>=17050772)){
#				if (ebony_min == 0){
#					ebony_min=i	
#				}
#				ebony_max=i	
#			}
#		}	
#	}
#	
#}

#highlight tan
#lines(c(tan_min, tan_max, tan_max, tan_min, tan_min), c(-0.35,-0.35,0.35,0.35,-0.35))
#text(tan_min+8, 0.4, labels="tan")
##highlight ebony
#lines(c(ebony_min, ebony_max, ebony_max, ebony_min, ebony_min), c(-0.35,-0.35,0.35,0.35,-0.35))
#text(ebony_min, 0.4, labels="ebony")
##highlight bab
#lines(c(bab_min, bab_max, bab_max, bab_min, bab_min), c(-0.35,-0.35,0.35,0.35,-0.35))
#text(bab_min+6, 0.4, labels="bab")
##highlight cip4
#lines(c(cip4_min, cip4_max, cip4_max, cip4_min, cip4_min), c(-0.35,-0.35,0.35,0.35,-0.35))	
#text(cip4_min, 0.4, labels="cip4")	#
	
#axis(1, at=c(0, ticks_labels, lastbase), lab=c("", chromosomes, ""), tck = 0)
#axis(1, at=ticks_diff, labels=c("","","","","",""))
##print sig/nonsig p-values in different colors
##print legend

##highlight bab, tan, ebony, cip4


#lightColSig = hsv(1,1,0.7)
#lightColNonSig = hsv(1,0,0.9)
#lightColNonDist = hsv(1,0,0) #hsv(0.6,0.8,0.8)

#earlyLateCol = hsv(1,0,0)
#darkColSig = hsv(1,1,0.4)
#darkColNonSig = hsv(1,0,0.6)
#darkColNonDist = hsv(1,0,0) #hsv(0.6,0.8,0.4)









##### plot p-value distribution for dark vs early late, light vs early late and dark vs light

#data_dark<-read.table("/Volumes/clc/tina/CMH_selection_trends/dark_CMH_test", sep = "\t", header=FALSE)
#data_light<-read.table("/Volumes/clc/tina/CMH_selection_trends/light_CMH_test", sep = "\t", header=FALSE)

#hist(data_light$V7, breaks=1000, main="Histogram of light vs. earlyLate CMH test p-value", xlab="CMH test p-value", ylab="Number of SNPs", ylim = c(0, 50000))
#hist(data_dark$V7, breaks=1000, main="Histogram of dark vs. earlyLate CMH test p-value", xlab="CMH test p-value", ylab="Number of SNPs", ylim = c(0, 50000))


#data_lightDark<-read.table("/Volumes/clc/tina/CMH_selection_trends/from_CMH_test_mq_20_bdtest_filtered_01_CMH_8_columns", sep = "\t", header=FALSE)
#hist(data_lightDark$V8, breaks=1000, main="Histogram of light vs. dark CMH test p-value", xlab="CMH test p-value", ylab="Number of SNPs",ylim=c(0,50000))








#lightFreq<-data.frame(lightPop1freq=data$lightPop1freq, lightPop2freq=data$lightPop2freq)
#trans_lightFreq_df <- data.frame(t(lightFreq))
#meansLight <- mean(trans_lightFreq_df)

#darkFreq <- data.frame(darkPop1freq=data$darkPop1freq, darkPop2freq=data$darkPop2freq)
#trans_darkFreq_df <- data.frame(t(darkFreq))
#meansDark <- mean(trans_darkFreq_df)

#data$meansLight = meansLight
#data$meansDark = meansDark
#data$increasedFreqInDark = abs(data$meansDark - data$earlyLatefreq)
#data$direction = 2*(data$meansDark > data$earlyLatefreq) -1
#data$increasedFreqInLight = data$direction * (data$meansLight - data$earlyLatefreq)

#sort according to chr/pos

#chromosomes=c("X", "2L", "2R", "3L", "3R", "4")

#data =  data[order( factor(data$chromosome, levels=chromosomes), data$position), ]

#ticks_labels = c()
#ticks_diff = c(0)
#lastbase = 0

#for (i in 1:length(chromosomes)){
#	chr = chromosomes[i]
#	length = length(subset(data, chromosome==chr)$position)
#	lastbase = length + lastbase
#	tick = lastbase - ceiling(length/2)
#	ticks_labels=c(ticks_labels, tick)
#	ticks_diff = c(ticks_diff, lastbase)
#}

#with( data, plot( seq(1:length(position)), seq(from=0, to=0, length.out=length(position)), xlim=c(0, length(data$position)+1), type = "l", ylim=c(-0.35,0.6), xaxt="n", xlab="chromosome", ylab="distance of light/dark freq from reference"))
