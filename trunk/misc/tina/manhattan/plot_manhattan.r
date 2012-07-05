args = commandArgs()

#usage: 
#
#gwas
#R --vanilla --slave --gwas --cutoff-positions ... --cutoff-labels ... --pval-space ... --input inFile < manhattan_plot.r  
#
#small region on one chromosome
#R --vanilla --slave --local --cutoff-positions ... --cutoff-labels ... --pval-space ... --input inFile < manhattan_plot.r


#
#   parse input parameters
#

options("show.error.messages" = TRUE)

counter = 1
arg =  args[1]

max = length(args)

gwas = 0
local_plot = 0
#chromosome=""
main=""
xlab=""
inFile = ""
cutoff = 0
positions = c()
pos_count=0
labels = c()
l_count=0
pvalSpace=150
inFile1=""
inFile2=""
cmp=0
at=c()
lab=c()

while (counter <= max){
	if (arg == "--global"){gwas=1;}
	if (arg == "--local"){local_plot=1;}
	if (arg == "--input"){
		counter=counter+1;
		arg = args[counter];
		inFile = arg; 
	}	
	
#	if (arg == "--chromosome"){
#		counter=counter+1;
#		arg = args[counter];
#		chromosome = arg; 
#	}
	if (arg == "--cutoff-positions"){
		counter=counter+1;
		arg = args[counter];
		positions = as.numeric(strsplit(arg,',')[[1]])
	}
	if (arg == "--cutoff-labels"){
		counter=counter+1;
		arg = args[counter];
		labels = as.character(strsplit(arg,',')[[1]])
	}
	if (arg == "--pval-space"){
		counter=counter+1;
		arg = args[counter];
		pvalSpace = as.numeric(arg); 	
	}

	if (arg == "--main"){
		counter=counter+1;
		arg = args[counter];
		main = arg;
	}
	if (arg == "--xlab"){
		counter=counter+1;
		arg = args[counter];
		xlab = arg;
	}
	
	
	if (arg == "--cmp"){cmp=1;}
	
	if (arg == "--input1"){
		counter=counter+1;
		arg = args[counter];
		inFile1 = arg; 
	}
	if (arg == "--input2"){
		counter=counter+1;
		arg = args[counter];
		inFile2 = arg; 
	}
	if (arg == "--at"){
		counter=counter+1;
		arg = args[counter];
		at = as.numeric(strsplit(arg,',')[[1]])
	}
	if (arg == "--lab"){
		counter=counter+1;
		arg = args[counter];
		lab = strsplit(arg,',')[[1]]
	}
	
	
	counter = counter+1;
	arg=args[counter];	
}

#
# some additional settings (colors)
#

print(labels)
print(positions)

l_count=length(labels)
pos_count=length(positions)
if (pos_count>0){cutoff=1}

for (i in 1:length(labels)){
	labels[i] = paste(strsplit(as.character(labels),'_')[[i]], collapse=" ")
}

print(labels)

col.names = c("CHR", "BP", "P")
lightCol = hsv(0,0,0.8)
darkCol = hsv(0,0,0.6)
cut1Col = hsv(0,0,0.2)
cut2Col = hsv(0,0,0.4)


#for highlighting genes 
#col.names = c("CHR", "BP", "P", "GENE_ID")
#hiCol = hsv(0.77, .86, 1)



#
# functions
#

plot.mht.png=function(
fileName, 
sep="\t",
col.names=col.names,	
cutoff_positions=NULL,
cutoff_labels=NULL,
cutoffsCol = c("gray70", "gray20"),
height.png = 4.5*72, 
width.png = 18*72,
mar = par()$mar,
mgp = par()$mgp,
ymax = "max",
xlab= "",
main="",
...){
	out <- preprocess.data(fileName, sep, col.names, cutoff_positions, cutoff_labels)
	if(is.null(out$GENE_ID)){
		data = data.frame(CHR=out$CHR, BP=out$BP, P=out$P)		
	}else{	
		data = data.frame(CHR=out$CHR, BP=out$BP, P=out$P, GENE_ID=out$GENE_ID)
	}
	
	cutoff_pos = out$cutoffs
			
	rm(out)
	
	outFile = paste(sep=".", fileName, "png")	
	png(outFile, height=height.png, width=width.png)
	print("plotting manhattan plot...")
	
	manhattan(
			  data=data, 
			  cutoff_positions = cutoff_pos,
			  cutoff_labels = cutoff_labels,
			  cutoffsCol = cutoffsCol,
			  mar=mar,
			  mgp=mgp,
			  ymax=ymax, 
			  xlab=xlab,
			  main=main,
			  ...)
	dev.off()	
}




plot.mht.small.png=function(
fileName, 
sep="\t",
col.names=col.names,	
cutoffsCol = c("gray70","gray20"),
cutoff_positions=NULL,
cutoff_labels=NULL,
height.png = 5.5*72, 
width.png = 5.5*72,
mar = par()$mar,
mgp = par()$mgp,
ymax = "max",
xlab="",
main="",
...){
	#print("plot.mht.small.png")
	out <- preprocess.data(fileName, sep, col.names, cutoff_positions, cutoff_labels)
	
	if(is.null(out$GENE_ID)){
		data = data.frame(CHR=out$CHR, BP=out$BP, P=out$P)		
	}else{	
		data = data.frame(CHR=out$CHR, BP=out$BP, P=out$P, GENE_ID=out$GENE_ID)
	}
	
	cutoff_pos=out$cutoffs 
	
#	print(cutoff_pos)
#	print(cutoff_lab)
#	print(cutoff_labels)	
	rm(out)
	
	outFile = paste(sep=".", fileName, "png")	
	

	
	png(outFile, height=height.png, width=width.png)
	print("plotting manhattan plot...")
	manhattan.small(
			  data=data, 
			  cutoff_positions = cutoff_pos,
			  cutoff_labels = cutoff_labels,
  			  cutoffsCol = cutoffsCol,
			  mar=mar,
			  mgp=mgp,
			  ymax=ymax, 
			  xlab=xlab,
			  main=main,	
			  ...)
	dev.off()	
}




preprocess.data = function(
fileName,
sep="\t",
col.names=c("CHR", "BP", "P"),
cutoff_positions=NULL,
cutoff_labels = NULL
){
	
	print("loading data...")
	data<-read.table(fileName,sep=sep, header=FALSE, col.names=col.names) 
	
	cutoffs=c()
	labels=c()
	
	if ((is.null(cutoff_labels)) || (length(cutoff_labels) == 0)){
		if (!is.null(cutoff_positions) && length(cutoff_positions)>0){
##cutoff lines defined, no labels 	

			print("calculating cutoff p-values...")
			data_sorted_pval<-data[order(data$P),]		

			for (i in 1:(length(cutoff_positions))){
				pos= positions[i]
				cutoffs = c(cutoffs, data_sorted_pval$P[pos]) 
			}
			cut= data.frame(cutoffs)
			rm(data_sorted_pval)					

			out=c(data, cut)
			
		}else{
#no cutoffs at all	
			cutoffs=c()
			cut = data.frame(cutoffs)		
			out=c(data)	
		}
	}else{
		if (length(cutoff_labels) == length(cutoff_positions)){
##cutoff lines and labels are defined			
			
			print("calculating cutoff p-values...")	
			data_sorted_pval<-data[order(data$P),]		

			for (i in 1:(length(cutoff_positions))){
				pos= positions[i]
				cutoffs = c(cutoffs, data_sorted_pval$P[pos]) 
			}

			cut= data.frame(cutoffs,cutoff_labels)
			rm(data_sorted_pval)			
			
			out=c(data, cut)
		}else{
		#write an error msg to stderr
			print(cutoff_labels)
			print(cutoff_positions)
			stop("lengths of cutoff_labels and cutoff_positions differs.") 	
		}	
	}	
	return(out)
}




# plot chromosomes without positions, only chromosome name
manhattan = function(
data = dataframe, 
colors=c("gray10", "gray50"), 
ymax="max", 
cex.x.axis=1, 
chromosomes=c("X", "2L", "2R", "3L", "3R", "4"), 
gap = 0, 
offsetForPvals=0,
pch = 1,
fdrPch =1,
hiPch = 19, 
cutoff_positions=NULL, 
cutoff_labels = NULL,
cutoffsXoffset = NULL,
cutoffsYoffset = NULL,
cutoffsCol=c("gray70", "gray20"),
hiGenes = NULL,
hiGenesLabels = NULL,
hiGenesXoffset = NULL,
hiGenesYoffset = NULL,
hiGenesCol= NULL,
mar = par()$mar, 
mgp = par()$mgp,
xlab="",
main="",
...){
	
    d=data
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your input data contains columns CHR, BP, and P")
	
    if (length(chromosomes)>=1) d=d[d$CHR %in% chromosomes, ]
	
    d=subset(na.omit(d[order( factor(d$CHR, levels=chromosomes), d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
	
    ticks=NULL
	
    colors <- rep(colors,length(chromosomes))[1:length(chromosomes)]
    cutoffsCol <- rep(cutoffsCol, length(cutoff_positions))[1:length(cutoff_positions)]  
    hiGenesCol <- rep(hiGenesCol, length(chromosomes))[1:length(chromosomes)]
    
    if (ymax=="max"){ylim<-ceiling(max(d$logp)+2.5)}else{ylim<-ymax}
	
    lastbase=0
	
	print ("calculating positions...")
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
	d$pos=d$BP
	ticks=floor(length(d$pos))/2+1
	lastbase = tail(subset(d,CHR==chromosomes[1])$BP, 1) + gap
    } else {
	for (i in 1:(length(chromosomes))){
	
		chromosome = chromosomes[i]
	
		if (length(row.names(d[d$CHR==chromosome, ]))!=0){
			if (i == 1){
				d[d$CHR==chromosome, ]$pos=d[d$CHR==chromosome, ]$BP       	     			
				lastbase = tail(subset(d,CHR==chromosome)$BP, 1) + gap
			}else{
				d[d$CHR==chromosome, ]$pos=d[d$CHR==chromosome, ]$BP + lastbase + gap	
				lastbase = lastbase + tail(subset(d, d$CHR==chromosomes[i] )$BP, 1) + gap		
			}
	
			chr_start = head(subset(d, CHR==chromosome)$BP ,1)			
			chr_end = tail(subset(d, CHR==chromosome)$BP ,1)			
			tick = lastbase - ceiling((chr_end-chr_start)/2)
			ticks=c(ticks, tick)
		}else{
			print(paste("your dataset does not contain any SNP from chromosome", chromosome))	
			chromosomes <- chromosomes[-i]
		}
	}
    }
	
	oldpar = par()
	par(mar=mar, mgp=mgp)
	
	print(numchroms)
	print(lastbase)
	print(offsetForPvals)
	
	if (xlab == ""){
	xlab="Chromosome";
	}
	
	print("plotting SNPs...")

	plot(c(), c(), ylim=c(0,ylim), xlim=c(0,lastbase+offsetForPvals), ylab=expression(-log[10](italic(p)-value)), xlab=xlab, main=main, xaxt="n", type="n", bty="n",...)
	axis(1, at=c(0, ticks, lastbase), lab=c("", chromosomes, ""), tck = 0,...)
	for (i in 1:length(chromosomes)) {
		with( d[(d$CHR==chromosomes[i]), ],points(pos, logp, col=colors[i], pch=pch,...))
	}
        
	
  	#highlight SNPs above cutoffs
	if (length(cutoff_positions)>0){
	    print("    highlighting cutoff levels...")    
	    l = length(d[1,])
	    
	    maxCut = cutoff_positions[1];
	    for (i in 1:length(cutoff_positions)){
	    	if (cutoff_positions[i]>maxCut){
	    		maxCut = cutoff_positions[i]
	    	}
	    }
	    
	    allHiLevel <-subset(d, logp >= -log10(maxCut), select=c(pos, logp))
	    
	    print(nrow(allHiLevel))
	    print(cutoff_positions)
	    
	   	for (i in 1:length(cutoff_positions)){
			inv = length(cutoff_positions) - i + 1
			
			hiLevelCol = cutoffsCol[i]
			
			print(i)
			print(cutoff_positions[i])
			print(-log10(cutoff_positions[i])) 
			
			if(inv==1){
				hiLevelData <- subset(allHiLevel, logp >= -log10(cutoff_positions[inv]), select=c(pos, logp))
				with(hiLevelData, points(pos,logp, col=hiLevelCol,pch=fdrPch,...))
				print(nrow(hiLevelData))
				rm(hiLevelData)
			}else{
				print(cutoff_positions[i-1])
								
				hiLevelData <- subset(allHiLevel, ((logp >= -log10(cutoff_positions[inv]))&(logp < -log10(cutoff_positions[inv-1]))) , select=c(pos, logp))
				
				with(hiLevelData, points(pos,logp, col=hiLevelCol,pch=fdrPch,...))
				print(nrow(hiLevelData))
				rm(hiLevelData)	
			}
		}	
	}
	
#	</->highlight genes    	
#	if (length(hiGenes)>0){
#    print("    highlighting specified genes...")    		
#    for (i in 1:length(hiGenes)){
	
#		hiGene = hiGenes[i]
	
#		hiGeneLabel = hiGenesLabels[i]
#		hiGeneXoffset = hiGenesXoffset[i]
#		hiGeneYoffset = hiGenesYoffset[i]
#		hiGeneCol = hiGenesCol[i]
	
#		hiGeneData <- subset(d, GENE_ID == hiGene, select=c(pos, logp))	
		
#		with(hiGeneData, points(pos, logp, col=hiGeneCol, pch=hiPch,...))	
#		max_value = max(hiGeneData$logp)
#		first_pos = head(hiGeneData$pos, 1)
	
#		if (length(hiGenesXoffset) == 0){
#			text(first_pos, max_value+hiGeneYoffset, pos=4, labels=hiGeneLabel, cex=2, font=2, offset=0,...)
#		}else{
#			text(hiGeneXoffset+first_pos, max_value+hiGeneYoffset, pos=4, labels=hiGeneLabel, font=2, offset=0,...)
#		}
#    }
#	}
	
    # print cutoff lines    
    if (!is.null(cutoff_positions) && (length(cutoff_positions)>0)){
	    print("    plotting cutoff lines...")  	
	    for (i in 1:(length(cutoff_positions))){
	
			cutoff = cutoff_positions[i]
			cutoffXoffset = cutoffsXoffset[i]
			cutoffYoffset = cutoffsYoffset[i]
			label=cutoff_labels[i]
	
			lines( x=c(0,lastbase+2.5*gap), y=c(-log10(cutoff), -log10(cutoff) ) , lty = 2, ...)
			text(x=lastbase+offsetForPvals/5+cutoffXoffset, y=-log10(cutoff)+cutoffYoffset, labels=label, pos=4, font=3, offset=0,...)
	
    	}
    }
	
	par(oldpar)
}
	


manhattan.small = function(
data = dataframe, 
color="gray10", 
ymax="max", 
cex.x.axis=1, 
chromosome="", 
gap = 0, 
offsetForPvals=0,
pch = 19,
fdrPch =19,
hiPch = 19, 
cutoff_positions=NULL, 
cutoff_labels = NULL,
cutoffsXoffset = NULL,
cutoffsYoffset = NULL,
cutoffsCol=c("gray70", "gray20"),
hiGenes = NULL,
hiGenesLabels = NULL,
hiGenesXoffset = NULL,
hiGenesYoffset = NULL,
hiGenesCol= NULL,
mar = par()$mar, 
mgp = par()$mgp,
main="",
xlab="",
at = NULL,	
lab = NULL,	
...){
	
	print(cutoff_positions)
	
    d=data
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your input data contains columns CHR, BP, and P")
	
    d=subset(na.omit(d[order(d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    
    if (ymax=="max"){ylim<-ceiling(max(d$logp)+2.5)}else{ylim<-ymax}

    cutoffsCol <- rep(cutoffsCol, length(cutoff_positions))[1:length(cutoff_positions)]  
    
	#d$pos=d$BP
	firstbase = head(d$BP, 1)
	lastbase = tail(d$BP, 1)
        
	if (is.null(at)){
		at = c(firstbase, lastbase)
	}
	if (is.null(lab)){
		lab = c(firstbase, lastbase)
	}
	
	oldpar = par()
	par(mar=mar, mgp=mgp)
	
	with(d, plot(BP, logp, ylim=c(0,ylim), xlim=c(firstbase,lastbase+offsetForPvals), ylab=expression(-log[10](italic(p)-value)), xlab="", xaxt="n", bty="n", col=color, pch=pch, ...))
	
	axis(1, at=at, lab=lab,...)
    	
  	#highlight SNPs above cutoffs
	if (length(cutoff_positions)>0){
	    print("    highlighting SNPs above cutoff levels...")    
    	for (i in 1:length(cutoff_positions)){
			hiLevelCol = cutoffsCol[i]
			hiLevelData <- subset(d, logp >= -log10(cutoff_positions[i]), select=c(BP, logp))
			if (i > 1){
				hiLevelData <- subset(hiLevelData, logp < -log10(cutoff_positions[i-1]), select=c(BP, logp))
			}
			with(hiLevelData, points(BP,logp, col=hiLevelCol,pch=fdrPch,...))
		}	
	}
		
    # print cutoff lines    
    if (!is.null(cutoff_positions) && (length(cutoff_positions)>0)){
	    for (i in 1:(length(cutoff_positions))){
	
			cutoff = cutoff_positions[i]
			cutoffXoffset = cutoffsXoffset[i]
			cutoffYoffset = cutoffsYoffset[i]
			label=cutoff_labels[i]
	
			lines( x=c(0,lastbase+2.5*gap), y=c(-log10(cutoff), -log10(cutoff) ) , lty = 2, ...)
			text(x=lastbase+cutoffXoffset+offsetForPvals/5, y=-log10(cutoff)+cutoffYoffset, labels=label, pos=4, font=3, offset=0,...)
    	}
    }
	par(oldpar)
}
	
manhattan.cmp.two.analyses = function(
data = dataframe, 
data_group=dataframe,
col="gray10", 
col_group="gray50",
ymax="max", 
#	cex.x.axis=1, 
pch = 1,
mar = par()$mar, 
mgp = par()$mgp,
at=NULL,
lab=NULL,
main="",
xlab = "",
...) {
	
	
    d=data
    
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
	
    d=subset(na.omit(d[order(d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)

    d_group = subset(na.omit(d_group[order(d_group$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
	d_group$logp = -log10(d_group$P)
	
	
    if (ymax=="max"){ylim<-ceiling(max(d$logp)+2.5)}else{ylim<-ymax}

	firstbase = min(head(d$BP,1), head(d_group$BP,1))
	lastbase = max(tail(d$BP,1), tail(d$BP,1))

	if (is.null(at)){
		at = c(firstbase, lastbase)
	}
	if (is.null(lab)){
		lab = c(firstbase, lastbase)
	}
	
	oldpar = par()
	par(mar=mar, mgp=mgp)
	
    with(d, plot(BP, logp, ylim=c(0,ylim), xlim=c(firstbase,lastbase), ylab=expression(-log[10](italic(p)-value)), xaxt="n",main = main, xlab= xlab, bty="n", pch=pch, col=col,...))		     
    axis(1, at=at, lab=lab,...)
	#bab
	#lastbase=1191093
	#firstbase=1029080
	#at =  c(firstbase, 1030000, 1040000, 1050000, 1060000,1070000, 1080000,1090000,1100000, 1110000, 1120000, 1130000, 1140000, 1150000, 1160000, 1170000, 1180000, 1190000, lastbase)
	#lab = c(firstbase, "1030k", "1040k", "1050k", "1060k", "1070k", "1080k","1090k","1100k", "1110k", "1120k", "1130k", "1140k", "1150k", "1160k", "1170k", "1180k", "1190k", lastbase)
	#main =  "Genes bab1 and bab2"
	#xlab = paste(sep="", "chromosome ",unique(d$CHR))
	d_group = data_group
    if (!("CHR" %in% names(d_group) & "BP" %in% names(d_group) & "P" %in% names(d_group))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    with(d_group, points(BP, logp, col=col_group, pch=pch,...))	
	
	par(oldpar)
}


















































########################
#
#	main
#
##################

if (gwas==1){
	write("gwas", stdout())
#	print(positions)
	
	if ( cutoff == 1 ){
		
		if ((l_count==0)&&(pos_count>0) || (l_count==pos_count)){
# gwas, print cutoff lines, print (do not print) labels of cutoff lines
			if (l_count == pos_count){
# print labels 				
				plot.mht.png(
						 inFile,
						 sep ="\t", 
						 col.names=col.names,
						 pch = 19,
						 fdrPch =19,
						 colors=c(darkCol, lightCol),  
						 chromosomes=c("X","2L", "2R", "3L", "3R", "4"), 
						 gap=1000000, 
						 offsetForPvals=pvalSpace*100000,
						 cutoffsXoffset=c(0,0), 
						 cutoffsYoffset=c(0,0),	 
						 cutoffsCol = c(cut1Col, cut2Col),
						 cutoff_positions = positions,
						 cutoff_labels = labels,
						 cex=2,
						 cex.axis=1.3,
						 cex.lab=1.4,
						 mar = c(5,5,1,1)+0.1,
						 mgp = c(3,1,0)
					)	
			}else{
#labels - no, cutoff colors - yes
				print(positions)
				plot.mht.png(
					 inFile,
					 sep ="\t", 
					 col.names=col.names,
					 pch = 19,
					 fdrPch =19,
					 colors=c(darkCol, lightCol),  
					 chromosomes=c("X","2L", "2R", "3L", "3R", "4"), 
					 gap=1000000, 
					 cutoff_positions = positions,
					 cutoffsXoffset=c(0,0), 
					 cutoffsYoffset=c(0,0),	 
					 cutoffsCol = c(cut1Col, cut2Col),
					 cex=2,
					 cex.axis=1.3,
					 cex.lab=1.4,
					 mar = c(5,5,1,1)+0.1,
					 mgp = c(3,1,0)
				)	
			}
		}
	}else{
#gwas no cutoff colors, no labels
		plot.mht.png(
			inFile,
			sep ="\t", 
			col.names=col.names,
			pch = 19,
			fdrPch =19,
			colors=c(darkCol, lightCol),  
			chromosomes=c("X","2L", "2R", "3L", "3R", "4"), 
			gap=1000000, 
			cutoffsXoffset=c(0,0), 
			cutoffsYoffset=c(0,0),	 
			cex=2,
			cex.axis=1.3,
			cex.lab=1.4,
			mar = c(5,5,1,1)+0.1,
			mgp = c(3,1,0)
		)		
	}
}

#warnings()


if (local_plot==1){
	if (cmp == 1){
#compare two different plots

		data = read.table(inFile1, header=FALSE, sep="\t");
		data_group = read.table(inFile2, header=FALSE, sep="\t");

		outFile = paste(sep=".", inFile1, inFile2, "png")	
		png(outFile, height=height.png, width=width.png)
		print("plotting manhattan plot...")
		manhattan.cmp.two.analyses(
								  data = data, 
								  data_group=data_group,
								  col="gray10", 
								  col_group="gray50",
								  ymax="max", 
								  pch = 1,
								  at=at,
								  lab=lab,
								  main=main,
								  xlab = xlab)
		
		dev.off()
			
	}else{	
	
	
	if ( cutoff == 1){
		if ((l_count==0)&&(pos_count>0) || (l_count==pos_count)){
			if (l_count == pos_count){
##cutoff colors - yes, labels - yes				
				plot.mht.small.png(
					inFile, 
					sep="\t",
					col.names=col.names,
					pch=19,
					fdrPch=19,
					color= darkCol,	
					chromosome=chromosome,  
					offsetForPvals=pvalSpace*1000,
					height.png = 5.5*72, 
					width.png = 5.5*72,
					gap = 10000,
					cutoffsCol = c(cut1Col, cut2Col),
					cutoff_positions = positions,
					cutoff_labels = labels,
					cutoffsXoffset=c(0,0), 
					cutoffsYoffset=c(0,0),	 
					mar = par()$mar,
					mgp = par()$mgp,
					ymax = "max"
				)	
			}else{
#cutoff colors - yes,labels - no 				
			plot.mht.small.png(
					inFile, 
					sep="\t",
					col.names=col.names,
					pch=19,
					fdrPch=19,
					color= darkCol,	
					chromosome=chromosome,  
   					offsetForPvals=pvalSpace*1000,
					height.png = 5.5*72, 
					width.png = 5.5*72,
					cutoffsCol = c(cut1Col, cut2Col),
					cutoff_positions = positions,
					cutoffsXoffset=c(0,0), 
					cutoffsYoffset=c(0,0),	 
					mar = par()$mar,
					mgp = par()$mgp,
					ymax = "max"
				)	
			}		
		}
	}else{
#	#no cutoffs
		plot.mht.small.png(
			inFile, 
			sep="\t",
			col.names=col.names,
			pch=19,
			fdrPch=19,
			color= darkCol,	
#			chromosome=chromosome,  
			height.png = 5.5*72, 
			width.png = 5.5*72,
			cutoffsXoffset=c(0,0), 
			cutoffsYoffset=c(0,0),	 
			mar = par()$mar,
			mgp = par()$mgp,
			ymax = "max"
		)	
	}
	}	
}




warnings()