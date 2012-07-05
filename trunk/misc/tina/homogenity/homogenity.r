args = commandArgs()

#usage: R --vanilla --slave --input inputName --output outputName --populations populations < homogenity.r > log 2> err

counter = 1
arg =  args[1]

m = length(args)

inFile = ""
outFile = ""
populations = list()

while (counter <=m){
	if (arg == "--input"){
		counter=counter+1;
		arg = args[counter];
		inFile = arg;
	}

	if (arg == "--output"){
		counter=counter+1;
		arg = args[counter];
		outFile = arg;		
	}
	
	if (arg == "--populations"){
		counter=counter+1;
		arg = args[counter];
		populations = strsplit(arg, ",")[[1]]
	}
	
	counter = counter+1;
	arg=args[counter];	
}

numberOfPairs = length(populations)

pairs=strsplit(populations,"-")

file.create(outFile)

woolf = function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  1 - pchisq(sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2), k - 1)
}

require(metafor)

data <- read.table(inFile, sep="\t", stringsAsFactors=FALSE)

woolfpval=0
bdpval = 0	

#write("", file = outFile, append=FALSE)

count = length(data$V1)
write(paste("need to be processed:",count,"data points.\n\n"), stderr())

write("processed until now:", stderr())
for (i in 1:count){	
		
	if ((i %% 50000 == 0)||(i == count)){
		out = i/count 
		write(paste(out*100, " %", collapse=""), stderr())
	}
	if(i==count){
		write("done.", stderr())
	}
	
	m = matrix(nrow=2*numberOfPairs+1, ncol=4)
	mi=1
	
	#split populations in right order, store in matrix, sum columns, find 2nd biggest value, select the columns, use the values in tests
	print(data[i,])
	
	for (k in 1:numberOfPairs){
		for (j in 1:2){			
			pop = strsplit(as.character(data[k, as.numeric(pairs[[k]][j])]), ":")[[1]]
			p = pop[1:4]
			m[mi,1:4]=as.numeric(p)
			mi=mi+1
		}
	}
	
	print(m)

	max=0
	second=0
		
	for (n in 1:4){ 

		s=0
		for (j in 1:(2*numberOfPairs)){		
			s = s+as.numeric(m[j,n])
		}
		m[2*numberOfPairs+1,n] = s			
	}

	print(m)
	
	for (n in 1:4){
		if (max<m[2*numberOfPairs+1,n]){
			second = max;
			max = m[2*numberOfPairs+1,n];
		}else{
			if(second< m[2*numberOfPairs+1,n]){
				second =  m[2*numberOfPairs+1,n];
			}
		}
	}
		
	selected = m[2*numberOfPairs+1,] >= second
	
	print(selected)
	
	print(length(selected[selected == TRUE]))
		  
	if (length(selected[selected == TRUE]) == 2){
	
		#check that there are two alleles with counts >= second 
		
		allele_counts = m[,selected]
		
		#print(allele_counts)
	
		for_woolf=vector(mode="numeric", length= 4*numberOfPairs)
		for_bd=matrix(nrow = 4, ncol = numberOfPairs)
	
		for (k in 1:numberOfPairs){
			index = 1+2*(k-1)
	
			for_woolf[1+4*(k-1)] = allele_counts[index,1]
			for_woolf[2+4*(k-1)] = allele_counts[index+1,1]
			for_woolf[3+4*(k-1)] = allele_counts[index,2]
			for_woolf[4+4*(k-1)] = allele_counts[index+1,2]
			for_bd[1,k] = allele_counts[index,1]
			for_bd[2,k] = allele_counts[index,2]
			for_bd[3,k] = allele_counts[index+1,1]
			for_bd[4,k] = allele_counts[index+1,2]
		}
	
	
		# run woolfs test, correction +1/2
		dim(for_woolf)<-c(2,2,numberOfPairs)
		woolfpval = woolf(for_woolf)
			
		# run BD test, correction +1/2
		res= rma.mh(ai=for_bd[1,], bi=for_bd[2,], ci=for_bd[3,], di=for_bd[4,], add=c(1/2,1/2), to=c("all", "all"))
		bdpval = res$BDp				
	
	}else{
		# for example if sums of counts are 41,20,20,3 then woolfpval, bdpval are set to NA 
		woolfpval = "NA"
		bdpval = "NA"
	}	
	write(paste(paste(paste(data[i,], collapse="\t"),woolfpval, sep="\t"), bdpval, sep="\t"), file = outFile, append=TRUE)
}
