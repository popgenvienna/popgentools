#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <getopt.h>
#include <math.h>
#include <iomanip> 


using namespace std;

int main(int argc, char *argv[]){
	
	//float tajimasPi(char *,int,int);
	//int min(int,int);
	//int factorial(int );
	//int binomialc(int,int);
	//double P(int,int,int);
	float Watterson_corr(int,int,int);
	float Tajima_corr(int, int, int);
	double WaCo(int);
	float Tajima_simp(int,int,int);
	float Watterson_simp(int,int,int);
	float Variance(int,int,int,float);
	float Variance_simp(int,int,int,float);
	float P(int k, int cov1, int n1);
	float Pi_b_split_correction(int,int,int,float,int,int,int,float,int);
	float Pi_b_split(int,int,int,float,int,int,int,float,float);
	float Pm(int,int,int,int);
	float Pi_b_split_simp(int,int,int,float,int,int,int,float,int);
	double* aj(int,int);
	double Pi_site(int,int,int);
	double Pjoint(int,int,int,int,int,int,int,int);
	double Pi_site_split(int,int,int,int,int,int,float,float);
	
	if(argc != 15){
		cout << "Gimme 1-MAF 1b-MAF_2 2-MAFaccepted 2b-MAFA_2 3-coverage 3b-coverage_2 4-poolSize 4b-poolsize_2 5-theta 6-MAFAtotal 7-a1 7b-a2 8-w1 8b-w2." << endl;

		}


	int i = atoi(argv[1]);    //allele count
	int i2 = atoi(argv[2]);    //allele count2
	int b = atoi(argv[3]);    //min MAF allowed
	int b2 = atoi(argv[4]);    //min MAF allowed2
	int cov = atoi(argv[5]);    // coverage 
	int cov2 = atoi(argv[6]);    // coverage2
	int n = atoi(argv[7]);	//sample size
	int n2 = atoi(argv[8]);	//sample size2
	float theta= atof(argv[9]);  //an estimate of theta
	int btot= atoi(argv[10]);  //MAFA of total data
	float a1= atoi(argv[11]);  //1st weight
	float a2= atoi(argv[12]);  //2nd weight
	float w1= atoi(argv[13]);  //1st frequency weight
	float w2= atoi(argv[14]);  //2nd frequency weight
	
	double Wat=0;
	double Watsimp=0;
	double* weights; 
	if((i>=b)&&(i<=cov-b)){
	Wat=1.0/Watterson_corr(cov,n,b);
	}else{Wat=0;}
	
	cout << Wat<< "      is Watterson estimate."<<endl ;

	double Taj=0;
	double Taj2=0;
	double Tajtot=0;
	double Tajsimp=0,site1=0,site2=0,sitetot=0,a1new=0,a2new=0,w1new=0,w2new=0;
	double Var=0;
	double D=0;
	double Fst=0,Tajtot_simp=0;
	if((i>=b)&&(i<=cov-b)){
	Taj=(i*(cov-i))/(Tajima_corr(cov,n,b));
	}else{Taj=0;}

	cout << Taj<< "      is Tajima estimate."<<endl ;
	
	Var=Variance(cov,n,b,theta);
	//cout << Var<< "      is variance estimate."<<endl ;
	Var=sqrt(Var);
	D=(Taj-Wat)/Var;
	cout << D<< "      is Tajima D estimate."<<endl ;
	
	Watsimp=Watterson_simp(cov,n,b);
	cout << Watsimp<< "      is Watterson simplified."<<endl ;
	Tajsimp=Tajima_simp(cov,n,b)*((cov*cov)-(i*i)-((cov-i)*(cov-i)))/(1.0 * cov*(cov-1));
	cout << Tajsimp<< "      is Tajima simplified."<<endl ;
	Var=Variance_simp(cov,n,b,theta);
	//cout << Var<< "      is variance simplified."<<endl ;
	Var=sqrt(Var);
	D=(Taj-Wat)/Var;
	cout << D<< "      is Tajima D simplified."<<endl ;
	
	if((i2>=b2)&&(i2<=cov2-b2)){
		Taj2=(i2*(cov2-i2))/(Tajima_corr(cov2,n2,b2));
	}else{Taj2=0;}
	
	cout << Taj2<< "      is the second pop. Tajima estimate."<<endl ;
	
	//cout << Pi_b_split(i,cov,n,w1,i2,cov2,n2,w2,btot)<< "      is the uncorrected total population Tajima estimate."<<endl ;
	
	Tajtot=Pi_b_split(i,cov,n,w1,i2,cov2,n2,w2,btot)/Pi_b_split_correction(i,cov,n,w1,i2,cov2,n2,w2,btot);
	
	cout << Tajtot<< "      is the total population Tajima estimate."<<endl ;
	
	Fst=1.0-(a1*Taj + a2*Taj2)/((a1+a2)*Tajtot);
	
	//cout << Fst << "      is the Fst estimate."<<endl ;
	
	Tajtot_simp=Pi_b_split(i,cov,n,w1,i2,cov2,n2,w2,btot)/Pi_b_split_simp(i,cov,n,w1,i2,cov2,n2,w2,btot);
	
	cout << Tajtot_simp << "      is the total population simplified Tajima estimate."<<endl ;
	
	cout << Fst << "      is the regional Fst estimate."<<endl ;
	
	weights=aj(cov,n);
	w1new=weights[0];
	a1new=weights[1];
	weights=aj(cov2,n2);
	w2new=weights[0];
	a2new=weights[1];
	
	cout <<"ENSI weights:  "<< w1new << " E[I_1] " <<a1new << " E[I_1*(I_1 - 1)/2] "<< w2new << " E[I_2] " <<a2new << " E[I_2*(I_2 - 1)/2] "<<endl ;
	
	Tajtot=Pi_b_split(i,cov,n,w1new,i2,cov2,n2,w2new,btot)/Pi_b_split_correction(i,cov,n,w1new,i2,cov2,n2,w2new,btot);
	Fst=1.0-(a1new*Taj + a2new*Taj2)/((a1new+a2new)*Tajtot);
	cout << Fst << "      is the Fst estimate with the ENSI weights."<<endl ;
	
	site1=Pi_site(i,cov,n);
	cout << site1 << "      is the taj estimate at a single locus."<<endl ;
	site2=Pi_site(i2,cov2,n2);
	cout << site2 << "      is the taj estimate for the second pop. at a single locus."<<endl ;
	sitetot=Pi_site_split(i, i2, cov,cov2, n,n2,w1,w2);
	cout << sitetot << "      is the taj estimate for the total pop. at a single locus."<<endl ;
	
	Fst=1.0-(a1*site1 + a2*site2)/((a1+a2)*sitetot);
	cout << Fst << "      is the Fst estimate for the total pop. at a single locus."<<endl ;
	
	
	sitetot=Pi_site_split(i, i2, cov,cov2, n,n2,w1new,w2new);
	Fst=1.0-(a1new*site1 + a2new*site2)/((a1new+a2new)*sitetot);
	cout << Fst << "      Fst estim. for the total pop. at a single locus with ENSI weights."<<endl ;
	
}









double* aj(int cov, int n){  //like in text 
	
	int l=0,s=0,k=0,k1=0,k2=0,sub1=0,top1=0,sub2=0,top2=0,top3=0,m1=0,m2=0,m=0,j=0;	
	double total=0,coefficient=0,p=0,term=0,sum=0, addendum=0,sum2=0,sum3=0,sum4=0,value=0,a=0,w=0;
	top1=min(cov,n);
	double PI[cov][top1];
	double results[2];
	PI[0][0]=1.0;
	for(m=1;m<=cov-1;m++){
		PI[m][0]=PI[m-1][0]/n;
	}
	for(s=1;s<=top1-1;s++){
		
		for(l=0;l<=cov-1;l++){
			sum=0;
			for(m=1;m<=l-s+1;m++){
				addendum=(PI[l-m][s-1]*(n-s))/(n);
				for(j=1;j<=m-1;j++){addendum=(addendum*(s+1))/n;}
				sum=sum+addendum;
			}
			PI[l][s]=sum;
		}
		
	
	}
	for(s=0;s<=top1-1;s++){
		w=w+PI[cov-1][s]*(s+1);
		a=a+PI[cov-1][s]*(s+1)*(s)/2;
		//cout << PI[cov-1][s]<<" probability "<<endl;
	}
	results[0]=w;
	results[1]=a;
	
	return results;
}




float Pi_b_split_simp(int i1, int cov1, int n1, float w1, int i2, int cov2, int n2, float w2, int b){  //like in text 
	
	float Pi_b_split(int,int,int,float,int,int,int,float,float);
	
	int l=0,s=0,k=0,k1=0,k2=0,n=0,sub1=0,top1=0,sub2=0,top2=0,top3=0,m1=0,m2=0,cov=0,m=0;	
	double total=0,coefficient=0,p=0,term=0,sum=0, addendum=0,sum2=0,sum3=0,sum4=0,value=0,pm=0,psplit=0;
	cov=cov1+cov2;
	for(m=b;m<=cov-b;m++){
		addendum=1.0/m;
		if(cov2<m){sub1=m-cov2;}else{sub1=0;}
		if(cov1<m){top1=cov1;}else{top1=m;}
		sum2=0;
		for(m1=sub1;m1<=top1;m1++){
			m2=m-m1;
			coefficient=1.0;
			for(s=1;s<=cov1;s++){
				coefficient=coefficient*s;
				coefficient=coefficient/(cov-cov1+s);
				if(s<=m1){coefficient=coefficient/s; coefficient=coefficient*(m-m1+s);}
				if(s<=(cov1-m1)){coefficient=coefficient/s; coefficient=coefficient*(cov-m-cov1+m1+s);}
			}
			coefficient=coefficient*Pi_b_split(m1,cov1,n1,w1,m2,cov2,n2,w2,b);
			sum2=sum2+coefficient;
		}
		sum=sum+addendum*sum2;
		//cout << sum<<endl;//
	}
	
	return sum;
}




float Pi_b_split_correction(int i1, int cov1, int n1, float w1, int i2, int cov2, int n2, float w2, int b){  //like in text 
	
	float Pi_b_split(int,int,int,float,int,int,int,float,float);
	float Pm(int,int,int,int);
	
	int l=0,s=0,k=0,k1=0,k2=0,n=0,sub1=0,top1=0,sub2=0,top2=0,top3=0,m1=0,m2=0;	
	double total=0,coefficient=0,p=0,term=0,sum=0, addendum=0,sum2=0,sum3=0,sum4=0,value=0,pm=0,psplit=0;
	n=n1+n2;
	for(k=1;k<=n-1;k++){
		addendum=1.0/k;
		if(n2<k){sub1=k-n2;}else{sub1=0;}
		if(n1<k){top1=n1;}else{top1=k;}
		sum2=0;
		for(k1=sub1;k1<=top1;k1++){
			k2=k-k1;
			coefficient=1.0;
			for(s=1;s<=n1;s++){
				coefficient=coefficient*s;
				coefficient=coefficient/(n-n1+s);
				if(s<=k1){coefficient=coefficient/s; coefficient=coefficient*(k-k1+s);}
				if(s<=(n1-k1)){coefficient=coefficient/s; coefficient=coefficient*(n-k-n1+k1+s);}
			}
			//cout << coefficient<< "coefficient"<<endl;
			sum3=0;
			for(m1=0;m1<=cov1;m1++){
				term=Pm(m1,cov1,n1,k1);
				//cout << term<< " Pm1 "<<endl;
				if(b-m1>0){sub2=b-m1;}else{sub2=0;}
				if(0<b-cov1+m1){top2=b-cov1+m1;}else{top2=0;}
				top2=cov2-top2;
				//cout << sub2<<" sub2 "<<top2<<" top2 "<<endl;
				//cout << k1<<" k1 "<<k2<<" k2 "<<endl;
				//cout << n1<<" n1 "<<n2<<" n2 "<<endl;
				sum4=0;
				for(m2=sub2;m2<=top2;m2++){
					pm=Pm(m2,cov2,n2,k2);
					psplit=Pi_b_split(m1,cov1,n1,w1,m2,cov2,n2,w2,b);
					//cout << psplit<<" psplit "<<pm<<" pm "<<m1<<" m1 "<<m2<<" m2 "<<endl;
					sum4=sum4+pm*psplit;
				}
				//cout << sum4<<" sum4 "<<m1<<" m1 "<<endl;
				term=term*sum4;
				
				sum3=sum3+term;
				//cout << sum3<<" sum3 "<<endl;
			//term here
			}
			sum2=sum2+coefficient*sum3;
			//cout << sum2<<"sum2"<<endl;
			
		}
		sum=sum+addendum*sum2;
		//cout << sum<<endl;//
	}

	return sum;
}




float Pi_b_split(int i1, int cov1, int n1, float w1, int i2, int cov2, int n2, float w2, float b){  //like in text 
	
	float Pm(int,int,int,int);
	
	int l=0,s=0;	
	double coefficient=0,p=0,term=0;
	if(i1+i2>=b){
		p=( (w1*i1)/cov1 + (w2*i2)/cov2 )/(w1+w2);
		term=2*p*(1 - p);
	}
	return term;
}




double Pjoint(int m1, int m2, int cov1, int cov2, int n1, int n2, int k1,int k2){
	float Pm(int,int,int,int);
	
	int j=0,s=0,n=0,top1=0,sub1=0,k=0;
	double sum=0,coefficient=0;
	
	if(n2<k){sub1=k-n2;}else{sub1=0;}
	if(n1<k){top1=n1;}else{top1=k;}
	sum=0;
	n=n1+n2;
	//for(k1=sub1;k1<=top1;k1++){
	k=k2+k1;
	coefficient=1.0;
	for(s=1;s<=n1;s++){
		coefficient=coefficient*s;
		coefficient=coefficient/(n-n1+s);
		if(s<=k1){coefficient=coefficient/s; coefficient=coefficient*(k-k1+s);}
		if(s<=(n1-k1)){coefficient=coefficient/s; coefficient=coefficient*(n-k-n1+k1+s);}
	}
	coefficient=coefficient*Pm(m1,cov1,n1,k1)*Pm(m2,cov2,n2,k2)/k;
	
	return coefficient;
}





double Pi_site_split(int m1, int m2, int cov1, int cov2, int n1, int n2, float w1, float w2){
	double Pjoint(int,int,int,int,int,int,int,int);
	
	int j=0, k=0,n=0,sub1=0,top1=0,k1=0,k2=0;
	double norm=0,result=0,p=0;
	n=n1+n2;
	for(k=1;k<=n-1;k++){
		if(n2<k){sub1=k-n2;}else{sub1=0;}
		if(n1<k){top1=n1;}else{top1=k;}
		for(k1=sub1;k1<=top1;k1++){
			k2=k-k1;
			p=Pjoint(m1,m2,cov1,cov2,n1,n2,k1,k2);
			result=result+(p*Pi_b_split(m1,cov1,n1,w1,m2,cov2,n2,w2,1));
			norm=norm+p;
		}
	}
	result=result/norm;
	return result;
	
}






double Pi_site(int m, int cov, int n){
	float Pm(int,int,int,int);
	
	int j=0, k=0;
	double norm=0,result=0,p=0;
	for(k=1;k<=n-1;k++){
		p=Pm(m,cov,n,k);
		result=result+(p*2*(n-k) )/(n*(n-1));
		norm=norm+p/k;
	
	}
	result=result/norm;
	return result;
	
}




float Pm(int m1, int cov1, int n1, int k1){
	int s=0;	
	double p=0,term=1.0;
	if(m1<cov1-m1){m1=cov1-m1;k1=n1-k1;}
	p=k1;
	p=p/n1;
	for(s=1;s<=cov1-m1;s++){
		term=term*(cov1-s+1);
		term=term/s;
		term=term*(1.0-p);
		term=term*p;
	}
	for(s=cov1-m1+1;s<=m1;s++){term=term*p;}
	return term;
}



float P(int k, int cov1, int n1){  //like in text but multiplied by WaCo(n1)

	int l=0,s=0;	
	double coefficient=0,p=0,term=0;
	for(l=1;l<=n1-1;l++){
		term=1.0/l;
		p=l;
		p=p/n1;
		for(s=1;s<=cov1-k;s++){
			term=term*(cov1-s+1);
			term=term/s;
			term=term*(1.0-p);
			if(s<=k){ term=term*p;}
		}
		for(s=cov1-k+1;s<=k;s++){term=term*p;}
		coefficient=coefficient+term;
	}
	return coefficient;
}





float Variance_simp(int cov1, int n1,int b1,float theta1){
	double WaCo(int);
	float Watterson_simp(int,int,int);
	float Tajima_simp(int, int, int);
	
	int k,l,s;
	double sum=0, T,W,diff,p,term,coefficient,Tc,Wc;
	
	Tc=Tajima_simp(cov1,n1,b1);
	W=Watterson_simp(cov1,n1,b1);
	for(k=b1;k<=cov1-b1;k++){	
		T=Tc*((cov1*cov1)-(k*k)-((cov1-k)*(cov1-k)))/(1.0 * cov1*(cov1-1));
		diff=T-W;
		diff=diff*diff;
		sum=sum+(diff/k);
	}
	sum=sum*theta1;
	return sum;
}





float Variance(int cov1, int n1,int b1,float theta1){
	double WaCo(int);
	float Watterson_corr(int,int,int);
	float Tajima_corr(int, int, int);

	int k,l,s;
	double sum=0, T,W,diff,p,term,coefficient,Tc,Wc;

	Tc=Tajima_corr(cov1,n1,b1);
	Wc=Watterson_corr(cov1,n1,b1);
	for(k=b1;k<=cov1-b1;k++){	
		T=(k*(cov1-k))/Tc;
		W=1.0/Wc;
		diff=T-W;
		diff=diff*diff;
		sum=sum+(diff*P(k,cov1,n1));
	}
	sum=sum*theta1;
	return sum;
}




double WaCo(int n1){  //sum(1/i)
	double coefficient,result=0;
	int i;
	for(i=1;i<=n1-1;i++){coefficient=1;coefficient=coefficient/i;result=result+coefficient;}
	return result;
	}
	




float Tajima_simp(int cov, int n, int b){
double result=1.0;
	result=(result*(cov-1))/(cov- (2*b) +1);

	return result;
	
}





float Watterson_simp(int cov,int n,int b){
double result=1.0;
	if(cov<2*b){return 0;}
	
result=result/(WaCo(cov-b+1)- WaCo(b));

return result ;
}


	
	
	
float Watterson_corr(int cov,int n,int b){  
	double coefficient,sum,addendum,factor;
	int i,j,m,h;
	
	if(cov<2*b){return 0;}
	
	sum=0;
	for(m=b;m<=cov-b;m++){
		sum=sum+P(m,cov,n);
	}
	return sum;
}





float Tajima_corr(int cov, int n, int b){
	double coefficient,sum,addendum,factor;
	int i,j,m,h;
	
	if(cov<2*b){return 0;}
	
	sum=0;
	for(m=b;m<=cov-b;m++){
		addendum=m*(cov-m);
		addendum=addendum*P(m,cov,n);
		sum=sum+addendum;
	}
	return sum;
}
