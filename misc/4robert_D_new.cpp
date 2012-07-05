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
	
	if(argc != 6){
		cout << "Gimme 1-MAF 2-MAFaccepted 3-coverage 4-poolSize 5-theta." << endl;

		}


	int i = atoi(argv[1]);    //allele count
	int b = atoi(argv[2]);    //min cov allowed
	int cov = atoi(argv[3]);    // coverage 
	int n = atoi(argv[4]);	//sample size
	float theta= atof(argv[5]);  //an estimate of theta
	
	double Wat=0;
	if((i>=b)&&(i<=cov-b)){
	Wat=1.0/Watterson_corr(cov,n,b);
	}else{Wat=0;}
	double Taj=0;
	double Var=0;
	double D=0;
	if((i>=b)&&(i<=cov-b)){
	Taj=(i*(cov-i))/(Tajima_corr(cov,n,b));
	}else{Taj=0;}
	cout << Wat<< "      is Watterson estimate."<<endl ;
	cout << Taj<< "      is Tajima estimate."<<endl ;
	
	Var=Variance(cov,n,b,theta);
	cout << Var<< "      is variance estimate."<<endl ;
	Var=sqrt(Var);
	D=(Taj-Wat)/Var;
	cout << D<< "      is Tajima D estimate."<<endl ;
	
	Wat=Watterson_simp(cov,n,b);
	cout << Wat<< "      is Watterson simplified."<<endl ;
	Taj=Tajima_simp(cov,n,b)*((cov*cov)-(i*i)-((cov-i)*(cov-i)))/(1.0 * cov*(cov-1));
	cout << Taj<< "      is Tajima simplified."<<endl ;
	Var=Variance_simp(cov,n,b,theta);
	cout << Var<< "      is variance simplified."<<endl ;
	Var=sqrt(Var);
	D=(Taj-Wat)/Var;
	cout << D<< "      is Tajima D simplified."<<endl ;
	
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




double WaCo(int n1){
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
