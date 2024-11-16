#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#define RAN() ((double)rand()/(double)(RAND_MAX))
#define PI M_PI
 // declaracao da ran1
   double ran1(long int *idum);
   // routine's declarations
   #define IA 16807
   #define IM 2147483647
   #define AM (1.0/IM)
   #define IQ 127773
   #define IR 2836
   #define NTAB 32
   #define NDIV (1+(IM-1)/NTAB)
   #define EPS 1.2e-7
   #define RNMX (1.0-EPS)
 // fim da declracao da ran1   
#define pmax 200000
#define nmax 6200

double _Complex C[nmax+1],CC2[nmax+1],y[nmax+1],gy,a;
double _Complex CC1[100];
double e[nmax+1],T[nmax+1],dt,t,PAR,psisum,WX,retorno,ext;
double fat,cent1,cent2,phi[nmax],mean,mean2,alpha,v[nmax+1],gg;
double aux;
long int amostra,SEED;
int n,l,l0,N,total_time,j,SS1,k,i,nm,ym,ym1,N2;
double tm,pr1[pmax],te1[pmax],norm1[pmax],borda1[pmax];
double module2(double _Complex r) {
    return creal(r)*creal(r)+cimag(r)*cimag(r);
}


   
   int main(int argc, char**argv){
   FILE *fil0, *fil1;
   char filename0[650], filename1[650];
   tm=600.;
   //########### dados de entrada #############//  
   gy=0.0-1.0*I;

   l0=8;
   dt=0.1;
   gg=0.1;
   total_time=300;

  N=5602;
  nm=300;
  SS1=1;

              printf("alpha:\n ");
        scanf("%lf",&alpha);
        aux=(double)N;
        aux=aux/2.;
        N2=(int)aux;
 printf("%i %i\n",N,N2);
   //################# dados de entrada ##############//

   //################# arquivos de saida ##############//
    sprintf(filename0,"PRN%ialpha%1.2lfM%iSEED%i.dat",N,alpha,nm,SS1);

    fil0=fopen(filename0,"w");
     sprintf(filename0,"CRTN%ialpha%1.2lfM%iSEED%i.dat",N,alpha,nm,SS1);

    fil1=fopen(filename0,"w");

   SEED=-SS1;
   //################# arquivos de saida ##############//

   //############ zerando as quantidades ############# //

  for(n=1;n<=pmax-1;n++){
      borda1[n]=norm1[n]=pr1[n]=te1[n]=0.;

  }
 for(ym=1;ym<=nm;ym++){
 gy=0.0-1.0*I;
   //N=302;
   l0=9;
   dt=0.1;
   gg=0.1;
      a=0.0+0.0*I;
      for(n=1;n<=N;n++){

            C[n]=y[n]=0.0+0.0*I;
	    CC2[n]=0.0+0.0*I;
	    e[n]=T[n]=0.0;

      }      
       a=(gy*dt);
       CC1[1]=a;

         for(l=2;l<=l0;l++){

	   CC1[l]=CC1[l-1]*a;
	 }
	   fat=1.;
        for(l=2;l<=l0;l++){

         fat=fat*(double)l;

	   CC1[l]=CC1[l]/fat;
	 }



   //############ zerando as quantidades ############# //
   
   
   //############ condicoes iniciais #################//
      C[N2+1]=y[N2+1]=1.0+0.0*I;
   //############# condicoes iniciais ################//
  
   //################## potencial ####################//

        mean=0.0;
     for(k=1;k<=N/2;k++)  {

         phi[k]=2.*PI*ran1(&SEED);
     }
//     printf("%lf\n",phi[k]);}
     for(i=1;i<=N;i++) {
       v[i]=0.0;
       for(k=1;k<=N/2;k++) {
         aux=((double)i)*((double)k)/((double)(N));
	 v[i]+=(exp(-0.5*alpha*log((double)k)))*cos(2.0*PI*aux+phi[k]);
      }
      mean+=v[i];
   }
   mean/=(double)N;
   for(i=1;i<=N;i++) v[i]-=mean;
   aux=0.0;
   for(i=1;i<=N;i++) aux+=v[i]*v[i];
   aux=sqrt(aux/(double)N);
   for(i=1;i<=N;i++) {
      v[i]*=1./aux;
   }

   for(i=1;i<=N;i++){
    if(i<(N2+1)){
     e[i]=v[i];
    }

    if (i==N2){
     e[i]=0.;
    }
    if(i>(N2)){
     e[i]=v[i-1];
    }



   }



   //################## potencial ####################//
    ym1=0;
   //############ time evolution operator #############//
   for(t=dt;t<=tm;t+=dt){
       ym1=ym1+1;
       for(l=1;l<=l0;l++){
	

          for(n=1;n<=N;n++){
         if (n<N2){
	      CC2[n]=e[n]*C[n]+C[n+1]+C[n-1];
	     }

	     if (n==N2){
         CC2[n]=e[n]*C[n]+gg*C[n+1]+C[n-1]+C[n+2];
         }

            if (n==(N2+1)){
         CC2[n]=e[n]*C[n]+gg*C[n-1];
         }
          if (n==(N2+2)){
         CC2[n]=e[n]*C[n]+C[n+1]+C[n-2];
         }


         if (n>(N2+2)){
	      CC2[n]=e[n]*C[n]+C[n+1]+C[n-1];
	     }






	   }


	   for(n=1;n<=N;n++){
		y[n]+=CC1[l]*CC2[n];
		 C[n]=CC2[n];
	      }
	     
	     
      }
   psisum=0.0;
        for(n=1;n<=N;n++){
	 C[n]=y[n];
     psisum=psisum+module2(y[n]);
          }
       
	 
	     

   PAR=0.0;
   retorno=module2(y[N2+1]);
   pr1[ym1]=pr1[ym1]+retorno/(double)nm;
   norm1[ym1]=norm1[ym1]+psisum/(double)nm;
   borda1[ym1]=borda1[ym1]+(module2(y[1])+module2(y[N]))/(double)nm;
   te1[ym1]=t;
   cent1=0.0;
   cent2=0.0;
   ext=0.0;


  // fprintf(fil0,"%12.6lf %20.6g %20.6g \n",t,retorno,fabs(1.0-psisum));
   } //loop of time
   
 }// loop de amostras
  for(n=1;n<=ym1;n++){
  fprintf(fil0,"%15.6lf %17.7g  %17.7g %17.7g \n",te1[n],pr1[n],fabs(1.-norm1[n]),borda1[n]);
  }

      fclose(fil0);

   for(n=40;n<=ym1;n=n+40){

      aux=0.;
      for(i=1;i<=(n-1);i=i+1){
          aux=aux+dt*(pr1[i]+pr1[i+1])*0.5;




      }

      fprintf(fil1,"%12.6lf %20.6g \n",te1[n],aux/te1[n]);



   }



   } //loop of main  
   
   
     double ran1(long int *idum){
      int j;
      long k;
      static long iy=0;
      static long iv[NTAB];
      double temp;
      if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
      }
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      j=iy/NDIV;
      iy=iv[j];
      iv[j] = *idum;
      if ((temp=AM*iy) > RNMX) return RNMX;
      else return temp;
}
 
   
