// Here our data structure in 'resid2.tmp' is made by 80 bytes:
  // 1 unsigned long (4 bytes)
  // 9 double        (8 bytes)
  // 1 unsigned long (4 bytes)
  // total bytes = 4 + 8*9 + 4 = 80 bytes
  // This is true in a 32bit machine
  // (On a 64 bit it's 88bytes 'cause ul becomes ull)
  // The data struct is as follows:
  // Length1  -->        unsigned long
  // TOA      -->        double
  // PhaseResid   -->    double    
  // TimeResid    -->    double 
  // OrbitalPhase -->    double 
  // ObsFreq      -->    double
  // FitWeight    -->    double
  // TimingUncertainty-->double
  // PrefitTimeResid --> double
  // Unused        -->   double
  // Length2   -->       unsigned long

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"

void fit_tempo(char *frequency, char *freqderiv);
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,double *b, double *siga, double *sigb, double *chi2, double *q);
double gammq(double a, double x);
//double rebin_phase(double number);
void errors(int *index_min, int *index_max);
//Global variables
char frequency[25],clean[1000], freqderiv[25];
double *var_doub,*res_phase,*res_time,*time_res,*err_time,*flux,*res,*sig,*time;
double period,chi2_array[100000],freq_array[100000],f1_array[100000],minchi2;
int    dof,num,num_f1;  //degrees of freedom.

int main()

{
  
  FILE *file,*output,*fitfile,*chisq,*fit_true,*bestvalue;
  unsigned long fileLen, vari[100000], varf[100000];
  int i,index_min, index_max;
  double truea,trueb,truesiga,truesigb;
  double T0,freq,truefreq,f1,truef1;
  char command1[1000],command2[1000];
  double a,b,siga,sigb,chi2,q;
  int mwt=1;
  
  sprintf(clean,"rm chisq.dat");
  system(clean);
  minchi2=1.e+8;
  bestvalue=fopen("bestvalue.dat","w");
 
  f1=-2e-13;

  dof=283;
  for(num_f1=0;num_f1<100;num_f1++)
   {
     freq=400.9752100;
      for(num=0;num<100;num++)
	{ 
    
	  file=fopen("resid2.tmp","rb");
	  output=fopen("output.asc","w");
	  chisq=fopen("chisq.dat","a");
	  
      if (!file)
	{
	  fprintf(stderr, "Unable to open file %s", file);
	  return;
	}
      sprintf(frequency,"%14.11lf",freq); 
      sprintf(freqderiv,"%e",f1);
      fit_tempo(frequency,freqderiv);
      
      
      //Get file length
      fseek(file, 0, SEEK_END);
      fileLen=ftell(file);
      fseek(file, 0, SEEK_SET);
      printf("File Length = %d\n",fileLen);
      
      //Allocate memory for the arrays
      var_doub =(double *) calloc(fileLen+1,sizeof(double));
      res_phase=(double *) calloc(fileLen+1,sizeof(double));
      res_time =(double *) calloc(fileLen+1,sizeof(double));
      time_res =(double *) calloc(fileLen+1,sizeof(double));
      err_time =(double *) calloc(fileLen+1,sizeof(double));
  
      // Store variables of interest in the proper arrays
      for(i=0; i<fileLen/80; i++)
	{
	  //vari & varf are just garbage arrays
	  fread(vari, sizeof(int), 1, file);
	  fread(var_doub,sizeof(double),9,file);
	  fread(varf,sizeof(int),1,file);
	  res_time[i]=var_doub[0];
	  res_phase[i]=var_doub[1];
	  time_res[i]=var_doub[2];
	  err_time[i]=var_doub[6];
	}
      
      // calculate T0 as the mean of the TOAs
      for(i=0; i<fileLen/80; i++)
	T0+=res_time[i];
      T0=T0/(fileLen/80);
      
  
      // Print output in the proper format
      period=time_res[0]/res_phase[0];  
      for(i=0; i<fileLen/80; i++)
	{
	  fprintf(output,"%lf %lf %lf\n", res_time[i],res_phase[i],err_time[i]*1e-6/period);
	}

      // close and free the pointers
      fclose(file);
      fclose(output);
      free(var_doub);
      free(res_phase);
      free(res_time);
      free(time_res);
      free(err_time);
      sprintf(command1,"paste flux.asc output.asc > fit");
      system(command1);
      sprintf(command2,"cat fit | awk '{print $1,$2,$4,$5}' > fit.asc");
      system(command2);

      //Start the fitting procedure 
      time =(double *) calloc(fileLen+1,sizeof(double));
      flux =(double *) calloc(fileLen+1,sizeof(double));
      sig  =(double *) calloc(fileLen+1,sizeof(double));
      res  =(double *) calloc(fileLen+1,sizeof(double));
      fitfile=fopen("fit.asc","r");
      for(i=0;i<dof;i++)
	fscanf(fitfile,"%lf %lf %lf %lf",&time[i],&flux[i],&res[i],&sig[i]);
      fclose(fitfile);
      fit(flux, res, dof, sig, mwt, &a, &b, &siga, &sigb, &chi2, &q);
      
      chi2_array[num]=chi2;
      fprintf(chisq,"%lf %lf %lf %lf %lf %lf %15.13lf %e\n",chi2,q,a,b,siga,sigb,freq,f1);
      if(chi2<minchi2)
	{
	  minchi2=chi2;
	  truefreq=freq;
	  truef1=f1;
	  truea=a;
	  trueb=b;
	  truesiga=siga;
	  truesigb=sigb;
	  fit_true=fopen("fit_true.dat","w");
	  for(i=0;i<fileLen/80;i++)
	    fprintf(fit_true,"%lf %lf %lf %lf\n",time[i],flux[i],res[i],sig[i]);
	  fclose(fit_true);
	}
      
      freq_array[num]=freq;
      f1_array[num_f1]=f1;
      free(time);
      free(flux);
      free(sig);
      free(res);
      fclose(chisq);
      freq+=2.0e-8;
	}
      f1+=5.e-14;
   }

 errors(&index_min, &index_max);
 printf("   Minimum chisqare=%lf       True frequency=%13.12lf   True F1=%e\n",minchi2,truefreq,truef1);
 fprintf(bestvalue,"#Minimum chi2    True frequency     True F1       A       B               sigA     sigB\n");
 fprintf(bestvalue,"     %lf  %15.13lf      %e     %lf   %lf    %lf     %lf\n",minchi2,truefreq,truef1,truea,trueb,truesiga,truesigb);
 fclose(bestvalue);
 return; 
  
}

void fit_tempo(char *frequency, char *freqderiv)
{
  char shellcom1[1000]="PSR       1808\nRAJ      18:08:27\nDECJ   -36:58:43\nPEPOCH      50915\nCLK        UTC(NIST)\nBINARY          ELL1\nA1       0.06281496\nPB       0.08390224\nTASC    50914.8784364\nEPS1    0\nEPS2      0\nEPHEM     DE405\nDM 0\n";

  char shellcom2[150],shellcom2b[150];
  char shellcom3[115];
  char timfile[]="1808.tim";
  FILE *input; 
  input=fopen("1808.par","w");
  
  strcat(shellcom2,"F0             ");
  strncat(shellcom2,frequency,16);
  fprintf(input,"%s",shellcom1);
  fprintf(input,"%s\n",shellcom2);
  strcat(shellcom2b,"F1             ");
  strncat(shellcom2b,freqderiv,16);
  fprintf(input,"%s\n",shellcom2b);
  fclose(input);
  sprintf(shellcom3,"tempo %s",timfile);
  system(shellcom3);

}


void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,double *b, double *siga, double *sigb, double *chi2, double *q)
{
	double gammq(double a, double x);
	int i;
	double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=0;i<ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=0;i<ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=0;i<ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=0;i<ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	*q=1.0;
	if (mwt == 0) {
		for (i=0;i<ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=0;i<ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}
#undef NRANSI

double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}




void errors(int *index_min, int *index_max) // calculate 68% confidence intervals (deltaChi2=1.00)
{

  int i,n;
  
  n=0;
  for(i=0;i<num;i++)
    {
      if(chi2_array[i]<minchi2+1.00)
	{
	  printf("Sup. Frequency= %13.11lf\n",freq_array[i-1]);
	  *index_min=i-1;
	  break;
	}
      n++;
    }
  for(i=n+1;i<num;i++)
    {
      if(chi2_array[i]>minchi2-1.00)
	{
	  printf("Inf. Frequency= %13.11lf\n",freq_array[i+1]);
	  *index_max=i+1;
	  break;
	}
    }
  
  n=0;
  for(i=0;i<num_f1; i++)
    {
      if(chi2_array[i]<minchi2+1.00)
	{
	  printf("Sup. Freq. Derivative= %e\n",f1_array[i-1]);
	  *index_min=i-1;
	  break;
	}
      n++;
    }
   for(i=n+1;i<num_f1;i++)
    {
      if(chi2_array[i]>minchi2-1.00)
	{
	  printf("Inf. Freq. Derivative= %e\n",f1_array[i+1]);
	  *index_max=i+1;
	  break;
	}
    }
  
 

  return;
}
