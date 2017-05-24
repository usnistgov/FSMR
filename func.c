/*   This file is part of Full Statistical Mode Reconstruction Project
 * 
 *   Copyright (C) 2017 Ivan Burenkov - All Rights Reserved
 *   You may use, distribute and modify this code under the
 *   terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////
struct data {
	size_t n;
	double * y;
	double * sigma;
  int *pp;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Thermal mode distribution
//mu - average phtotn number [0.0 ... MAX_DOUBLE]
//*P - distrbution array
//nmax - size of disstribution array (max phtotn number-1)
inline int Pt(double mu, double *P, int nmax){
	int k;
	for(k=0;k<nmax;k++){
	P[k]=1/(1 + mu)*pow((mu/(1 + mu)),k);
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Poisson mode distribution
//mu - average phtotn number [0.0 ... MAX_DOUBLE]
//*P - distrbution array
//nmax - size of disstribution array (max phtotn number-1)
inline int Pp(double mu, double *P, int nmax){
	int k;
  P[0]=exp(-mu);
	for(k=1;k<nmax;k++){
	P[k]=P[k-1]*mu/k;
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Single photon mode distribution
//mu - average phtotn number [0.0 ... 1.0]
//*P - distrbution array
//nmax - size of disstribution array (max phtotn number-1)
inline int Ps(double mu, double *P, int nmax){
	int k;
  P[0]=1.0-mu;
  P[1]=mu;	
  for(k=2;k<nmax;k++){
	P[k]=0.0;
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Vacuum mode distribution
//*P - distrbution array
//nmax - size of disstribution array (max phtotn number-1)
inline int P0(double *P, int nmax){
	int k;
  P[0]=1.0;
  for(k=1;k<nmax;k++){
	P[k]=0.0;
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int Max(int a, int b){
	if(a>=b){
		return a;
	} else {
		return b;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int Min(int a, int b){
	if(a<=b){
		return a;
	} else {
		return b;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add another mode distribution PM to existiong mode distribution P
//*P - distrbution array
//*PM - distrbution array
//nmax - size of disstribution array (max phtotn number-1)
inline int PAdd(double *P, double *PM, int nmax){
  int nn,i;
  int nm=nmax;
  double * curP=(double *)malloc(nm*sizeof(double)); 
	for (nn=0;nn<nm;nn++){
    curP[nn]=0.0;
    for(i=0;i<=nn;i++){
      curP[nn]+=P[i]*PM[nn-i];
    }
  }    
  for (nn=0;nn<nm;nn++){
    P[nn]=curP[nn];
  }
  free(curP);
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Loss factors
//eta - efficiency [0.0 ... 1.0] (1-Loss)
//*Loss - loss factors array
//nm - size of disstribution array (max phtotn number-1)
//nmax - (max phtotn number-1) state to contribute in distribution due to losses 
inline int Losses(double eta, double *Loss, int nm, int nmax){
	int n,ipk,q,p;
#ifdef _OPENMP
#pragma omp parallel for private(n,ipk,q,p)
#endif
	for(n=0;n<nm;n++){
		for(ipk=n;ipk<2*nmax;ipk++){
			Loss[n*2*nmax+ipk]=1.0;
			for(q=0;q<ipk-n;q++){
				Loss[n*2*nmax+ipk] *= (1.0-eta);
			}
			for(p=0;p<n;p++){
				Loss[n*2*nmax+ipk] *= eta*(ipk-n+p+1.0)/(p+1.0);
			}
			//printf("%d\t%d\t%.20g\n",n,ipk,LossL[n*2*nmax+ipk]);
		}
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generate JPD
//*x vector of parameters of lenght (\Sum pp[i]) + 2 [Number of modes + 2 efficiency losses]  
//pp[9] matrix of mode numbers vs arms (Conjugated,Signal,Idler) vs type (Thermal,Single Photon,Poisson) 
//nn - size of JPD matrix nn*nn
//*z - JPD array of length nn*nn
inline int gen_jpd_z (double * x_init, int * pp, int nn, double * z){

//Total number of fitting (generating) parameters
  int p = pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7]+pp[8]+2;
//Number of states to contubute in JPD due to the losses
	int nmax=3*nn;
//Counters
  int i,ii,n,m,j,k,s,q,r,ipk;
//Array of parameters
  double *mu=(double *)malloc((p-2)*sizeof(double));
  for(j=0;j<p-2;j++){
    mu[j]=x_init[j];
  }
  double etaL=x_init[p-2];
  double etaR=x_init[p-1];
//Raw probability distributions for All represented modes
  double *PM = (double *)malloc(((pp[0]+pp[1]+pp[2])*nmax+(p-2-pp[0]-pp[1]-pp[2])*nn)*sizeof(double));
//Correlated and Signal/Idler PDs
  double *P_c = (double *)malloc(nmax*sizeof(double));
  double *P_s = (double *)malloc(nn*sizeof(double));
  double *P_i = (double *)malloc(nn*sizeof(double));
  double *LossL = (double *)malloc(2*nmax*nn*sizeof(double));
	double *LossR = (double *)malloc(2*nmax*nn*sizeof(double));
//Initialize PDs
  P0(P_c,nmax);
  P0(P_s,nn);
  P0(P_i,nn);

  int pc=0;
//Add correlated modes
  for(j=0;j<pp[0];j++){
    //Thermal
    Pt(mu[j+pc],PM+(j+pc)*nmax,nmax);
    PAdd(P_c,PM+(j+pc)*nmax,nmax);
  }
  pc+=pp[0];
  for(j=0;j<pp[1];j++){
    //Single photon
    Ps(mu[pc+j],PM+(pc+j)*nmax,nmax);
    PAdd(P_c,PM+(j+pc)*nmax,nmax);
  }
  pc+=pp[1];
  for(j=0;j<pp[2];j++){
    //Poissonian
    Pp(mu[pc+j],PM+(pc+j)*nmax,nmax);
    PAdd(P_c,PM+(j+pc)*nmax,nmax);
  }
  pc+=pp[2];
//Add uncorrelated signal modes
  for(j=0;j<pp[3];j++){
    //Thermal
    Pt(mu[j+pc],PM+(j+pc)*nn,nn);
    PAdd(P_s,PM+(j+pc)*nn,nn);
  }
  pc+=pp[3];
  for(j=0;j<pp[4];j++){
    //Single photon
    Ps(mu[pc+j],PM+(pc+j)*nn,nn);
    PAdd(P_s,PM+(j+pc)*nn,nn);
  }
  pc+=pp[4];
  for(j=0;j<pp[5];j++){
    //Poissonian
    Pp(mu[pc+j],PM+(pc+j)*nn,nn);
    PAdd(P_s,PM+(j+pc)*nn,nn);
  }
  pc+=pp[5];
//Add uncorrelated idler modes
  for(j=0;j<pp[6];j++){
    //Thermal
    Pt(mu[j+pc],PM+(j+pc)*nn,nn);
    PAdd(P_i,PM+(j+pc)*nn,nn);
  }
  pc+=pp[6];
  for(j=0;j<pp[7];j++){
    //Single photon
    Ps(mu[pc+j],PM+(pc+j)*nn,nn);
    PAdd(P_i,PM+(j+pc)*nn,nn);
  }
  pc+=pp[7];
  for(j=0;j<pp[8];j++){
    //Poissonian
    Pp(mu[pc+j],PM+(pc+j)*nn,nn);
    PAdd(P_i,PM+(j+pc)*nn,nn);
  }
  
//Calculate Loss Factors
  Losses(etaL,LossL,nn,nmax);
	Losses(etaR,LossR,nn,nmax);

  double sum;
//Conjugated part of JPD
  double *P_cl = (double *)malloc(nn*nn*sizeof(double));

//JPD for conjugated part, including losses
#ifdef _OPENMP
#pragma omp parallel for private(n,m,k)
#endif
	for(n=0;n<nn;n++){
    for(m=0;m<nn;m++){
      P_cl[n*nn+m]=0.0;
      for(k=Max(m,n);k<nmax;k++){
        P_cl[n*nn+m]+=P_c[k]*LossL[n*2*nmax+k]*LossR[m*2*nmax+k];
			}
    }
  }
//Total JPD, losses are not considered for the background
#ifdef _OPENMP
#pragma omp parallel for private(n,m,k,j,i,ii,sum)
#endif
	for(n=0;n<nn;n++){
    for(m=0;m<nn;m++){
      sum=0.0;
      for(j=0;j<=m;j++){
			  for(i=0;i<=n;i++){
          sum+=P_s[i]*P_i[j]*P_cl[(n-i)*nn+m-j];
				}
			}
			ii=n*nn+m;
 			z[ii]=sum;
    }
	}

  free(PM);
  free(P_c);
  free(P_cl);
  free(P_i);
  free(P_s);
  free(mu);
  free(LossL);
  free(LossR);

	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate JPD for parameters *x and returns fit errors *f for *data for GSL function gsl_multifit_fdfsolver
inline int mode_reconstruction_f (const gsl_vector * x, void *data, gsl_vector * f){
  size_t nnsq = ((struct data *)data)->n;
  int nn = (int)sqrt(nnsq);
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *) data)->sigma;
  int *pp = ((struct data *)data)->pp;
//Total number of fitting (generating) parameters
  int p = pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7]+pp[8]+2;
//Array of parameters
  double *x_init=(double *)malloc(p*sizeof(double));
  int n,m,i,j;
  for(j=0;j<p;j++){
    x_init[j]=gsl_vector_get (x, j);
  }
//JPD array
  double *z = (double *)malloc(nn*nn*sizeof(double));
//Generage JPD
  gen_jpd_z(x_init, pp, nn, z);
//Calculated fit errors
	for(n=0;n<nn;n++){
    for(m=0;m<nn;m++){
      i=n*nn+m;
      //gsl_vector_set (f, i, (z[i] - y[i])/sigma[i]); //Linear scoring function
			gsl_vector_set (f, i, (pow(z[i],0.5) - pow(y[i],0.5))/sigma[i]);
    }
	}
  free(z);
 	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generate RPD
//*x vector of parameters of lenght (\Sum pp[i]) [Number of modes]  
//pp[3] vector of mode numbers vs type (Thermal,Single Photon,Poisson) 
//nn - length of RPD vector
//*z - RPD array of length nn
inline int gen_rpd_z (double * x_init, int * pp, int nn, double * z){

//Total number of fitting (generating) parameters
  int p = pp[0]+pp[1]+pp[2];
//Number of states to contubute in JPD due to the losses
	int nm=nn;
//Counters
  int i,ii,n,m,j,k,s,q,r,ipk;
//Array of parameters
  double *mu=(double *)malloc(p*sizeof(double));
  for(j=0;j<p;j++){
    mu[j]=x_init[j];
  }
//Raw probability distributions for All represented modes
  double *PM = (double *)malloc(p*nm*sizeof(double));
//Correlated and Left Right PDs
  double *Pc = z;
//Initialize PDs
  P0(Pc,nm);
  int pc=0;
//Add modes
  for(j=0;j<pp[0];j++){
    //Thermal
    Pt(mu[j+pc],PM+(j+pc)*nm,nm);
    PAdd(Pc,PM+(j+pc)*nm,nm);
  }
  pc+=pp[0];
  for(j=0;j<pp[1];j++){
    //Single photon
    Ps(mu[pc+j],PM+(pc+j)*nm,nm);
    PAdd(Pc,PM+(j+pc)*nm,nm);
  }
  pc+=pp[1];
  for(j=0;j<pp[2];j++){
    //Poissonian
    Pp(mu[pc+j],PM+(pc+j)*nm,nm);
    PAdd(Pc,PM+(j+pc)*nm,nm);
  }
  free(PM);
  free(mu);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate RPD for parameters *x and returns fit errors *f for *data for GSL function gsl_multifit_fdfsolver
inline int mode_reconstruction_RPD_f (const gsl_vector * x, void *data, gsl_vector * f){
  size_t n = ((struct data *)data)->n;
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *) data)->sigma;
	size_t i;
  int *pp = ((struct data *)data)->pp; //array of number of modes {N_Th,N_Sp,N_P}
//Total number of fitting (generating) parameters
  int p = pp[0]+pp[1]+pp[2];
//RPD length	
	int nm=n;
//Counters
	int ii,nn,m,j,k,s,q,r,ipk;
//Array of parameters
  double *mu=(double *)malloc(p*sizeof(double));
  for(j=0;j<p;j++){
    mu[j]=gsl_vector_get (x, j);
  }
//RPD array
  double *z = (double *)malloc(nm*sizeof(double));
//Generage JPD
  gen_rpd_z(mu, pp, nm, z);
//Calculated fit errors
	for(nn=0;nn<nm;nn++){
		ii=nn;
		//gsl_vector_set (f, ii, (z[nn] - y[ii])/sigma[ii]);
		//gsl_vector_set (f, ii, (pow(z[nn],0.2) - pow(y[ii],0.2))/sigma[ii]);
      gsl_vector_set (f, ii, (pow(z[nn],0.5) - pow(y[ii],0.5))/sigma[ii]);
	}
  free(mu);
	return GSL_SUCCESS;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
