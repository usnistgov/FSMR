/*   This file is part of Full Statistical Mode Reconstruction Project
 * 
 *   Original FSMR algorithm is created by Ivan Burenkov
 *   Nonlinear fitting is implemented with use of GSL library under GPL license.
 *
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <dirent.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#define thread_num 8
#else
#define thread_num 1
#endif

//All necessary functions for fitting
#include "func.c"

void print_state (size_t iter, gsl_multifit_fdfsolver * s, int p);

int main (int argc, char *argv[]){

//Set desired number of parallel CPU threads here (j or j+1 is recomended for j cores of CPU)
#ifdef _OPENMP
  omp_set_num_threads(thread_num);
#endif

//Timing
  clock_t begin, end;
  double time_spent;

  struct timeval start, stop;
  gettimeofday(&start, NULL);

  begin = clock();

// Number of modes of each type in each arm {c_th,c_p,c_sp,i_th,i_p,i_sp,s_th,s_p,s_sp}, where c is for conjugated, i is for idler, s is for signal
	int pp[9]={atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9])};

  if(pp[2]>1||pp[5]>1||pp[8]>1){
    printf("WARNING! Number of Poisson modes is greater than 1! Please, set the number of Poisson modes to 0 or 1 for better reconstruction accuracy.\n");
  }

//total number of fit parameters (number of modes +2 efficiencies (1-loss))
  const size_t p = pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7]+pp[8]+2;


	char fName[200];
  strncpy(fName, argv[9+p+1], 200);
	
//Read data from file (single column data table)
	FILE *ifp, *ofp;
	char *mode = "r";
	ifp = fopen(fName, mode);

	if (ifp == NULL) {
		fprintf(stderr, "Can't open input file in.list!\n");
		exit(1);
	}

	double data;
	int ndata=0;
  double  mindata=1.0;
  double norm=0.0;
//Find minimum non-zero value to define accuracy
  while (fscanf(ifp, "%lf", &data) != EOF) {
    if((data<mindata)&&(data!=0.0)){
    mindata=data;
  }
		ndata++;
    norm+=data;
	}
  fclose(ifp);
//Number of trials calculated from minimum non-zero value, assuming it corresponds to single measurement event
  printf("mindata = %0.15g\n Total numeber of trials = %g\n Norm = %g\n",mindata,1.0/mindata/norm,norm);
  double * y = (double *)malloc(ndata*sizeof(double));
  double * sigma = (double *)malloc(ndata*sizeof(double));
  ifp = fopen(fName, mode);
//Additional gsl parameters
  const gsl_rng_type * type;
  gsl_rng * r; 
  gsl_rng_env_setup();
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
//Initializing data
  ndata=0;
	while (fscanf(ifp, "%lf", &data) != EOF) {
		y[ndata] = data;
		sigma[ndata] = sqrt(mindata/data);
		ndata++;
	}
	fclose(ifp);

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = ndata; //total number of JPD values
  
//Output table header
  printf ("Iteration\t");
  for(i=0;i<p-2;i++){
    printf("mu_%d\t\t",i);
  }
  printf("eta_s\t\teta_i\t\t");
	printf("Error\n");

	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	
	struct data d = { n, y, sigma, pp};
	gsl_multifit_function_fdf f;
//fitting parameters space
	double * x_init = (double *)malloc(p*sizeof(double));
  for(i=0;i<p;i++){
    x_init[i]=atof(argv[10+i]);;
  }

	gsl_vector_view x = gsl_vector_view_array (x_init, p);
//Set fitting algorithm
	f.f = &mode_reconstruction_f; //Reference to the function used to calculate fit errors for the model
	f.df = 0; //&mode_reconstruction_df; -|- model Jacobian (NOT USED, Jacobian is internally computed using finite difference approximations of the function f, see GSL reference chapter 38.3)
	f.fdf = 0; //&mode_reconstruction_fdf; -|-
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	print_state (iter, s,p);
//Iterator
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		printf ("status = %s\n", gsl_strerror (status));
  	print_state (iter, s,p);
		if (status)
		break;
		status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 100);

//Results output	
	gsl_multifit_covar (s->J, mindata, covar);
	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));
		printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
    for(i=0;i<p-2;i++){
      if(i<pp[0]){
        printf("mu%d (Thermal Conjugated) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((pp[0]<=i)&&(i<pp[0]+pp[1])){
        printf("mu%d (Single Conjugated) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      } 
      if((pp[0]+pp[1]<=i)&&(i<pp[0]+pp[1]+pp[2])){
        printf("mu%d (Poisson Conjugated) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((i>=pp[0]+pp[1]+pp[2])&&(i<pp[0]+pp[1]+pp[2]+pp[3])){
        printf("mu%d (Thermal Signal) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((i>=pp[0]+pp[1]+pp[2]+pp[3])&&(i<pp[0]+pp[1]+pp[2]+pp[3]+pp[4])){
        printf("mu%d (Single Signal) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      } 
      if((i>=pp[0]+pp[1]+pp[2]+pp[3]+pp[4])&&(i<pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5])){
        printf("mu%d (Poisson Signal) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((i>=pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5])&&(i<pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6])){
        printf("mu%d (Thermal Idler) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((i>=pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6])&&(i<pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7])){
        printf("mu%d (Single Idler) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      } 
      if(pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7]<=i){
        printf("mu%d (Poisson Idler) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
    }
      printf("eta_Signal = % 15.10f +/-% 15.10f\n",FIT(p-2),c*ERR(p-2));
      printf("eta_Idler = % 15.10f +/-% 15.10f\n",FIT(p-1),c*ERR(p-1));
	}
//Get and show calculation times
  gettimeofday(&stop, NULL);

  double delta = ((stop.tv_sec  - start.tv_sec) * 1000000u + 
         stop.tv_usec - start.tv_usec) / 1.e6;

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC ;

	printf ("status = %s in %lf sec, CPU time = %lf sec\n", gsl_strerror (status),delta,time_spent);

//Print in file data fit
  char * cptr;
//Check if reconstruction DIR exists
  DIR* dir = opendir("recs");
  if (dir)
  {
    closedir(dir);
    sprintf(fName,"./recs/fit_JPD%d_%d%d%d_%d%d%d_%d%d%d",(int)sqrt(n),atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
  }
  else
  {
    sprintf(fName,"./fit_JPD%d_%d%d%d_%d%d%d_%d%d%d",(int)sqrt(n),atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
  }
  
  for(i=0;i<p-2;i++){
    cptr=fName+strlen(fName);
  	sprintf(cptr,"_mu%i_%.5f",i+1,FIT(i));
  }
  cptr=fName+strlen(fName);
  	sprintf(cptr,"_eta_l_%.5f",FIT(p-2));
  cptr=fName+strlen(fName);
  	sprintf(cptr,"_eta_r_%.5f",FIT(p-1));
  cptr=fName+strlen(fName);
 	sprintf(cptr,".dat");
 	ofp = fopen(fName, "w+t");

  for(i=0;i<ndata;i++){
			fprintf(ofp,"%lf\n",pow(gsl_vector_get(s->f, i)*sigma[i]+pow(y[i],0.5),2));
	}
	fclose(ofp);	


	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
  free(y);
  free(sigma);
  free(x_init);

	return 0;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s, int p){
	int i;
  printf ("iter: %3zu x =",iter);
  for(i=0;i<p;i++){
    printf("% 15.10f ",	gsl_vector_get (s->x, i));
  }
  printf("|f(x)| = %g\n",gsl_blas_dnrm2 (s->f));
}
