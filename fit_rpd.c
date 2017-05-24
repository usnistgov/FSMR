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
#endif

#include "func.c"

//#define N 60

void print_state (size_t iter, gsl_multifit_fdfsolver * s, int p);

int main (int argc, char *argv[]){

  #ifdef _OPENMP
  omp_set_num_threads(thread_num);
  #endif

  clock_t begin, end;
  double time_spent;

  begin = clock();

	char fName[2000];
	strncpy(fName, argv[4+atoi(argv[1])+atoi(argv[2])+atoi(argv[3])], 2000);
	
//Read data from file
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
		sigma[ndata] = sqrt(mindata*mindata/data);
		ndata++;
	}

	fclose(ifp);

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = ndata;
  int pp[3] ={atoi(argv[1]),atoi(argv[2]),atoi(argv[3])};
  if(pp[2]>1){
   printf("WARNING! Number of Poisson modes is greater than 1! Please, set the number of Poisson modes to 0 or 1 for better reconstruction accuracy.\n");
  }
  const size_t p = pp[0]+pp[1]+pp[2];
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	
	struct data d = { n, y, sigma, pp};
	gsl_multifit_function_fdf f;
	double * x_init = (double *)malloc(p*sizeof(double));
  for(i=0;i<p;i++){
    x_init[i]=atof(argv[4+i]);
  }
	gsl_vector_view x = gsl_vector_view_array (x_init, p);

  f.f = &mode_reconstruction_RPD_f;
  f.df = 0;//&expb_df;
	f.fdf = 0;//&expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	print_state (iter, s,p);

	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		printf ("status = %s\n", gsl_strerror (status));
		print_state (iter, s,p);
		if (status)
		break;
		status = gsl_multifit_test_delta (s->dx, s->x, 1e-14, 1e-14);
	}
	
	while (status == GSL_CONTINUE && iter < 50000);
	
	gsl_multifit_covar (s->J, mindata, covar);
	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	{
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));
		printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
    for(i=0;i<p;i++){
      if(i<pp[0]){
        printf("mu%d (Thermal) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }
      if((pp[0]<=i)&&(i<pp[0]+pp[1])){
        printf("mu%d (Single) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      } 
      if(pp[0]+pp[1]<=i){
        printf("mu%d (Poisson) = % 15.10f +/-% 15.10f\n",i,FIT(i),c*ERR(i));
      }

    }
	}

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	printf ("status = %s in %lf sec\n", gsl_strerror (status),time_spent);
	
  char * cptr;
  DIR* dir = opendir("recs");
  if (dir)
  {
    closedir(dir);
    sprintf(fName,"./recs/fit_RPD_%d%d%d",atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
  }
  else
  {
    sprintf(fName,"./fit_RPD_%d%d%d",atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
  }
  for(i=0;(i<p)&&(i<12);i++){
    cptr=fName+strlen(fName);
  	sprintf(cptr,"_mu%i_%.5f",i+1,FIT(i));
  }
  cptr=fName+strlen(fName);
 	sprintf(cptr,".dat");
	ofp = fopen(fName, "w+t");
  for(i=0;i<ndata;i++){
  //fprintf(ofp,"%lf\n",pow(gsl_vector_get(s->f, i)*sigma[i]+pow(y[i],0.2),5));
    fprintf(ofp,"%lf\n",pow(gsl_vector_get(s->f, i)*sigma[i]+pow(y[i],0.5),2));
	}
	fclose(ofp);

	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
  free(y);
  free(sigma);

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

