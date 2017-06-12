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
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#ifdef _OPENMP
# include <omp.h>
#define thread_num 8
#endif

//All necessary functions for JPD generation
#include "func.c"

int main (int argc, char *argv[]){
//Counters
  int n,m,i,j;
//Read input parameters
//Array of number of modes {NcTh,NcSp,NcP,NlTh...}
  int pp[9]={atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9])};
//Total number of modes + 2 efficiency losses
  int p = pp[0]+pp[1]+pp[2]+pp[3]+pp[4]+pp[5]+pp[6]+pp[7]+pp[8]+2;
  if(pp[2]>1||pp[5]>1||pp[8]>1){
    printf("WARNING! Number of Poisson modes is greater than 1! Please, set the number of Poisson modes to 0 or 1 for better reconstruction accuracy.\n");
  }
//JPD size nn*nn
  int nn=atoi(argv[10]);
//Number of states to contubute in JPD due to the losses
  int nmax=3*nn;
//Vector of mode parameters
  double * x_init = (double *)malloc(p*sizeof(double));
  for(i=0;i<p;i++){
    x_init[i]=atof(argv[11+i]);
  }
  //gsl_vector_view x = gsl_vector_view_array (x_init, p);
//Output JPD array of length nn*nn
  double *z = (double *)malloc(nn*nn*sizeof(double));

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generage JPD
  gen_jpd_z(x_init, pp, nn, z);
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//Create filename incuding parameters
  char fName[200];
  char * cptr;
  DIR* dir = opendir("data");
  if (dir)
  {
    closedir(dir);
    sprintf(fName,"./data/gen_JPD%d_%d%d%d_%d%d%d_%d%d%d",nn,atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
  }
  else
  {
    sprintf(fName,"./gen_JPD%d_%d%d%d_%d%d%d_%d%d%d",nn,atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
  }

  for(i=0;i<p-2;i++){
    cptr=fName+strlen(fName);
  	sprintf(cptr,"_mu%i_%.5f",i+1,x_init[i]);
  }
  cptr=fName+strlen(fName);
 	sprintf(cptr,"_etaL_%.5f",x_init[p-2]);
  cptr=fName+strlen(fName);
 	sprintf(cptr,"_etaR_%.5f",x_init[p-1]);
  cptr=fName+strlen(fName);
 	sprintf(cptr,".dat");
  FILE *ofp;
	ofp = fopen(fName, "w+t");
////////////////////
//Write JPD in file
	for(n=0;n<nn;n++){
    for(m=0;m<nn;m++){
      i=n*nn+m;
			fprintf(ofp,"%.16lf\n",z[i]);
    }
	}

  free(z);
  free(x_init);

return 0;
}
