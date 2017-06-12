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
#include <math.h>

int main (int argc, char *argv[]){

	char fName[200];
  char sfName[200];
  char ifName[200];
	strncpy(fName, argv[1], 200);
	
//Read data from file
	FILE *ifp, *ofp;
	char *mode = "r";
	ifp = fopen(fName, mode);

	if (ifp == NULL) {
		fprintf(stderr, "Can't open input file in.list!\n");
		exit(1);
	}

	double data;
	int i,j,ndata=0;
  while (fscanf(ifp, "%lf", &data) != EOF) {
		ndata++;
	}
  fclose(ifp);

  double * y = (double *)malloc(ndata*sizeof(double));

  ifp = fopen(fName, mode);
  ndata=0;
	while (fscanf(ifp, "%lf", &data) != EOF) {
		y[ndata] = data;
		ndata++;
	}
	fclose(ifp);
  double sum;
  int n=(int)sqrt(ndata);
  //Calculate and write in file Signal RPD
  printf("JPD data file: %s \n",fName);
  strncpy(sfName, fName, 200);
  sfName[strlen(sfName)]='s';
  printf("Signal RPD data file: %s \n",sfName);
  ofp = fopen(sfName, "w+t");
	for(i=0;i<n;i++){
    sum=0.0;
  	for(j=0;j<n;j++){
      sum+=y[n*i+j];
    }
		fprintf(ofp,"%.16lf\n",sum);
	}
	fclose(ofp);
  //Calculate and write in file Idler RPD
  strncpy(ifName, fName, 200);
  ifName[strlen(ifName)]='i';
  printf("Idler RPD data file: %s \n",ifName);
  ofp = fopen(ifName, "w+t");
	for(i=0;i<n;i++){
    sum=0.0;
  	for(j=0;j<n;j++){
      sum+=y[n*j+i];
    }
		fprintf(ofp,"%.16lf\n",sum);
	}
	fclose(ofp);


 return 0;
}
