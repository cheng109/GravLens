/**/
/* xamoeba_test.c  -- C code is from NUMERICAL RECIPES IN C --
 *                       Additional code written by Keith Poole
 *                       September - October 2009
 *
   Uses Nedler-Mead downhill simplex method to find a minimum/maximum
   of a multi-dimensional function
                 |->amotry->func 
Main -> amoeba --|
                 |->func

                 Amoeba initializes the simplex and calls func for each
                 point in the simplex.  The simplex starts at the usual
                 triangluar definition plus the origin -- ndim +1.

                 func is the user defined function being
                 minimized/maximized;

                 Amotry -- (from description in NUMERICAL RECIPES)
                 Extrapolates by a factor, fac, through the face of the
                 simplex across from a high point, tries, and replaces
                 the high point if the new point is better.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NRANSI

#define MP 26
#define NP 25
#define nrowX 1000
#define ncolX 25
#define FTOL 0.0000001
#define TINY 1.0e-10
#define NMAX 100000
#define NR_END 1
#define FREE_ARG char*

void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk);
double func(double *);
double amotry(double **p, double y[], double psum[], int ndim,
	     double (*funk)(double []), int ihi, double fac);
double keithrules(double x[]);

double func(double x[])
{
	return keithrules(x);
}
static double *Y,*X;
FILE *fp;
FILE *jp;

int main(void)
{
	int i,nfunc,j,ndim=NP;
	double *x,*y,**p;

	X = (double *) malloc (nrowX*ncolX*sizeof(double));
	Y = (double *) malloc (nrowX*sizeof(double));
	x = (double *) malloc ((NP+1)*sizeof(double));
	y = (double *) malloc ((MP+1)*sizeof(double));
/* Dynamic Allocation of a Matrix: */
/* First, allocate pointers to rows */
/* Second, allocate rows and set pointers to them -- note that p[0] is
 * a pointer to p[0][0:n] -- that is a row -- so the syntax below is
 * p[0][0:nrow*ncol] */
/* Third, form the matrix from step two -- the pointer to the kth row,
 * p[k] is set equal to the memory
 * location of the number of nrows-1 times the number of columns plus
 * the memory location of p[0] */
	p= (double **) malloc((MP+1)*sizeof(double*));
	p[0]=(double *) malloc((NP+1)*(MP+1)*sizeof(double));
	for (i=0;i<MP+1;i++)p[i]=p[0]+i*(NP+1);
//
	jp = fopen("c:/docs_c_summer_course/data_AMOEBA.txt","w");
	if((fp = fopen("c:/docs_c_summer_course/data_ols.txt","r"))==NULL)
	{
		printf("\nUnable to open file OLS_DATA.TXT: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	else {

		fprintf(jp," Y and X = \n");
		for(i=0;i<nrowX;i++)
		{
			fscanf(fp,"%lf",&Y[i]);
			for(j=0;j<ncolX;j++)
			{
				fscanf(fp,"%lf",&X[i+j*nrowX]);
			}
			fprintf(jp,"%10d %12.6f", i,Y[i]);
			for(j=0;j<ncolX;j++)
			{
				fprintf(jp,"%12.6f",X[i+j*nrowX]);
			}
			fprintf(jp,"\n");
		}
	}
/* if i == (j+1) is true then the value of x[j] = 1.0, otherwise =0.0 */
	for (i=1;i<=MP;i++) {
		for (j=1;j<=NP;j++)
			x[j]=p[i][j]=(i == (j+1) ? 1.0 : 0.0);
		y[i]=func(x);
	}
	amoeba(p,y,ndim,FTOL,func,&nfunc);
	printf("\nNumber of function evaluations: %3d\n",nfunc);
	fprintf(jp,"\nNumber of function evaluations: %3d\n",nfunc);
	printf("Vertices of final 3-d simplex and\n");
	fprintf(jp,"Vertices of final 3-d simplex and\n");
	printf("function values at the vertices:\n\n");
	fprintf(jp,"function values at the vertices:\n\n");
	printf("%3s %10s %12s %12s %14s\n\n",
	       "i","x[i]","y[i]","z[i]","function");
	for (i=1;i<=MP;i++) {
		printf("%3d ",i);
		fprintf(jp,"%3d ",i);
		for (j=1;j<=NP;j++) {
			fprintf(jp,"%12.6f ",p[i][j]);
		}
		printf("%12.6f\n",y[i]);
		fprintf(jp,"%12.6f\n",y[i]);
	}
	free(X);
	free(Y);
	free(x);
	free(y);
	free(p);
	return 0;
}

void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk)
{
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;

	psum = (double *) malloc ((ndim+1)*sizeof(double));
	*nfunk=0;
			for (j=1;j<=ndim;j++) {
		for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];
		psum[j]=sum;}
			for (;;) {
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
		if (rtol < ftol) {
			swap= y[1];
			y[1]=y[ilo];
			y[ilo]=swap;
			for (i=1;i<=ndim;i++)
			{
				swap=p[1][i];
				p[1][i]=p[ilo][i];
				p[ilo][i]=swap;
			}
					break;
		}
		if (*nfunk >= NMAX) 
		{
			printf("NMAX exceeded");
			exit(EXIT_FAILURE);
		}
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
						for (j=1;j<=ndim;j++) {
					for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];
					psum[j]=sum;}
			}
		} else --(*nfunk);
	}
	free(psum);
}


double amotry(double **p, double y[], double psum[], int ndim,
	     double (*funk)(double []), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry = (double *) malloc ((ndim+1)*sizeof(double));
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free(ptry);
	return ytry;
}

double keithrules(double x[])
{
	int i, j;
	double sum=0;
	double sumsquared=0;
	for(i=0;i<nrowX;i++)
	{
		sum=0.0;
		for(j=0;j<ncolX;j++)
		{
			sum=sum+X[i+j*nrowX]*x[j+1];
		}
		sumsquared = sumsquared + (Y[i]-sum)*(Y[i]-sum);
	}
	printf("%lf\n",sumsquared);
	return sumsquared;

}