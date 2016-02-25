/**************************************
*                                     *
* minfo.c                             *
*                                     *
***************************************
*                                     *
* Calcula la informaci�n mutua entre  *
* dos variables aleatorias            *
*                                     *
***************************************
* (c) 2002 Pedro Ortega               *
* peortega@dcc.uchile.cl			  *
*									  *	
* Modificacion Michel Tesmer		  *
* mtesmer@ing.uchile.cl				  *
* 04/2003							  *
*									  *									
* Modificacion algoritmo de Fraser	  *
* para elementos repetidos en una	  *
* caracteristica continua			  *	
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "fifoqueues.h"
#include "minfomodr.h"
#include "memoria.h"


/**************************************************************/

/* Este heapsort fue modificado para que cree */ 
/* simult�neamente un arreglo de �ndices de   */
/* permutaciones                              */

void hpsortfloat(int n, double ra[], int idx[])
{
	int i,ir,j,l;
	int ridx;
	double rra; if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra  = ra[--l];
			ridx = idx[l];
		} else {
			rra  = ra[ir];
			ridx = idx[ir];
			ra[ir]  = ra[1];
			idx[ir] = idx[1];
			if (--ir == 1) {
				ra[1]  = rra;
				idx[1] = ridx;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i] = ra[j];
				idx[i]= idx[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]  = rra;
		idx[i] = ridx;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

/**************************************************************/

/* Este heapsort fue modificado para que cree */ 
/* simult�neamente un arreglo de �ndices de   */
/* permutaciones                              */

void hpsortint(int n, int ra[], int idx[])
{
	int i,ir,j,l;
	int ridx;
	int rra; if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra  = ra[--l];
			ridx = idx[l];
		} else {
			rra  = ra[ir];
			ridx = idx[ir];
			ra[ir]  = ra[1];
			idx[ir] = idx[1];
			if (--ir == 1) {
				ra[1]  = rra;
				idx[1] = ridx;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i] = ra[j];
				idx[i]= idx[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]  = rra;
		idx[i] = ridx;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

/**************************************************************/

/* Escribe dentro de f[] la funci�n f[i] = y[x[i]] */

void connFeat(int x[], int y[], int f[], int n)
{
	int *idx;
	int i;

	idx = (int *) malloc(sizeof(int) * (n+1));
	for (i=1; i<=n; ++i)
	{
		idx[i] = x[i];
		f[i]   = y[i];
	}
	hpsortint(n, idx, f);
	free(idx);
}

/**************************************************************/

/* Escribe en qx[] un ranking de los valores  */
/* del vector x[] en orden ascendente, que va */
/* de 0 hasta n-1                             */

void  quantify(double x[], int qx[], int n)
{
	double *wx;
	int i;
	int *ix;

	wx = (double *) malloc(sizeof(double) * (n+1));
	ix = (int *) malloc(sizeof(int) * (n+1));

	for (i=1;i<=n; ++i)
	{
		wx[i] = x[i];
		qx[i] = (i - 1);
		ix[i] = (i - 1);
	}

	hpsortfloat(n, wx, ix);
	hpsortint(n, ix, qx);

	free(wx);
	free(ix);

	/* mantener el n�mero de la caracter�stica */
	qx[0] = (int) x[0];

}

/************ Modificado Mtesmer Alg. Fraser Elementos repetidos ***********/
/* Escribe en qx[] un ranking de los valores  */
/* del vector x[] en orden ascendente, que va */
/* de 0 hasta n-1                             */

void  quantifyRepeat(const double x[], int qx[], int n, int feature_index)
{
	double *wx;
	int i, idx;
	int *ix;

	wx = (double *) malloc(sizeof(double) * (n+1));
	ix = (int *) malloc(sizeof(int) * (n+1));

	for (i=1; i <= n; ++i)
	{
		wx[i] = x[i-1];
		ix[i] = (i - 1);
	}
	hpsortfloat(n, wx, ix);
	
	/* Elementos repetidos */
	idx = 0;
	for (i=1; i <= (n-1); ++i)
	{
		if (wx[i] == wx[i+1])
		{
			qx[i] = idx;
  		}
  		else
  		{
  			qx[i] = idx;
  			idx++;
    	}
	}
	qx[n] = idx;
			
	hpsortint(n, ix, qx);

	free(wx);
	free(ix);

	/* mantener el n�mero de la caracter�stica */
	qx[0] = feature_index; //(int) x[0];
}

/*****************************************************************************/

int compara(const void *x, const void *y)
{
	float *a, *b;

	a = (float *) x;
	b = (float *) y;

	if(*a < *b)
		return -1;
	else if(*a > *b)
		return 1;
	else 
		return 0;
}

/**************************************************************/

int repeat(int *x, int lx, int ux)
{
	int i, *aux, nelem, Naux;

	Naux = ux - lx + 1;

	aux = vectorInt(Naux);
	for(i=lx; i <= ux; i++)
		aux[i-lx] = x[i];

	qsort(aux, Naux, sizeof(int), compara);
	
	nelem = 0;
	for(i=0; i < (Naux-1); i++)
		if( aux[i] != aux[i+1] )
			nelem++;
	nelem++;
	
	free(aux);
	
	return nelem;
}

/*****************************************************************************/

int last(int *x, int ne)
{
	int i, *aux, nelem;

	aux = vectorInt(ne);
	for(i=1; i <= ne; i++)
		aux[i-1] = x[i];

	qsort(aux, ne, sizeof(int), compara);
	nelem = aux[ne-1];
	
	free(aux);
	return nelem;
}

/*****************************************************************************/

int first(int *x, int ne)
{
	int i, *aux, nelem;

	aux = vectorInt(ne);
	for(i=1; i <= ne; i++)
		aux[i-1] = x[i];

	qsort(aux, ne, sizeof(int), compara);
	nelem = aux[0];

	free(aux);
	return nelem;
}
/*****************************************************************************/

int *Elements(int *feature, int n, int *nelem)
{
	int i, idx, *aux, *elements;

	aux = vectorInt(n);
	for(i=1; i <= n; i++)
		aux[i-1] = feature[i];

	qsort(aux, n, sizeof(int), compara);

	*nelem = 0;
	for(i=0; i < (n-1); i++)
		if( aux[i] != aux[i+1] )
			(*nelem)++;
	(*nelem)++;

	/* Elementos distintos en feature */
	elements = vectorInt(*nelem);
	idx = 0;
	for(i=0; i < (n-1); i++)
	{
		if( aux[i] != aux[i+1] )
		{
			elements[idx] = aux[i];
			idx++;
		}
	}
	elements[idx] = aux[n-1];

	free(aux);
	return elements;
}

/**************************************************************/

/**************************************************************/
/* PARTE BIDIMENSIONAL                                        */
/**************************************************************/

void sDiv2D(int buf[], int x[], int lx, int ux, int ly, int uy)
{
	int mx, my, i;

	mx = (int) ceil((lx + ux) / 2);
	my = (int) ceil((ly + uy) / 2);

	for (i=0; i < 4; i++)
 		buf[i] = 0;

	for (i=lx; i < mx; i++)
	{
		if (ly <= x[i] && x[i] < my)
			buf[0]++;
		else if(my <= x[i] && x[i] <= uy)
			buf[2]++;
	}
	for (i=mx; i <= ux; i++)
	{
		if (ly <= x[i] && x[i] < my)
			buf[1]++;
		else if(my <= x[i] && x[i] <= uy)
			buf[3]++;
	}
}

/**************************************************************/

float fraser2D(int a[], int x[], int lx, int ux, int ly, int uy)
{
	int* b[4];
	int  N;
	int  i, j, mx, my;
	double chi2_3, chi2_15;
	float mInfo;

	/* Crear el espacio para almacenar los coeficientes b[i] */
	b[0] = (int *) malloc(sizeof(int) * 4);
	b[1] = (int *) malloc(sizeof(int) * 4);
	b[2] = (int *) malloc(sizeof(int) * 4);
	b[3] = (int *) malloc(sizeof(int) * 4);

	/* Calcular los coeficientes bi */
	mx = (int) ceil((lx + ux) / 2);
	my = (int) ceil((ly + uy) / 2);
	
	sDiv2D(b[0], x, lx, mx-1, ly, my-1);
	sDiv2D(b[1], x, mx, ux,   ly, my-1);
	sDiv2D(b[2], x, lx, mx-1, my, uy);
	sDiv2D(b[3], x, mx, ux,   my, uy);

	/* Calcular los Chi-Cuadrado's */
	chi2_3 = chi2_15 = 0.0;
	N = a[0] + a[1] + a[2] + a[3];

	for (i=0; i<4; i++)
		chi2_3 += (a[i] - 0.25 * N) * (a[i] - 0.25 * N);
	chi2_3  *= 1.777778/N;

	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			chi2_15 += (float) (b[i][j] - 0.0625 * N) * (b[i][j] - 0.0625 * N);
	chi2_15 *= 1.137778/N;

	/* Hacer los tests y retornar el valor correspondiente */
	if (N <= MIN_PUNTOS || (chi2_3 < TEST_CHI2_3 && chi2_15 < TEST_CHI2_15))
	{
		/* no hay subestructura */
		mInfo = (float)((N==0)? 0.0 : N * (log(N)*TOLOG2));
	}
	else
	{
		/* hay subestructura */
		mInfo = (float) 2.0 * N
					+ fraser2D(b[0], x, lx, mx-1, ly, my-1)
					+ fraser2D(b[1], x, mx, ux,   ly, my-1)
					+ fraser2D(b[2], x, lx, mx-1, my, uy)
					+ fraser2D(b[3], x, mx, ux,   my, uy);
	}
	free(b[0]);
	free(b[1]);
	free(b[2]);
	free(b[3]);
	return mInfo;
}

/**************************************************************/

/* mutualInfo calcula la informaci�n mutua entre */
/* el vector x[] e y[] para la regi�n delimitada */
/* por (lx,ly)-(ux,uy)                           */
float mutualInfoFF(int x[], int lx, int ux, int ly, int uy)
{
	int *a;
	int N;
	float mInfo;

	/* Pedir espacio para coefientes */
	a = (int *) malloc(sizeof(int) * 4);

	/* Calcular los coeficientes a0 */
	sDiv2D(a, x, lx, ux, ly, uy);

	/* N�mero de ejemplos */
	N = a[0] + a[1] + a[2] + a[3];
	mInfo = (N == 0)? 0.0f: (float)(fraser2D(a, x, lx, ux, ly, uy) / N
 									- (log(1.0 * N)*TOLOG2));
	free(a);
	return mInfo;
}

/**************************************************************/
/* PARTE UNIDIMENSIONAL                                       */
/**************************************************************/

void sDiv1D(int buf[], int x[], int lx, int ux, int clase)
{
	int mx, i;

	mx = (int) ceil((lx + ux) / 2);

	for (i=0; i<2; i++)
 		buf[i] = 0;

	for (i=lx; i < mx; i++)
		if (x[i] == clase)
			buf[0]++;

	for (i=mx; i <= ux; i++)
		if (x[i] == clase)
			buf[1]++;
}

/**************************************************************/

float fraser1D(int a[], int x[], int lx, int ux, int clase)
{
	int* b[2];
	int  N;
	int  i, j, mx;
	double chi2_1, chi2_3;
	float mInfo;

	/* Crear el espacio para almacenar los coeficientes b[i] */
	b[0] = (int *) malloc(sizeof(int) * 2);
	b[1] = (int *) malloc(sizeof(int) * 2);

	/* Calcular los coeficientes bi */
	mx = (int) ceil((lx + ux) / 2);
	
	sDiv1D(b[0], x, lx, mx-1, clase);
	sDiv1D(b[1], x, mx, ux,   clase);

	/* Calcular los Chi-Cuadrado's */
	chi2_1 = chi2_3 = 0.0;
	N = a[0] + a[1];

	for (i=0; i<2; i++)
		chi2_1 += (a[i] - 0.5 * N) * (a[i] - 0.5 * N);
	chi2_1 *= 4.0/N;

	for (i=0; i<2; i++)
		for (j=0; j<2; j++)
			chi2_3 += (b[i][j] - 0.25 * N) * (b[i][j] - 0.25 * N);
	chi2_3 *= 1.77778/N;

	/* Hacer los tests y retornar el valor correspondiente */
	if (N <= MIN_PUNTOS || (chi2_1 < TEST_CHI2_1 && chi2_3 < TEST_CHI2_3))
	{
		/* no hay subestructura */
		mInfo = (float)((N==0)? 0.0 : N * (log(N)*TOLOG2));
	}
	else
	{
		/* hay subestructura */
		mInfo = N + fraser1D(b[0], x, lx, mx-1, clase)
			      + fraser1D(b[1], x, mx, ux,   clase);
	}
	free(b[0]);
	free(b[1]);
	return mInfo;
}

/**************************************************************/

/* mutualInfo calcula la informaci�n mutua entre */
/* el vector x[] e y[] para la regi�n delimitada */
/* por (lx,ly)-(ux,uy)                           */
float mutualInfoCF(int x[], int lx, int ux, int clases)
{
	int *a;
	int N;
	float mInfo;
	int i;

	/* Pedir espacio para coefientes */
	a = (int *) malloc(sizeof(int) * 2);

	/* Calcular la sumatoria de integrales */
	mInfo = 0.0;
	for (i=1; i<=clases; ++i)
	{
		int cCount;

		/* Calcular los coeficientes a0 */
		sDiv1D(a, x, lx, ux, i);

		/* N�mero de ejemplos */
		cCount = a[0] + a[1];
		N = ux - lx;

		/* Acumular la integral */
		mInfo += (cCount == 0)? 0.0f: (float)(- cCount * log(cCount) *TOLOG2
  											  + fraser1D(a, x, lx, ux, i));
	}
	mInfo /= N;
	free(a);
	return mInfo;
}

/**************************************************************/
/* PARTE ENTROPIA                                             */
/**************************************************************/

/* Modificacion MTT Calculo Entropia */
/* Calculo Entropia mediante Fraser [ H(f) = I(f;f) ] */
/******************************************************************************/

double entropyFraser(int *feature, int n)
{
	int i, max, min;
        int *feature_aux, *f;
	double entropy;
	
	feature_aux = (int *) malloc(sizeof(int) * (n+1));
	f = (int *) malloc(sizeof(int) * (n+1));
	
	for(i=0; i <= n; i++) {
		feature_aux[i] = feature[i];
	}
	connFeat(feature, feature_aux, f, n);
	max = last(f, n);
	min = first(f, n);
	if (max==min) { //Entropy is zero if the variable is fixed
		return 0.0;
	}
	entropy = mutualInfoFF(f, 1, n, 0, max);

	free(feature_aux);
	free(f);
	
	return entropy;
}

/******************************************************************************/


/* Calcula la entropia de una variable aleatoria mediante histograma */
/* Recibe como entrada un arreglo que representa */
/* f(x), pero transformada a enteros             */
/******************************************************************************/

float entropy(float f[], int n)
{
	double mi, max, min, delta;
	int  *h;
	int   ncubes, i;

	/* inicializar las variables */	
	mi = 0.0;
	ncubes = n / (MIN_PUNTOS / 2) + 1;
	h  = (int *) malloc(sizeof(int) * (ncubes+1));
	
	/* buscar el maximo para */
	/* hacer el histograma   */
	min = max = f[1];
	for (i=1; i <= n; ++i)
	{
		if (f[i] > max)
			max = f[i];
		if (f[i] < min)
			min = f[i];
	}
	delta = (max - min) / ncubes;
	
	/* armar el histograma */
	for (i=0; i < ncubes; ++i)
		h[i] = 0;
		
	for (i=1; i <= n; ++i)
		h[(int) ((f[i]-min) / delta)]++;
	
	/* calcular la entropia */
	for (i=0; i < ncubes; ++i)
		if (h[i] != 0)
			mi -= (double)h[i]/n * log((double)h[i]/n) * TOLOG2;
	
	return (float) mi;
}

/******************************************************************************/



/* Modificacion MTT */
/***********************************************************************************/
/********* Calculo de MI mediante Histogramas (Caracteristicas Discretas) **********/
/***********************************************************************************/

/***********************************************************************************/

float discreteMeasureFeature(float *feature, int nrow)
//int discreteMeasureFeature(float *feature, int nrow, float criterion)
{
	int i, j, *flag;
	float measure;

	flag = vectorInt(nrow+1);
	measure = 0.0f;

	for(i=1; i <= nrow; i++)
		flag[i] = 0;
		
	for(i=1; i <= (nrow-1); i++)
	{
		for(j=(i+1); j <= nrow; j++)
		{
			if(flag[i] == 1)
				break;
			
			if( (flag[j] == 0) && (feature[i] == feature[j]) )
			{
				measure++;
				flag[j] = 1;
			}
		}
	}
	measure /= nrow;
	
	free(flag);
	return measure;
	/*
	if (measure < criterion)
		return 1;		/* Caracteristica Continua */
	/*else
		return 0;		/* Caracteristica Discreta */
}

/***********************************************************************************/
/*
float *discreteMeasureBase(float **base, int nrow, int ncol)
//int *discreteMeasureBase(float **base, int nrow, int ncol, float criterion)
{
	int i, j;
// 	int *indicator;
 	float *indicator;
	float *feature;

	feature = vectorFloat(nrow+1);
	indicator = vectorFloat(ncol+1);
	//indicator = vectorInt(ncol);

	for(i=0; i < ncol ; i++)
	{
		for(j=0; j < nrow; j++)
			feature[i] = base[j][i];
		
		indicator[i] = discreteMeasureFeature(feature, nrow); 	
	}
	
	free(feature);
	return indicator;
}
*/
/************************************************************************************/


/***********************************************************************************/

float entropyDiscrete(int *f, int n)
{
	int i, cont, *hist, *aux;
 	float pbb, entropy;

	hist = vectorInt(n);
	aux  = vectorInt(n);
	for(i=1; i <= n; i++)
	{
		hist[i-1] = 0;
		aux[i-1]  = f[i];
	}

	cont = 0;
	entropy = 0.0f;

	qsort(aux, n, sizeof(int), compara);
	for(i=0; i < (n-1); i++)
	{
		if(aux[i] != aux[i+1])
		{
			hist[cont]++;
			cont++;
		}
		else
		{
			hist[cont]++;
		}
	}
	hist[cont]++;

	for(i=0; i <= cont; i++)
	{
		pbb = (float)hist[i] / n;
		entropy -= (float)(pbb * log(pbb) * TOLOG2);
	}

	free(hist);
	free(aux);
	return ((float) entropy);
}

/***********************************************************************************/

float **histMatrix(int *feature1, int *feature2, int n, int *row, int *col)
{
	int i, j, *f1, *f2, *elem1, *elem2, idx1, idx2;
	float **pbbmatrix;

	f1 = vectorInt(n+1);
	f2 = vectorInt(n+1);
	*row = *col = idx1 = idx2 = 0;
	
	for(i=1; i <= n; i++)
	{
		f1[i] = feature1[i];
		f2[i] = feature2[i];
	}

	/* Determinando N� columnas de matriz de pbb */
	hpsortint(n, f2, f1);
	for(i=1; i <= (n-1); i++)
		if( f2[i] != f2[i+1] )
			(*col)++;
	(*col)++;

	/* Elementos distintos en feature 2 */
	elem2 = vectorInt(*col);
	for(i=1; i <= (n-1); i++)
	{
		if( f2[i] != f2[i+1] )
		{
			elem2[idx2] = f2[i];
			idx2++;
		}
	}
	elem2[idx2] = f2[n];

	/* Determinando N� filas de matriz de pbb */
	hpsortint(n, f1, f2);
	for(i=1; i <= (n-1); i++)
		if( f1[i] != f1[i+1] )
			(*row)++;
	(*row)++;

	/* Elementos distintos en feature 1 */
	elem1 = vectorInt(*row);
	for(i=1; i <= (n-1); i++)
	{
		if( f1[i] != f1[i+1] )
		{
			elem1[idx1] = f1[i];
			idx1++;
		}
	}
	elem1[idx1] = f1[n];

	/* Inicializacion de matriz de probabilidades */
	pbbmatrix = matrizFloat(*row, *col);
	for(i=0; i < *row; i++)
		for(j=0; j < *col; j++)
			pbbmatrix[i][j] = 0.0f; 

	/* Calculo de probabilidades */
	for(i=1; i <= n; i++)
	{
		for(j=0; j < *row; j++)
		{
			if(f1[i] == elem1[j])
			{
				idx1 = j;
				break;
			}
		}
		for(j=0; j < *col; j++)
		{
			if(f2[i] == elem2[j])
			{
				idx2 = j;
				break;
			}	
		}
		pbbmatrix[idx1][idx2]++;
	}

	for(i=0; i < *row; i++)
		for(j=0; j < *col; j++)
			pbbmatrix[i][j] /= n;

	free(f1);
	free(f2);
	free(elem1);
	free(elem2);

	return pbbmatrix;
}

/***********************************************************************************/

float mutualInfoDD(int *feature1, int *feature2, int n)
{
	int i, j, k, l, row, col;
	float minfo, Pmarg1, Pmarg2, **pbbmatrix;

	minfo = 0.0f;
	pbbmatrix = histMatrix(feature1, feature2, n, &row, &col);

	for(i=0; i < row; i++)
	{
		for(j=0; j < col; j++)
		{
			if(pbbmatrix[i][j] != 0.0f)
			{
				Pmarg1 = Pmarg2 = 0.0f;

				/* Pbb Marginal de la 1ra componente */
				for(k=0; k < col; k++)
					Pmarg1 += pbbmatrix[i][k];
				
				/* Pbb Marginal de la 2da componente */
				for(l=0; l < row; l++)
					Pmarg2 += pbbmatrix[l][j];

				/* Calculo de MI */
				minfo += (float)(pbbmatrix[i][j] * 
								log(pbbmatrix[i][j] / (Pmarg1 * Pmarg2)) * TOLOG2);
			}
		}
	}

	return minfo;
}

/***********************************************************************************/

float MinfoDF(int x[], int lx, int ux, int *FeatElements, int nelem)
{
	int *a;
	int i, N;
	float mInfo;

	/* Pedir espacio para coefientes */
	a = vectorInt(2);

	/* Calcular la sumatoria de integrales */
	mInfo = 0.0;
	for(i=0; i < nelem; i++)
	{
		int cCount;

		/* Calcular los coeficientes a0 */
		sDiv1D(a, x, lx, ux, FeatElements[i]);

		/* N�mero de ejemplos */
		cCount = a[0] + a[1];
		N = ux - lx;

		/* Acumular la integral */
		mInfo += (cCount == 0)? 0.0f: (float)(- cCount * log(cCount) * TOLOG2
  											  + fraser1D(a, x, lx, ux, FeatElements[i]));
	}
	mInfo /= N;
	
	free(a);
	return mInfo;
}

/***********************************************************************************/

float mutualInfoDF(int *disc_feature, int *cont_feature, int n)
{
	int *f, nelem, *elements;
	float minfo;

	f = vectorInt(n+1);

	elements = Elements(disc_feature, n, &nelem);
	connFeat(cont_feature, disc_feature, f, n);
	minfo = MinfoDF(f, 1, n, elements, nelem);

	free(f);
	return minfo;
}

/***********************************************************************************/
