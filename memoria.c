
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************************************************************/ 

int *vectorInt(int size)
{
	int *a;

	a = (int *) malloc(size * sizeof(int));
	return a;
}
/************************************************************************/ 

float *vectorFloat(int size)
{
	float *a;

	a = (float *) malloc(size * sizeof(float));
	return a;
}
/************************************************************************/ 

int **matrizInt(int fila, int colu)
{
	int i, **a;
	
	a = (int **) malloc(fila * sizeof(int *));
	for (i=0; i<fila; ++i)
		a[i] = (int *) malloc(colu * sizeof(int));

	return a;
}

/************************************************************************/

float **matrizFloat(int fila, int colu)
{
	int i;
	float **a;

	a = (float **) malloc(fila * sizeof(float *));
	for (i=0; i<fila; ++i)
		a[i] = (float *) malloc(colu * sizeof(float));

	return a;
}

/************************************************************************/ 

void borrarMatriz(void **matriz, int fila)
{
	int i;
	
	for (i=0; i < fila; ++i)
		free(matriz[i]);

	free(matriz);
}

/************************************************************************/ 

int ***cuboInt(int fila, int colu, int planos)
{
	int i, j;
	int ***a;

	a = (int ***) malloc(fila * sizeof(int **));
	for (i=0; i < fila; ++i)
	{
		a[i] = (int **) malloc(colu * sizeof(int *));
		for (j=0; j < colu; ++j)
			a[i][j] = (int *) malloc(planos * sizeof(int));
	}

	return a;
}
/************************************************************************/ 

float ***cuboFloat(int fila, int colu, int planos)
{
	int i, j;
	float ***a;

	a = (float ***) malloc(fila * sizeof(float **));
	for (i=0; i < fila; ++i)
	{
		a[i] = (float **) malloc(colu * sizeof(float *));
		for (j=0; j < colu; ++j)
			a[i][j] = (float *) malloc(planos * sizeof(float));
	}

	return a;
}
/************************************************************************/

void borrarCubo(void ***cubo, int fila, int colu)
{
	int i, j;
	
	for (i=0; i < fila; ++i)
	{
		for(j=0; j < colu; j++)
		{
			free(cubo[i][j]);
		}
		free(cubo[i]);
	}

	free(cubo);
}

/************************************************************************/
