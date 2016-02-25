/***********************************************************************
*                             "memoria.h"	   						   *
*	Funciones de reserva de memoria para arreglos multidimensionales   *	
***********************************************************************/

/* Asigna memoria a un vector de int */
int *vectorInt(int size);

/* Asigna memoria a un vector de float */
float *vectorFloat(int size);

/* Asigna memoria a una matriz de int */
int **matrizInt(int fila, int colu);

/* Asigna memoria a una matriz de float */
float **matrizFloat(int fila, int colu);

/* Borra una matriz */
void borrarMatriz(void **matriz, int fila);

/* Asigna memoria a un cubo de int */
int ***cuboInt(int fila, int colu, int planos);

/* Asigna memoria a un cubo de float */
float ***cuboFloat(int fila, int colu, int planos);

/* Borra una matriz de 3 dimensiones */
void borrarCubo(void ***cubo, int fila, int colu);
