/**************************************
*                                     *
* minfo.h                             *
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
#ifndef MINFOMODR_H
#define MINFOMODR_H

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define TOLOG2       1.4426f
#define MIN_PUNTOS   16

/* Test Chi Cuadrado 50% */
//#define TEST_CHI2_1  0.4549f
//#define TEST_CHI2_3  0.7887f
//#define TEST_CHI2_15 0.9559f

/* Test Chi Cuadrado 20% Original */
#define TEST_CHI2_1  1.6424f
#define TEST_CHI2_3  1.5472f
#define TEST_CHI2_15 1.2874f

/* Test Chi Cuadrado 15% */
//#define TEST_CHI2_1  2.1739f
//#define TEST_CHI2_3  1.8155f
//#define TEST_CHI2_15 1.3872f

/* Test Chi Cuadrado 10% */
//#define TEST_CHI2_1  2.7055f
//#define TEST_CHI2_3  2.0838f
//#define TEST_CHI2_15 1.4871f

/* Test Chi Cuadrado 5% */
//#define TEST_CHI2_1  3.8415f
//#define TEST_CHI2_3  2.6049f
//#define TEST_CHI2_15 1.6664f



/* Macro para comparar n�meros en punto flotante */
#define COMP(x,y) ( fabs((x) - (y)) < (10 * FLT_EPSILON) )

/** Funciones para Calcular MI mediante Algoritmo de Fraser (Caracteristicas Continuas) **/

/*****************************************************************************************/
/* FUNCIONES DE ORDENAMIENTO Y DE UTILIDAD ***********************************************/
/*****************************************************************************************/

void hpsortfloat(int n, double ra[], int idx[]);
void hpsortint(int n, int ra[], int idx[]);


void  connFeat(int x[], int y[], int f[], int n);
void  quantify(double x[], int qx[], int n);

/************ Modificado Mtesmer Alg. Fraser Elementos repetidos ***********/
void  quantifyRepeat(const double x[], int qx[], int n, int feature_index);


/* Funcion de comparacion para ordenacion ascendente por "qsort" */
int compara(const void *x, const void *y);

/*************** Calculo de Elementos Repetidos en un vector *****************/
int repeat(int *x, int lx, int ux);

/****************** Entrega el maximo elemento de un vector *******************/
/************* Ultimo elemento del vector ordenado de menor a mayor ***********/
int last(int *x, int ne);

/* Entrega en un vector los elementos distintos de un vector y el N� total
de elementos distintos */
int *Elements(int *feature, int n, int *nelem);
/*****************************************************************************/


/*****************************************************************************************/
/* ALGORITMO DE FRASER BIDIMENSIONAL, PARA DOS CARACTERISTICAS ***************************/
/*****************************************************************************************/

void  sDiv2D(int buf[], int x[], int lx, int ux, int ly, int uy);
float fraser2D(int a[], int x[], int lx, int ux, int ly, int uy);
float mutualInfoFF(int x[], int lx, int ux, int ly, int uy);


/*****************************************************************************************/
/* ALGORITMO DE FRASER DE CLASE-CARACTERISTICA UNIDIMENSIONAL ****************************/
/*****************************************************************************************/

void  sDiv1D(int buf[], int x[], int lx, int ux, int clase);
float fraser1D(int a[], int x[], int lx, int ux, int clase);
float mutualInfoCF(int x[], int lx, int ux, int clases);

/*****************************************************************************************/
/* CALCULO DE ENTROPIA *******************************************************************/
/*****************************************************************************************/

/* Calculo de Entropia mediante Fraser [ H(f) = I(f;f) ] */
double entropyFraser(int *feature, int n);

/* Calculo de Entropia mediante histograma */
float entropy(float f[], int n);

/*****************************************************************************************/


/* Modificacion MTT */
/*****************************************************************************************/
/****** Funciones para Calcular MI mediante Histogramas (Caracteristicas Discretas) ******/
/*****************************************************************************************/

/* Determina si una caracteristica tiene muchos valores repetidos */
/* Entrega (0 -> discreta) (1 -> continua) */
float discreteMeasureFeature(float *feature, int nrow);
//int discreteMeasureFeature(float *feature, int nrow, float criterion);

/* Determina cuan discretas o continuas son las columnas de una matriz */
//float *discreteMeasureBase(float **base, int nrow, int ncol);
//int *discreteMeasureBase(float **base, int nrow, int ncol, float criterion);

/* Calculo de Entropia */
float entropyDiscrete(int *f, int n);

/* Calculo de Matriz de Probabilidades para 2 vectores de elementos discretos */
float **histMatrix(int *feature1, int *feature2, int n, int *row, int *col);

/* Informacion Mutua Caracteristica Discreta - Caracteristica Discreta */
float mutualInfoDD(int *feature1, int *feature2, int n);

/* Calcula la MI entre una caracteristica discreta y otra continua mediante
el algoritmo de Fraser en 1D. Se debe ordenar previamente ambas caracteristicas */
float MinfoDF(int x[], int lx, int ux, int *FeatElements, int nelem);

/* Informacion Mutua Caracteristica Discreta - Caracteristica Continua */
/* La 1ra caracteristica debe ser discreta y la 2da continua */
float mutualInfoDF(int *disc_feature, int *cont_feature, int n);

#endif
