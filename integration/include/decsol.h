/***************************************************************************
                                decsol.h
                             -------------------
    modified for C by:   : Blake Ashby
    last modified        : Nov 15, 2002
    email                : bmashby@stanford.edu

Helper routines for StiffIntegratorT class. Modified from original Fortran
code provided by:

	E. Hairer & G. Wanner
	Universite de Geneve, dept. de Mathematiques
	CH-1211 GENEVE 4, SWITZERLAND
	E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

 ***************************************************************************/

#ifndef _DECSOL_H_
#define _DECSOL_H_

// template<class T>
// inline const T& max(const T& a, const T& b)
// { return a > b ? a : b; }

// template<class T>
// inline const T& min(const T& a, const T& b)
// { return a < b ? a : b; }

// Matrix Triangularization by Gaussian Elimination
int dec(const int n, NB_TYPE **A, int *ip);

// Solution of linear system A*x = b
void sol(const int n, NB_TYPE **A, NB_TYPE *b, int *ip);

// Matrix Triangularization by Gaussian Elimination of a Hessenberg
// matrix with lower bandwidth lb
int dech(const int n, NB_TYPE **A, int lb, int *ip);

// Solution of linear system A*x = b -- Hessenberg matrix
void solh(const int n, NB_TYPE **A, int lb, NB_TYPE *b, int *ip);

// Matrix Triangularization by Gaussian Elimination for complex matrices
int decc(const int n, NB_TYPE **AR, NB_TYPE **AI, int *ip);

// Solution of linear system A*x = b -- complex matrices
void solc(const int n, NB_TYPE **AR, NB_TYPE **AI, NB_TYPE *br,
	NB_TYPE *bi, int *ip);

// Matrix Triangularization by Gaussian Elimination -- Hessenberg, complex
// matrices
int dechc(const int n, NB_TYPE **AR, NB_TYPE **AI, int lb, int *ip);

// Solution of linear system A*x = b -- Hessenberg, complex matrices
void solhc(const int n, NB_TYPE **AR, NB_TYPE **AI, int lb,
	NB_TYPE *br, NB_TYPE *bi, int *ip);

//Matrix Triangularization by Gaussian Elimination -- banded matrix
int decb(const int n, NB_TYPE **A, int ml, int mu, int *ip);

// Solution of linear system A*x = b -- banded matrix
void solb(const int n, NB_TYPE **A, int ml, int mu, NB_TYPE *b, int *ip);

//Matrix Triangularization by Gaussian Elimination -- banded, complex matrices
int decbc(const int n, NB_TYPE **AR, NB_TYPE **AI, int ml, int mu, int *ip);

// Solution of linear system A*x = b -- banded, complex matrices
void solbc(const int n, NB_TYPE **AR, NB_TYPE **AI, int ml, int mu,
	NB_TYPE *br, NB_TYPE *bi, int *ip);

// reduces a submatrix to upper Hessenberg form
void elmhes(const int n, int low, int igh, NB_TYPE **A, int *inter);

#endif /* _DECSOL_H_ */
