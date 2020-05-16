/*****************************************************************
 *  Name       : unit.h                                           *
 *  Verwendung : Schnittstelle zu Praktikumsumgebung (MAPRA),     *
 *               Iterative Loesungsverfahren                      *
 *  Autor      : Y. Zhang, IGPM RWTH Aachen                       *
 *  Datum      : Feb 2014                                         *
 *****************************************************************/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

class Matrix;
class Vector;

extern const int num_examples;
/**
   Eingabe: int ex_id, Matrix a, Vektro x0, Vector b, double tol, int maxiter.
   Ausgabe: Matrix a, Vector x0, Vector b, double tol, int maxiter.
 */
void getExample(int ex_id, Matrix &a, Vector &x0, Vector &b, double &tol, int &maxiter);

/**
 Methode: 0 -> Jacobi, 1-> GS, 2->CG
 */
void checkSolution(Vector &x, int iterations, int method);
