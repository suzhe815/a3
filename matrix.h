#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <vector>
#include <string>
#include "vector.h"


class Vector;

class Matrix
{
private:
	size_t Zeil, Spalt;
	double* Mat;


public:
	explicit Matrix(size_t z = 1, size_t s = 1);  // Konstruktor mit Laengenangabe. Typ-Umwandlung der Eingabedaten kann nur explizit duerchgefuehrt werden
	~Matrix() {};                    // Destruktor
	Matrix(const Matrix&);  // Kopierkonstruktor

	double& operator()(size_t, size_t);       // Zugriff auf Einträge des Vektors
	double operator()(size_t,size_t) const;  // Zugriff auf falls Vektor const

	Matrix& operator=(const Matrix&);   // Zuweisung
	Matrix& operator+=(const Matrix&);  // Zuweisungen mit arithm. Operation
	Matrix& operator-=(const Matrix&);
	Matrix& operator*=(double);
	Matrix& operator/=(double);

	Matrix& redim(size_t,size_t);  // neue Laenge festlegen
	size_t Zeilen() const { return Zeil; };   // zeilenanzahl
	size_t Spalten() const { return Spalt; }; // spaltenanzahl
	size_t Index(size_t, size_t);     // berechnet den index des eintrages (z,s) im vektor
	size_t Index(size_t, size_t) const;  // nur Ablesen des indexes

	double norm2() const;    // Euklidische Norm des Vektors
	double normMax() const;  // Maximum-Norm des Vektors

	static void matError(const char str[]);  // Fehlermeldung ausgeben

	friend Matrix operator+(const Matrix&, const Matrix&);  // Addition
	friend Matrix operator-(const Matrix&, const Matrix&);  // Subtraktion
	friend Matrix operator-(const Matrix&);                 // Vorzeichen

	friend Matrix operator*(const Matrix&, const Matrix&);  // matrixpordukt
	friend Matrix operator*(double, const Matrix&);         // Vielfache
	friend Matrix operator*(const Matrix&, double);
	friend Matrix operator/(const Matrix&, double);

	friend bool operator==(const Matrix&, const Matrix&);  // Vergleich
	friend bool operator!=(const Matrix&, const Matrix&);

	friend std::istream& operator>>(std::istream&, Matrix&);        // Eingabe
	friend std::ostream& operator<<(std::ostream&, const Matrix&);  // Ausgabe

	friend Vector operator*(const Matrix&, const Vector&);  // Matrix-Vektor-
	friend Vector operator*(const Vector&, const Matrix&);  // Multiplikation


};



#endif 

