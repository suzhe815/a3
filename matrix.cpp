#include "matrix.h"
#include <iomanip>
#include <cmath>
#include <cstdlib>
#define NDEBUG

using namespace std;
/*
Matrix::Matrix(size_t z, size_t s)
{
#ifndef NDEBUG
	if (z <= 0 || s <= 0)
		matError("Nur Matrixen mit positiver Dimension!");
#endif

	this->redim(z, s);
	this->Zeil = z;
	this->Spalt = s;
}

//---- - Indexberechner---- -
// auf den ersten Eintrag der Matrix (Element linke Ecke oben) 
// greift man mit Mat(0,0) zu

// ----- Kopierkonstruktor -----
Matrix::Matrix(const Matrix& x) :Mat(x.Mat) {
	Zeil = x.Zeilen();
	Spalt = x.Spalten();
}*/

Matrix::Matrix(size_t z, size_t s)
{
#ifndef NDEBUG
	if (z <= 0 || s <= 0)
		matError("Nur Matrixen mit positiver Dimension!");
#endif
	Zeil = z;
	Spalt = s;

	Mat = new (std::nothrow) double[z * s];
	if (Mat == nullptr)
	{
		matError("nicht genuegend speicher!");

	}

	for (size_t i = 0; i< z*s; i++)
	{
		Mat[i] = 0;
	
}
}
//---- - Indexberechner---- -
// auf den ersten Eintrag der Matrix (Element linke Ecke oben) 
// greift man mit Mat(0,0) zu

// ----- Kopierkonstruktor -----

// ----- Kopierkonstruktor -----

Matrix::Matrix(const Matrix& x) : Mat(x.Mat)
{
	Zeil = x.Zeilen();
	Spalt = x.Spalten();
	Mat = new (std::nothrow) double[Zeil*Spalt];
	if (Mat == nullptr) {
		matError("Nicht genuegend Speicher!");
	}
	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			size_t d = Index(z, s);
			Mat[d] = x(z, s);
		}
	}
} 

size_t Matrix::Index(size_t z, size_t s) 
{
#ifndef NDEBUG
	if (z < 0 || s < 0)
		matError("Nicht gueltiger Index2");
#endif

	return (z * Spalt + s);
}

size_t Matrix::Index(size_t z, size_t s) const
{
#ifndef NDEBUG
	if (z < 0 || s < 0)
		matError("Nicht gueltiger Index2");
#endif

	return (z * Spalt + s);
}

// ----- Schreib- und Lesezugriff auf Matrixelemente -----

double& Matrix::operator()(size_t z, size_t s) {
#ifndef NDEBUG
	if (z < 0 || s < 0)
		matError("Nicht gueltiger Index3");
#endif


	return Mat[Index(z,s)];
}

// ----- Lesezugriff auf Elemente konstanter Matrizen -----

double Matrix::operator()(size_t z, size_t s) const {
#ifndef NDEBUG
	if (z < 0 || s < 0)
		matError("Nicht gueltiger Index4");
#endif
	
	double d = Mat[Index(z, s)];
	return d;
}

// =====================
//      Zuweisungen
// =====================

// ----- Zuweisungsoperator "=" ----

Matrix& Matrix::operator=(const Matrix& x) {
#ifndef NDEBUG
	if (Zeil != x.Zeilen() || Spalt != x.Spalten())
		matError("Inkompatible Dimensionen fuer 'Matrix = Matrix'!");
#endif

	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			(*this)(z, s) = x(z, s);
		}
	}

	return *this;
}

// ----- Zuweisungsoperator mit Addition "+=" ----

Matrix& Matrix::operator+=(const Matrix& x) {

#ifndef NDEBUG
	if ((Zeil != x.Zeilen()) || (Spalt != x.Spalten()))
		matError("Inkompatible Dimensionen fuer 'Matrix += Matrix'!");
#endif
	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			(*this)(z, s) += x(z, s);
		}
	}
	return *this;
}

// ----- Zuweisungsoperator mit Subtraktion "-=" ----

Matrix& Matrix::operator-=(const Matrix& x) {
#ifndef NDEBUG
	if ((Zeil != x.Zeilen()) || (Spalt != x.Spalten()))
		matError("Inkompatible Dimensionen fuer 'Matrix -= Matrix'!");
#endif
	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			(*this)(z, s) -= x(z, s);
		}
	}
	return *this;
}

// ----- Zuweisungsoperator mit Multiplikation "*=" ----

Matrix& Matrix::operator*=(double c) {
	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			(*this)(z, s) *= c;
		}
	}
	return *this;
}

// ----- Zuweisungsoperator mit Divsion "/=" ----

Matrix& Matrix::operator/=(double c) {
#ifndef NDEBUG
	if (abs(c) < epsilon)
	{
		vecError("Kann nicht durch 0 dividieren!");
	}
#endif

	for (size_t z = 0; z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			(*this)(z, s) /= c;
		}
	}
	return *this;
}

// ==============================
//      Matrixlaenge aendern
// ==============================

// ----- Matrixlaenge aendern -----

Matrix& Matrix::redim(size_t z, size_t s) {

#ifndef NDEBUG
	if (z <= 0 || s <= 0)
		matError("Nur Matrixen mit positiver Dimension!");
#endif

	delete[] Mat;
	Mat = new(nothrow) double[z*s];
	this->Zeil = z;
	this->Spalt = s;
	
	for (size_t i = 0; i < z*s; i++)
	{
		Mat[i] = 0;
	}
	return (*this);

}

// ======================
//      Matrixnormen
// ======================

// ----- Frobeniusnorm -----

double Matrix::norm2() const {

	double sum = 0;
	for (size_t z = 0;z < Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			sum += pow((*this)(z, s), 2);
		}
	}

	sum = sqrt(sum);
	return sum;
}

// ----- Maximumsnorm/Zeilensummennorm -----

double Matrix::normMax() const {

	double max = 0;
	double tmp = 0;
	for (size_t z = 0; z<Zeil; z++)
	{
		for (size_t s = 0; s < Spalt; s++)
		{
			double tmp = tmp + abs((*this)(z,s));
		}

		if (tmp > max)
		{
			max = tmp;
		}
	}

	return max;
}

// ==================================
//      arithmetische Operatoren
// ==================================

// ----- Addition "+" -----
//Unterschied zu =+ : hier wird die vorgegebene matrix nicht geaendert, sondern
// nur eine kopie zurueckgegeben.

Matrix operator+(const Matrix& x, const Matrix& y) {
#ifndef NDEBUG
	if (x.Zeilen() != y.Zeilen() || x.Spalten() != y.Spalten())
		Matrix::matError("Inkompatible Dimensionen fuer 'Matrix + Matrix'!");
#endif

	Matrix z = x;
	return z += y;
}

// ----- Subtraktion "-" -----

Matrix operator-(const Matrix& x, const Matrix& y) {
#ifndef NDEBUG
if (x.Zeilen() != y.Zeilen() || x.Spalten() != y.Spalten())
Matrix::matError("Inkompatible Dimensionen fuer 'Matrix + Matrix'!");
#endif

Matrix z = x;
return z -= y;
}

// ----- Vorzeichen wechseln "-" -----

Matrix operator-(const Matrix& x) {

	Matrix z = x;
	z *= -1;
	return z;
}

// ----- Matirxpordukt "*" -----

Matrix operator*(const Matrix& x, const Matrix& y) {
#ifndef NDEBUG
	if (x.Spalten() != y.Zeilen())
		Matrix::matError("Inkompatible Dimensionen fuer 'Matrix * Matrix'!");
#endif

	// speichere x*y auf matrix "Produkt"
	Matrix Produkt(x.Zeilen(), y.Spalten());

	double eintrag = 0;

	for (size_t z = 0; z < x.Zeilen();z++)
	{
		// Berechne den (z,i)-Eintrag fuer eine Spalte i
		for (size_t i = 0; i < y.Spalten(); i++)
		{
			// Berechne fuer alle Spalten
			for (size_t s = 0; s < x.Spalten(); s++)
			{
				eintrag = eintrag + x(z, s) * y(s, i);
			}

			Produkt(z, i) = eintrag;
			eintrag = 0;
		}
	}

	return Produkt;
}

// ----- Multiplikation Skalar*Matrix "*" -----

Matrix operator*(double c, const Matrix& x) {

	Matrix m = x;
	for (size_t z = 0;z < x.Zeilen(); z++)
	{
		for (size_t s = 0; s < x.Spalten(); s++)
		{
			m(z, s) = c * x(z, s);
		}
	}

	return m;
}

// ----- Multiplikation Matrix*Skalar "*" -----

Matrix operator*(const Matrix& x, double c) {

	Matrix m = x;
	for (size_t z = 0;z < x.Zeilen(); z++)
	{
		for (size_t s = 0; s < x.Spalten(); s++)
		{
			m(z, s) = c * x(z, s);
		}
	}

	return m;
}

// ----- Division Vektor/Skalar "/" -----

Matrix operator/(const Matrix& x, double c) {

	Matrix m = x;
	for (size_t z = 0;z < x.Zeilen(); z++)
	{
		for (size_t s = 0; s < x.Spalten(); s++)
		{
			m(z, s) = x(z, s) / c;
		}
	}

	return m;
}

//Matrix *Vektor
Vector operator * (const Matrix& x, const Vector& v)
{

#ifndef NDEBUG
if (x.Spalten() != v.length())
Matrix::matError("Inkompatible Dimensionen fuer 'Matrix * Vektor'!");
#endif

    Vector Produkt(x.Zeilen());
	double eintrag = 0;
	for (size_t z = 0; z < x.Zeilen();z++)
	{
		for (size_t s = 0; s < x.Spalten(); s++)
		{
			eintrag = eintrag + v(s) * x(z, s);
		}
		Produkt(z) = eintrag;
		eintrag = 0;
	}

	return Produkt;
}

//Vektor*Matrix  (Zeilenvektor * Martrix)
Vector operator * (const Vector& v,const Matrix& x)
{

#ifndef NDEBUG
	if (x.Zeilen() != v.length())
		Matrix::matError("Inkompatible Dimensionen fuer 'Vektor * Matrix'!");
#endif

	Vector Produkt(x.Spalten());
	double eintrag = 0;
	for (size_t z = 0; z < x.Spalten();z++)
	{
		for (size_t s = 0; s < v.length; s++)
		{
			eintrag = eintrag + v(s) * x(s, z);
		}
		Produkt(z) = eintrag;
		eintrag = 0;
	}

	return Produkt;
}



// ==============================
//      Vergleichsoperatoren
// ==============================



// ----- Test auf Gleichheit "==" -----

bool operator==(const Matrix& x, const Matrix& y) {
	if (x.Zeilen()!=y.Zeilen() || x.Spalten() !=y.Spalten()) {
		return false;
	}

	for (size_t z = 0; z < x.Zeilen(); z++) {
		for (size_t s = 0; s < x.Spalten(); s++)
		{
			if (x(z, s) != y(z, s))
			{
				return false;
			}
		}
	}

	return true;
}

// ----- Test auf Ungleichheit "!=" -----

bool operator!=(const Matrix& x, const Matrix& y) {

	return !(x == y);
}

// ==========================
//      Ein- und Ausgabe
// ==========================

// ----- Ausgabe "<<" -----

std::ostream& operator<<(std::ostream& s, const Matrix& x) {
	s << std::setiosflags(std::ios::right);
	s << "# Zeilen: " << x.Zeilen();
	s << "\n";
	s << "# Spalten: " << x.Spalten();
	s << "\n";
	for (size_t z = 0; z < x.Zeilen(); z++)
	{
		for (size_t m = 0; m < x.Spalten(); m++)
		{
			s << "  " << x(z, m);
		}
		s << "\n";
	}

	return s << std::endl;
}

// ----- Eingabe ">>" -----

std::istream& operator>>(std::istream& s, Matrix& x) {
	std::cout << std::setiosflags(std::ios::right);
	for (size_t z = 0; z < x.Zeilen(); z++)
	{
		for (size_t m = 0; m < x.Spalten(); m++)
		{
			s >> x(z, m);
		}
	}
	return s;
}

// ==========================
//      Fehlerbehandlung
// ==========================

// ----- Ausgabe der Fehlermeldung "str" -----

void Matrix::matError(const char str[]) {
	std::cerr << "\nMatrixfehler: " << str << '\n' << std::endl;
	exit(1);
}

