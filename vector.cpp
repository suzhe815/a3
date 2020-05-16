/************************************************
 *  Name       : vector.cpp                     *
 *  Verwendung : Vektorklasse                   *
 *  Autor      :                                *
 *  Datum      :                                *
 ***********************************************/

#include "vector.h"

#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <new>
#define epsilon 1e-9
#define NDEBUG

using std::nothrow;
using std::cout;

// =======================
//      Konstruktoren
// =======================

// ----- Konstruktor -----

Vector::Vector(size_t len) {
#ifndef NDEBUG
    if (len <= 0) {
        vecError("Nur Vektoren mit positiver Laenge!");
    }
#endif

    length = len;

    elems = new (std::nothrow) double[length];
    if (elems == nullptr) {
        vecError("Nicht genuegend Speicher!");
    }

    for (size_t i = 0; i < length; i++) {
        (*this)(i) = 0;
    }
}

// ----- Kopierkonstruktor -----

Vector::Vector(const Vector &x) {
    length = x.length;

    elems = new (std::nothrow) double[length];
    if (elems == nullptr) {
        vecError("Nicht genuegend Speicher!");
    }

    for (size_t i = 0; i < length; i++) {
        (*this)(i) = x(i);
    }
}

// ===========================================
//      Vektorelement schreiben und lesen
// ===========================================

// ----- Schreib- und Lesezugriff auf Vektorelemente -----

double &Vector::operator()(size_t i) {
#ifndef NDEBUG
    if (i >= length) vecError("Ungueltiger Index!");
#endif

    return elems[i];
}

// ----- Lesezugriff auf Elemente konstanter Vektoren -----

double Vector::operator()(size_t i) const {
#ifndef NDEBUG
    if (i >= length) vecError("Ungueltiger Index!");
#endif

    return elems[i];
}

// =====================
//      Zuweisungen
// =====================

// ----- Zuweisungsoperator "=" ----

Vector &Vector::operator=(const Vector &x) {
#ifndef NDEBUG
    if (length != x.length) {
        vecError("Inkompatible Dimensionen fuer 'Vektor = Vektor'!");
    }
#endif

    for (size_t i = 0; i < length; i++) {
        (*this)(i) = x(i);
    }

    return *this;
}

// ----- Zuweisungsoperator mit Addition "+=" ----

Vector &Vector::operator+=(const Vector &x) {
    
#ifndef NDEBUG
	if (length != x.length)
		vecError("Inkompatible Dimensionen fuer 'Vektor += Vektor'!");
#endif
	for (size_t i = 0; i < length; i++)
	{
		(*this)(i) += x(i);
	}
	return *this;
}

// ----- Zuweisungsoperator mit Subtraktion "-=" ----

Vector &Vector::operator-=(const Vector &x) {
#ifndef NDEBUG
	if (length != x.length)
		vecError("Inkompatible Dimensionen fuer 'Vektor -= Vektor'!");
#endif
	for (size_t i = 0; i < length; i++)
	{
		(*this)(i) -= x(i);
	}
	return *this;
}

// ----- Zuweisungsoperator mit Multiplikation "*=" ----

Vector& Vector::operator*=(double c) {
	for (size_t i = 0; i < length;i++)
	{
			(*this)(i) *= c;
    }

	return *this;
}

// ----- Zuweisungsoperator mit Divsion "/=" ----

Vector &Vector::operator/=(double c) {
#ifndef NDEBUG
	if (abs(c) < epsilon)
	{
		vecError("Kann nicht durch 0 dividieren!");
	}
#endif // !NDEBUG

	for (size_t i = 0; i < length;i++)
	{
		(*this)(i) /= c;
	}

	return *this;
}

// ==============================
//      Vektorlaenge aendern
// ==============================

// ----- Vektorlaenge aendern -----

Vector &Vector::redim(size_t l) {
	delete[] elems;
	elems = new(nothrow) double[l];
	length = l;
#ifndef NDEBUG
	if (elems == 0)
	{
		vecError("nicht genuegende plaetze!");
	}
#endif

	for (size_t i = 0; i < length; i++)
	{
		elems[i] = 0;
	}
	return *this;
}

// ======================
//      Vektornormen
// ======================

// ----- Euklidische Norm -----

double Vector::norm2() const {
#ifndef NDEBUG
	if (elems == 0)
	{
		vecError("vektor existiert nicht!");
	}
#endif

	double sum = 0;
	for (size_t i = 0; i < length; i++)
	{
		sum += elems[i] * elems[i];
	}
	return sqrt(sum);
}

// ----- Maximum-Norm -----

double Vector::normMax() const {
#ifndef NDEBUG
	if (elems == 0)
	{
		vecError("vektor existiert nicht!");
	}
#endif

	double max = 0;
	for (size_t i = 0; i < length; i++)
	{
		double tmp = abs(elems[i]);
		if (max < tmp)
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
//Unterschied zu =+ : hier wird der vorgegebene vector nicht geaendert, sondern
// nur eine kopie zurueckgegeben.

Vector operator+(const Vector &x, const Vector &y) {
#ifndef NDEBUG
    if (x.length != y.length) {
        Vector::vecError("Inkompatible Dimensionen fuer 'Vektor + Vektor'!");
    }
#endif

    Vector z = x;
    return z += y;
}

// ----- Subtraktion "-" -----

Vector operator-(const Vector &x, const Vector &y) {
#ifndef NDEBUG
	if (x.length != y.length) {
		Vector::vecError("Inkompatible Dimensionen fuer 'Vektor + Vektor'!");
	}
#endif

	Vector z = x;
	return z -= y;
}

// ----- Vorzeichen wechseln "-" -----

Vector operator-(const Vector &x) {
#ifndef NDEBUG
	if (x.length == 0) {
		Vector::vecError("vektor existiert nicht!");
	}
#endif
	Vector z = x;
	z *= -1;
	return z;
}

// ----- Skalarprodukt "*" -----

double operator*(const Vector &x, const Vector &y) {
#ifndef NDEBUG
	if (x.length != y.length) {
		Vector::vecError("Inkompatible Dimensionen fuer 'Vektor + Vektor'!");
	}
#endif
	double sum = 0;
	for (size_t i = 0; i < x.length; i++)
	{
		sum += x.elems[i] * y.elems[i];
	}

	return sum;
}

// ----- Multiplikation Skalar*Vektor "*" -----

Vector operator*(double c, const Vector &x) {
#ifndef NDEBUG
	if (x.length == 0) {
		Vector::vecError("vektor existiert nicht!");
	}
#endif
	Vector z = x;
	for (size_t i = 0; i < z.length; i++)
	{
		z.elems[i] = c * z.elems[i];
	}

	return z;
}

// ----- Multiplikation Vektor*Skalar "*" -----

Vector operator*(const Vector &x, double c) {
#ifndef NDEBUG
	if (x.length == 0) {
		Vector::vecError("vektor existiert nicht!");
	}
#endif
	Vector z = x;
	for (size_t i = 0; i < z.length; i++)
	{
		z.elems[i] = c * z.elems[i];
	}

	return z;
}

// ----- Division Vektor/Skalar "/" -----

Vector operator/(const Vector &x, double c) {
#ifndef NDEBUG
	if (x.length < epsilon) {
		Vector::vecError("vektor existiert nicht!");
	}
#endif

#ifndef NDEBUG
	if (abs(c) < epsilon)
	{
		Vector::vecError("Kann nicht durch 0 dividieren!");
	}
#endif 
	Vector z = x;
	for (size_t i = 0; i < z.length;i++)
	{
		z.elems[i] /= c;
	}

	return z;
}

// ==============================
//      Vergleichsoperatoren
// ==============================

// ----- Test auf Gleichheit "==" -----

bool operator==(const Vector &x, const Vector &y) {
    if (x.length != y.length) {
        return false;
    }

    for (size_t i = 0; i < x.length; i++) {
        if (x(i) != y(i)) {
            return false;
        }
    }

    return true;
}

// ----- Test auf Ungleichheit "!=" -----

bool operator!=(const Vector &x, const Vector &y) {

	return !(x == y);
}

// ==========================
//      Ein- und Ausgabe
// ==========================

// ----- Ausgabe "<<" -----

std::ostream &operator<<(std::ostream &s, const Vector &x) {
    s << std::setiosflags(std::ios::right);
    for (size_t i = 0; i < x.length; i++) {
        s << "\n(" << std::setw(4) << i << ") " << x(i);
    }

    return s << std::endl;
}

// ----- Eingabe ">>" -----

std::istream &operator>>(std::istream &s, Vector &x) {
    std::cout << std::setiosflags(std::ios::right);
    for (size_t i = 0; i < x.length; i++) {
        std::cout << "\n(" << std::setw(4) << i << ") ";
        s >> x(i);
    }
    return s;
}

// ==========================
//      Fehlerbehandlung
// ==========================

// ----- Ausgabe der Fehlermeldung "str" -----

void Vector::vecError(const char str[]) {
    std::cerr << "\nVektorfehler: " << str << '\n' << std::endl;
    exit(1);
}
