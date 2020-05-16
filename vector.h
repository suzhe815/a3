/************************************************
 *  Name        : vector.h                      *
 *  Verwendung  : Header zu Vektorklasse        *
 *  Autor       : V. Reichelt, IGPM RWTH Aachen *
 *  Datum       : Nov '96 - Mai '00             *
 *  Aktualisiert: K. Brix, Apr '07              *
 *                F. Gruber, Apr '16            *
 ***********************************************/

#ifndef _VECTOR_H  // vector.h nicht doppelt benutzen
#define _VECTOR_H

#include <iostream>

class Matrix;  // fuer friend Matrix * Vektor etc.

class Vector {
   private:
    double* elems;  // Zeiger auf Feld fuer Vektorelemente
    size_t length;  // Vektorlaenge

   public:
    explicit Vector(size_t len = 1);  // Konstruktor mit Laengenangabe
    ~Vector() {
        if (elems) delete[] elems;
    }                       // Destruktor
    Vector(const Vector&);  // Kopierkonstruktor

    double& operator()(size_t);       // Zugriff auf Einträge des Vektors
    double operator()(size_t) const;  // Zugriff auf falls Vektor const

    Vector& operator=(const Vector&);   // Zuweisung
    Vector& operator+=(const Vector&);  // Zuweisungen mit arithm. Operation
    Vector& operator-=(const Vector&);
    Vector& operator*=(double);
    Vector& operator/=(double);

    Vector& redim(size_t);  // neue Laenge festlegen
    size_t getLength() const {
        return length;
    }                        // Laenge
    double norm2() const;    // Euklidische Norm des Vektors
    double normMax() const;  // Maximum-Norm des Vektors

    static void vecError(const char str[]);  // Fehlermeldung ausgeben

    friend Vector operator+(const Vector&, const Vector&);  // Addition
    friend Vector operator-(const Vector&, const Vector&);  // Subtraktion
    friend Vector operator-(const Vector&);                 // Vorzeichen

    friend double operator*(const Vector&, const Vector&);  // Skalarprodukt
    friend Vector operator*(double, const Vector&);         // Vielfache
    friend Vector operator*(const Vector&, double);
    friend Vector operator/(const Vector&, double);

    friend bool operator==(const Vector&, const Vector&);  // Vergleich
    friend bool operator!=(const Vector&, const Vector&);

    friend std::istream& operator>>(std::istream&, Vector&);        // Eingabe
    friend std::ostream& operator<<(std::ostream&, const Vector&);  // Ausgabe

    friend Vector operator*(const Matrix&, const Vector&);  // Matrix-Vektor-
    friend Vector operator*(const Vector&, const Matrix&);  // Multiplikation
};

#endif
