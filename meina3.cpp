#include <iostream>
#include <ostream>
#include <sstream>
#include "vector.h"
#include "matrix.h"
#include "unit.h"
#include <fstream>

//Jakobiverfahren
Vector jakobi(Matrix& A, Vector& x0, Vector& b, double& tol, int& maxiter,  int& iterationenJ, Vector& rResiduumJ)
{
	std::cout << "Jakobi:";
	// speichere die benoetigten It.schritte auf iterationenJ
	// k = Anzahl der Iterationen
	int k = 0;

	size_t n = x0.getLength();

	// speichere die Zwischenergebnisse auf xk und am Ende gib xk als die annaehrende Loesung zurueck.
	Vector xk = x0;

	// //Berechnung des Residuums
	Vector residuum = b - A * xk;

	// relative Fehler von norm(b - A * xk)/ norm(b)
	// speichere den rel.Fehler jedes Schrittes jeweils auf der 0,..k,... Stelle von rResiduumJ
	rResiduumJ(k) = residuum.norm2() / b.norm2();

	// Jakobi-Verfahren im Buch
	// solange das Residuum groesser der uebergebenen Toleranz ist 
	// und die maximale Anzahl von Iterationen nicht erreicht ist, wird 
	// mit gegebenen Algorithmus weiter iteriert
	while (k < maxiter && rResiduumJ(k) >= tol)
	{
		//Da x_k akutualisiert werden wird und man die x_k in letztem Schritt verwenden muss,
		// definieren man eine Kopie xk_tmp in jeder Iteration zur Verfuegung zu stellen
		Vector xk_tmp = xk;
		for (size_t i = 0; i < n; i++)
		{
			double sum = 0;
			for (size_t j = 0; j < n;j++)
			{
				if (j!= i)
					sum = sum + A(i, j) * xk_tmp(j);
			}
			xk(i) = (1 / A(i, i)) * (b(i) - sum);
		}

		//bis jetzt ist x_k+1 berechnet
		k = k + 1;
		residuum = b - A * xk;

		rResiduumJ(k) = residuum.norm2() / b.norm2();
	}

	iterationenJ = k;

	return xk;
}


// GS - Methode
Vector GaussSeidel(Matrix& A, Vector& x0, Vector& b, double& tol, int& maxiter, int& iterationenGS, Vector& rResiduumGS)
{
	std::cout << "Gauss - Seidel:";
	// speichere die benoetigten It.schritte auf iterationenJ
	// k = Anzahl der Iterationen
	int k = 0;

	size_t n = x0.getLength();

	// speichere die Zwischenergebnisse auf xk und am Ende gib xk als die annaehrende Loesung zurueck.
	Vector xk = x0;
	Vector tmp_xk = xk;
	// //Berechnung des Residuums
	Vector residuum = b - A * xk;

	// relative Fehler von norm(b - A * xk)/ norm(b)
	// speichere den rel.Fehler jedes Schrittes jeweils auf der 0,..k,... Stelle von rResiduumJ
	rResiduumGS(k) = residuum.norm2() / b.norm2();

	double sum1 = 0;
	double sum2 = 0;

	while (rResiduumGS(k) >= tol && k < maxiter)
	{
		
		// in sum2 muss man das x_k in letztem schritt verwenden 
		
		for (size_t i = 0; i < n; i++)
		{
			// in jedem schritt muessen sum1 und sum2  wieder als 0 eingesetzt werden
			sum1 = 0;
			sum2 = 0;
			// Unterscheide den Fall x_k+1(1) und den Rest von x_k+1
			// Beachte dazu, dass im Programm A(0,0) = "A(1,1) schreibweise im Buch"
			// Berechne zunaechst die zweite Summation in GS Verfahren:
			
			// die erste Summantion tritt nur fuer die 2-te bis n-te Komponente von x_k+1 auf
			if (i > 0)
			{
				
				for (size_t j = 0; j <= i - 1;j++)
				{
					sum1 = sum1 + A(i, j) * xk(j);
				}
			}

			for (size_t j = i+1 ; j < n; j++)
				{
					sum2 = sum2 + A(i, j) * tmp_xk(j);
				}
			xk(i) = (1 / A(i, i)) * (b(i) - sum1 - sum2);
			
		}
		tmp_xk = xk;
		k = k + 1;
		residuum = b - A * xk;
		rResiduumGS(k) = residuum.norm2() / b.norm2();




	}
	iterationenGS = k;

	return xk;
}




//CG - Methode
Vector CG(Matrix& A, Vector& x0, Vector& b, double& tol, int& maxiter, int& iterationenCG, Vector& rResiduumCG)
{
	std::cout << "CG Verfahren:";
	int k = 0;

	// speichere die Zwischenergebnisse auf xk und am Ende gib xk als die annaehrende Loesung zurueck.
	Vector xk = x0;;

	// //Berechnung des Residuums
	Vector residuum = b - A * xk;

	// relative Fehler von norm(b - A * xk)/ norm(b)
	rResiduumCG(k) = residuum.norm2() / b.norm2();

	// speichere die Zwischenergebnisse.
	Vector pk = residuum;
	Vector qk = A * pk;
	Vector rk = residuum;
	double gammak = residuum * residuum;
	double alphak = gammak / (qk * pk);
	double tmp_gammak = gammak;
	//CG Verfahren wie im Ue3
	while (rResiduumCG(k) >= tol && k < maxiter)
	{
		xk = xk + alphak * pk;
		rk = rk - alphak * qk;
		
		gammak = rk * rk;
		pk = rk + (gammak / tmp_gammak) * pk;
		tmp_gammak = gammak;
		residuum = b - A * xk;
		k = k + 1;
		rResiduumCG(k) = residuum.norm2() / b.norm2();
		// jetzt muessen noch qk und alphak akutalisiert werden
		qk = A * pk;
		alphak = gammak / (qk * pk);
		
	}

	iterationenCG = k;

	return xk;
}



// der Benutzer gibt mit der ersten Zahl an welches Beispiel er benutzen moechte und mit der zweiten Zahl welches Iterations-Verfahren
// 0:Jacobi, 1: GaussSeidel, 2: CG
// auf argc sind die anzahl der parametern(Programmname als das erste argument gezaehlt)
int main(int argc, char* argv[])
{
	std::istringstream IStr(argv[1]);
	int i;
	IStr >> i;

	std::istringstream IStr1(argv[2]);
	int j;
	IStr1 >> j;

	if (i > num_examples || i < 0 || j>2 || j < 0)
	{
		std::cout << "es gibt nur " << num_examples << "beispiele und drei verfahren: 0: Jakobi, 1: GS, 2:CG\n" ;
	}
	
	else
	{
		Matrix A;
		Vector b;
		Vector x0;
		double tol = 0;
		int maxiter = 0;

		// Iterationsvariablen werden an Algorithmen Ã¼bergeben und in diesen verÃ¤ndert, sodass nach Durchlauf der Algorithmen die Anzahl der Iterationen in diesen gespeichert ist
		int iterationenJ = 0, iterationenGS = 0, iterationenCG = 0;
		

		//Jakobi
		if (j == 0)
		{
			getExample(i, A, x0, b, tol, maxiter);
			Vector rResiduumJ(maxiter + 1);

			// speichere das endliche Ergebnis von x_k auf j
			Vector j = jakobi(A, x0, b, tol, maxiter, iterationenJ, rResiduumJ);

			std::ofstream Aufgabe3J("jakobi.txt", std::ofstream::out);

			if (Aufgabe3J.is_open())
			{
				for (int i = 0; i <= iterationenJ; i++)
				{
					Aufgabe3J << rResiduumJ(i) << "\n";
				}
			}
			Aufgabe3J.close();
			checkSolution(j, iterationenJ, 0);


		}

		else if (j == 1)
		{
			getExample(i, A, x0, b, tol, maxiter);
			Vector rResiduumGS(maxiter + 1);

			// speichere das endliche Ergebnis von x_k auf j
			Vector gs = GaussSeidel(A, x0, b, tol, maxiter, iterationenGS, rResiduumGS);

			std::ofstream Aufgabe3GS("GS.txt", std::ofstream::out);

			if (Aufgabe3GS.is_open())
			{
				for (int i = 0; i <= iterationenGS; i++)
				{
					Aufgabe3GS << rResiduumGS(i) << "\n";
				}
			}
			Aufgabe3GS.close();
			checkSolution(gs, iterationenGS, 1);


		}

		else if (j == 2)
		{
			getExample(i, A, x0, b, tol, maxiter);
			Vector rResiduumCG(maxiter + 1);

			// speichere das endliche Ergebnis von x_k auf j
			Vector cg = CG(A, x0, b, tol, maxiter, iterationenCG, rResiduumCG);

			std::ofstream Aufgabe3CG("CG.txt", std::ofstream::out);

			if (Aufgabe3CG.is_open())
			{
				for (int i = 0; i <= iterationenCG; i++)
				{
					Aufgabe3CG << rResiduumCG(i) << "\n";
				}
			}
			Aufgabe3CG.close();
			checkSolution(cg, iterationenCG, 2);


		}
		return 0;
	}
}
