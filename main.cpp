/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include <fstream>
#include <iomanip>
#include <filesystem>   
 
using namespace std;

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1()
{
	cout << "------------ lab 1 ------------" << endl;
// 	random_device rd;
// 	mt19937 gen(rd());
// 	uniform_real_distribution<double> dist(-100.0, 100.0);
//
// 	double x0 = dist(gen);  //randomowy punkt pocz
// 	double d = 0.5;
// 	double ALPHAS[] = { 1.1, 1.5, 2.0 };
// 	int Nmax_exp = 1000;
// 	double eps = 0.001;
// 	matrix ud1 = NAN, ud2 = NAN;
// 	// --- Lagrange ---
// 	double gammaL = 1e-8;          // czułość na zmianę d(i)
// 	int    NmaxL = 1000;          // maks. wywołań f dla tej metody
//
// 	/*
// 	cout << "Losowy punkt startowy: " << x0 << endl;
//
// 	double* przedzial = expansion(ff1L, x0, d, ALPHAS[1], Nmax, ud1, ud2);
//
// 	cout << "Przedzial ekspansji: [" << przedzial[0] << ", " << przedzial[1] << "]" << endl;
//
// 	solution wynik = fib(ff1L, przedzial[0], przedzial[1], eps, ud1, ud2);
// 	cout << "Minimum funkcji Fibonacciego:\n";
// 	cout << "x* = " << m2d(wynik.x) << endl;
// 	cout << "f(x*) = " << m2d(wynik.y) << endl;
// 	cout << "f_calls = " << solution::f_calls << endl;
// 	cout << "Exit flag = " << wynik.flag << endl;
//
// 	solution wynikL = lag(ff1L, przedzial[0], przedzial[1], eps, gammaL, NmaxL, ud1, ud2);
// 	cout << "\nMinimum (Lagrange):\n";
// 	cout << "x* = " << m2d(wynikL.x) << endl;
// 	cout << "f(x*) = " << m2d(wynikL.y) << endl;
// 	cout << "f_calls = " << solution::f_calls << endl;
// 	cout << "Exit flag = " << wynikL.flag << endl;
// 	*/
//
// 	// CSV – nowy plik na każde uruchomienie
// 	ofstream log("wyniki_lab1.csv", ios::out | ios::trunc);
// 	if (!log.is_open()) {
// 		cerr << "Nie moge otworzyc pliku wyniki_lab1.csv" << endl;
// 		return;
// 	}
// 	log << fixed << setprecision(6);
// 	log << "alpha,lp,x0,a,b,"
// 		"f_calls_exp,"
// 		"fib_x,fib_y,f_calls_fib,"
// 		"lag_x,lag_y,f_calls_lag\n";
// 	log.flush();
//
// 	// === pętle: dla każdej alfy zrób 100 uruchomień ===
// 	for (int a_id = 0; a_id < 3; a_id++)
// 	{
// 		double alpha = ALPHAS[a_id];
// 		cout << "\n===== 100 biegow dla alpha = " << alpha << " =====" << endl;
//
// 		int lp = 0;
//
// 		for (int k = 0; k < 100; k++)
// 		{
// 			double x0 = dist(gen);
// 			cout << "Losowy punkt startowy: " << x0 << endl;
//
// 			// --- ekspansja (czyści licznik przed fazą) ---
// 			solution::clear_calls();
// 			double* przedzial = expansion(ff1L, x0, d, alpha, Nmax_exp, ud1, ud2);
// 			int f_calls_exp = solution::f_calls;
//
// 			cout << "Przedzial ekspansji: [" << przedzial[0] << ", " << przedzial[1] << "]" << endl;
//
// 			// --- Fibonacci ---
// 			solution wynik = fib(ff1L, przedzial[0], przedzial[1], eps, ud1, ud2);
// 			cout << "Minimum funkcji Fibonacciego:\n";
// 			cout << "x* = " << m2d(wynik.x) << endl;
// 			cout << "f(x*) = " << m2d(wynik.y) << endl;
// 			cout << "f_calls = " << solution::f_calls << endl;
// 			cout << "Exit flag = " << wynik.flag << endl;
// 			int f_calls_fib = solution::f_calls;
//
// 			// --- Lagrange ---
// 			solution wynikL = lag(ff1L, przedzial[0], przedzial[1], eps, gammaL, NmaxL, ud1, ud2);
// 			cout << "\nMinimum (Lagrange):\n";
// 			cout << "x* = " << m2d(wynikL.x) << endl;
// 			cout << "f(x*) = " << m2d(wynikL.y) << endl;
// 			cout << "f_calls = " << solution::f_calls << endl;
// 			cout << "Exit flag = " << wynikL.flag << endl;
// 			int f_calls_lag = solution::f_calls;
//
// 			// --- zapis jednego "bloku" do CSV (jeden wiersz) ---
// 			lp++;
// 			log << alpha << "," << lp << "," << x0 << ","
// 				<< przedzial[0] << "," << przedzial[1] << ","
// 				<< f_calls_exp << ","
// 				<< m2d(wynik.x) << "," << m2d(wynik.y) << "," << f_calls_fib << ","
// 				<< m2d(wynikL.x) << "," << m2d(wynikL.y) << "," << f_calls_lag
// 				<< "\n";
// 			log.flush();
//
// 			delete[] przedzial;
// 			cout << "---------------------------------------------\n";
// 		}
// 	}
//
// 	log.close();
// 	cout << "\nZapisano 3x100 wierszy do wyniki_lab1.csv" << endl;

	matrix ud1 = NAN, ud2 = NAN;
	double a = -100;
	double b = 100;
	double eps = 0.001;

	fib_no_exp_plot(ff1L, a, b, eps, ud1, ud2);
}



void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}



//ad lab2
/*
ofstream log("wyniki_lab1.csv", ios::out | ios::trunc);
if (!log.is_open()) {
	cerr << "Nie moge otworzyc pliku wyniki_lab1.csv" << endl;
	return;
}
log << fixed << setprecision(6);
log << "alpha,lp,x0,a,b,"
<< "f_calls_exp,"
<< "fib_x,fib_y,f_calls_fib,"
<< "lag_x,lag_y,f_calls_lag\n";
log.flush();

for (int a_id = 0; a_id < 3; a_id++)
{
	double alpha = ALPHAS[a_id];
	cout << "\n===== 100 biegow dla alpha = " << alpha << " =====" << endl;

	int lp = 0;

	for (int k = 0; k < 100; k++)
	{
		double x0 = dist(gen);

		// --- EKSPANSJA ---
		solution::clear_calls();
		double* przedzial = expansion(ff1L, x0, d, alpha, Nmax_exp, ud1, ud2);
		int f_calls_exp = solution::f_calls;

		double a = przedzial[0];
		double b = przedzial[1];

		cout << "[" << (k + 1) << "] x0=" << x0
			<< "  przedzial=[" << a << ", " << b << "]" << endl;

		// --- FIBONACCI ---
		solution wynikF = fib(ff1L, a, b, eps, ud1, ud2);
		int f_calls_fib = solution::f_calls;

		// --- LAGRANGE ---
		solution wynikL = lag(ff1L, a, b, eps, gammaL, NmaxL, ud1, ud2);
		int f_calls_lag = solution::f_calls;
		if (wynikL.flag == 0) {
			// drobne "rozruszanie" – przesuwamy c wewnętrznie przez mikrozwężenie [a,b]
			double wiggle = 1e-9 * (b - a);
			double a2 = a + wiggle;
			double b2 = b - wiggle;
			if (a2 < b2) {
				wynikL = lag(ff1L, a2, b2, eps, 10.0 * gammaL, NmaxL, ud1, ud2);
				f_calls_lag = solution::f_calls;
			}
		}

		// --- LOG ---
		lp++;
		log << alpha << "," << lp << "," << x0 << ","
			<< a << "," << b << ","
			<< f_calls_exp << ","
			<< m2d(wynikF.x) << "," << m2d(wynikF.y) << "," << f_calls_fib << ","
			<< m2d(wynikL.x) << "," << m2d(wynikL.y) << "," << f_calls_lag
			<< "\n";
		log.flush();

		delete[] przedzial;
	}
}

log.close();
cout << "\nZapisano 3x100 wierszy do wyniki_lab1.csv" << endl;*/