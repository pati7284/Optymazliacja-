#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}


matrix ff1L(matrix x, matrix ud1, matrix ud2)
{
	//----------- lab 1 funkcja celu
	//f(x) = -cos(0.1x) * exp(-(0.1x - 2pi)^2) + 0.002 * (0.1x)^2
	const double PI = 3.14159265358979323846;
	double t = 0.1 * x(0);
	matrix y;
	y = -cos(t) * exp(-pow(t - 2.0 * PI, 2)) + 0.002 * (t * t);
	return y;
}


matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {

	// TRZYMAMY SI

	// tu dla zbiornika A
	double P_A = 2.0 ; // m^3
	double V_A_in = 5.0; // to do main
	double T_A = 95.0;


	// tu dla zbiornika B
	double P_B = 1.0;
	double V_B_in = 1.0; // to do main
	double T_B_in = 20.0;

	double F_B_in =0.010; // m^3/s
	double D_B = 0.00365665; //m^2


	double a = 0.98 ; //wspolczynnik odpowiadajacy za lepkosc cieczy
	double b = 0.64; // wspolaczynik odpowiadajacy za zwezenie strumieni cieczy
	double g = 9.81; //przysieszenie



	double F_A_out , F_B_out; // gdzie dV_A/dt = - F_A_out
	double D_A = m2d(ud1); // to ma byc przekazane


	if(Y(0)> 0.0) {

		F_A_out = a*b*D_A*sqrt(2.0*g*(Y(0)/ P_A));
	}
	else {
		F_A_out = 0.0;
	}

	if(Y(1)> 0.0) {
		F_B_out = a*D_B*sqrt(2.0*g*(Y(1) / P_B));
	}
	else {
		F_B_out = 0.0;
	}

	// rownania rozniczkowe
	matrix dY = matrix(3,1);


	dY(0) = -F_A_out;
	dY(1) = F_A_out - F_B_out + F_B_in;
	dY(2) = (F_B_in / Y(1)) * (T_B_in - Y(2)) + (F_A_out / Y(1)) * (T_A - Y(2));

	return dY;




}