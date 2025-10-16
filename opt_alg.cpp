#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

//helper do zliczania l wykonan funkcji celu
static double f_eval(matrix(*ff)(matrix, matrix, matrix), double x, matrix ud1, matrix ud2)
{
	solution tmp(matrix(1, 1, x));
	tmp.fit_fun(ff, ud1, ud2);     // tu rośnie solution::f_calls
	return m2d(tmp.y);
}



double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 }; //przedzial
		//Tu wpisz kod funkcji
		int i = 0;

		// obliczenie f(x0) i f(x1)
		double fx0 = f_eval(ff, x0, ud1, ud2);
		double x1 = x0 + d;
		double fx1 = f_eval(ff, x1, ud1, ud2);

		//sprawdzenie kieruknu
		if (fx1 == fx0)
		{
			p[0] = min(x0, x1); // tu mamy przyklad [x(-1); x(1)]
			p[1] = max(x0, x1);
			return p;
		}

		// jesli idziemy w prawo i rosnie, odwracamy
		if (fx1 > fx0)
		{
			d = -d;
			x1 = x0 + d;
			fx1 = f_eval(ff, x1, ud1, ud2);

			// jesli nadal rosnie zwracamy przedzial miedzy x0
			if (fx1 >= fx0)
			{
				p[0] = min(x1, x0 - d);
				p[1] = max(x1, x0 - d);
				return p;
			}
		}


		double xi_prev = x0;
		double xi = x1;
		double fx_prev = fx0;
		double fx = fx1;

		do
		{
			if (solution::f_calls > Nmax)
			{
				throw string("expansion(): przekroczono maksymalną liczbę wywołań funkcji celu.");
			}

			i++;
			double xi_next = x0 + pow(alpha, i) * d;
			double fx_next = fx_next = f_eval(ff, xi_next, ud1, ud2);

			// jak przestaje malec
			if (fx_next >= fx)
			{
				if (d > 0)
				{
					p[0] = xi_prev;
					p[1] = xi_next;
				}
				else
				{
					p[0] = xi_next;
					p[1] = xi_prev;
				}
				return p;
			}

			// przesuniecie okna
			xi_prev = xi;
			fx_prev = fx;
			xi = xi_next;
			fx = fx_next;

		} while (true);

	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt; // finalnie ma byc Xopt(x*, y*, flaga itd.)
        solution::clear_calls(); // wyczyszczenie licznika

        // --------- tu tworzymy ciag Fibonacciego - vector, zeby latwo na indeksy sie odnosic
        vector<double> F;
        F.push_back(0);
        F.push_back(1);
        while (F.back() < (b - a) / epsilon)
        {
            int n = F.size();
            F.push_back(F[n - 1] + F[n - 2]);
        }

        int n = F.size() - 1;  // ostatni indeks

        // --------- punkty startowe C i D -----------
        solution C, D;
        C.x = matrix(1, 1, b - (F[n - 2] / F[n]) * (b - a));
        D.x = matrix(1, 1, a + b - m2d(C.x));   // d(0) = a + b - c(0)

        C.fit_fun(ff, ud1, ud2);
        D.fit_fun(ff, ud1, ud2);

        // --------- petla Fibonacciego -----------
        for (int i = 0; i <= n - 3; i++)
        {
            if (m2d(C.y) < m2d(D.y))
            {
                // minimum po lewej stronie
                b = m2d(D.x);
                D = C;  // przesuwamy d w miejsce c

                C.x = matrix(1, 1, b - (F[n - i - 3] / F[n - i - 2]) * (b - a));
                C.fit_fun(ff, ud1, ud2); // nowe obliczenie f(c)
            }
            else
            {
                // minimum po prawej stronie
                a = m2d(C.x);
                C = D;  // przesuwamy c w miejsce d

                D.x = matrix(1, 1, a + (F[n - i - 3] / F[n - i - 2]) * (b - a));
                D.fit_fun(ff, ud1, ud2); // nowe obliczenie f(d)
            }
        }

        // --------- wynik koncowy -----------
        Xopt.x = matrix(1, 1, (m2d(C.x) + m2d(D.x)) / 2.0); // wspolrzedna x*
        Xopt.fit_fun(ff, ud1, ud2); // obliczenie y*
        Xopt.flag = 1; // ustawienie flagi sukcesu
        Xopt.ud = matrix(1, 1, solution::f_calls); // user data = liczba wywolan funkcji celu

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution fib(...):\n" + ex_info);
    }
}


void fib_no_exp_plot(matrix(*ff)(matrix, matrix, matrix), double a, double b, double eps, matrix ud1, matrix ud2)
{
    try
    {
        vector<double> interval_lengths;

        // ----------- tworzenie ciagu Fib -----------
        vector<double> F;
        F.push_back(0);
        F.push_back(1);
        while (F.back() < (b - a) / eps)
        {
            int n = F.size();
            F.push_back(F[n - 1] + F[n - 2]);
        }

        int n = F.size() - 1; // ostatni indeks

        // ----------- punkty poczatkowe c(0), d(0) -----------
        double a_i = a;
        double b_i = b;
        double c = b_i - (F[n - 2] / F[n]) * (b_i - a_i);
        double d = a_i + b_i - c;

        double f_c = m2d(ff(matrix(1, 1, c), ud1, ud2));
        double f_d = m2d(ff(matrix(1, 1, d), ud1, ud2));

        // zapisanie dlugosci przedzialu (na starcie = 200)
        interval_lengths.push_back(b_i - a_i);

        for (int k = n - 1; k >= 2; k--)
        {
            if (f_c < f_d)
            {
                // minimum jest w lewej części
                b_i = d;
                d = c;
                f_d = f_c;
                c = b_i - (F[k - 3] / F[k - 2]) * (b_i - a_i);
                f_c = m2d(ff(matrix(1, 1, c), ud1, ud2));
            }
            else
            {
                // minimum jest w prawej części
                a_i = c;
                c = d;
                f_c = f_d;
                d = a_i + (F[k - 3] / F[k - 2]) * (b_i - a_i);
                f_d = m2d(ff(matrix(1, 1, d), ud1, ud2));
            }

            // zapisujemy dlugosc przedzialu po tej iteracji
            interval_lengths.push_back(b_i - a_i);
        }

        // ----------- zapis danych do pliku CSV -----------
        ofstream file("fib_interval_lengths.csv");
        file << "Iteracja,Dlugosc_przedzialu\n";
        for (int i = 0; i < interval_lengths.size(); i++)
        {
            file << i << "," << interval_lengths[i] << "\n";
        }
        file.close();

        cout << "Dane zapisano do pliku fib_interval_lengths.csv" << endl;
    }
    catch (string ex_info)
    {
        throw("fib_no_exp_plot(...): " + ex_info);
    }
}

solution lag(matrix(*ff)(matrix, matrix, matrix),
	double a, double b, double epsilon, double gamma, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution::clear_calls();

		double ai = a, bi = b, ci = 0.5 * (a + b);

		solution SA(matrix(1, 1, ai)), SB(matrix(1, 1, bi)), SC(matrix(1, 1, ci));
		SA.fit_fun(ff, ud1, ud2);
		SB.fit_fun(ff, ud1, ud2);
		SC.fit_fun(ff, ud1, ud2);

		int i = 0;
		double d_prev = std::numeric_limits<double>::quiet_NaN();
		double di = std::numeric_limits<double>::quiet_NaN();

		while (true)
		{
			const double a2 = ai * ai, b2 = bi * bi, c2 = ci * ci;
			const double fA = m2d(SA.y), fB = m2d(SB.y), fC = m2d(SC.y);
			const double l = fA * (b2 - c2) + fB * (c2 - a2) + fC * (a2 - b2);
			const double m = fA * (bi - ci) + fB * (ci - ai) + fC * (ai - bi);

			// --- osłona numeryczna dla m≈0
			double tau = 1e-15 * (std::fabs(fA) + std::fabs(fB) + std::fabs(fC) + 1.0)
				* (std::fabs(ai) + std::fabs(bi) + std::fabs(ci) + 1.0);

			d_prev = di;
			if (!std::isfinite(m) || std::fabs(m) <= tau) {
				// zamiast "return error" z linii 6–8 -> wykonaj bezpieczny krok (połówkowanie)
				di = 0.5 * (ai + bi);
			}
			else {
				di = 0.5 * (l / m);
			}

			// --- jeśli d wyszło poza (a,b), też zrób bezpieczny krok zamiast error (linie 31–33)
			if (!(ai < di && di < bi)) {
				di = 0.5 * (ai + bi);
			}

			// wybór podprzedziału (linie 10–29) – bez zmian
			if (ai < di && di < ci) {
				solution SD(matrix(1, 1, di)); SD.fit_fun(ff, ud1, ud2);
				if (m2d(SD.y) < fC) { bi = ci; SB = SC; ci = di; SC = SD; }
				else { ai = di; SA = SD; }
			}
			else if (ci < di && di < bi) {
				solution SD(matrix(1, 1, di)); SD.fit_fun(ff, ud1, ud2);
				if (m2d(SD.y) < fC) { ai = ci; SA = SC; ci = di; SC = SD; }
				else { bi = di; SB = SD; }
			}
			else {
				// rzadki przypadek równości wskutek zaokrągleń -> minimalny skurcz
				di = 0.5 * (ai + bi);
				solution SD(matrix(1, 1, di)); SD.fit_fun(ff, ud1, ud2);
				if (di <= ci) { ai = di; SA = SD; }
				else { bi = di; SB = SD; }
			}

			++i;
			if (solution::f_calls > Nmax) { Xopt.x = NAN; Xopt.y = NAN; Xopt.flag = 0; return Xopt; }

			const bool stop_len = (bi - ai) < epsilon;
			const bool stop_step = (i > 1 && std::fabs(di - d_prev) < gamma);

			if (stop_len || stop_step) {
				Xopt.x = matrix(1, 1, di);
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				Xopt.ud = matrix(1, 1, solution::f_calls);
				return Xopt;
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


/*solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution::clear_calls();              // spójnie z fib()

		// startowy punkt wewnętrzny (brak 'c' w sygnaturze → bierzemy środek)
		double ai = a, bi = b, ci = 0.5 * (a + b);

		// a(0), b(0), c(0) – utwórz "solution", żeby fit_fun() zliczał f_calls
		solution SA(matrix(1, 1, ai)), SB(matrix(1, 1, bi)), SC(matrix(1, 1, ci));
		SA.fit_fun(ff, ud1, ud2);
		SB.fit_fun(ff, ud1, ud2);
		SC.fit_fun(ff, ud1, ud2);

		int i = 0;
		double d_prev = std::numeric_limits<double>::quiet_NaN();
		double di = std::numeric_limits<double>::quiet_NaN();

		while (true)
		{
			// l, m zgodnie z pseudokodem
			const double a2 = ai * ai, b2 = bi * bi, c2 = ci * ci;
			const double fA = m2d(SA.y), fB = m2d(SB.y), fC = m2d(SC.y);
			const double l = fA * (b2 - c2) + fB * (c2 - a2) + fC * (a2 - b2);
			const double m = fA * (bi - ci) + fB * (ci - ai) + fC * (ai - bi);

			// --- osłona numeryczna dla m≈0
			double tau = 1e-15 * (std::fabs(fA) + std::fabs(fB) + std::fabs(fC) + 1.0)
				* (std::fabs(ai) + std::fabs(bi) + std::fabs(ci) + 1.0);

			if (!std::isfinite(m) || m <= 0.0) {         // m ≤ 0 → error
				Xopt.x = NAN; Xopt.y = NAN; Xopt.flag = 0;
				return Xopt;
			}
			d_prev = di;
			di = 0.5 * (l / m);
			
d_prev = di;
if (!std::isfinite(m) || std::fabs(m) <= tau) {
	// zamiast "return error" z linii 6–8 -> wykonaj bezpieczny krok (połówkowanie)
	di = 0.5 * (ai + bi);
}
else {
	di = 0.5 * (l / m);
}

// --- jeśli d wyszło poza (a,b), też zrób bezpieczny krok zamiast error (linie 31–33)
if (!(ai < di && di < bi)) {
	di = 0.5 * (ai + bi);
}

// wybór podprzedziału i aktualizacja (linie 10–34 pseudokodu)
if (ai < di && di < ci) {
	solution SD(matrix(1, 1, di)); SD.fit_fun(ff, ud1, ud2);
	if (m2d(SD.y) < fC) {
		// a(i+1)=a(i); c(i+1)=d(i); b(i+1)=c(i)
		bi = ci;  SB = SC;
		ci = di;  SC = SD;
		// ai, SA bez zmian
	}
	else {
		// a(i+1)=d(i); c(i+1)=c(i); b(i+1)=b(i)
		ai = di;  SA = SD;
		// ci, SC; bi, SB bez zmian
	}
}
else if (ci < di && di < bi) {
	solution SD(matrix(1, 1, di)); SD.fit_fun(ff, ud1, ud2);
	if (m2d(SD.y) < fC) {
		// a(i+1)=c(i); c(i+1)=d(i); b(i+1)=b(i)
		ai = ci;  SA = SC;
		ci = di;  SC = SD;
	}
	else {
		// a(i+1)=a(i); c(i+1)=c(i); b(i+1)=d(i)
		bi = di;  SB = SD;
	}
}
else {
	// d(i) poza (a,c)∪(c,b) → error
	Xopt.x = NAN; Xopt.y = NAN; Xopt.flag = 0;
	return Xopt;
}

++i;
if (solution::f_calls > Nmax) {               // limit wywołań
	Xopt.x = NAN; Xopt.y = NAN; Xopt.flag = 0;
	return Xopt;
}

// kryteria stopu: b(i)-a(i) < eps LUB |d(i)-d(i-1)| < gamma
const bool stop_len = (bi - ai) < epsilon;
const bool stop_step = (i > 1 && std::fabs(di - d_prev) < gamma);
if (stop_len || stop_step) {
	Xopt.x = matrix(1, 1, di);
	Xopt.fit_fun(ff, ud1, ud2);               // f(x*)
	Xopt.flag = 1;
	Xopt.ud = matrix(1, 1, solution::f_calls); // w ud: liczba f-calls (jak we fib)
	return Xopt;
}
		}
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}*/

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
