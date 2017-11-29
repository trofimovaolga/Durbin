#define _USE_MATH_DEFINES
#include <cstdio> 
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <fstream>
#include <random>
#include <complex>

using namespace std;

complex <double> im(0, 1);
complex <double> M, int_h, m, ex1, ex2, exz1, exz2;
double c = 2, myu_s = 0.0, myu_a = 0.06, myu = myu_a + myu_s;
double delta = 1.0 / 5000.0; //R = N * delta 
int N = 1000000, n = 20;
double z = 0.5;
double a = 0.01, b = 5.0, h = (b - a) / double(n), l1 = 0.8;

void halt(int s) {
	system("pause");
	exit(s);
}

void slau(complex <double> eq[2][3], complex <double>& c1, complex <double>& c2) {
	complex <double> delta = eq[0][0] * eq[1][1] - eq[1][0] * eq[0][1];
	complex <double> delta_c1 = eq[1][1] * eq[0][2] - eq[0][1] * eq[1][2];
	complex <double> delta_c2 = eq[0][0] * eq[1][2] - eq[1][0] * eq[0][2];

	if (delta == complex<double>(0, 0)) {
		cout << "delta 0" << endl << eq[0][0] << ' ' << eq[1][1] << endl << eq[1][0] * eq[0][1] << endl;
		halt(0);
	}
	if (delta != delta) {
		cout << "delta nan" << endl << eq[0][0] << ' ' << eq[1][1] << ' ' << eq[1][0] << ' ' << eq[0][1] << endl;
		halt(0);
	}

	c1 = delta_c1 / delta;
	c2 = delta_c2 / delta;

	if (c1 != c1) {
		cout << "c1 nan" << endl << delta_c1 << ' ' << delta << ' ' << delta_c1 / delta << endl;
		cout << "determ" << endl << eq[0][1] << ' ' << eq[1][2] << ' ' << eq[1][1] << ' ' << eq[0][2] << endl;
		halt(0);
	}
	if (c2 != c2) {
		cout << "c2 nan" << endl << delta_c2 << ' ' << delta << ' ' << delta_c2 / delta << endl;
		halt(0);
	}
	//cout << "c " << c1 << "\t" << c2 << endl;
	//system("pause");
}

complex <double> lapl(double z, complex <double> lam) {
	M = complex <double>(3.0 * (myu + lam / c) * (myu_a + lam / c));
	if (M == complex <double>(0, 0)) {
		cout << "M = 0" << endl << lam << endl;
		halt(0);
	}
	int_h = 1.0 / (lam + log(2));
	m = (0.5 / (myu + lam / c));
	ex1 = exp(sqrt(M));
	ex2 = exp(-sqrt(M));
	if (ex1 == complex <double>(0, 0)) {
		cout << "ex1 = 0" << endl;
		halt(0);
	}

	complex <double> eq[2][3] = {
		{ 1.0 - m * sqrt(M),		1.0 + m * sqrt(M),			int_h },
		{ ex1 + m * sqrt(M) * ex1,	ex2 - m * sqrt(M) * ex2,	0 }
	};

	complex <double> c1, c2;
	slau(eq, c1, c2);

	exz1 = pow(ex1, z);
	exz2 = pow(ex2, z);
	return c1 * exz1 + c2 * exz2;
	//return (double(1) / lam);
}

complex <double> fun(double z, complex <double> lam) {
	return lam / (lam*lam + 1.); //cos(w*t) w=1
	//return 1. / lam;
}

int main() {
	complex <double> *f_0 = new complex <double>[n]();
	double *t = new double[n + 1]();
	t[0] = a;

	complex <double> ikd, ikd1, ikdn, f1, f2, fk, fn;
	ikdn = im * double(N) * delta;
	double q = delta / (2.0 * M_PI);

	/*for (int k = 1; k < (N - 1); k++) {
		f_0[0] = fun(z, (l1 + im * double(k) * delta)) + fun(z, (l1 + im * double(k) * delta));
		//f_0[0] = lapl(z, (l1 + im * double(k) * delta)) + lapl(z, (l1 + im * double(k) * delta));
	}*/

	for (int i = 0; i < 5; i++) {                                          //|w|<= w2
		if (i < (n - 1))
			t[i + 1] = t[i] + h;

		double w = t[i] * delta / 2;
		complex <double> D1 = exp(-w*im);
		complex <double> D2 = exp(w*im);

		fk = (0, 0);
		for (int k = 1; k < (N - 1); k++) {
			ikd = im * double(k) * delta;
			ikd1 = im * double(k + 1) * delta;
			f1 = fun(z, (l1 + ikd)) * D1;
			f2 = fun(z, (l1 + ikd1)) * D2;

			f1 = lapl(z, (l1 + ikd)) * D1;
			f2 = lapl(z, (l1 + ikd1)) * D2;
			//f1 = fun(z, (l1 + ikd)) * exp(ikd * t[i]);
			//f2 = fun(z, (l1 + ikd1)) * exp(ikd1 * t[i]);
			//f1 = lapl(z, (l1 + ikd)) * exp(ikd * t[i]);
			//f2 = lapl(z, (l1 + ikd1)) * exp(ikd1 * t[i]);
			//fk += (f1 + f2);
			fk += (f1 + f2) * exp(im * t[i] * delta * (double(k) + 0.5));
		}
		f_0[i] = exp(l1*t[i]) * fk;
	}

	for (int i = 5; i < n; i++) {											//|w|> w2
		if (i < (n - 1))
			t[i + 1] = t[i] + h;

		double w = t[i] * delta / 2;
		complex <double> D1 = sin(w) / w + (w*cos(w) - sin(w))*im / (w*w);
		complex <double> D2 = sin(w) / w - (w*cos(w) - sin(w))*im / (w*w);
		if (D1 != D1) {
			cout << "D1 nan" << endl << sin(w) << ' ' << (w*cos(w) - sin(w))*im << ' ' << t[i] << ' ' << i << ' ' << t[i - 1] << endl;
			halt(0);
		}

		fk = (0, 0);
		for (int k = 1; k < (N - 1); k++) {
			ikd = im * double(k) * delta;
			ikd1 = im * double(k + 1) * delta;
			//f1 = fun(z, (l1 + ikd)) * D1;
			//f2 = fun(z, (l1 + ikd1)) * D2;
			f1 = lapl(z, (l1 + ikd)) * D1;
			f2 = lapl(z, (l1 + ikd1)) * D2;
			fk += (f1 + f2) * exp( im * t[i] * delta * (double(k) + 0.5) );
		}
		f_0[i] = exp(l1*t[i]) * fk;
		//cout << q*f_0[i] << endl;
	}

	ofstream fin("output.txt");
	fin << "[";
	
	/*for (int t = 0; t < (n - 1); t++) {
	fin << q * f_0[t] << "," << endl;
	}
	fin << q * f_0[n - 1] << "];";
	*/
	
	for (int i = 0; i < (n - 1); i++) {
		fin << q * f_0[i] << "," << endl;
		//fin << cos(t[i]) << "     t =  " << t[i] << endl;
		//fin << q * f_0[t] - complex <double>(1.0) << "," << endl;
		//fin << q * f_0[t] << "," << endl;
	}
	fin <<  q * f_0[n - 1] << "];";

	//fin << q * f_0[n - 1] << "];";
	fin.close();
	system("pause");
	delete[] t, f_0;
	return 0;
}