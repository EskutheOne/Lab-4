#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "winbgi2.h"
#include "rk4.h"

double euler(double y, double t, double h, double (*pf)(double, double));
double fun(double t, double y);
double anali(double y, double t, double t0);
void scan(double* x);
double round(double num);		//funkcja zaokraglajaca ktora natywnie nie wystepuje w c tylko w c++ dlatego tutaj dodana
double l = 1;
double p = 6;		//najwieksza potega liczby N kotra wyznacza ilosc krokow

void main()
{

	double h, N, Ee = 0, Er = 0;		//zmienne obliczane
	double t, ye, yr;		//zmienne uzywane w programie
	double tm, t0, y0;		//zmienne podane przez uzytkownika
	double p1 = 0, p2 = 0, p3 = 0;			//zmienne pomocnicze
	double *yee, *yrr;		//tablice przechowujace wartosci obliczone numerycznie w poprzednim kroku
	double rzade, rzadrk;		//zmienne rzedu obu metod

	FILE* g;
	g = fopen("wynik.txt", "w");

	if (!g)		//oblsuga pliku
	{

		printf("Wystapil blad przy otwieraniu pliku!");
		exit(1);

	}

	fprintf(g, "L.p.\tLiczba krokow N\tDlugosc kroku h\tBlad met. Eulera\tBlad met. RK4\n");		//naglowek do pliku


	printf("Podaj gorna krawedz przedzialu [tm]: \n");
	scan(&tm);
	
	printf("Podaj warunek poczatkowy [t_0, y_0]: \n");
	do
	{
		scan(&t0);
		if (t0 >= tm) printf("t0 nie moze byc wieksze niz tm! \n");
	} while (t0 >= tm);
	scan(&y0);

	printf("Podaj wartosc lambda: \n");
	scan(&l);


	yee = (double*)malloc((p+1) * sizeof(double));
	yrr = (double*)malloc((p+1) * sizeof(double));

	if (yee == NULL && yrr == NULL)
	{

		printf("Wystapil blad alokacji :( \n");
		exit(1);

	}


	graphics(600,600);			//oblsuga grafiki
	scale(0.0, 0.0, (tm-t0), 2.2);
	title("Blad metod:	Eulera [Czerwony]\tRK4 [Zolty]", "", "");

	
	for (int j = 0; j < (p + 1); j++)
	{

		int k = 0;		//zmienna pomocnicza
		double Ee_max = -10000, Er_max = -10000;		//zmienne pomocnicze

		N = pow(2.0, j);
		h = (tm - t0) / N;
		t = t0;
		ye = y0;
		yr = y0;
		while ((tm - t) > h / 2.0)
		{
			t += h;

			yr = rk4(t, yr, h, fun);
			printf(" Met. RK4 t[%i] = %f , y[%i] = %f , N = %f , h = %f \n", k, t, k, yr, N, h);
			Er = (fabs((yr - anali(y0, t, t0))) / fabs(anali(y0, t, t0)));		//blad metody RK4
			printf("Blad met. RK4 = %lf \n\n", Er);
			setcolor(0.7);		  //zolty
			point(h, Er);
			if (Er_max < Er) Er_max = Er;

			ye = euler(ye, t, h, fun);
			printf(" Metoda Eulera t[%i] = %f , y[%i] = %f , N = %f , h = %f \n", k, t, k, ye, N, h);
			Ee = (fabs((ye - anali(y0, t, t0))) / fabs(anali(y0, t, t0)));		//blad metody Eulera
			printf("Blad met. Eulera = %lf \n\n", Ee);
			setcolor(1);		//czerwony
			point(h, Ee);
			if (Ee_max < Ee) Ee_max = Ee;

			k++;

		}
	
		if (j>0)
		{

			rzade = round(log(fabs((ye - anali(y0, t, t0)))) / log(fabs((yee[j-1] - anali(y0, t, t0)))));
			rzadrk = round(log(fabs((yr - anali(y0, t, t0)))) / log(fabs((yrr[j-1] - anali(y0, t, t0)))));		//szacowanie rzêdu zbie¿noœci
			printf("\n RZAD ZBIEZONSCI E = %lf RZAD ZBIEZONOSCI RK = %lf \n\n", rzade, rzadrk);

		}

		yee[j] = ye;
		yrr[j] = yr;

		if (j > 0)		//rysowanie wykresu	
		{

			line(p1, p2, h, Ee_max);
			setcolor(0.7);
			line(p1, p3, h, Er_max);

		}
		
		p1 = h;
		p2 = Ee_max;
		p3 = Er_max;

		fprintf(g, "%d\t%lf\t%lf\t%lf\t%lf\n", j, N, h, Ee, Er);   //zapisuje dane do pliku na koncu kroku

	}
	
	free(yee);
	free(yrr);
	fclose(g);
	wait();

}


void scan(double *x)		//funkcja do wczytywania i sprawdzania zgodnosci inputu
{
	double temp;
	while (scanf("%lf", &temp) != 1)
	{
		printf("Nieprawid³owy format danych! \n");
		int n;
		while ((n = getchar()) != EOF && n != '\n');
	}
	*x = temp;
}

double fun(double t, double y)	 //funkcja liczaca prawa strone rownania, wywolywana jest z dwiema zmiennymi poniewaz potrzebuje tego funkcja rk4
{

	return l*y;

}

double euler(double y, double t, double h, double (*pf)(double, double))
{

	return y += h * pf(t, y);

}

double anali(double y0, double t, double t0)	 //analityczne rozwiazanie
{

	return y0 * exp(l * (t - t0));

}

double round(double num)
{
	return (num >= 0) ? (int)(num + 0.5) : (int)(num - 0.5);
}