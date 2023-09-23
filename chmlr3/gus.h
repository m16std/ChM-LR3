#pragma once
#include "Header.h"

double a1 = -4.5, b1 = -1.8, c1 = 0.26, d1 = 4.61, e1 = 8.86;
double a2 = 7.21, b2 = 4.15, c2 = -5.7, d2 = -20, e2 = -28;
double a3 = -2.9, b3 = -3.8, c3 = -2.9, d3 = 8.49, e3 = 5.32;
double a4 = -1.6, b4 = 1.01, c4 = -0.3, d4 = -9.9, e4 = -6.9;
double a5 = -3.4, b5 = -3, c5 = -0.9, d5 = 2.34, e5 = -2.8;


double x[5];

/*
-12, 9.65, -3.4, -4.6, -24
-9.1, 10.5, 0.53, -3.9, -9.8
-8, 5.51, 2.86, -3.7, -9
4.99, -1.7, 2.17, 2.87, 13.2
5, -2.8, 0.56, 0.35, 13.5
*/

/*
double matrix2[5][5] =
{
-5.496032, -4.153618, 11.60026057, 1.557029, -3.164633,
-2.225157, -0.193178, -19.97299538, -2.907853, 6.0289,
-3.626436, -6.880966, -6.903825654, -6.496388, -1.037427,
3.141419, 0.927751, 9.797346915, -1.489755, -2.778278,
-4.149118, -9.196849, 0.470804217, -8.367168, -7.91721
};
for (int i = 0; i < 5; i++)
	for (int j = 0; j < 5; j++)
		A->r->at(i)->c->at(j) = matrix2[j][i];
*/

void Gaussian(Matrix* A, Matrix* X) //С выбором главного элемента
{

	double a[5][6];
	int i, j, m, k;
	double aa, bb;

	for (k = 0; k < 5; k++)
	{
		for (j = 0; j < 5; j++)
			a[k][j] = A->r->at(k)->c->at(j);
		a[k][5] = 0;
	}

	for (k = 0; k < 5; k++)
	{
		aa = abs(a[k][k]);
		i = k;
		for (m = k + 1; m < 5; m++)
			if (abs(a[m][k]) > aa)
			{
				i = m;
				aa = abs(a[m][k]);
			}

		if (i != k)
		{
			for (j = k; j < 6; j++)
			{
				bb = a[k][j];
				a[k][j] = a[i][j];
				a[i][j] = bb;
			}
		}
		aa = a[k][k];
		a[k][k] = 1;

		for (j = k + 1; j < 6; j++)
			a[k][j] = a[k][j] / aa;

		for (i = k + 1; i < 5; i++)
		{
			bb = a[i][k];
			a[i][k] = 0;
			if (bb != 0)
				for (j = k + 1; j < 6; j++)
					a[i][j] = a[i][j] - bb * a[k][j];
		}
	}

	X->r->at(4)->c->at(0) = 1;

	for (i = 3; i >= 0; i--)
	{
		X->r->at(i)->c->at(0) = 0;
		aa = a[i][5];
		for (j = 4; j > i; j--)
			aa = aa - a[i][j] * X->r->at(j)->c->at(0);
		X->r->at(i)->c->at(0) = aa;
	}
}

void SetValuesToMatrix(Matrix* A)
{
	A->r->at(0)->c->at(0) = a1;
	A->r->at(0)->c->at(1) = b1;
	A->r->at(0)->c->at(2) = c1;
	A->r->at(0)->c->at(3) = d1;
	A->r->at(0)->c->at(4) = e1;

	A->r->at(1)->c->at(0) = a2;
	A->r->at(1)->c->at(1) = b2;
	A->r->at(1)->c->at(2) = c2;
	A->r->at(1)->c->at(3) = d2;
	A->r->at(1)->c->at(4) = e2;

	A->r->at(2)->c->at(0) = a3;
	A->r->at(2)->c->at(1) = b3;
	A->r->at(2)->c->at(2) = c3;
	A->r->at(2)->c->at(3) = d3;
	A->r->at(2)->c->at(4) = e3;

	A->r->at(3)->c->at(0) = a4;
	A->r->at(3)->c->at(1) = b4;
	A->r->at(3)->c->at(2) = c4;
	A->r->at(3)->c->at(3) = d4;
	A->r->at(3)->c->at(4) = e4;

	A->r->at(4)->c->at(0) = a5;
	A->r->at(4)->c->at(1) = b5;
	A->r->at(4)->c->at(2) = c5;
	A->r->at(4)->c->at(3) = d5;
	A->r->at(4)->c->at(4) = e5;
}

void PrintMatrix(Matrix* A)
{
	std::cout << "\n";
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
			std::cout << "\t" << Element(A, i, j);
		std::cout << "\n";
	}
}

double VectorLength(double* A)
{
	double sum = 0;
	for (int i = 0; i < 5; i++)
	{
		sum += pow(A[i], 2);
	}
	sum = sqrt(sum);
	return sum;
}

double VectorMax(double* A)
{
	double max = -1000;
	for (int i = 0; i < 5; i++)
	{
		if (abs(A[i]) > max)
			max = abs(A[i]);
	}
	return max;
}

double VectorVectorMultiply(double* A, double* B)
{
	double sum = 0;
	for (int i = 0; i < 5; i++)
	{
		sum += A[i] * B[i];
	}
	return sum;
}

void VectorValueMultiply(double* X, double* A, double B)
{
	for (int i = 0; i < 5; i++)
	{
		X[i] = A[i] * B;
	}
}

void VectorMatrixMultiply(double* C, Matrix* B, double* A)
{

	for (int i = 0; i < 5; i++)
	{
		double sum = 0;
		for (int j = 0; j < 5; j++)
		{
			sum += A[j] * Element(B, i, j);
		}
		C[i] = sum;
	}
}

void Gaussian(Matrix* A, double* B, double* x)
{
	double a[5][6];
	int i, j, m, k;
	double aa, bb;

	for (k = 0; k < 5; k++)
	{
		for (j = 0; j < 5; j++)
			a[k][j] = A->r->at(k)->c->at(j);
		a[k][5] = B[k];
	}

	for (k = 0; k < 5; k++)
	{
		aa = abs(a[k][k]);
		i = k;
		for (m = k + 1; m < 5; m++)
			if (abs(a[m][k]) > aa)
			{
				i = m;
				aa = abs(a[m][k]);
			}

		if (i != k)
		{
			for (j = k; j < 6; j++)
			{
				bb = a[k][j];
				a[k][j] = a[i][j];
				a[i][j] = bb;
			}
		}
		aa = a[k][k];
		a[k][k] = 1;

		for (j = k + 1; j < 6; j++)
			a[k][j] = a[k][j] / aa;

		for (i = k + 1; i < 5; i++)
		{
			bb = a[i][k];
			a[i][k] = 0;
			if (bb != 0)
				for (j = k + 1; j < 6; j++)
					a[i][j] = a[i][j] - bb * a[k][j];
		}
	}

	for (i = 4; i >= 0; i--)
	{
		x[i] = 0;
		aa = a[i][5];
		for (j = 4; j > i; j--)
			aa = aa - a[i][j] * x[j];
		x[i] = aa;
	}
}