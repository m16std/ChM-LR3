#include "gus.h" //мой заголовок со вспомогательными функциями
#include "poly345.h" //решалка уравнений 5го порялка
#include <iomanip>

void Le_Verrier()
{
	std::cout << "\n\nLe Verrier Algorithm: " << "\n";

	Matrix* A, * A2, * A3, * A4, * A5;
	double n, determinant, S1, S2, S3, S4, S5, P1, P2, P3, P4, P5, x[5];

	n = 5;
	A = CreateSquareMatrix(n);
	SetValuesToMatrix(A);
	PrintMatrix(A);
	determinant = Determinant(A);

	A2 = CreateSquareMatrix(n);
	A3 = CreateSquareMatrix(n);
	A4 = CreateSquareMatrix(n);
	A5 = CreateSquareMatrix(n);

	Multiply(A2, A, A);
	Multiply(A3, A2, A);
	Multiply(A4, A3, A);
	Multiply(A5, A4, A);

	S1 = Trace(A);
	S2 = Trace(A2);
	S3 = Trace(A3);
	S4 = Trace(A4);
	S5 = Trace(A5);

	P1 = S1;
	P2 = (S2 - P1 * S1) / 2.0;
	P3 = (S3 - P1 * S2 - P2 * S1) / 3.0;
	P4 = (S4 - P1 * S3 - P2 * S2 - P3 * S1) / 4.0;
	P5 = (S5 - P1 * S4 - P2 * S3 - P3 * S2 - P4 * S1) / 5.0;

	//0 = -1 * (a ^ 5 - P1 * a ^ 4 - P2 * a ^ 3 - P3 * a ^ 2 - P4 * a - P5)

	SolveP5(x, -P1, -P2, -P3, -P4, -P5);

	std::cout << "\n\tEigenvalues: " << "\n";
	for (int i = 0; i < n; i++)
		std::cout << "\n\t" << x[i];

	FreeMatrix(A);
	FreeMatrix(A2);
	FreeMatrix(A3);
	FreeMatrix(A4);
	FreeMatrix(A5);
}

void Fadeev() //Без выбора главного элемента
{
	std::cout << "\n\nFadeev Method: " << "\n";

	Matrix* A[6], * B[6], * X[6], * A0, * E, * E_sample;
	double n, determinant, P[6], x[6];

	n = 5;

	A0 = CreateSquareMatrix(n);
	SetValuesToMatrix(A0);
	determinant = Determinant(A0);

	for (int i = 0; i <= n; i++)
	{
		A[i] = CreateSquareMatrix(n);
		B[i] = CreateSquareMatrix(n);
		X[i] = CreateMatrix(n, 1);
	}

	E = CreateSquareMatrix(n);
	E_sample = CreateIdentityMatrix(n);

	Assign(E, E_sample);
	Assign(A[1], A0);
	P[1] = Trace(A[1]);
	ScalarMultiply(E, P[1]);
	Subtract(A[1], E);
	Assign(B[1], A[1]);

	for (int i = 2; i <= n; i++)
	{
		Assign(E, E_sample);
		Multiply(A[i], A0, B[i - 1]);
		P[i] = Trace(A[i]) / (double)i;
		ScalarMultiply(E, P[i]);
		Subtract(A[i], E);
		Assign(B[i], A[i]);
	}

	std::cout << "\n\n\tMatrices A1, A2 .. An: " << "\n\n";

	for (int k = 1; k <= n; k++)
	{
		PrintMatrix(A[k]);
	}

	std::cout << "\n\n\tMatrices B1, B2 .. Bn: " << "\n\n";

	for (int k = 1; k <= n; k++)
	{
		PrintMatrix(B[k]);
	}

	std::cout << "\n\n\tEigenvalues: " << "\n\n";

	SolveP5(x, -P[1], -P[2], -P[3], -P[4], -P[5]);

	for (int i = 0; i < n; i++)
	{
		std::cout << "\n\t" << x[i];
	}

	std::cout << "\n\n\n\tInverse matrix: " << "\n\n";
	std::cout.precision(3);

	ScalarDivide(B[4], P[5]);

	PrintMatrix(B[4]);

	std::cout << "\n\n\tInverse matrix multiply original matrix: " << "\n\n";  //поверка на обратность

	Multiply(B[5], A0, B[4]);

	PrintMatrix(B[5]);

	std::cout.precision(5);
	for (int i = 1; i <= n; i++)
	{
		Assign(E, E_sample);        //берём единичную матрицу
		Assign(A[i], A0);           //записываем в ячеку массиваа оригинальную матрицу
		ScalarMultiply(E, x[i - 1]);  //умножаем единичную на собственное число
		Subtract(A[i], E);          //вычитаем из A
		Gaussian(A[i], X[i]);       //решаем СЛАУ
	}

	std::cout << "\n\n\tVectors X1, X2 .. Xn: " << "\n\n\n";

	for (int k = 1; k <= n; k++)
	{
		for (int j = 0; j < n; j++)
			std::cout << "\t" << Element(X[k], j, 0) << "\n";
		std::cout << "\n";
	}

	for (int i = 0; i <= n; i++)
	{
		FreeMatrix(A[i]);
		FreeMatrix(B[i]);
		FreeMatrix(X[i]);
	}
}

void Danilevsky()
{
	std::cout << "\n\nDanilevsky Method: " << "\n";

	Matrix* A, * S, * S2, * S3, * inverseM, * M;
	double x[5];
	int n = 5;
	A = CreateSquareMatrix(n);

	S2 = CreateSquareMatrix(n);
	S3 = CreateSquareMatrix(n);

	SetValuesToMatrix(A);

	for (int i = n - 2; i >= 0; i--)
	{
		inverseM = CreateIdentityMatrix(n);
		M = CreateIdentityMatrix(n);

		for (int k = 0; k < n; k++)
		{
			inverseM->r->at(i)->c->at(k) = A->r->at(i + 1)->c->at(k);
			if (i == k)
				M->r->at(i)->c->at(k) = 1 / A->r->at(i + 1)->c->at(i);
			else
				M->r->at(i)->c->at(k) = -A->r->at(i + 1)->c->at(k) / A->r->at(i + 1)->c->at(i);
		}

		std::cout << "\n\tinverseM: " << "\n";
		PrintMatrix(inverseM);
		std::cout << "\n\tM: " << "\n";
		PrintMatrix(M);

		Multiply(S2, inverseM, A);
		Multiply(S3, S2, M);

		std::cout << "\n\tA: " << "\n";
		PrintMatrix(S3);

		Assign(A, S3);
	}

	SolveP5(x, -Element(A, 0, 0), -Element(A, 0, 1), -Element(A, 0, 2), -Element(A, 0, 3), -Element(A, 0, 4));

	std::cout << "\n\tEigenvalues: " << "\n";
	for (int i = 0; i < n; i++)
		std::cout << "\n\t" << x[i];
}

void SimpleIteration()
{
	std::cout << "\n\nSimple Iteration Method: " << "\n\n";

	const int n = 5;
	const int N = 50;

	double X[N][n], Y[N][n], Lamda[N], e = 0.0001;

	Matrix* A;

	A = CreateSquareMatrix(n);
	SetValuesToMatrix(A);

	for (int k = 0; k < N - 2; k++)
	{
		X[0][k] = 2;
	}

	for (int k = 0; k < N - 2; k++)
	{
		VectorMatrixMultiply(Y[k + 1], A, X[k]);
		Lamda[k + 1] = VectorVectorMultiply(Y[k + 1], X[k]);
		VectorValueMultiply(X[k + 1], Y[k + 1], 1 / VectorLength(Y[k + 1]));
		std::cout << "\t" << Lamda[k + 1] << "\n";
		if (abs(Lamda[k + 1] - Lamda[k]) < e) break;
	}
}

void StraightIteration()
{
	std::cout << "\n\Straight Iteration Method: " << "\n\n";

	const int n = 5;
	const int N = 50;

	double X[N][n], Y[N][n], Lamda[N], c[N], a[N], e = 0.0001;

	Matrix* A;

	A = CreateSquareMatrix(n);
	SetValuesToMatrix(A);

	for (int k = 0; k < 5; k++)
	{
		X[0][k] = 1;
	}

	for (int k = 1; k < N - 2; k++)
	{
		a[k] = VectorMax(X[k - 1]);
		VectorValueMultiply(X[k - 1], X[k - 1], 1 / a[k]);
		VectorMatrixMultiply(X[k], A, X[k - 1]);

		for (int i = 0; i < n; i++)
		{
			c[i] = X[k][i] / X[k - 1][i];
		}

		Lamda[k] = (c[0] + c[1] + c[2] + c[3] + c[4]) / 5.0;

		std::cout << "\t" << Lamda[k] << "\n";

		if (k > 1 && abs(Lamda[k - 1] - Lamda[k]) < e) break;
	}
}

void ReverseIteration()
{
	std::cout << "\n\Reverse Iteration Method: " << "\n\n";

	const int n = 5;
	const int N = 50;

	double X[N][n], Y[N][n], Lamda[N], c[N], a[N], e = 0.0001;

	Matrix* A;

	A = CreateSquareMatrix(n);
	SetValuesToMatrix(A);

	for (int k = 0; k < n; k++)
	{
		X[0][k] = 1;
	}

	for (int k = 1; k < N - 2; k++)
	{
		a[k] = VectorMax(X[k - 1]);
		VectorValueMultiply(X[k - 1], X[k - 1], 1 / a[k]);

		Gaussian(A, X[k - 1], X[k]);

		Lamda[k] = 1/a[k];

		std::cout << "\n\t" << Lamda[k];

		if (k > 1 && abs(Lamda[k - 1] - Lamda[k]) < e) break;
	}
}

void Householder()
{
	std::cout << "\n\nHouseholder Method: " << "\n";

	Matrix* A = CreateMatrix(5, 5);
	SetValuesToMatrix(A);

	Matrix* QR, * R, * Q, * r, * q,* N, * ra, * rq;
	vector <double> Rdiag = { 0,0,0,0,0,0 };
	int i, j, k, count=0, flag;
	double m, n, s, nrm, E = 0.00001, Lamda[5], e;

	for (int t = 0; t < 500; t++) {
		QR = CreateCopyOfMatrix(A);
		m = NumberOfRows(A);
		n = NumberOfColumns(A);
		r = CreateSquareMatrix(n);
		q = CreateSquareMatrix(n);
		for (k = 0; k < n; k++) {
			/* Compute 2-norm of k-th column
			without under/overflow. */
			nrm = 0;
			for (i = k; i < m; i++) {
				nrm = Hypothenuse(nrm, Element(QR, i, k));
			}

			if (nrm != 0) {
				/* Form k-th Householder vector. */
				if (Element(QR, k, k) < 0) {
					nrm = -nrm;
				}
				for (i = k; i < m; i++) {
					QR->r->at(i)->c->at(k) = Element(QR, i, k) / nrm;
				}
				QR->r->at(k)->c->at(k) = Element(QR, k, k) + 1;

				/* Apply transformation to remaining columns. */
				for (j = k + 1; j < n; j++) {
					s = 0.0;
					for (i = k; i < m; i++) {
						s = s + Element(QR, i, k) * Element(QR, i, j);
					}
					s = -s / Element(QR, k, k);
					for (i = k; i < m; i++) {
						QR->r->at(i)->c->at(j) = Element(QR, i, j) + s * Element(QR, i, k);
					}
				}
			}
			Rdiag.at(k) = -nrm;
		}

		R = CreateSquareMatrix(n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i < j) {
					R->r->at(i)->c->at(j) = Element(QR, i, j);
				}
				else if (i == j) {
					R->r->at(i)->c->at(j) = Rdiag.at(i);
				}
				else {
					R->r->at(i)->c->at(j) = 0.0;
				}
			}
		}
		Assign(r, R);

		Q = CreateMatrix(m, n);
		for (k = n - 1; k >= 0; k--) {
			for (i = 0; i < m; i++) {
				Q->r->at(i)->c->at(k) = 0;
			}
			Q->r->at(k)->c->at(k) = 1;
			for (j = k; j < n; j++) {
				if (Element(QR, k, k) != 0) {
					s = 0;
					for (i = k; i < m; i++) {
						s = s + Element(QR, i, k) * Element(Q, i, j);
					}
					s = -s / Element(QR, k, k);
					for (i = k; i < m; i++) {
						Q->r->at(i)->c->at(j) = Element(Q, i, j) + s * Element(QR, i, k);
					}
				}
			}
		}
		Assign(q, Q);

		/* Adjust for positive R. */
		n = NumberOfRows(r);

		N = CreateIdentityMatrix(n);

		for (i = 0; i < n; i++) {
			e = Element(r, i, i);

			if (e < 0.0) {
				N->r->at(i)->c->at(i) = -1;
			}
		}

		ra = MultiplyToNew(N, r);
		Assign(r, ra);
		rq = MultiplyToNew(q, N);
		Assign(q, rq);
		Multiply(A, r, q);
		count++;

		flag = 0;
		for (i = 0; i < n; i++)
		{
			if (t > 1 && abs(Element(A, i, i) - Lamda[i]) > E) flag = 1;
		}

		if (t > 1 && !flag) break;

		for (i = 0; i < n; i++)
		{
			Lamda[i] = Element(A, i, i);
		}
	}
	std::cout << "\n\tA:\n";
	PrintMatrix(A);
	std::cout << "\n\tEigenvalues: " << "\n";
	for (int i = 0; i < n; i++)
		std::cout << "\n\t" << Lamda[i];
	std::cout << "\n\n\n\t" << count;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	std::cout.precision(5);

	Le_Verrier();
	Fadeev();
	Danilevsky();
	SimpleIteration();
	StraightIteration();
	ReverseIteration();
	Householder();
	std::cout << "\n\n\n\tWork done\n";
}