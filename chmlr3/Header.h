// Downloaded from https://repo.progsbase.com - Code Developed Using progsbase.

#include <cmath>
#include <cstring>
#include <vector>
#include <cwchar>
#include <iostream>

using namespace std;

#define toVector(s) (new vector<wchar_t> ((s), (s) + wcslen(s)))

struct Matrix;

struct MatrixArrayReference;

struct MatrixReference;

struct MatrixRow;

struct ComplexMatrix;

struct ComplexMatrixArrayReference;

struct ComplexMatrixReference;

struct ComplexMatrixRow;

struct BooleanArrayReference;

struct BooleanReference;

struct CharacterReference;

struct NumberArrayReference;

struct NumberReference;

struct StringArrayReference;

struct StringReference;

struct cComplexNumber;

struct cComplexNumberArrayReference;

struct cComplexNumberReference;

struct cPolarComplexNumber;

struct pComplexPolynomial;

struct LinkedListNodeStrings;

struct LinkedListStrings;

struct LinkedListNodeNumbers;

struct LinkedListNumbers;

struct LinkedListCharacters;

struct LinkedListNodeCharacters;

struct DynamicArrayNumbers;

struct Matrix {
	vector<MatrixRow*>* r;
};

struct MatrixArrayReference {
	vector<Matrix*>* matrices;
};

struct MatrixReference {
	Matrix* matrix;
};

struct MatrixRow {
	vector<double>* c;
};

struct ComplexMatrix {
	vector<ComplexMatrixRow*>* r;
};

struct ComplexMatrixArrayReference {
	vector<ComplexMatrix*>* matrices;
};

struct ComplexMatrixReference {
	ComplexMatrix* matrix;
};

struct ComplexMatrixRow {
	vector<cComplexNumber*>* c;
};

struct BooleanArrayReference {
	vector<bool>* booleanArray;
};

struct BooleanReference {
	bool booleanValue;
};

struct CharacterReference {
	wchar_t characterValue;
};

struct NumberArrayReference {
	vector<double>* numberArray;
};

struct NumberReference {
	double numberValue;
};

struct StringArrayReference {
	vector<StringReference*>* stringArray;
};

struct StringReference {
	vector<wchar_t>* string;
};

struct cComplexNumber {
	double re;
	double im;
};

struct cComplexNumberArrayReference {
	vector<cComplexNumber*>* complexNumbers;
};

struct cComplexNumberReference {
	cComplexNumber* complexNumbers;
};

struct cPolarComplexNumber {
	double r;
	double phi;
};

struct pComplexPolynomial {
	vector<cComplexNumber*>* cs;
};

struct LinkedListNodeStrings {
	bool end;
	vector<wchar_t>* value;
	LinkedListNodeStrings* next;
};

struct LinkedListStrings {
	LinkedListNodeStrings* first;
	LinkedListNodeStrings* last;
};

struct LinkedListNodeNumbers {
	LinkedListNodeNumbers* next;
	bool end;
	double value;
};

struct LinkedListNumbers {
	LinkedListNodeNumbers* first;
	LinkedListNodeNumbers* last;
};

struct LinkedListCharacters {
	LinkedListNodeCharacters* first;
	LinkedListNodeCharacters* last;
};

struct LinkedListNodeCharacters {
	bool end;
	wchar_t value;
	LinkedListNodeCharacters* next;
};

struct DynamicArrayNumbers {
	vector<double>* array;
	double length;
};

void Add(Matrix* a, Matrix* b);
void Assign(Matrix* A, Matrix* B);
void Resize(Matrix* A, double r, double c);
void Subtract(Matrix* a, Matrix* b);
Matrix* SubtractToNew(Matrix* a, Matrix* b);
void ScalarMultiply(Matrix* A, double b);
void ScalarDivide(Matrix* A, double b);
void ElementWisePower(Matrix* A, double p);
Matrix* ScalarMultiplyToNew(Matrix* A, double b);
Matrix* MultiplyToNew(Matrix* a, Matrix* b);
void Multiply(Matrix* x, Matrix* a, Matrix* b);
Matrix* CreateSquareMatrix(double d);
Matrix* CreateMatrix(double rows, double cols);
Matrix* CreateIdentityMatrix(double d);
void Transpose(Matrix* a);
void TransposeAssign(Matrix* t, Matrix* a);
Matrix* TransposeToNew(Matrix* a);
void CofactorOfMatrix(Matrix* mat, Matrix* temp, double p, double q, double n);
double DeterminantOfSubmatrix(Matrix* mat, double n);
double Determinant(Matrix* m);
void Adjoint(Matrix* A, Matrix* adj);
bool Inverse(Matrix* A, Matrix* inverseResult);
bool InverseUsingAdjoint(Matrix* A, Matrix* inverseResult);
bool InverseUsingLUDecomposition(Matrix* A, Matrix* inverseResult);
bool LUDecomposition(Matrix* A, Matrix* L, Matrix* U);
bool IsSymmetric(Matrix* A);
bool IsSquare(Matrix* A);
bool Cholesky(Matrix* A, Matrix* L);
void Clear(Matrix* a);
void Fill(Matrix* a, double value);
double Element(Matrix* matrix, double m, double n);
double Trace(Matrix* a);
Matrix* ColumnCombineMatricesToNew(Matrix* A, Matrix* B);
double NumberOfRows(Matrix* A);
double NumberOfColumns(Matrix* A);
vector<double>* CharacteristicPolynomial(Matrix* A);
void CharacteristicPolynomialWithInverse(Matrix* A, Matrix* AInverse, NumberArrayReference* cp, NumberReference* determinant);
void FaddeevLeVerrierAlgorithm(Matrix* A, Matrix* AInverse, NumberArrayReference* cp, NumberReference* determinant);
Matrix* InverseUsingCharacteristicPolynomial(Matrix* A);
bool Eigenvalues(Matrix* A, NumberArrayReference* eigenValuesReference);
bool EigenvaluesUsingQRAlgorithm(Matrix* A, NumberArrayReference* eigenValuesReference, double precision, double maxIterations);
bool EigenvaluesUsingLaguerreIterations(Matrix* A, NumberArrayReference* eigenValuesReference);
void GaussianElimination(Matrix* A);
Matrix* GaussianEliminationToNew(Matrix* A);
Matrix* CreateCopyOfMatrix(Matrix* A);
void SwapRows(Matrix* A, double to, double from);
void UnnormalizeVector(vector<double>* numberArray);
bool InversePowerMethod(Matrix* A, double eigenvalue, double maxIterations, NumberArrayReference* eigenvector);
bool Eigenvectors(Matrix* A, MatrixArrayReference* eigenVectorsReference);
bool Eigenpairs(Matrix* A, NumberArrayReference* eigenValuesReference, MatrixArrayReference* eigenVectorsReference);
bool EigenpairsUsingQRAlgorithmAndInversePowerMethod(Matrix* M, NumberArrayReference* eigenValuesReference, MatrixArrayReference* eigenVectorsReference, double precision, double maxIterations);
bool CheckEigenpairPrecision(Matrix* a, double lambda, Matrix* e, double precision);
bool EigenvectorsLaguerreIterationsAndGaussianEliminations(Matrix* A, MatrixArrayReference* eigenVectorsReference);
bool RowIsZero(Matrix* X, double r);
void FreeMatrix(Matrix* X);
void FreeMatrixRows(vector<MatrixRow*>* r);
Matrix* CreateDiagonalMatrixFromArray(vector<double>* array);
Matrix* CreateMatrixFromRowCopies(vector<double>* row, double times);
void ExtractDiagonal(Matrix* X, vector<double>* diag);
vector<double>* ExtractDiagonalToNew(Matrix* X);
bool MatrixEqualsEpsilon(Matrix* a, Matrix* b, double epsilon);
Matrix* Minor(Matrix* x, double row, double column);
void QRDecomposition(Matrix* m, Matrix* Q, Matrix* R);
void HouseholderTriangularizationAlgorithm(Matrix* m, Matrix* Qout, Matrix* Rout);
void HouseholderMethod(Matrix* A, Matrix* q, Matrix* r);
double Hypothenuse(double a, double b);
double Norm(Matrix* a);
Matrix* ExtractSubMatrix(Matrix* M, double r1, double r2, double c1, double c2);
bool QRAlgorithm(Matrix* M, Matrix* R, Matrix* A, Matrix* Q, double precision, double maxIterations);
bool InvertUpperTriangularMatrix(Matrix* A, Matrix* inverse);
bool InvertLowerTriangularMatrix(Matrix* A, Matrix* inverse);
bool ParseMatrixFromString(MatrixReference* aref, vector<wchar_t>* matrixString, StringReference* errorMessage);
void FreeStringReferenceArray(vector<StringReference*>* stringReferencesArray);
vector<wchar_t>* MatrixToString(Matrix* matrix, double digitsAfterPoint);
double RoundToDigits(double element, double digitsAfterPoint);
vector<wchar_t>* MatrixArrayToString(vector<Matrix*>* matrices, double digitsAfterPoint);
void RoundMatrixElementsToDigits(Matrix* a, double digits);
bool SingularValueDecomposition(Matrix* Ap, MatrixReference* URef, MatrixReference* SigmaRef, MatrixReference* VRef);

ComplexMatrix* CreateComplexMatrix(double rows, double cols);
ComplexMatrix* CreateComplexMatrixFromMatrix(Matrix* a);
Matrix* CreateReMatrixFromComplexMatrix(ComplexMatrix* a);
Matrix* CreateImMatrixFromComplexMatrix(ComplexMatrix* a);
double NumberOfRowsComplex(ComplexMatrix* A);
double NumberOfColumnsComplex(ComplexMatrix* A);
cComplexNumber* IndexComplex(ComplexMatrix* a, double m, double n);
void AddComplex(ComplexMatrix* a, ComplexMatrix* b);
void SubtractComplex(ComplexMatrix* a, ComplexMatrix* b);
ComplexMatrix* SubtractComplexToNew(ComplexMatrix* a, ComplexMatrix* b);
void MultiplyComplex(ComplexMatrix* x, ComplexMatrix* a, ComplexMatrix* b);
ComplexMatrix* MultiplyComplexToNew(ComplexMatrix* a, ComplexMatrix* b);
void Conjugate(ComplexMatrix* a);
void AssignComplexMatrix(ComplexMatrix* A, ComplexMatrix* B);
void ScalarMultiplyComplex(ComplexMatrix* A, cComplexNumber* b);
ComplexMatrix* ScalarMultiplyComplexToNew(ComplexMatrix* A, cComplexNumber* b);
void ScalarDivideComplex(ComplexMatrix* A, cComplexNumber* b);
void ElementWisePowerComplex(ComplexMatrix* A, double p);
ComplexMatrix* CreateComplexIdentityMatrix(double d);
ComplexMatrix* CreateSquareComplexMatrix(double d);
void ClearComplex(ComplexMatrix* a);
void FillComplex(ComplexMatrix* a, double re, double im);
cComplexNumber* TraceComplex(ComplexMatrix* a);
void CofactorOfComplexMatrix(ComplexMatrix* mat, ComplexMatrix* temp, double p, double q, double n);
cComplexNumber* DeterminantOfComplexSubmatrix(ComplexMatrix* mat, double n);
void DeleteComplexMatrix(ComplexMatrix* X);
cComplexNumber* DeterminantComplex(ComplexMatrix* m);
void AdjointComplex(ComplexMatrix* A, ComplexMatrix* adj);
bool InverseComplex(ComplexMatrix* A, ComplexMatrix* inverseResult);
bool ComplexMatrixEqualsEpsilon(ComplexMatrix* b, ComplexMatrix* f, double epsilon);
ComplexMatrix* MinorComplex(ComplexMatrix* x, double row, double column);
void AssignComplex(ComplexMatrix* A, ComplexMatrix* B);
ComplexMatrix* CreateCopyOfComplexMatrix(ComplexMatrix* A);
bool TransposeComplex(ComplexMatrix* a);
bool ConjugateTransposeComplex(ComplexMatrix* a);
bool IsSquareComplexMatrix(ComplexMatrix* A);
void TransposeComplexAssign(ComplexMatrix* t, ComplexMatrix* a);
ComplexMatrix* TransposeComplexToNew(ComplexMatrix* a);
ComplexMatrix* ExtractComplexSubMatrix(ComplexMatrix* M, double r1, double r2, double c1, double c2);
double NormComplex(ComplexMatrix* a);
void ComplexCharacteristicPolynomial(ComplexMatrix* A, pComplexPolynomial* p);
void ComplexCharacteristicPolynomialWithInverse(ComplexMatrix* A, ComplexMatrix* AInverse, pComplexPolynomial* p, cComplexNumber* determinant);
void FaddeevLeVerrierAlgorithmComplex(ComplexMatrix* A, ComplexMatrix* AInverse, pComplexPolynomial* p, cComplexNumber* determinant);
bool EigenvaluesComplex(ComplexMatrix* A, cComplexNumberArrayReference* eigenValuesReference);
bool EigenvectorsComplex(ComplexMatrix* A, ComplexMatrixArrayReference* eigenVectorsReference);
bool InversePowerMethodComplex(ComplexMatrix* A, cComplexNumber* eigenvalue, double maxIterations, cComplexNumberArrayReference* eigenvector);
bool EigenpairsComplex(ComplexMatrix* M, cComplexNumberArrayReference* eigenValuesReference, ComplexMatrixArrayReference* eigenVectorsReference);
bool ComplexEigenpairsUsingDurandKernerAndInversePowerMethod(ComplexMatrix* M, cComplexNumberArrayReference* eigenValuesReference, ComplexMatrixArrayReference* eigenVectorsReference, double precision, double maxIterations);
bool CheckComplexEigenpairPrecision(ComplexMatrix* a, cComplexNumber* lambda, ComplexMatrix* e, double precision);

void testDeterminant(NumberReference* failures);
void testInverse(NumberReference* failures);
Matrix* createAnswerInverse4x4_1();
Matrix* createExample4x4_1();
void testInverseAlt2(NumberReference* failures);
void TestCholesky(NumberReference* failures);
void testCharacteristicPolynomial(NumberReference* failures);
void testCharacteristicPolynomial2(NumberReference* failures);
void testCharacteristicPolynomial3(NumberReference* failures);
void testEigenValues(NumberReference* failures);
void testEigenValues2(NumberReference* failures);
void testGaussianEliminationInto(NumberReference* failures);
void testEigenpairs(NumberReference* failures);
void testQRDecomposition(NumberReference* failures);
void testQRDecomposition2(NumberReference* failures);
void testQRDecomposition3(NumberReference* failures);
void testEigenValues3QR(NumberReference* failures);
void SetsEqualEpsilon(vector<double>* a, vector<double>* b, double epsilon, NumberReference* failures);
void testLUDecomposition(NumberReference* failures);
void TestSingularValueDecomposition1(NumberReference* failures);
void CheckSingularValueDecomposition(Matrix* a, MatrixReference* URef, MatrixReference* SRef, MatrixReference* VRef, bool success, NumberReference* failures);
void TestSingularValueDecomposition2(NumberReference* failures);
void TestSingularValueDecomposition3(NumberReference* failures);
void TestSingularValueDecomposition4(NumberReference* failures);
void TestSingularValueDecomposition5(NumberReference* failures);

Matrix* CreateMatrixWithComplexEigenvalues();
Matrix* CreateMatrixWithRepeatedEigenvalues();

double test();

void testComplexMatrices(NumberReference* failures);
void testComplexDeterminant(NumberReference* failures);
ComplexMatrix* CreateExampleComplexMatrix2x2();
ComplexMatrix* CreateExampleComplexMatrix2x2Inverse();
void testComplexInverse(NumberReference* failures);
void testComplexCharacteristicPolynomial(NumberReference* failures);
void testComplexCharacteristicPolynomial2(NumberReference* failures);
void testComplexEigenValues(NumberReference* failures);
void testComplexEigenValues2(NumberReference* failures);
ComplexMatrix* CreateComplex3x3Matrix();
void testComplexInversePowerMethod(NumberReference* failures);
void testComplexEigenpairs(NumberReference* failures);

bool FindRoots(vector<double>* p, NumberArrayReference* rootsReference);
bool LaguerresMethodWithRepeatedDivision(vector<double>* p, double maxIterations, double precision, double guess, NumberArrayReference* rootsReference);
bool LaguerresMethod(vector<double>* p, double guess, double maxIterations, double precision, NumberReference* rootReference);
bool DurandKernerMethod(vector<double>* p, double precision, double maxIterations, NumberArrayReference* rootsReference);

bool FindRootsComplex(pComplexPolynomial* p, cComplexNumberArrayReference* rootsReference);
bool DurandKernerMethodComplex(pComplexPolynomial* p, double precision, double maxIterations, cComplexNumberArrayReference* rootsReference);

BooleanReference* CreateBooleanReference(bool value);
BooleanArrayReference* CreateBooleanArrayReference(vector<bool>* value);
BooleanArrayReference* CreateBooleanArrayReferenceLengthValue(double length, bool value);
void FreeBooleanArrayReference(BooleanArrayReference* booleanArrayReference);
CharacterReference* CreateCharacterReference(wchar_t value);
NumberReference* CreateNumberReference(double value);
NumberArrayReference* CreateNumberArrayReference(vector<double>* value);
NumberArrayReference* CreateNumberArrayReferenceLengthValue(double length, double value);
void FreeNumberArrayReference(NumberArrayReference* numberArrayReference);
StringReference* CreateStringReference(vector<wchar_t>* value);
StringReference* CreateStringReferenceLengthValue(double length, wchar_t value);
void FreeStringReference(StringReference* stringReference);
StringArrayReference* CreateStringArrayReference(vector<StringReference*>* strings);
StringArrayReference* CreateStringArrayReferenceLengthValue(double length, vector<wchar_t>* value);
void FreeStringArrayReference(StringArrayReference* stringArrayReference);

double mathNegate(double x);
double mathPositive(double x);
double mathFactorial(double x);
double mathRound(double x);
double mathBankersRound(double x);
double mathCeil(double x);
double mathFloor(double x);
double mathTruncate(double x);
double mathAbsolute(double x);
double mathLogarithm(double x);
double mathNaturalLogarithm(double x);
double mathSin(double x);
double mathCos(double x);
double mathTan(double x);
double mathAsin(double x);
double mathAcos(double x);
double mathAtan(double x);
double mathAtan2(double y, double x);
double mathSquareroot(double x);
double mathExp(double x);
bool mathDivisibleBy(double a, double b);
double mathCombinations(double n, double k);
double mathPermutations(double n, double k);
bool mathEpsilonCompare(double a, double b, double epsilon);
double mathGreatestCommonDivisor(double a, double b);
double mathGCDWithSubtraction(double a, double b);
bool mathIsInteger(double a);
bool mathGreatestCommonDivisorWithCheck(double a, double b, NumberReference* gcdReference);
double mathLeastCommonMultiple(double a, double b);
double mathSign(double a);
double mathMax(double a, double b);
double mathMin(double a, double b);
double mathPower(double a, double b);
double mathGamma(double x);
double mathLogGamma(double x);
double mathLanczosApproximation(double z);
double mathBeta(double x, double y);
double mathSinh(double x);
double mathCosh(double x);
double mathTanh(double x);
double mathCot(double x);
double mathSec(double x);
double mathCsc(double x);
double mathCoth(double x);
double mathSech(double x);
double mathCsch(double x);
double mathError(double x);
double mathErrorInverse(double x);
double mathFallingFactorial(double x, double n);
double mathRisingFactorial(double x, double n);
double mathHypergeometric(double a, double b, double c, double z, double maxIterations, double precision);
double mathHypergeometricDirect(double a, double b, double c, double z, double maxIterations, double precision);
double mathBernouilliNumber(double n);
double mathAkiyamaTanigawaAlgorithm(double n);

cComplexNumber* cCreateComplexNumber(double re, double im);
cPolarComplexNumber* cCreatePolarComplexNumber(double r, double phi);
void cAdd(cComplexNumber* z1, cComplexNumber* z2);
cComplexNumber* cAddToNew(cComplexNumber* z1, cComplexNumber* z2);
void cSub(cComplexNumber* z1, cComplexNumber* z2);
cComplexNumber* cSubToNew(cComplexNumber* z1, cComplexNumber* z2);
void cMul(cComplexNumber* z1, cComplexNumber* z2);
cComplexNumber* cMulToNew(cComplexNumber* z1, cComplexNumber* z2);
void cDiv(cComplexNumber* z1, cComplexNumber* z2);
cComplexNumber* cDivToNew(cComplexNumber* z1, cComplexNumber* z2);
void cConjugate(cComplexNumber* z);
cComplexNumber* cConjugateToNew(cComplexNumber* z);
double cAbs(cComplexNumber* z);
double cArg(cComplexNumber* z);
cPolarComplexNumber* cCreatePolarFromComplexNumber(cComplexNumber* z);
cComplexNumber* cCreateComplexFromPolar(cPolarComplexNumber* p);
double cRe(cComplexNumber* z);
double cIm(cComplexNumber* z);
void cAddPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
cPolarComplexNumber* cAddPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
void cSubPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
cPolarComplexNumber* cSubPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
void cMulPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
cPolarComplexNumber* cMulPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
void cDivPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
cPolarComplexNumber* cDivPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2);
void cSquareRoot(cComplexNumber* z);
void cPowerPolar(cPolarComplexNumber* p, double n);
cComplexNumber* cPowerToNew(cComplexNumber* z, double n);
void cPower(cComplexNumber* z, double n);
void cNegate(cComplexNumber* z);
void cAssignComplexByValues(cComplexNumber* s, double re, double im);
void cAssignComplex(cComplexNumber* a, cComplexNumber* b);
bool cEpsilonCompareComplex(cComplexNumber* a, cComplexNumber* b, double epsilon);
void cExpComplex(cComplexNumber* x);
void cSineComplex(cComplexNumber* x);
void cCosineComplex(cComplexNumber* x);
vector<wchar_t>* cComplexToString(cComplexNumber* a);

vector<wchar_t>* pPolynomialToTextDirect(vector<double>* p, vector<wchar_t>* x);
void pLinkedListCharactersAddString(LinkedListCharacters* ll, vector<wchar_t>* str);
void pGenerateCommonRenderSpecification(vector<double>* p, BooleanArrayReference* showCoefficient, StringReference* sign, NumberArrayReference* coefficient, BooleanArrayReference* showPower, BooleanArrayReference* showX);
vector<wchar_t>* pPolynomialToText(vector<double>* p, vector<wchar_t>* x);
vector<wchar_t>* pComplexPolynomialToTextDirect(pComplexPolynomial* p, vector<wchar_t>* x);

void pAdd(vector<double>* a, vector<double>* b);
void pSubtract(vector<double>* a, vector<double>* b);
void pMultiply(vector<double>* c, vector<double>* a, vector<double>* b);
void pDivide(vector<double>* q, vector<double>* r, vector<double>* n, vector<double>* d);
bool pIsZero(vector<double>* a);
void pAssign(vector<double>* a, vector<double>* b);
vector<double>* pCreatePolynomial(double deg);
void pFill(vector<double>* p, double value);
double pDegree(vector<double>* A);
double pLead(vector<double>* A);
double pEvaluate(vector<double>* A, double x);
double pEvaluateWithHornersMethod(vector<double>* A, double x);
double pEvaluateWithPowers(vector<double>* A, double x);
double pEvaluateDerivative(vector<double>* A, double x, double n);
void pDerivative(vector<double>* A);

void pAddComplex(pComplexPolynomial* a, pComplexPolynomial* b);
void pSubtractComplex(pComplexPolynomial* a, pComplexPolynomial* b);
bool pIsZeroComplex(pComplexPolynomial* a);
void pAssignComplex(pComplexPolynomial* a, pComplexPolynomial* b);
pComplexPolynomial* pCreateComplexPolynomial(double deg);
void pFillComplex(pComplexPolynomial* p, double re, double im);
double pDegreeComplex(pComplexPolynomial* A);
cComplexNumber* pLeadComplex(pComplexPolynomial* A);
cComplexNumber* pEvaluateComplex(pComplexPolynomial* A, cComplexNumber* x);

double pTotalNumberOfRoots(vector<double>* p);

void WriteStringToStingStream(vector<wchar_t>* stream, NumberReference* index, vector<wchar_t>* src);
void WriteCharacterToStingStream(vector<wchar_t>* stream, NumberReference* index, wchar_t src);
void WriteBooleanToStingStream(vector<wchar_t>* stream, NumberReference* index, bool src);

bool SubstringWithCheck(vector<wchar_t>* string, double from, double to, StringReference* stringReference);
vector<wchar_t>* Substring(vector<wchar_t>* string, double from, double to);
vector<wchar_t>* AppendString(vector<wchar_t>* s1, vector<wchar_t>* s2);
vector<wchar_t>* ConcatenateString(vector<wchar_t>* s1, vector<wchar_t>* s2);
vector<wchar_t>* AppendCharacter(vector<wchar_t>* string, wchar_t c);
vector<wchar_t>* ConcatenateCharacter(vector<wchar_t>* string, wchar_t c);
vector<StringReference*>* SplitByCharacter(vector<wchar_t>* toSplit, wchar_t splitBy);
bool IndexOfCharacter(vector<wchar_t>* string, wchar_t character, NumberReference* indexReference);
bool SubstringEqualsWithCheck(vector<wchar_t>* string, double from, vector<wchar_t>* substring, BooleanReference* equalsReference);
bool SubstringEquals(vector<wchar_t>* string, double from, vector<wchar_t>* substring);
bool IndexOfString(vector<wchar_t>* string, vector<wchar_t>* substring, NumberReference* indexReference);
bool ContainsCharacter(vector<wchar_t>* string, wchar_t character);
bool ContainsString(vector<wchar_t>* string, vector<wchar_t>* substring);
void ToUpperCase(vector<wchar_t>* string);
void ToLowerCase(vector<wchar_t>* string);
bool EqualsIgnoreCase(vector<wchar_t>* a, vector<wchar_t>* b);
vector<wchar_t>* ReplaceString(vector<wchar_t>* string, vector<wchar_t>* toReplace, vector<wchar_t>* replaceWith);
vector<wchar_t>* ReplaceCharacter(vector<wchar_t>* string, wchar_t toReplace, wchar_t replaceWith);
vector<wchar_t>* Trim(vector<wchar_t>* string);
bool StartsWith(vector<wchar_t>* string, vector<wchar_t>* start);
bool EndsWith(vector<wchar_t>* string, vector<wchar_t>* end);
vector<StringReference*>* SplitByString(vector<wchar_t>* toSplit, vector<wchar_t>* splitBy);
bool StringIsBefore(vector<wchar_t>* a, vector<wchar_t>* b);

vector<wchar_t>* nCreateStringScientificNotationDecimalFromNumber(double decimal);
vector<wchar_t>* nCreateStringDecimalFromNumber(double decimal);
bool nCreateStringFromNumberWithCheck(double decimal, double base, StringReference* stringReference);
double nGetMaximumDigitsForBase(double base);
double nGetFirstDigitPosition(double decimal, double base);
bool nGetSingleDigitCharacterFromNumberWithCheck(double c, double base, CharacterReference* characterReference);
vector<wchar_t>* nGetDigitCharacterTable();

bool nCreateNumberFromDecimalStringWithCheck(vector<wchar_t>* string, NumberReference* decimalReference, StringReference* errorMessage);
double nCreateNumberFromDecimalString(vector<wchar_t>* string);
bool nCreateNumberFromStringWithCheck(vector<wchar_t>* string, double base, NumberReference* numberReference, StringReference* errorMessage);
double nCreateNumberFromParts(double base, bool numberIsPositive, vector<double>* beforePoint, vector<double>* afterPoint, bool exponentIsPositive, vector<double>* exponent);
bool nExtractPartsFromNumberString(vector<wchar_t>* n, double base, BooleanReference* numberIsPositive, NumberArrayReference* beforePoint, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages);
bool nExtractPartsFromNumberStringFromSign(vector<wchar_t>* n, double base, double i, NumberArrayReference* beforePoint, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages);
bool nExtractPartsFromNumberStringFromPointOrExponent(vector<wchar_t>* n, double base, double i, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages);
bool nExtractPartsFromNumberStringFromExponent(vector<wchar_t>* n, double base, double i, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages);
double nGetNumberFromNumberCharacterForBase(wchar_t c, double base);
bool nCharacterIsNumberCharacterInBase(wchar_t c, double base);
vector<double>* nStringToNumberArray(vector<wchar_t>* str);
bool nStringToNumberArrayWithCheck(vector<wchar_t>* str, NumberArrayReference* numberArrayReference, StringReference* errorMessage);

void AssertFalse(bool b, NumberReference* failures);
void AssertTrue(bool b, NumberReference* failures);
void AssertEquals(double a, double b, NumberReference* failures);
void AssertBooleansEqual(bool a, bool b, NumberReference* failures);
void AssertCharactersEqual(wchar_t a, wchar_t b, NumberReference* failures);
void AssertStringEquals(vector<wchar_t>* a, vector<wchar_t>* b, NumberReference* failures);
void AssertNumberArraysEqual(vector<double>* a, vector<double>* b, NumberReference* failures);
void AssertBooleanArraysEqual(vector<bool>* a, vector<bool>* b, NumberReference* failures);
void AssertStringArraysEqual(vector<StringReference*>* a, vector<StringReference*>* b, NumberReference* failures);

vector<double>* AddNumber(vector<double>* list, double a);
void AddNumberRef(NumberArrayReference* list, double i);
vector<double>* RemoveNumber(vector<double>* list, double n);
double GetNumberRef(NumberArrayReference* list, double i);
void RemoveNumberRef(NumberArrayReference* list, double i);

vector<StringReference*>* AddString(vector<StringReference*>* list, StringReference* a);
void AddStringRef(StringArrayReference* list, StringReference* i);
vector<StringReference*>* RemoveString(vector<StringReference*>* list, double n);
StringReference* GetStringRef(StringArrayReference* list, double i);
void RemoveStringRef(StringArrayReference* list, double i);

vector<bool>* AddBoolean(vector<bool>* list, bool a);
void AddBooleanRef(BooleanArrayReference* list, bool i);
vector<bool>* RemoveBoolean(vector<bool>* list, double n);
bool GetBooleanRef(BooleanArrayReference* list, double i);
void RemoveDecimalRef(BooleanArrayReference* list, double i);


LinkedListStrings* CreateLinkedListString();
void LinkedListAddString(LinkedListStrings* ll, vector<wchar_t>* value);
vector<StringReference*>* LinkedListStringsToArray(LinkedListStrings* ll);
double LinkedListStringsLength(LinkedListStrings* ll);
void FreeLinkedListString(LinkedListStrings* ll);


LinkedListNumbers* CreateLinkedListNumbers();
vector<LinkedListNumbers*>* CreateLinkedListNumbersArray(double length);
void LinkedListAddNumber(LinkedListNumbers* ll, double value);
double LinkedListNumbersLength(LinkedListNumbers* ll);
double LinkedListNumbersIndex(LinkedListNumbers* ll, double index);
void LinkedListInsertNumber(LinkedListNumbers* ll, double index, double value);
void LinkedListSet(LinkedListNumbers* ll, double index, double value);
void LinkedListRemoveNumber(LinkedListNumbers* ll, double index);
void FreeLinkedListNumbers(LinkedListNumbers* ll);
void FreeLinkedListNumbersArray(vector<LinkedListNumbers*>* lls);
vector<double>* LinkedListNumbersToArray(LinkedListNumbers* ll);
LinkedListNumbers* ArrayToLinkedListNumbers(vector<double>* array);
bool LinkedListNumbersEqual(LinkedListNumbers* a, LinkedListNumbers* b);

LinkedListCharacters* CreateLinkedListCharacter();
void LinkedListAddCharacter(LinkedListCharacters* ll, wchar_t value);
vector<wchar_t>* LinkedListCharactersToArray(LinkedListCharacters* ll);
double LinkedListCharactersLength(LinkedListCharacters* ll);
void FreeLinkedListCharacter(LinkedListCharacters* ll);



DynamicArrayNumbers* CreateDynamicArrayNumbers();
DynamicArrayNumbers* CreateDynamicArrayNumbersWithInitialCapacity(double capacity);
void DynamicArrayAddNumber(DynamicArrayNumbers* da, double value);
void DynamicArrayNumbersIncreaseSize(DynamicArrayNumbers* da);
bool DynamicArrayNumbersDecreaseSizeNecessary(DynamicArrayNumbers* da);
void DynamicArrayNumbersDecreaseSize(DynamicArrayNumbers* da);
double DynamicArrayNumbersIndex(DynamicArrayNumbers* da, double index);
double DynamicArrayNumbersLength(DynamicArrayNumbers* da);
void DynamicArrayInsertNumber(DynamicArrayNumbers* da, double index, double value);
void DynamicArraySet(DynamicArrayNumbers* da, double index, double value);
void DynamicArrayRemoveNumber(DynamicArrayNumbers* da, double index);
void FreeDynamicArrayNumbers(DynamicArrayNumbers* da);
vector<double>* DynamicArrayNumbersToArray(DynamicArrayNumbers* da);
DynamicArrayNumbers* ArrayToDynamicArrayNumbersWithOptimalSize(vector<double>* array);
DynamicArrayNumbers* ArrayToDynamicArrayNumbers(vector<double>* array);
bool DynamicArrayNumbersEqual(DynamicArrayNumbers* a, DynamicArrayNumbers* b);
LinkedListNumbers* DynamicArrayNumbersToLinkedList(DynamicArrayNumbers* da);
DynamicArrayNumbers* LinkedListToDynamicArrayNumbers(LinkedListNumbers* ll);

vector<wchar_t>* AddCharacter(vector<wchar_t>* list, wchar_t a);
void AddCharacterRef(StringReference* list, wchar_t i);
vector<wchar_t>* RemoveCharacter(vector<wchar_t>* list, double n);
wchar_t GetCharacterRef(StringReference* list, double i);
void RemoveCharacterRef(StringReference* list, double i);

wchar_t charToLowerCase(wchar_t character);
wchar_t charToUpperCase(wchar_t character);
bool charIsUpperCase(wchar_t character);
bool charIsLowerCase(wchar_t character);
bool charIsLetter(wchar_t character);
bool charIsNumber(wchar_t character);
bool charIsWhiteSpace(wchar_t character);
bool charIsSymbol(wchar_t character);
bool charCharacterIsBefore(wchar_t a, wchar_t b);

vector<double>* StringToNumberArray(vector<wchar_t>* string);
vector<wchar_t>* NumberArrayToString(vector<double>* array);
bool NumberArraysEqual(vector<double>* a, vector<double>* b);
bool BooleanArraysEqual(vector<bool>* a, vector<bool>* b);
bool StringsEqual(vector<wchar_t>* a, vector<wchar_t>* b);
void FillNumberArray(vector<double>* a, double value);
void FillString(vector<wchar_t>* a, wchar_t value);
void FillBooleanArray(vector<bool>* a, bool value);
bool FillNumberArrayRange(vector<double>* a, double value, double from, double to);
bool FillBooleanArrayRange(vector<bool>* a, bool value, double from, double to);
bool FillStringRange(vector<wchar_t>* a, wchar_t value, double from, double to);
vector<double>* CopyNumberArray(vector<double>* a);
vector<bool>* CopyBooleanArray(vector<bool>* a);
vector<wchar_t>* CopyString(vector<wchar_t>* a);
bool CopyNumberArrayRange(vector<double>* a, double from, double to, NumberArrayReference* copyReference);
bool CopyBooleanArrayRange(vector<bool>* a, double from, double to, BooleanArrayReference* copyReference);
bool CopyStringRange(vector<wchar_t>* a, double from, double to, StringReference* copyReference);
bool IsLastElement(double length, double index);
vector<double>* CreateNumberArray(double length, double value);
vector<bool>* CreateBooleanArray(double length, bool value);
vector<wchar_t>* CreateString(double length, wchar_t value);
void SwapElementsOfNumberArray(vector<double>* A, double ai, double bi);
void SwapElementsOfStringArray(StringArrayReference* A, double ai, double bi);
void ReverseNumberArray(vector<double>* array);

void Add(Matrix* a, Matrix* b) {
	double m, n;
	double r, c;

	r = NumberOfRows(a);
	c = NumberOfColumns(a);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			a->r->at(m)->c->at(n) = Element(a, m, n) + Element(b, m, n);
		}
	}
}
void Assign(Matrix* A, Matrix* B) {
	double m, n;
	double r, c;

	r = NumberOfRows(A);
	c = NumberOfColumns(A);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			A->r->at(m)->c->at(n) = Element(B, m, n);
		}
	}
}
void Resize(Matrix* A, double r, double c) {
	double m, n, ar, ac;
	Matrix* C;

	C = CreateMatrix(r, c);

	ar = NumberOfRows(A);
	ac = NumberOfColumns(A);

	for (m = 0.0; m < fmin(r, ar); m = m + 1.0) {
		for (n = 0.0; n < fmin(c, ac); n = n + 1.0) {
			C->r->at(m)->c->at(n) = Element(A, m, n);
		}
	}

	FreeMatrixRows(A->r);
	A->r = C->r;
}
void Subtract(Matrix* a, Matrix* b) {
	double m, n;
	double r, c;

	r = NumberOfRows(a);
	c = NumberOfColumns(a);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			a->r->at(m)->c->at(n) = Element(a, m, n) - Element(b, m, n);
		}
	}
}
Matrix* SubtractToNew(Matrix* a, Matrix* b) {
	Matrix* X;

	X = CreateCopyOfMatrix(a);
	Subtract(X, b);

	return X;
}
void ScalarMultiply(Matrix* A, double b) {
	double m, n;
	double r, c;

	r = NumberOfRows(A);
	c = NumberOfColumns(A);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			A->r->at(m)->c->at(n) = b * A->r->at(m)->c->at(n);
		}
	}
}
void ScalarDivide(Matrix* A, double b) {
	double m, n;
	double r, c;

	r = NumberOfRows(A);
	c = NumberOfColumns(A);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			A->r->at(m)->c->at(n) = Element(A, m, n) / b;
		}
	}
}
void ElementWisePower(Matrix* A, double p) {
	double m, n;
	double r, c;

	r = NumberOfRows(A);
	c = NumberOfColumns(A);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			A->r->at(m)->c->at(n) = pow(A->r->at(m)->c->at(n), p);
		}
	}
}
Matrix* ScalarMultiplyToNew(Matrix* A, double b) {
	Matrix* matrix;

	matrix = CreateCopyOfMatrix(A);
	ScalarMultiply(matrix, b);

	return matrix;
}
Matrix* MultiplyToNew(Matrix* a, Matrix* b) {
	double rows, cols;
	Matrix* x;

	rows = NumberOfRows(a);
	cols = NumberOfColumns(b);
	x = CreateMatrix(rows, cols);
	Multiply(x, a, b);

	return x;
}
void Multiply(Matrix* x, Matrix* a, Matrix* b) {
	double m, n;
	double rows, cols, d;
	double i, s;

	rows = NumberOfRows(a);
	cols = NumberOfColumns(b);
	d = NumberOfColumns(a);

	for (m = 0.0; m < rows; m = m + 1.0) {
		for (n = 0.0; n < cols; n = n + 1.0) {
			s = 0.0;

			for (i = 0.0; i < d; i = i + 1.0) {
				s = s + a->r->at(m)->c->at(i) * b->r->at(i)->c->at(n);
			}

			x->r->at(m)->c->at(n) = s;
		}
	}
}
Matrix* CreateSquareMatrix(double d) {
	double m, n;
	Matrix* matrix;

	matrix = new Matrix();
	matrix->r = new vector<MatrixRow*>(d);
	for (m = 0.0; m < d; m = m + 1.0) {
		matrix->r->at(m) = new MatrixRow();
		matrix->r->at(m)->c = new vector<double>(d);
		for (n = 0.0; n < d; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = 0.0;
		}
	}

	return matrix;
}
Matrix* CreateMatrix(double rows, double cols) {
	double m, n;
	Matrix* matrix;

	matrix = new Matrix();
	matrix->r = new vector<MatrixRow*>(rows);
	for (m = 0.0; m < rows; m = m + 1.0) {
		matrix->r->at(m) = new MatrixRow();
		matrix->r->at(m)->c = new vector<double>(cols);
		for (n = 0.0; n < cols; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = 0.0;
		}
	}

	return matrix;
}
Matrix* CreateIdentityMatrix(double d) {
	double m;
	Matrix* matrix;

	matrix = CreateSquareMatrix(d);
	Fill(matrix, 0.0);

	for (m = 0.0; m < d; m = m + 1.0) {
		matrix->r->at(m)->c->at(m) = 1.0;
	}

	return matrix;
}
void Transpose(Matrix* a) {
	Matrix* ap;

	ap = TransposeToNew(a);

	FreeMatrixRows(a->r);
	a->r = ap->r;
}
void TransposeAssign(Matrix* t, Matrix* a) {
	double m, n;
	double rows, cols;

	cols = NumberOfRows(a);
	rows = NumberOfColumns(a);

	for (m = 0.0; m < cols; m = m + 1.0) {
		for (n = 0.0; n < rows; n = n + 1.0) {
			t->r->at(n)->c->at(m) = a->r->at(m)->c->at(n);
		}
	}
}
Matrix* TransposeToNew(Matrix* a) {
	double m, n;
	double rows, cols;
	Matrix* c;

	cols = NumberOfRows(a);
	rows = NumberOfColumns(a);

	c = CreateMatrix(rows, cols);

	for (m = 0.0; m < cols; m = m + 1.0) {
		for (n = 0.0; n < rows; n = n + 1.0) {
			c->r->at(n)->c->at(m) = a->r->at(m)->c->at(n);
		}
	}

	return c;
}
void CofactorOfMatrix(Matrix* mat, Matrix* temp, double p, double q, double n) {
	double i, j;
	double row, col;

	i = 0.0;
	j = 0.0;

	for (row = 0.0; row < n; row = row + 1.0) {
		for (col = 0.0; col < n; col = col + 1.0) {
			if (row != p && col != q) {
				temp->r->at(i)->c->at(j) = mat->r->at(row)->c->at(col);
				j = j + 1.0;

				if (j == n - 1.0) {
					j = 0.0;
					i = i + 1.0;
				}
			}
		}
	}
}
double DeterminantOfSubmatrix(Matrix* mat, double n) {
	double D, f, sign;
	Matrix* temp;

	D = 0.0;

	if (n == 1.0) {
		D = mat->r->at(0)->c->at(0);
	}
	else {
		temp = CreateSquareMatrix(n);

		sign = 1.0;

		for (f = 0.0; f < n; f = f + 1.0) {
			CofactorOfMatrix(mat, temp, 0.0, f, n);
			D = D + sign * mat->r->at(0)->c->at(f) * DeterminantOfSubmatrix(temp, n - 1.0);
			sign = -sign;
		}

		FreeMatrix(temp);
	}

	return D;
}
double Determinant(Matrix* m) {
	double D, n;

	n = NumberOfRows(m);
	D = DeterminantOfSubmatrix(m, n);

	return D;
}
void Adjoint(Matrix* A, Matrix* adj) {
	double n, sign;
	Matrix* cofactors;
	double i, j;

	n = A->r->size();

	if (n == 1.0) {
		adj->r->at(0)->c->at(0) = 1.0;
	}
	else {
		cofactors = CreateSquareMatrix(n);

		for (i = 0.0; i < n; i = i + 1.0) {
			for (j = 0.0; j < n; j = j + 1.0) {
				CofactorOfMatrix(A, cofactors, i, j, n);

				if (fmod(i + j, 2.0) == 0.0) {
					sign = 1.0;
				}
				else {
					sign = -1.0;
				}

				adj->r->at(j)->c->at(i) = sign * DeterminantOfSubmatrix(cofactors, n - 1.0);
			}
		}

		FreeMatrix(cofactors);
	}
}
bool Inverse(Matrix* A, Matrix* inverseResult) {
	return InverseUsingLUDecomposition(A, inverseResult);
}
bool InverseUsingAdjoint(Matrix* A, Matrix* inverseResult) {
	bool success;
	Matrix* adj;
	double n, i, j;
	double det;

	if (NumberOfColumns(A) == NumberOfRows(A)) {
		n = NumberOfColumns(A);

		det = Determinant(A);
		if (det != 0.0) {
			adj = CreateSquareMatrix(n);
			Adjoint(A, adj);

			for (i = 0.0; i < n; i = i + 1.0) {
				for (j = 0.0; j < n; j = j + 1.0) {
					inverseResult->r->at(i)->c->at(j) = adj->r->at(i)->c->at(j) / det;
				}
			}

			success = true;
			FreeMatrix(adj);
		}
		else {
			success = false;
		}
	}
	else {
		success = false;
	}

	return success;
}
bool InverseUsingLUDecomposition(Matrix* A, Matrix* inverseResult) {
	bool success;
	Matrix* l, * u, * li, * ui;

	l = CreateCopyOfMatrix(A);
	u = CreateCopyOfMatrix(A);
	li = CreateCopyOfMatrix(A);
	ui = CreateCopyOfMatrix(A);
	inverseResult->r = CreateCopyOfMatrix(A)->r;

	success = LUDecomposition(A, l, u);
	if (success) {
		success = InvertLowerTriangularMatrix(l, li);
		if (success) {
			success = InvertUpperTriangularMatrix(u, ui);
			if (success) {
				Multiply(inverseResult, ui, li);
			}
		}
	}

	FreeMatrix(l);
	FreeMatrix(u);
	FreeMatrix(li);
	FreeMatrix(ui);

	return success;
}
bool LUDecomposition(Matrix* A, Matrix* L, Matrix* U) {
	double n, i, j, k, sum;
	bool success;

	n = NumberOfRows(A);

	L->r = CreateSquareMatrix(n)->r;
	U->r = CreateSquareMatrix(n)->r;

	if (IsSquare(A)) {
		success = true;

		for (i = 0.0; i < n && success; i = i + 1.0) {
			for (k = i; k < n; k = k + 1.0) {
				sum = 0.0;
				for (j = 0.0; j < i; j = j + 1.0) {
					sum = sum + (Element(L, i, j) * Element(U, j, k));
				}

				U->r->at(i)->c->at(k) = Element(A, i, k) - sum;
			}

			for (k = i; k < n && success; k = k + 1.0) {
				if (i == k) {
					L->r->at(i)->c->at(i) = 1.0;
				}
				else {
					sum = 0.0;
					for (j = 0.0; j < i; j = j + 1.0) {
						sum = sum + (Element(L, k, j) * Element(U, j, i));
					}

					if (Element(U, i, i) == 0.0) {
						success = false;
					}
					else {
						L->r->at(k)->c->at(i) = (Element(A, k, i) - sum) / Element(U, i, i);
					}
				}
			}
		}
	}
	else {
		success = false;
	}

	return success;
}
bool IsSymmetric(Matrix* A) {
	double N;
	double i, j;
	bool is, done;

	N = NumberOfRows(A);

	done = false;
	is = true;
	for (i = 0.0; i < N && !done; i = i + 1.0) {
		for (j = 0.0; j < i && !done; j = j + 1.0) {
			if (A->r->at(i)->c->at(j) != A->r->at(j)->c->at(i)) {
				is = false;
				done = true;
			}
		}
	}

	return is;
}
bool IsSquare(Matrix* A) {
	bool is;

	if (NumberOfRows(A) == NumberOfColumns(A)) {
		is = true;
	}
	else {
		is = false;
	}

	return is;
}
bool Cholesky(Matrix* A, Matrix* L) {
	bool success;
	double N;
	double i, j, k, s;

	Clear(L);

	if (IsSquare(A) && IsSymmetric(A)) {
		success = true;

		N = NumberOfRows(A);

		for (i = 0.0; i < N && success; i = i + 1.0) {
			for (j = 0.0; j <= i && success; j = j + 1.0) {
				s = 0.0;
				for (k = 0.0; k < j; k = k + 1.0) {
					s = s + L->r->at(i)->c->at(k) * L->r->at(j)->c->at(k);
				}
				if (i == j) {
					L->r->at(i)->c->at(i) = sqrt(A->r->at(i)->c->at(i) - s);
				}
				else {
					L->r->at(i)->c->at(j) = 1.0 / L->r->at(j)->c->at(j) * (A->r->at(i)->c->at(j) - s);
				}
			}
			if (L->r->at(i)->c->at(i) <= 0.0) {
				success = false;
			}
		}

		success = true;
	}
	else {
		success = false;
	}

	return success;
}
void Clear(Matrix* a) {
	Fill(a, 0.0);
}
void Fill(Matrix* a, double value) {
	double m, n;

	for (m = 0.0; m < NumberOfRows(a); m = m + 1.0) {
		for (n = 0.0; n < NumberOfColumns(a); n = n + 1.0) {
			a->r->at(m)->c->at(n) = value;
		}
	}
}
double Element(Matrix* matrix, double m, double n) {
	return matrix->r->at(m)->c->at(n);
}
double Trace(Matrix* a) {
	double m;
	double d, tr;

	tr = 0.0;

	d = a->r->size();
	for (m = 0.0; m < d; m = m + 1.0) {
		tr = tr + a->r->at(m)->c->at(m);
	}

	return tr;
}
Matrix* ColumnCombineMatricesToNew(Matrix* A, Matrix* B) {
	Matrix* X;
	double m, n;

	X = CreateMatrix(NumberOfRows(A), NumberOfColumns(A) + NumberOfColumns(B));

	for (m = 0.0; m < NumberOfRows(A); m = m + 1.0) {
		for (n = 0.0; n < NumberOfColumns(A); n = n + 1.0) {
			X->r->at(m)->c->at(n) = A->r->at(m)->c->at(n);
		}
	}

	for (m = 0.0; m < NumberOfRows(B); m = m + 1.0) {
		for (n = 0.0; n < NumberOfColumns(B); n = n + 1.0) {
			X->r->at(m)->c->at(NumberOfColumns(A) + n) = B->r->at(m)->c->at(n);
		}
	}

	return X;
}
double NumberOfRows(Matrix* A) {
	return A->r->size();
}
double NumberOfColumns(Matrix* A) {
	return A->r->at(0)->c->size();
}
vector<double>* CharacteristicPolynomial(Matrix* A) {
	Matrix* dummy;
	NumberArrayReference* coeffs;
	NumberReference* determinant;

	dummy = CreateSquareMatrix(NumberOfRows(A));

	coeffs = new NumberArrayReference();
	determinant = new NumberReference();
	CharacteristicPolynomialWithInverse(A, dummy, coeffs, determinant);

	FreeMatrix(dummy);

	return coeffs->numberArray;
}
void CharacteristicPolynomialWithInverse(Matrix* A, Matrix* AInverse, NumberArrayReference* cp, NumberReference* determinant) {
	FaddeevLeVerrierAlgorithm(A, AInverse, cp, determinant);
}
void FaddeevLeVerrierAlgorithm(Matrix* A, Matrix* AInverse, NumberArrayReference* cp, NumberReference* determinant) {
	vector<double>* p;
	Matrix* Mk, * Mkm1, * t1, * I;
	double n, k;

	n = NumberOfRows(A);
	p = new vector<double>(n + 1.0);
	p->at(n) = 1.0;
	Mkm1 = CreateSquareMatrix(n);
	Fill(Mkm1, 0.0);
	I = CreateIdentityMatrix(n);
	Mk = CreateSquareMatrix(n);
	t1 = CreateSquareMatrix(n);

	for (k = 1.0; k <= n; k = k + 1.0) {
		/* M_k = A * M_(k-1) + c_(n-k+1) * I */
		Multiply(Mk, A, Mkm1);
		Assign(t1, I);
		ScalarMultiply(t1, p->at(n - k + 1.0));
		Add(Mk, t1);

		/* c_(n-k) = -1/k * trace(A * M_k) */
		Multiply(t1, A, Mk);
		p->at(n - k) = -1.0 / k * Trace(t1);

		/* done */
		Assign(Mkm1, Mk);

		if (k == n) {
			Assign(AInverse, Mk);
			determinant->numberValue = -p->at(0);
			if (p->at(0) == 0.0) {
			}
			else {
				ScalarDivide(AInverse, determinant->numberValue);
			}
		}
	}

	FreeMatrix(Mkm1);
	FreeMatrix(I);
	FreeMatrix(Mk);
	FreeMatrix(t1);

	cp->numberArray = p;
}
Matrix* InverseUsingCharacteristicPolynomial(Matrix* A) {
	Matrix* inverse;
	NumberArrayReference* coeffs;
	NumberReference* determinant;

	inverse = CreateSquareMatrix(NumberOfRows(A));
	coeffs = new NumberArrayReference();
	determinant = new NumberReference();
	CharacteristicPolynomialWithInverse(A, inverse, coeffs, determinant);
	delete coeffs->numberArray;
	delete coeffs;

	return inverse;
}
bool Eigenvalues(Matrix* A, NumberArrayReference* eigenValuesReference) {
	MatrixArrayReference* eigenVectorsReference;
	bool success;
	double i;

	eigenVectorsReference = new MatrixArrayReference();
	success = Eigenpairs(A, eigenValuesReference, eigenVectorsReference);
	if (success) {
		for (i = 0.0; i < eigenVectorsReference->matrices->size(); i = i + 1.0) {
			FreeMatrix(eigenVectorsReference->matrices->at(i));
		}
		delete eigenVectorsReference->matrices;
		delete eigenVectorsReference;
	}

	return success;
}
bool EigenvaluesUsingQRAlgorithm(Matrix* A, NumberArrayReference* eigenValuesReference, double precision, double maxIterations) {
	Matrix* x, * q, * r;
	bool success;
	double i, n, v, v1, v2, ev, found;
	vector<double>* cp;

	n = NumberOfRows(A);
	x = CreateSquareMatrix(n);
	q = CreateSquareMatrix(n);
	r = CreateSquareMatrix(n);
	eigenValuesReference->numberArray = new vector<double>(n);
	success = QRAlgorithm(A, r, x, q, precision, maxIterations);
	found = 0.0;
	if (success) {
		ExtractDiagonal(x, eigenValuesReference->numberArray);

		/* find the correct sign of the eigenvalue. */
		cp = CharacteristicPolynomial(A);
		for (i = 0.0; i < n; i = i + 1.0) {
			ev = eigenValuesReference->numberArray->at(i);

			v1 = pEvaluate(cp, ev);
			v2 = pEvaluate(cp, -ev);

			if (abs(v2) < abs(v1)) {
				eigenValuesReference->numberArray->at(i) = -ev;
				v = v2;
			}
			else {
				v = v1;
			}

			if (abs(v) < precision * pow(10.0, 4.0)) {
				found = found + 1.0;
			}
		}

		FreeMatrix(x);
		FreeMatrix(q);
		FreeMatrix(r);
	}

	if (found != n) {
		success = false;
	}

	return success;
}
bool EigenvaluesUsingLaguerreIterations(Matrix* A, NumberArrayReference* eigenValuesReference) {
	vector<double>* p;
	bool success;

	p = CharacteristicPolynomial(A);
	success = FindRoots(p, eigenValuesReference);

	return success;
}
void GaussianElimination(Matrix* A) {
	double h, k, m, n, maxElement, i, j, max, maxCandidate, f;

	m = NumberOfRows(A);
	n = NumberOfColumns(A);

	h = 0.0;
	k = 0.0;
	for (; h < m && k < n; ) {
		maxElement = h;
		max = 0.0;
		for (i = h; i < m; i = i + 1.0) {
			maxCandidate = abs(Element(A, i, k));
			if (max < maxCandidate) {
				maxElement = i;
				max = maxCandidate;
			}
		}
		if (A->r->at(maxElement)->c->at(k) == 0.0) {
			k = k + 1.0;
		}
		else {
			SwapRows(A, h, maxElement);
			for (i = h + 1.0; i < m; i = i + 1.0) {
				f = Element(A, i, k) / Element(A, h, k);
				A->r->at(i)->c->at(k) = 0.0;
				for (j = k + 1.0; j < n; j = j + 1.0) {
					A->r->at(i)->c->at(j) = Element(A, i, j) - Element(A, h, j) * f;
				}
			}
			h = h + 1.0;
			k = k + 1.0;
		}
	}
}
Matrix* GaussianEliminationToNew(Matrix* A) {
	Matrix* X;

	X = CreateCopyOfMatrix(A);
	GaussianElimination(X);

	return X;
}
Matrix* CreateCopyOfMatrix(Matrix* A) {
	Matrix* X;

	X = CreateMatrix(NumberOfRows(A), NumberOfColumns(A));
	Assign(X, A);

	return X;
}
void SwapRows(Matrix* A, double to, double from) {
	double n;
	double c, t;

	c = NumberOfRows(A);
	for (n = 0.0; n < c; n = n + 1.0) {
		t = A->r->at(to)->c->at(n);
		A->r->at(to)->c->at(n) = A->r->at(from)->c->at(n);
		A->r->at(from)->c->at(n) = t;
	}
}
void UnnormalizeVector(vector<double>* numberArray) {
	double i, m;
	bool mSet;

	mSet = false;
	m = 0.0;

	for (i = 0.0; i < numberArray->size(); i = i + 1.0) {
		if (numberArray->at(i) - mathTruncate(numberArray->at(i)) < 0.001) {
			if (!mSet) {
				m = abs(numberArray->at(i));
				mSet = true;
			}
			else {
				m = fmin(m, abs(numberArray->at(i)));
			}
		}
	}

	if (mSet) {
		for (i = 0.0; i < numberArray->size(); i = i + 1.0) {
			numberArray->at(i) = numberArray->at(i) / m;
		}
	}
}
bool InversePowerMethod(Matrix* A, double eigenvalue, double maxIterations, NumberArrayReference* eigenvector) {
	Matrix* x, * y, * z, * b, * t;
	double n, i, c;
	bool singular;

	n = NumberOfRows(A);

	x = CreateIdentityMatrix(n);
	ScalarMultiply(x, eigenvalue);
	y = SubtractToNew(A, x);
	z = CreateSquareMatrix(n);
	singular = !Inverse(y, z);
	if (singular) {
		/* Try again with more erroneous eigenvalue estimate. */
		x = CreateIdentityMatrix(n);
		ScalarMultiply(x, eigenvalue * 1.01);
		y = SubtractToNew(A, x);
		z = CreateSquareMatrix(n);
		singular = !Inverse(y, z);
	}

	if (!singular) {
		b = CreateMatrix(n, 1.0);

		for (i = 0.0; i < n; i = i + 1.0) {
			b->r->at(i)->c->at(0) = 1.0;
		}

		for (i = 0.0; i < maxIterations; i = i + 1.0) {
			t = MultiplyToNew(z, b);
			c = Norm(t);
			ScalarDivide(t, c);
			Assign(b, t);
		}

		eigenvector->numberArray = new vector<double>(n);
		for (i = 0.0; i < n; i = i + 1.0) {
			eigenvector->numberArray->at(i) = b->r->at(i)->c->at(0);
		}
	}

	return  !singular;
}
bool Eigenvectors(Matrix* A, MatrixArrayReference* eigenVectorsReference) {
	NumberArrayReference* evsReference;
	bool success;

	evsReference = new NumberArrayReference();
	success = Eigenpairs(A, evsReference, eigenVectorsReference);
	if (success) {
		delete evsReference->numberArray;
		delete evsReference;
	}

	return success;
}
bool Eigenpairs(Matrix* A, NumberArrayReference* eigenValuesReference, MatrixArrayReference* eigenVectorsReference) {
	return EigenpairsUsingQRAlgorithmAndInversePowerMethod(A, eigenValuesReference, eigenVectorsReference, 0.00000000001, 100.0);
}
bool EigenpairsUsingQRAlgorithmAndInversePowerMethod(Matrix* M, NumberArrayReference* eigenValuesReference, MatrixArrayReference* eigenVectorsReference, double precision, double maxIterations) {
	NumberArrayReference* evecReference;
	bool done, inverseSuccess;
	double i, j, k, N, v1, v2, eigenValue, withinPrecision;
	Matrix* A, * Q, * R, * eigenVector;
	vector<double>* cp;

	N = NumberOfRows(M);

	A = CreateCopyOfMatrix(M);
	Q = CreateCopyOfMatrix(M);
	R = CreateCopyOfMatrix(M);

	done = false;
	eigenVectorsReference->matrices = new vector<Matrix*>(N);
	evecReference = new NumberArrayReference();
	eigenValuesReference->numberArray = new vector<double>(N);
	cp = CharacteristicPolynomial(M);

	for (j = 0.0; j < N; j = j + 1.0) {
		eigenVectorsReference->matrices->at(j) = CreateMatrix(N, 1.0);
	}

	for (i = 0.0; i < maxIterations && !done; i = i + 1.0) {
		QRDecomposition(A, Q, R);
		Multiply(A, R, Q);

		/* Check */
		withinPrecision = 0.0;
		ExtractDiagonal(R, eigenValuesReference->numberArray);

		for (j = 0.0; j < N; j = j + 1.0) {
			/* Find the correct sign of the eigenvalue. */
			eigenValue = eigenValuesReference->numberArray->at(j);
			v1 = pEvaluate(cp, eigenValue);
			v2 = pEvaluate(cp, -eigenValue);
			if (abs(v2) < abs(v1)) {
				eigenValuesReference->numberArray->at(j) = -eigenValue;
				eigenValue = -eigenValue;
			}

			/* Calculate the eigenvector corresponding to the eigenvalue. */
			inverseSuccess = InversePowerMethod(M, eigenValue, i + 1.0, evecReference);
			if (inverseSuccess) {
				for (k = 0.0; k < N; k = k + 1.0) {
					eigenVectorsReference->matrices->at(j)->r->at(k)->c->at(0) = evecReference->numberArray->at(k);
				}

				/* Check eigenpair agains precision. */
				eigenVector = eigenVectorsReference->matrices->at(j);

				if (CheckEigenpairPrecision(M, eigenValue, eigenVector, precision)) {
					withinPrecision = withinPrecision + 1.0;
				}
			}
		}

		if (withinPrecision == N) {
			done = true;
		}
	}

	FreeMatrix(A);
	FreeMatrix(Q);
	FreeMatrix(R);
	delete evecReference;
	delete cp;

	return done;
}
bool CheckEigenpairPrecision(Matrix* a, double lambda, Matrix* e, double precision) {
	Matrix* vec1, * vec2;
	bool equal;

	vec1 = MultiplyToNew(a, e);
	vec2 = ScalarMultiplyToNew(e, lambda);

	equal = MatrixEqualsEpsilon(vec1, vec2, precision);

	return equal;
}
bool EigenvectorsLaguerreIterationsAndGaussianEliminations(Matrix* A, MatrixArrayReference* eigenVectorsReference) {
	Matrix* Id, * t1, * B, * v;
	vector<Matrix*>* eigenVectorsResult;
	bool success;
	NumberArrayReference* eigenValuesReference;
	double i, lambda, j, N, x, k;
	vector<double>* ev;

	N = NumberOfRows(A);

	eigenValuesReference = new NumberArrayReference();
	success = Eigenvalues(A, eigenValuesReference);

	eigenVectorsResult = new vector<Matrix*>(N);

	if (success) {
		ev = eigenValuesReference->numberArray;

		Id = CreateIdentityMatrix(N);
		t1 = CreateSquareMatrix(N);
		B = CreateSquareMatrix(N);

		for (j = 0.0; j < ev->size() && success; j = j + 1.0) {
			lambda = ev->at(j);

			/* B = A - lambda * Id */
			Assign(t1, Id);
			ScalarMultiply(t1, lambda);

			Assign(B, A);
			Subtract(B, t1);

			GaussianElimination(B);

			v = CreateMatrix(N, 1.0);
			v->r->at(N - 1.0)->c->at(0) = 1.0;
			for (i = N - 2.0; i >= 0.0 && success; i = i - 1.0) {
				if (!RowIsZero(B, i)) {
					x = 0.0;

					for (k = N - 1.0; k > i; k = k - 1.0) {
						x = x - Element(B, i, k) * Element(v, k, 0.0);
					}

					v->r->at(i)->c->at(0) = x / Element(B, i, i);
				}
				else {
					success = false;
				}
			}

			eigenVectorsResult->at(j) = v;
		}

		FreeMatrix(t1);
		FreeMatrix(B);
		FreeMatrix(Id);

		eigenVectorsReference->matrices = eigenVectorsResult;
	}

	return success;
}
bool RowIsZero(Matrix* X, double r) {
	bool isZero;
	double columns, i;

	isZero = true;

	columns = NumberOfColumns(X);
	for (i = 0.0; i < columns && isZero; i = i + 1.0) {
		if (Element(X, r, i) != 0.0) {
			isZero = false;
		}
	}

	return isZero;
}
void FreeMatrix(Matrix* X) {
	FreeMatrixRows(X->r);
	delete X;
}
void FreeMatrixRows(vector<MatrixRow*>* r) {
	double m, rows;

	rows = r->size();
	for (m = 0.0; m < rows; m = m + 1.0) {
		delete r->at(m)->c;
		delete r->at(m);
	}

	delete r;
}
Matrix* CreateDiagonalMatrixFromArray(vector<double>* array) {
	double m;
	Matrix* matrix;

	matrix = CreateSquareMatrix(array->size());
	Fill(matrix, 0.0);

	for (m = 0.0; m < array->size(); m = m + 1.0) {
		matrix->r->at(m)->c->at(m) = array->at(m);
	}

	return matrix;
}
Matrix* CreateMatrixFromRowCopies(vector<double>* row, double times) {
	double m, n;
	Matrix* matrix;

	matrix = CreateMatrix(times, row->size());

	for (m = 0.0; m < times; m = m + 1.0) {
		for (n = 0.0; n < row->size(); n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = row->at(n);
		}
	}

	return matrix;
}
void ExtractDiagonal(Matrix* X, vector<double>* diag) {
	double n, i;

	n = NumberOfRows(X);

	for (i = 0.0; i < n; i = i + 1.0) {
		diag->at(i) = X->r->at(i)->c->at(i);
	}
}
vector<double>* ExtractDiagonalToNew(Matrix* X) {
	vector<double>* diag;
	double n, i;

	n = NumberOfRows(X);
	diag = new vector<double>(n);

	for (i = 0.0; i < n; i = i + 1.0) {
		diag->at(i) = X->r->at(i)->c->at(i);
	}

	return diag;
}
bool MatrixEqualsEpsilon(Matrix* a, Matrix* b, double epsilon) {
	double x, y, columns, rows;
	bool equals;

	equals = true;

	if (NumberOfRows(a) == NumberOfRows(b) && NumberOfColumns(a) == NumberOfColumns(b)) {
		columns = NumberOfColumns(a);
		rows = NumberOfRows(a);

		for (x = 0.0; x < rows; x = x + 1.0) {
			for (y = 0.0; y < columns; y = y + 1.0) {
				equals = equals && mathEpsilonCompare(Element(a, x, y), Element(b, x, y), epsilon);
			}
		}
	}
	else {
		equals = false;
	}

	return equals;
}
Matrix* Minor(Matrix* x, double row, double column) {
	Matrix* theMinor;
	double cols, rows, i, j, m, n;

	rows = NumberOfRows(x) - 1.0;
	cols = NumberOfColumns(x) - 1.0;

	theMinor = CreateMatrix(rows, cols);

	for (i = 0.0; i < rows; i = i + 1.0) {
		if (i < row) {
			m = i;
		}
		else {
			m = i + 1.0;
		}

		for (j = 0.0; i != row && j < cols; j = j + 1.0) {
			if (j != column) {

				if (j < column) {
					n = j;
				}
				else {
					n = j + 1.0;
				}

				theMinor->r->at(m)->c->at(n) = x->r->at(i)->c->at(j);
			}
		}
	}

	return theMinor;
}
void QRDecomposition(Matrix* m, Matrix* Q, Matrix* R) {
	HouseholderMethod(m, Q, R);
}
void HouseholderTriangularizationAlgorithm(Matrix* m, Matrix* Qout, Matrix* Rout) {
	Matrix* P, * AA, * A, * PP, * v, * vt, * t;
	double i, j, rows, cols, r, s;

	rows = NumberOfRows(m);
	cols = NumberOfColumns(m);

	P = CreateIdentityMatrix(rows);
	A = CreateCopyOfMatrix(m);
	AA = CreateMatrix(rows, cols);

	t = CreateMatrix(rows, cols);
	PP = CreateIdentityMatrix(rows);
	vt = CreateMatrix(1.0, rows);

	for (j = 0.0; j < cols; j = j + 1.0) {
		v = ExtractSubMatrix(A, 0.0, rows - 1.0, j, j);
		if (j > 0.0) {
			for (i = 0.0; i < j; i = i + 1.0) {
				v->r->at(i)->c->at(0) = 0.0;
			}
		}

		s = mathSign(v->r->at(j)->c->at(0));
		if (s == 0.0) {
			s = 1.0;
		}
		v->r->at(j)->c->at(0) = v->r->at(j)->c->at(0) + Norm(v) * s;
		r = -2.0 / (Norm(v) * Norm(v));
		Assign(AA, A);

		TransposeAssign(vt, v);
		Multiply(t, vt, A);
		Assign(A, t);
		Multiply(t, v, A);
		Assign(A, t);
		ScalarMultiply(A, r);
		Assign(t, AA);
		Add(A, AA);

		Assign(PP, P);

		Multiply(t, vt, P);
		Assign(P, t);
		Multiply(t, v, P);
		Assign(P, t);
		ScalarMultiply(P, r);
		Assign(t, AA);
		Add(P, PP);
	}

	Assign(Rout, A);
	Assign(Qout, P);
	Transpose(Qout);
}
void HouseholderMethod(Matrix* A, Matrix* q, Matrix* r) {
	Matrix* QR, * R, * Q;
	vector<double>* Rdiag;
	double m, n;
	double i, j, k;
	double s, nrm;
	double e;
	Matrix* N, * ra, * rq;

	/* Initialize. */
	QR = CreateCopyOfMatrix(A);
	m = NumberOfRows(A);
	n = NumberOfColumns(A);
	Rdiag = new vector<double>(n);

	/* Main loop. */
	for (k = 0.0; k < n; k = k + 1.0) {
		/* Compute 2-norm of k-th column without under/overflow. */
		nrm = 0.0;
		for (i = k; i < m; i = i + 1.0) {
			nrm = Hypothenuse(nrm, Element(QR, i, k));
		}

		if (nrm != 0.0) {
			/* Form k-th Householder vector. */
			if (Element(QR, k, k) < 0.0) {
				nrm = -nrm;
			}
			for (i = k; i < m; i = i + 1.0) {
				QR->r->at(i)->c->at(k) = Element(QR, i, k) / nrm;
			}
			QR->r->at(k)->c->at(k) = Element(QR, k, k) + 1.0;

			/* Apply transformation to remaining columns. */
			for (j = k + 1.0; j < n; j = j + 1.0) {
				s = 0.0;
				for (i = k; i < m; i = i + 1.0) {
					s = s + Element(QR, i, k) * Element(QR, i, j);
				}
				s = -s / Element(QR, k, k);
				for (i = k; i < m; i = i + 1.0) {
					QR->r->at(i)->c->at(j) = Element(QR, i, j) + s * Element(QR, i, k);
				}
			}
		}
		Rdiag->at(k) = -nrm;
	}

	/* Compute R */
	R = CreateSquareMatrix(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		for (j = 0.0; j < n; j = j + 1.0) {
			if (i < j) {
				R->r->at(i)->c->at(j) = Element(QR, i, j);
			}
			else if (i == j) {
				R->r->at(i)->c->at(j) = Rdiag->at(i);
			}
			else {
				R->r->at(i)->c->at(j) = 0.0;
			}
		}
	}
	Assign(r, R);

	/* Compute Q */
	Q = CreateMatrix(m, n);
	for (k = n - 1.0; k >= 0.0; k = k - 1.0) {
		for (i = 0.0; i < m; i = i + 1.0) {
			Q->r->at(i)->c->at(k) = 0.0;
		}
		Q->r->at(k)->c->at(k) = 1.0;
		for (j = k; j < n; j = j + 1.0) {
			if (Element(QR, k, k) != 0.0) {
				s = 0.0;
				for (i = k; i < m; i = i + 1.0) {
					s = s + Element(QR, i, k) * Element(Q, i, j);
				}
				s = -s / Element(QR, k, k);
				for (i = k; i < m; i = i + 1.0) {
					Q->r->at(i)->c->at(j) = Element(Q, i, j) + s * Element(QR, i, k);
				}
			}
		}
	}
	Assign(q, Q);

	/* Adjust for positive R. */
	n = NumberOfRows(r);

	N = CreateIdentityMatrix(n);

	for (i = 0.0; i < n; i = i + 1.0) {
		e = Element(r, i, i);

		if (e < 0.0) {
			N->r->at(i)->c->at(i) = -1.0;
		}
	}

	ra = MultiplyToNew(N, r);
	Assign(r, ra);
	rq = MultiplyToNew(q, N);
	Assign(q, rq);

	FreeMatrix(ra);
	FreeMatrix(rq);
	FreeMatrix(Q);
	FreeMatrix(R);
	FreeMatrix(QR);
}
double Hypothenuse(double a, double b) {
	return sqrt(pow(a, 2.0) + pow(b, 2.0));
}
double Norm(Matrix* a) {
	double l, i, j, rows, cols;

	l = 0.0;

	rows = NumberOfRows(a);
	cols = NumberOfColumns(a);

	for (i = 0.0; i < rows; i = i + 1.0) {
		for (j = 0.0; j < cols; j = j + 1.0) {
			l = l + a->r->at(i)->c->at(j) * a->r->at(i)->c->at(j);
		}
	}
	l = sqrt(l);

	return l;
}
Matrix* ExtractSubMatrix(Matrix* M, double r1, double r2, double c1, double c2) {
	Matrix* A;
	double i, j;

	A = CreateMatrix(r2 - r1 + 1.0, c2 - c1 + 1.0);

	for (i = r1; i <= r2; i = i + 1.0) {
		for (j = c1; j <= c2; j = j + 1.0) {
			A->r->at(i - r1)->c->at(j - c1) = M->r->at(i)->c->at(j);
		}
	}

	return A;
}
bool QRAlgorithm(Matrix* M, Matrix* R, Matrix* A, Matrix* Q, double precision, double maxIterations) {
	double n, i, j, v;
	vector<double>* previous;
	double withinPrecision;
	bool done, previousSet;

	Assign(A, M);
	n = NumberOfRows(M);
	previous = new vector<double>(n);
	previousSet = false;

	done = false;
	for (i = 0.0; i < maxIterations && !done; i = i + 1.0) {
		QRDecomposition(A, Q, R);
		Multiply(A, R, Q);

		/* Check precision. */
		if (previousSet) {
			withinPrecision = 0.0;
			for (j = 0.0; j < n; j = j + 1.0) {
				v = previous->at(j) - Element(A, j, j);
				if (abs(v) < precision || v == 0.0) {
					withinPrecision = withinPrecision + 1.0;
				}
			}
			if (withinPrecision == n) {
				done = true;
			}
		}

		for (j = 0.0; j < n; j = j + 1.0) {
			previous->at(j) = Element(A, j, j);
		}
		previousSet = true;
	}

	return done;
}
bool InvertUpperTriangularMatrix(Matrix* A, Matrix* inverse) {
	double sum, i, j, k, n;
	bool success;

	inverse->r = CreateCopyOfMatrix(A)->r;
	n = NumberOfRows(inverse);
	success = true;

	for (i = n - 1.0; i >= 0.0 && success; i = i - 1.0) {
		if (Element(inverse, i, i) == 0.0) {
			success = false;
		}
		else {
			inverse->r->at(i)->c->at(i) = 1.0 / Element(inverse, i, i);
			for (j = i - 1.0; j >= 0.0 && success; j = j - 1.0) {
				sum = 0.0;
				for (k = i; k > j; k = k - 1.0) {
					sum = sum - Element(inverse, j, k) * Element(inverse, k, i);
				}
				if (Element(inverse, j, j) == 0.0) {
					success = false;
				}
				else {
					inverse->r->at(j)->c->at(i) = sum / Element(inverse, j, j);
				}
			}
		}
	}

	return success;
}
bool InvertLowerTriangularMatrix(Matrix* A, Matrix* inverse) {
	double sum, i, j, k, n;
	bool success;

	inverse->r = CreateCopyOfMatrix(A)->r;
	n = NumberOfRows(inverse);
	success = true;

	for (i = 0.0; i < n && success; i = i + 1.0) {
		if (Element(inverse, i, i) == 0.0) {
			success = false;
		}
		else {
			inverse->r->at(i)->c->at(i) = 1.0 / Element(inverse, i, i);
			for (j = i + 1.0; j < n; j = j + 1.0) {
				sum = 0.0;
				for (k = i; k < j && success; k = k + 1.0) {
					sum = sum - Element(inverse, j, k) * Element(inverse, k, i);
				}
				if (Element(inverse, j, j) == 0.0) {
					success = false;
				}
				else {
					inverse->r->at(j)->c->at(i) = sum / Element(inverse, j, j);
				}
			}
		}
	}

	return success;
}
bool ParseMatrixFromString(MatrixReference* aref, vector<wchar_t>* matrixString, StringReference* errorMessage) {
	bool success;
	vector<StringReference*>* lines;
	double rows, cols, i;
	vector<double>* row;
	vector<wchar_t>* replaced, * trimmed;

	replaced = ReplaceString(matrixString, toVector(L"\r"), toVector(L""));
	trimmed = Trim(replaced);
	lines = SplitByCharacter(trimmed, '\n');

	delete replaced;
	delete trimmed;

	success = true;

	rows = lines->size();
	if (rows == 0.0) {
		aref->matrix = CreateMatrix(0.0, 0.0);
	}
	else {
		row = nStringToNumberArray(lines->at(0)->string);
		cols = row->size();
		delete row;

		aref->matrix = CreateMatrix(rows, cols);

		for (i = 0.0; i < rows && success; i = i + 1.0) {
			delete aref->matrix->r->at(i)->c;
			aref->matrix->r->at(i)->c = nStringToNumberArray(lines->at(i)->string);

			if (aref->matrix->r->at(i)->c->size() != cols) {
				success = false;
				errorMessage->string = toVector(L"All rows must have the same number of columns.");
			}
		}
	}

	FreeStringReferenceArray(lines);

	return success;
}
void FreeStringReferenceArray(vector<StringReference*>* stringReferencesArray) {
	double i;
	for (i = 0.0; i < stringReferencesArray->size(); i = i + 1.0) {
		delete stringReferencesArray->at(i);
	}
	delete stringReferencesArray;
}
vector<wchar_t>* MatrixToString(Matrix* matrix, double digitsAfterPoint) {
	vector<wchar_t>* s1, * s2;
	double n, m, element;

	s1 = new vector<wchar_t>(0.0);

	for (n = 0.0; n < NumberOfRows(matrix); n = n + 1.0) {
		for (m = 0.0; m < NumberOfColumns(matrix); m = m + 1.0) {
			element = Element(matrix, n, m);
			element = RoundToDigits(element, digitsAfterPoint);
			s2 = AppendString(s1, nCreateStringDecimalFromNumber(element));
			delete s1;
			s1 = s2;
			if (m + 1.0 != NumberOfColumns(matrix)) {
				s2 = AppendString(s1, toVector(L", "));
				delete s1;
				s1 = s2;
			}
		}
		s2 = AppendString(s1, toVector(L"\n"));
		delete s1;
		s1 = s2;
	}

	return s1;
}
double RoundToDigits(double element, double digitsAfterPoint) {
	return mathRound(element * pow(10.0, digitsAfterPoint)) / pow(10.0, digitsAfterPoint);
}
vector<wchar_t>* MatrixArrayToString(vector<Matrix*>* matrices, double digitsAfterPoint) {
	vector<wchar_t>* s1, * s2;
	double i;

	s1 = new vector<wchar_t>(0.0);

	for (i = 0.0; i < matrices->size(); i = i + 1.0) {
		s2 = AppendString(s1, MatrixToString(matrices->at(i), digitsAfterPoint));
		delete s1;
		s1 = s2;

		s2 = AppendString(s1, toVector(L"\n"));
		delete s1;
		s1 = s2;
	}

	return s1;
}
void RoundMatrixElementsToDigits(Matrix* a, double digits) {
	double m, n;

	for (m = 0.0; m < NumberOfRows(a); m = m + 1.0) {
		for (n = 0.0; n < NumberOfColumns(a); n = n + 1.0) {
			a->r->at(m)->c->at(n) = RoundToDigits(Element(a, m, n), digits);
			if (a->r->at(m)->c->at(n) == -0.0) {
				a->r->at(m)->c->at(n) = 0.0;
			}
		}
	}
}
bool SingularValueDecomposition(Matrix* Ap, MatrixReference* URef, MatrixReference* SigmaRef, MatrixReference* VRef) {
	Matrix* A, * U, * V;
	double m, n, nu, nct, nrt, i, j, k, t, pp, iter, eps, tiny, kase, f, cs, sn, ks, size, orgm, orgn;
	double scale, sp, spm1, epm1, sk, ek, b, c, shift, g, p;
	bool done;
	vector<double>* s, * e, * work;

	/* Square matrix, adjust results correspondingly. */
	orgm = NumberOfRows(Ap);
	orgn = NumberOfColumns(Ap);
	size = fmax(orgm, orgn);
	A = CreateCopyOfMatrix(Ap);
	Resize(A, size, size);

	/* Initialize. */
	m = size;
	n = size;

	/* Compute */
	nu = fmin(m, n);
	s = new vector<double>(fmin(m + 1.0, n));
	U = CreateMatrix(m, nu);
	V = CreateSquareMatrix(n);
	e = new vector<double>(n);
	work = new vector<double>(m);

	/* Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e. */
	nct = fmin(m - 1.0, n);
	nrt = fmax(0.0, fmin(n - 2.0, m));
	for (k = 0.0; k < fmax(nct, nrt); k = k + 1.0) {
		if (k < nct) {

			/* Compute the transformation for the k-th column and place the k-th diagonal in s[k]. */
			/* Compute 2-norm of k-th column without under/overflow. */
			s->at(k) = 0.0;
			for (i = k; i < m; i = i + 1.0) {
				s->at(k) = Hypothenuse(s->at(k), Element(A, i, k));
			}
			if (s->at(k) != 0.0) {
				if (Element(A, k, k) < 0.0) {
					s->at(k) = -s->at(k);
				}
				for (i = k; i < m; i = i + 1.0) {
					A->r->at(i)->c->at(k) = Element(A, i, k) / s->at(k);
				}
				A->r->at(k)->c->at(k) = Element(A, k, k) + 1.0;
			}
			s->at(k) = -s->at(k);
		}
		for (j = k + 1.0; j < n; j = j + 1.0) {
			if ((k < nct) && (s->at(k) != 0.0)) {

				/* Apply the transformation. */
				t = 0.0;
				for (i = k; i < m; i = i + 1.0) {
					t = t + Element(A, i, k) * Element(A, i, j);
				}
				t = -t / Element(A, k, k);
				for (i = k; i < m; i = i + 1.0) {
					A->r->at(i)->c->at(j) = Element(A, i, j) + t * Element(A, i, k);
				}
			}

			/* Place the k-th row of A into e for the subsequent calculation of the row transformation. */
			e->at(j) = Element(A, k, j);
		}
		if (k < nct) {

			/* Place the transformation in U for subsequent back */
			/* multiplication. */
			for (i = k; i < m; i = i + 1.0) {
				U->r->at(i)->c->at(k) = Element(A, i, k);
			}
		}
		if (k < nrt) {
			/* Compute the k-th row transformation and place the k-th super-diagonal in e[k]. */
			/* Compute 2-norm without under/overflow. */
			e->at(k) = 0.0;
			for (i = k + 1.0; i < n; i = i + 1.0) {
				e->at(k) = Hypothenuse(e->at(k), e->at(i));
			}
			if (e->at(k) != 0.0) {
				if (e->at(k + 1.0) < 0.0) {
					e->at(k) = -e->at(k);
				}
				for (i = k + 1.0; i < n; i = i + 1.0) {
					e->at(i) = e->at(i) / e->at(k);
				}
				e->at(k + 1.0) = e->at(k + 1.0) + 1.0;
			}
			e->at(k) = -e->at(k);
			if ((k + 1.0 < m) && (e->at(k) != 0.0)) {

				/* Apply the transformation. */
				for (i = k + 1.0; i < m; i = i + 1.0) {
					work->at(i) = 0.0;
				}
				for (j = k + 1.0; j < n; j = j + 1.0) {
					for (i = k + 1.0; i < m; i = i + 1.0) {
						work->at(i) = work->at(i) + e->at(j) * Element(A, i, j);
					}
				}
				for (j = k + 1.0; j < n; j = j + 1.0) {
					t = -e->at(j) / e->at(k + 1.0);
					for (i = k + 1.0; i < m; i = i + 1.0) {
						A->r->at(i)->c->at(j) = Element(A, i, j) + t * work->at(i);
					}
				}
			}

			/* Place the transformation in V for subsequent back multiplication. */
			for (i = k + 1.0; i < n; i = i + 1.0) {
				V->r->at(i)->c->at(k) = e->at(i);
			}
		}
	}

	/* Set up the final bidiagonal matrix or order p. */
	p = fmin(n, m + 1.0);
	if (nct < n) {
		s->at(nct) = Element(A, nct, nct);
	}
	if (m < p) {
		s->at(p - 1.0) = 0.0;
	}
	if (nrt + 1.0 < p) {
		e->at(nrt) = Element(A, nrt, p - 1.0);
	}
	e->at(p - 1.0) = 0.0;

	/* Generate U. */
	for (j = nct; j < nu; j = j + 1.0) {
		for (i = 0.0; i < m; i = i + 1.0) {
			U->r->at(i)->c->at(j) = 0.0;
		}
		U->r->at(j)->c->at(j) = 1.0;
	}
	for (k = nct - 1.0; k >= 0.0; k = k - 1.0) {
		if (s->at(k) != 0.0) {
			for (j = k + 1.0; j < nu; j = j + 1.0) {
				t = 0.0;
				for (i = k; i < m; i = i + 1.0) {
					t = t + Element(U, i, k) * Element(U, i, j);
				}
				t = -t / Element(U, k, k);
				for (i = k; i < m; i = i + 1.0) {
					U->r->at(i)->c->at(j) = Element(U, i, j) + t * Element(U, i, k);
				}
			}
			for (i = k; i < m; i = i + 1.0) {
				U->r->at(i)->c->at(k) = -Element(U, i, k);
			}
			U->r->at(k)->c->at(k) = 1.0 + Element(U, k, k);
			for (i = 0.0; i < k - 1.0; i = i + 1.0) {
				U->r->at(i)->c->at(k) = 0.0;
			}
		}
		else {
			for (i = 0.0; i < m; i = i + 1.0) {
				U->r->at(i)->c->at(k) = 0.0;
			}
			U->r->at(k)->c->at(k) = 1.0;
		}
	}

	/* Generate V. */
	for (k = n - 1.0; k >= 0.0; k = k - 1.0) {
		if ((k < nrt) && (e->at(k) != 0.0)) {
			for (j = k + 1.0; j < nu; j = j + 1.0) {
				t = 0.0;
				for (i = k + 1.0; i < n; i = i + 1.0) {
					t = t + Element(V, i, k) * Element(V, i, j);
				}
				t = -t / Element(V, k + 1.0, k);
				for (i = k + 1.0; i < n; i = i + 1.0) {
					V->r->at(i)->c->at(j) = Element(V, i, j) + t * Element(V, i, k);
				}
			}
		}
		for (i = 0.0; i < n; i = i + 1.0) {
			V->r->at(i)->c->at(k) = 0.0;
		}
		V->r->at(k)->c->at(k) = 1.0;
	}

	/* Main iteration loop for the singular values. */
	pp = p - 1.0;
	iter = 0.0;
	eps = pow(2.0, -52.0);
	tiny = pow(2.0, -966.0);
	for (; p > 0.0; ) {
		/* Here is where a test for too many iterations would go. */
		/* This section of the program inspects for negligible elements in the s and e arrays. */
		/* On completion the variables kase and k are set as follows. */
		/* kase = 1, if s(p) and e[k-1] are negligible and k<p */
		/* kase = 2, if s(k) is negligible and k<p */
		/* kase = 3, if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step). */
		/* kase = 4, if e(p-1) is negligible (convergence). */
		done = false;
		for (k = p - 2.0; k > -1.0 && !done; ) {
			if (abs(e->at(k)) <= tiny + eps * (abs(s->at(k)) + abs(s->at(k + 1.0)))) {
				e->at(k) = 0.0;
				done = true;
			}
			else {
				k = k - 1.0;
			}
		}
		if (k == p - 2.0) {
			kase = 4.0;
		}
		else {
			done = false;
			for (ks = p - 1.0; ks > k && !done; ) {
				if (ks != p) {
					t = abs(e->at(ks));
				}
				else {
					t = 0.0;
				}

				if (ks != k + 1.0) {
					t = t + abs(e->at(ks - 1.0));
				}

				if (abs(s->at(ks)) <= tiny + eps * t) {
					s->at(ks) = 0.0;
					done = true;
				}
				else {
					ks = ks - 1.0;
				}
			}
			if (ks == k) {
				kase = 3.0;
			}
			else if (ks == p - 1.0) {
				kase = 1.0;
			}
			else {
				kase = 2.0;
				k = ks;
			}
		}
		k = k + 1.0;

		/* Perform the task indicated by kase. */
		if (kase == 1.0) {
			/* Deflate negligible s(p). */
			f = e->at(p - 2.0);
			e->at(p - 2.0) = 0.0;
			for (j = p - 2.0; j >= k; j = j - 1.0) {
				t = Hypothenuse(s->at(j), f);
				cs = s->at(j) / t;
				sn = f / t;
				s->at(j) = t;
				if (j != k) {
					f = -sn * e->at(j - 1.0);
					e->at(j - 1.0) = cs * e->at(j - 1.0);
				}

				for (i = 0.0; i < n; i = i + 1.0) {
					t = cs * Element(V, i, j) + sn * Element(V, i, p - 1.0);
					V->r->at(i)->c->at(p - 1.0) = -sn * Element(V, i, j) + cs * Element(V, i, p - 1.0);
					V->r->at(i)->c->at(j) = t;
				}
			}
		}
		else if (kase == 2.0) {
			/* Split at negligible s(k). */
			f = e->at(k - 1.0);
			e->at(k - 1.0) = 0.0;
			for (j = k; j < p; j = j + 1.0) {
				t = Hypothenuse(s->at(j), f);
				cs = s->at(j) / t;
				sn = f / t;
				s->at(j) = t;
				f = -sn * e->at(j);
				e->at(j) = cs * e->at(j);

				for (i = 0.0; i < m; i = i + 1.0) {
					t = cs * Element(U, i, j) + sn * Element(U, i, k + 1.0);
					U->r->at(i)->c->at(k - 1.0) = -sn * Element(U, i, j) + cs * Element(U, i, k + 1.0);
					U->r->at(i)->c->at(j) = t;
				}
			}
		}
		else if (kase == 3.0) {
			/* Perform one qr step. */
			/* Calculate the shift. */
			scale = fmax(fmax(fmax(fmax(abs(s->at(p - 1.0)), abs(s->at(p - 2.0))), abs(e->at(p - 2.0))), abs(s->at(k))), abs(e->at(k)));
			sp = s->at(p - 1.0) / scale;
			spm1 = s->at(p - 2.0) / scale;
			epm1 = e->at(p - 2.0) / scale;
			sk = s->at(k) / scale;
			ek = e->at(k) / scale;
			b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
			c = (sp * epm1) * (sp * epm1);
			shift = 0.0;
			if ((b != 0.0) || (c != 0.0)) {
				shift = sqrt(b * b + c);
				if (b < 0.0) {
					shift = -shift;
				}
				shift = c / (b + shift);
			}
			f = (sk + sp) * (sk - sp) + shift;
			g = sk * ek;

			/* Chase zeros. */
			for (j = k; j < p - 1.0; j = j + 1.0) {
				t = Hypothenuse(f, g);
				cs = f / t;
				sn = g / t;
				if (j != k) {
					e->at(j - 1.0) = t;
				}
				f = cs * s->at(j) + sn * e->at(j);
				e->at(j) = cs * e->at(j) - sn * s->at(j);
				g = sn * s->at(j + 1.0);
				s->at(j + 1.0) = cs * s->at(j + 1.0);
				for (i = 0.0; i < n; i = i + 1.0) {
					t = cs * Element(V, i, j) + sn * Element(V, i, j + 1.0);
					V->r->at(i)->c->at(j + 1.0) = -sn * Element(V, i, j) + cs * Element(V, i, j + 1.0);
					V->r->at(i)->c->at(j) = t;
				}
				t = Hypothenuse(f, g);
				cs = f / t;
				sn = g / t;
				s->at(j) = t;
				f = cs * e->at(j) + sn * s->at(j + 1.0);
				s->at(j + 1.0) = -sn * e->at(j) + cs * s->at(j + 1.0);
				g = sn * e->at(j + 1.0);
				e->at(j + 1.0) = cs * e->at(j + 1.0);
				if (j < m - 1.0) {
					for (i = 0.0; i < m; i = i + 1.0) {
						t = cs * Element(U, i, j) + sn * Element(U, i, j + 1.0);
						U->r->at(i)->c->at(j + 1.0) = -sn * Element(U, i, j) + cs * Element(U, i, j + 1.0);
						U->r->at(i)->c->at(j) = t;
					}
				}
			}
			e->at(p - 2.0) = f;
			iter = iter + 1.0;
		}
		else if (kase == 4.0) {
			/* Make the singular values positive. */
			if (s->at(k) <= 0.0) {
				if (s->at(k) < 0.0) {
					s->at(k) = -s->at(k);
				}
				else {
					s->at(k) = 0.0;
				}

				for (i = 0.0; i <= pp; i = i + 1.0) {
					V->r->at(i)->c->at(k) = -Element(V, i, k);
				}
			}

			/* Order the singular values. */
			for (; k < pp && s->at(k) < s->at(k + 1.0); ) {
				t = s->at(k);
				s->at(k) = s->at(k + 1.0);
				s->at(k + 1.0) = t;
				if (k < n - 1.0) {
					for (i = 0.0; i < n; i = i + 1.0) {
						t = Element(V, i, k + 1.0);
						V->r->at(i)->c->at(k + 1.0) = Element(V, i, k);
						V->r->at(i)->c->at(k) = t;
					}
				}
				if (k < m - 1.0) {
					for (i = 0.0; i < m; i = i + 1.0) {
						t = Element(U, i, k + 1.0);
						U->r->at(i)->c->at(k + 1.0) = Element(U, i, k);
						U->r->at(i)->c->at(k) = t;
					}
				}
				k = k + 1.0;
			}
			iter = 0.0;
			p = p - 1.0;
		}
	}

	Resize(U, orgm, orgm);
	Resize(V, orgn, orgn);

	URef->matrix = U;
	VRef->matrix = V;
	SigmaRef->matrix = CreateMatrix(orgm, orgn);
	for (i = 0.0; i < fmin(orgm, orgn); i = i + 1.0) {
		SigmaRef->matrix->r->at(i)->c->at(i) = s->at(i);
	}

	return true;
}
ComplexMatrix* CreateComplexMatrix(double rows, double cols) {
	double m, n;
	ComplexMatrix* matrix;

	matrix = new ComplexMatrix();
	matrix->r = new vector<ComplexMatrixRow*>(rows);
	for (m = 0.0; m < rows; m = m + 1.0) {
		matrix->r->at(m) = new ComplexMatrixRow();
		matrix->r->at(m)->c = new vector<cComplexNumber*>(cols);
		for (n = 0.0; n < cols; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = cCreateComplexNumber(0.0, 0.0);
		}
	}

	return matrix;
}
ComplexMatrix* CreateComplexMatrixFromMatrix(Matrix* a) {
	double m, n, rows, cols;
	ComplexMatrix* matrix;

	rows = NumberOfRows(a);
	cols = NumberOfColumns(a);

	matrix = new ComplexMatrix();
	matrix->r = new vector<ComplexMatrixRow*>(rows);
	for (m = 0.0; m < rows; m = m + 1.0) {
		matrix->r->at(m) = new ComplexMatrixRow();
		matrix->r->at(m)->c = new vector<cComplexNumber*>(cols);
		for (n = 0.0; n < cols; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = cCreateComplexNumber(a->r->at(m)->c->at(n), 0.0);
		}
	}

	return matrix;
}
Matrix* CreateReMatrixFromComplexMatrix(ComplexMatrix* a) {
	double m, n, rows, cols;
	Matrix* matrix;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfColumnsComplex(a);

	matrix = new Matrix();
	matrix->r = new vector<MatrixRow*>(rows);
	for (m = 0.0; m < rows; m = m + 1.0) {
		matrix->r->at(m) = new MatrixRow();
		matrix->r->at(m)->c = new vector<double>(cols);
		for (n = 0.0; n < cols; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = IndexComplex(a, m, n)->re;
		}
	}

	return matrix;
}
Matrix* CreateImMatrixFromComplexMatrix(ComplexMatrix* a) {
	double m, n, rows, cols;
	Matrix* matrix;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfColumnsComplex(a);

	matrix = new Matrix();
	matrix->r = new vector<MatrixRow*>(rows);
	for (m = 0.0; m < rows; m = m + 1.0) {
		matrix->r->at(m) = new MatrixRow();
		matrix->r->at(m)->c = new vector<double>(cols);
		for (n = 0.0; n < cols; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = IndexComplex(a, m, n)->im;
		}
	}

	return matrix;
}
double NumberOfRowsComplex(ComplexMatrix* A) {
	return A->r->size();
}
double NumberOfColumnsComplex(ComplexMatrix* A) {
	return A->r->at(0)->c->size();
}
cComplexNumber* IndexComplex(ComplexMatrix* a, double m, double n) {
	return a->r->at(m)->c->at(n);
}
void AddComplex(ComplexMatrix* a, ComplexMatrix* b) {
	double m, n;
	double d;

	d = NumberOfRowsComplex(a);

	for (m = 0.0; m < d; m = m + 1.0) {
		for (n = 0.0; n < d; n = n + 1.0) {
			cAdd(IndexComplex(a, m, n), IndexComplex(b, m, n));
		}
	}
}
void SubtractComplex(ComplexMatrix* a, ComplexMatrix* b) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(a);
	c = NumberOfColumnsComplex(a);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cSub(IndexComplex(a, m, n), IndexComplex(b, m, n));
		}
	}
}
ComplexMatrix* SubtractComplexToNew(ComplexMatrix* a, ComplexMatrix* b) {
	ComplexMatrix* X;

	X = CreateCopyOfComplexMatrix(a);
	SubtractComplex(X, b);

	return X;
}
void MultiplyComplex(ComplexMatrix* x, ComplexMatrix* a, ComplexMatrix* b) {
	double m, n;
	double rows, cols, d;
	double i;
	cComplexNumber* s, * t;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfColumnsComplex(b);
	d = NumberOfColumnsComplex(a);
	t = cCreateComplexNumber(0.0, 0.0);

	for (m = 0.0; m < rows; m = m + 1.0) {
		for (n = 0.0; n < cols; n = n + 1.0) {
			s = cCreateComplexNumber(0.0, 0.0);

			for (i = 0.0; i < d; i = i + 1.0) {
				cAssignComplex(t, s);
				cAssignComplex(s, IndexComplex(a, m, i));
				cMul(s, IndexComplex(b, i, n));
				cAdd(s, t);
			}

			x->r->at(m)->c->at(n) = s;
		}
	}
}
ComplexMatrix* MultiplyComplexToNew(ComplexMatrix* a, ComplexMatrix* b) {
	double rows, cols;
	ComplexMatrix* x;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfColumnsComplex(b);
	x = CreateComplexMatrix(rows, cols);
	MultiplyComplex(x, a, b);

	return x;
}
void Conjugate(ComplexMatrix* a) {
	double m, n;
	double rows, cols;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfRowsComplex(a);

	for (m = 0.0; m < rows; m = m + 1.0) {
		for (n = 0.0; n < cols; n = n + 1.0) {
			cConjugate(IndexComplex(a, m, n));
		}
	}
}
void AssignComplexMatrix(ComplexMatrix* A, ComplexMatrix* B) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(A);
	c = NumberOfColumnsComplex(A);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cAssignComplex(IndexComplex(A, m, n), IndexComplex(B, m, n));
		}
	}
}
void ScalarMultiplyComplex(ComplexMatrix* A, cComplexNumber* b) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(A);
	c = NumberOfColumnsComplex(A);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cMul(IndexComplex(A, m, n), b);
		}
	}
}
ComplexMatrix* ScalarMultiplyComplexToNew(ComplexMatrix* A, cComplexNumber* b) {
	ComplexMatrix* matrix;

	matrix = CreateCopyOfComplexMatrix(A);
	ScalarMultiplyComplex(matrix, b);

	return matrix;
}
void ScalarDivideComplex(ComplexMatrix* A, cComplexNumber* b) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(A);
	c = NumberOfColumnsComplex(A);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cDiv(IndexComplex(A, m, n), b);
		}
	}
}
void ElementWisePowerComplex(ComplexMatrix* A, double p) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(A);
	c = NumberOfColumnsComplex(A);

	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cPower(IndexComplex(A, m, n), p);
		}
	}
}
ComplexMatrix* CreateComplexIdentityMatrix(double d) {
	double m;
	ComplexMatrix* matrix;

	matrix = CreateSquareComplexMatrix(d);
	FillComplex(matrix, 0.0, 0.0);

	for (m = 0.0; m < d; m = m + 1.0) {
		IndexComplex(matrix, m, m)->re = 1.0;
	}

	return matrix;
}
ComplexMatrix* CreateSquareComplexMatrix(double d) {
	double m, n;
	ComplexMatrix* matrix;

	matrix = new ComplexMatrix();
	matrix->r = new vector<ComplexMatrixRow*>(d);
	for (m = 0.0; m < d; m = m + 1.0) {
		matrix->r->at(m) = new ComplexMatrixRow();
		matrix->r->at(m)->c = new vector<cComplexNumber*>(d);
		for (n = 0.0; n < d; n = n + 1.0) {
			matrix->r->at(m)->c->at(n) = cCreateComplexNumber(0.0, 0.0);
		}
	}

	return matrix;
}
void ClearComplex(ComplexMatrix* a) {
	FillComplex(a, 0.0, 0.0);
}
void FillComplex(ComplexMatrix* a, double re, double im) {
	double m, n;

	for (m = 0.0; m < NumberOfRowsComplex(a); m = m + 1.0) {
		for (n = 0.0; n < NumberOfColumnsComplex(a); n = n + 1.0) {
			IndexComplex(a, m, n)->re = re;
			IndexComplex(a, m, n)->im = im;
		}
	}
}
cComplexNumber* TraceComplex(ComplexMatrix* a) {
	double m;
	double d;
	cComplexNumber* tr;

	tr = cCreateComplexNumber(0.0, 0.0);

	d = a->r->size();
	for (m = 0.0; m < d; m = m + 1.0) {
		cAdd(tr, IndexComplex(a, m, m));
	}

	return tr;
}
void CofactorOfComplexMatrix(ComplexMatrix* mat, ComplexMatrix* temp, double p, double q, double n) {
	double i, j;
	double row, col;

	i = 0.0;
	j = 0.0;

	for (row = 0.0; row < n; row = row + 1.0) {
		for (col = 0.0; col < n; col = col + 1.0) {
			if (row != p && col != q) {
				cAssignComplex(IndexComplex(temp, i, j), IndexComplex(mat, row, col));
				j = j + 1.0;

				if (j == n - 1.0) {
					j = 0.0;
					i = i + 1.0;
				}
			}
		}
	}
}
cComplexNumber* DeterminantOfComplexSubmatrix(ComplexMatrix* mat, double n) {
	double f, sign;
	cComplexNumber* D, * t;
	ComplexMatrix* temp;

	D = cCreateComplexNumber(0.0, 0.0);
	t = cCreateComplexNumber(0.0, 0.0);

	if (n == 1.0) {
		D = mat->r->at(0)->c->at(0);
	}
	else {
		temp = CreateSquareComplexMatrix(n);

		sign = 1.0;

		for (f = 0.0; f < n; f = f + 1.0) {
			CofactorOfComplexMatrix(mat, temp, 0.0, f, n);
			cAssignComplexByValues(t, sign, 0.0);
			cMul(t, IndexComplex(mat, 0.0, f));
			cMul(t, DeterminantOfComplexSubmatrix(temp, n - 1.0));
			cAdd(D, t);
			sign = -sign;
		}

		DeleteComplexMatrix(temp);
	}

	return D;
}
void DeleteComplexMatrix(ComplexMatrix* X) {
	double m, n, rows, cols;

	rows = NumberOfRowsComplex(X);
	cols = NumberOfColumnsComplex(X);
	for (m = 0.0; m < rows; m = m + 1.0) {
		for (n = 0.0; n < cols; n = n + 1.0) {
			delete X->r->at(m)->c->at(n);
		}
		delete X->r->at(m)->c;
		delete X->r->at(m);
	}

	delete X->r;
	delete X;
}
cComplexNumber* DeterminantComplex(ComplexMatrix* m) {
	double n;
	cComplexNumber* D;

	n = NumberOfRowsComplex(m);
	D = DeterminantOfComplexSubmatrix(m, n);

	return D;
}
void AdjointComplex(ComplexMatrix* A, ComplexMatrix* adj) {
	double n;
	ComplexMatrix* cofactors;
	double i, j;
	cComplexNumber* t, * sign;

	n = A->r->size();
	t = cCreateComplexNumber(0.0, 0.0);
	sign = cCreateComplexNumber(0.0, 0.0);

	if (n == 1.0) {
		cAssignComplexByValues(IndexComplex(adj, 0.0, 0.0), 1.0, 0.0);
	}
	else {
		cofactors = CreateSquareComplexMatrix(n);

		for (i = 0.0; i < n; i = i + 1.0) {
			for (j = 0.0; j < n; j = j + 1.0) {
				CofactorOfComplexMatrix(A, cofactors, i, j, n);

				if (fmod(i + j, 2.0) == 0.0) {
					cAssignComplexByValues(sign, 1.0, 0.0);
				}
				else {
					cAssignComplexByValues(sign, -1.0, 0.0);
				}

				cAssignComplex(t, sign);
				cMul(t, DeterminantOfComplexSubmatrix(cofactors, n - 1.0));
				cAssignComplex(IndexComplex(adj, j, i), t);
			}
		}

		DeleteComplexMatrix(cofactors);
	}
}
bool InverseComplex(ComplexMatrix* A, ComplexMatrix* inverseResult) {
	bool success;
	ComplexMatrix* adj;
	double n, i, j;
	cComplexNumber* det, * t;

	t = cCreateComplexNumber(0.0, 0.0);

	if (NumberOfColumnsComplex(A) == NumberOfRowsComplex(A)) {
		n = NumberOfColumnsComplex(A);

		det = DeterminantComplex(A);
		if (det->re != 0.0 || det->im != 0.0) {
			adj = CreateSquareComplexMatrix(n);
			AdjointComplex(A, adj);

			for (i = 0.0; i < n; i = i + 1.0) {
				for (j = 0.0; j < n; j = j + 1.0) {
					cAssignComplex(t, IndexComplex(adj, i, j));
					cDiv(t, det);
					cAssignComplex(IndexComplex(inverseResult, i, j), t);
				}
			}

			success = true;
			DeleteComplexMatrix(adj);
		}
		else {
			success = false;
		}
	}
	else {
		success = false;
	}

	return success;
}
bool ComplexMatrixEqualsEpsilon(ComplexMatrix* b, ComplexMatrix* f, double epsilon) {
	double x, y, columns, rows;
	bool equals;

	equals = true;

	if (NumberOfRowsComplex(b) == NumberOfRowsComplex(f) && NumberOfColumnsComplex(b) == NumberOfColumnsComplex(f)) {
		columns = NumberOfColumnsComplex(b);
		rows = NumberOfRowsComplex(b);

		for (x = 0.0; x < rows; x = x + 1.0) {
			for (y = 0.0; y < columns; y = y + 1.0) {
				equals = equals && cEpsilonCompareComplex(IndexComplex(b, x, y), IndexComplex(f, x, y), epsilon);
			}
		}
	}
	else {
		equals = false;
	}

	return equals;
}
ComplexMatrix* MinorComplex(ComplexMatrix* x, double row, double column) {
	ComplexMatrix* minor;
	double cols, rows, i, j, m, n;

	rows = NumberOfRowsComplex(x) - 1.0;
	cols = NumberOfColumnsComplex(x) - 1.0;

	minor = CreateComplexMatrix(rows, cols);

	for (i = 0.0; i < rows; i = i + 1.0) {
		if (i < row) {
			m = i;
		}
		else {
			m = i + 1.0;
		}

		for (j = 0.0; i != row && j < cols; j = j + 1.0) {
			if (j != column) {

				if (j < column) {
					n = j;
				}
				else {
					n = j + 1.0;
				}

				cAssignComplex(IndexComplex(minor, m, n), IndexComplex(x, i, j));
			}
		}
	}

	return minor;
}
void AssignComplex(ComplexMatrix* A, ComplexMatrix* B) {
	double m, n;
	double r, c;

	r = NumberOfRowsComplex(A);
	c = NumberOfColumnsComplex(A);
	for (m = 0.0; m < r; m = m + 1.0) {
		for (n = 0.0; n < c; n = n + 1.0) {
			cAssignComplex(IndexComplex(A, m, n), IndexComplex(B, m, n));
		}
	}
}
ComplexMatrix* CreateCopyOfComplexMatrix(ComplexMatrix* A) {
	ComplexMatrix* X;

	X = CreateComplexMatrix(NumberOfRowsComplex(A), NumberOfColumnsComplex(A));
	AssignComplex(X, A);

	return X;
}
bool TransposeComplex(ComplexMatrix* a) {
	double m, n;
	double rows;
	bool square;
	cComplexNumber* tmp;

	tmp = cCreateComplexNumber(0.0, 0.0);

	square = IsSquareComplexMatrix(a);
	if (square) {
		rows = NumberOfColumnsComplex(a);

		for (m = 0.0; m < rows; m = m + 1.0) {
			for (n = 0.0; n < m; n = n + 1.0) {
				cAssignComplex(tmp, IndexComplex(a, n, m));
				cAssignComplex(IndexComplex(a, n, m), IndexComplex(a, m, n));
				cAssignComplex(IndexComplex(a, m, n), tmp);
			}
		}
	}

	return square;
}
bool ConjugateTransposeComplex(ComplexMatrix* a) {
	double m, n;
	double rows;
	bool square;
	cComplexNumber* tmp;

	tmp = cCreateComplexNumber(0.0, 0.0);

	square = IsSquareComplexMatrix(a);
	if (square) {
		rows = NumberOfColumnsComplex(a);

		for (m = 0.0; m < rows; m = m + 1.0) {
			for (n = 0.0; n < m; n = n + 1.0) {
				cAssignComplex(tmp, IndexComplex(a, n, m));
				cAssignComplex(IndexComplex(a, n, m), IndexComplex(a, m, n));
				cAssignComplex(IndexComplex(a, m, n), tmp);
				cConjugate(IndexComplex(a, m, n));
			}
		}
	}

	return square;
}
bool IsSquareComplexMatrix(ComplexMatrix* A) {
	bool is;

	if (NumberOfRowsComplex(A) == NumberOfColumnsComplex(A)) {
		is = true;
	}
	else {
		is = false;
	}

	return is;
}
void TransposeComplexAssign(ComplexMatrix* t, ComplexMatrix* a) {
	double m, n;
	double rows, cols;

	cols = NumberOfRowsComplex(a);
	rows = NumberOfColumnsComplex(a);

	for (m = 0.0; m < cols; m = m + 1.0) {
		for (n = 0.0; n < rows; n = n + 1.0) {
			cAssignComplex(IndexComplex(t, n, m), IndexComplex(a, m, n));
		}
	}
}
ComplexMatrix* TransposeComplexToNew(ComplexMatrix* a) {
	double m, n;
	double rows, cols;
	ComplexMatrix* c;

	cols = NumberOfRowsComplex(a);
	rows = NumberOfColumnsComplex(a);

	c = CreateComplexMatrix(rows, cols);

	for (m = 0.0; m < cols; m = m + 1.0) {
		for (n = 0.0; n < rows; n = n + 1.0) {
			cAssignComplex(IndexComplex(c, n, m), IndexComplex(a, m, n));
		}
	}

	return c;
}
ComplexMatrix* ExtractComplexSubMatrix(ComplexMatrix* M, double r1, double r2, double c1, double c2) {
	ComplexMatrix* A;
	double i, j;

	A = CreateComplexMatrix(r2 - r1 + 1.0, c2 - c1 + 1.0);

	for (i = r1; i <= r2; i = i + 1.0) {
		for (j = c1; j <= c2; j = j + 1.0) {
			cAssignComplex(IndexComplex(A, i - r1, j - c1), IndexComplex(M, i, j));
		}
	}

	return A;
}
double NormComplex(ComplexMatrix* a) {
	double l, i, j, rows, cols;
	cComplexNumber* cComplexNumber;

	l = 0.0;

	rows = NumberOfRowsComplex(a);
	cols = NumberOfColumnsComplex(a);

	for (i = 0.0; i < rows; i = i + 1.0) {
		for (j = 0.0; j < cols; j = j + 1.0) {
			cComplexNumber = IndexComplex(a, i, j);
			l = l + pow(cComplexNumber->re, 2.0) + pow(cComplexNumber->im, 2.0);
		}
	}
	l = sqrt(l);

	return l;
}
void ComplexCharacteristicPolynomial(ComplexMatrix* A, pComplexPolynomial* p) {
	ComplexMatrix* cp;
	cComplexNumber* determinant;

	cp = CreateSquareComplexMatrix(NumberOfRowsComplex(A));
	determinant = new cComplexNumber();

	ComplexCharacteristicPolynomialWithInverse(A, cp, p, determinant);

	DeleteComplexMatrix(cp);
}
void ComplexCharacteristicPolynomialWithInverse(ComplexMatrix* A, ComplexMatrix* AInverse, pComplexPolynomial* p, cComplexNumber* determinant) {
	FaddeevLeVerrierAlgorithmComplex(A, AInverse, p, determinant);
}
void FaddeevLeVerrierAlgorithmComplex(ComplexMatrix* A, ComplexMatrix* AInverse, pComplexPolynomial* p, cComplexNumber* determinant) {
	ComplexMatrix* Mk, * Mkm1, * t1, * Id;
	cComplexNumber* t, * t2, * n1;
	double i, n, k;

	n1 = cCreateComplexNumber(-1.0, 0.0);
	n = NumberOfRowsComplex(A);
	p->cs = new vector<cComplexNumber*>(n + 1.0);
	for (i = 0.0; i < n + 1.0; i = i + 1.0) {
		p->cs->at(i) = new cComplexNumber();
	}
	cAssignComplexByValues(p->cs->at(n), 1.0, 0.0);
	Mkm1 = CreateSquareComplexMatrix(n);
	FillComplex(Mkm1, 0.0, 0.0);
	Id = CreateComplexIdentityMatrix(n);
	Mk = CreateSquareComplexMatrix(n);
	t1 = CreateSquareComplexMatrix(n);

	for (k = 1.0; k <= n; k = k + 1.0) {
		MultiplyComplex(Mk, A, Mkm1);
		AssignComplex(t1, Id);
		ScalarMultiplyComplex(t1, p->cs->at(n - k + 1.0));
		AddComplex(Mk, t1);

		MultiplyComplex(t1, A, Mk);
		t = TraceComplex(t1);
		t2 = cCreateComplexNumber(-1.0 / k, 0.0);
		cMul(t, t2);
		cAssignComplex(p->cs->at(n - k), t);

		/* done */
		AssignComplex(Mkm1, Mk);

		if (k == n) {
			AssignComplex(AInverse, Mk);
			cAssignComplex(t, p->cs->at(0));
			cMul(t, n1);
			cAssignComplex(determinant, t);
			if (t->re == 0.0 && t->im == 0.0) {
			}
			else {
				ScalarDivideComplex(AInverse, t);
			}
		}
	}

	DeleteComplexMatrix(Mkm1);
	DeleteComplexMatrix(Id);
	DeleteComplexMatrix(Mk);
	DeleteComplexMatrix(t1);
}
bool EigenvaluesComplex(ComplexMatrix* A, cComplexNumberArrayReference* eigenValuesReference) {
	ComplexMatrixArrayReference* eigenVectorsReference;
	bool success;
	double i;

	eigenVectorsReference = new ComplexMatrixArrayReference();
	success = EigenpairsComplex(A, eigenValuesReference, eigenVectorsReference);
	if (success) {
		for (i = 0.0; i < eigenVectorsReference->matrices->size(); i = i + 1.0) {
			DeleteComplexMatrix(eigenVectorsReference->matrices->at(i));
		}
		delete eigenVectorsReference->matrices;
		delete eigenVectorsReference;
	}

	return success;
}
bool EigenvectorsComplex(ComplexMatrix* A, ComplexMatrixArrayReference* eigenVectorsReference) {
	cComplexNumberArrayReference* evsReference;
	bool success;
	double i;

	evsReference = new cComplexNumberArrayReference();
	success = EigenpairsComplex(A, evsReference, eigenVectorsReference);
	if (success) {
		for (i = 0.0; i < evsReference->complexNumbers->size(); i = i + 1.0) {
			delete evsReference->complexNumbers->at(i);
		}
		delete evsReference->complexNumbers;
		delete evsReference;
	}

	return success;
}
bool InversePowerMethodComplex(ComplexMatrix* A, cComplexNumber* eigenvalue, double maxIterations, cComplexNumberArrayReference* eigenvector) {
	ComplexMatrix* t1, * t2, * t3, * t4, * b;
	double n, i, c;
	bool isSingular;
	cComplexNumber* c101, * k, * cc;

	n = NumberOfRowsComplex(A);

	t2 = CreateComplexIdentityMatrix(n);
	ScalarMultiplyComplex(t2, eigenvalue);
	t3 = SubtractComplexToNew(A, t2);
	t4 = CreateSquareComplexMatrix(n);
	isSingular = !InverseComplex(t3, t4);
	cc = cCreateComplexNumber(0.0, 0.0);
	t1 = CreateComplexMatrix(n, 1.0);

	if (isSingular) {
		DeleteComplexMatrix(t2);
		DeleteComplexMatrix(t3);
		DeleteComplexMatrix(t4);

		c101 = cCreateComplexNumber(1.01, 0.0);
		/* Try again with more erroneous eigenvalue estimate. */
		t2 = CreateComplexIdentityMatrix(n);
		k = cMulToNew(eigenvalue, c101);
		ScalarMultiplyComplex(t2, k);
		t3 = SubtractComplexToNew(A, t2);
		t4 = CreateSquareComplexMatrix(n);
		isSingular = !InverseComplex(t3, t4);
		delete c101;
	}

	if (!isSingular) {
		b = CreateComplexMatrix(n, 1.0);

		for (i = 0.0; i < n; i = i + 1.0) {
			cAssignComplexByValues(b->r->at(i)->c->at(0), 1.0, 1.0);
		}

		for (i = 0.0; i < maxIterations; i = i + 1.0) {
			MultiplyComplex(t1, t4, b);
			c = NormComplex(t1);
			cAssignComplexByValues(cc, c, 0.0);
			ScalarDivideComplex(t1, cc);
			AssignComplex(b, t1);
		}

		eigenvector->complexNumbers = new vector<cComplexNumber*>(n);
		for (i = 0.0; i < n; i = i + 1.0) {
			eigenvector->complexNumbers->at(i) = b->r->at(i)->c->at(0);
		}
	}

	DeleteComplexMatrix(t1);
	DeleteComplexMatrix(t2);
	DeleteComplexMatrix(t3);
	DeleteComplexMatrix(t4);
	delete cc;

	return  !isSingular;
}
bool EigenpairsComplex(ComplexMatrix* M, cComplexNumberArrayReference* eigenValuesReference, ComplexMatrixArrayReference* eigenVectorsReference) {
	return ComplexEigenpairsUsingDurandKernerAndInversePowerMethod(M, eigenValuesReference, eigenVectorsReference, 0.000001, 100.0);
}
bool ComplexEigenpairsUsingDurandKernerAndInversePowerMethod(ComplexMatrix* M, cComplexNumberArrayReference* eigenValuesReference, ComplexMatrixArrayReference* eigenVectorsReference, double precision, double maxIterations) {
	bool success, inverseSuccess;
	double n, i, j, k, withinPrecision;
	cComplexNumber* t1, * t2, * t3, * xn1, * eigenValue;
	vector<cComplexNumber*>* rs, * rsPrev;
	pComplexPolynomial* p;
	cComplexNumberArrayReference* evecReference;
	ComplexMatrix* eigenVector;

	evecReference = new cComplexNumberArrayReference();

	p = new pComplexPolynomial();
	ComplexCharacteristicPolynomial(M, p);

	n = p->cs->size() - 1.0;
	rs = new vector<cComplexNumber*>(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		rs->at(i) = cCreateComplexNumber(0.0, 0.0);
	}
	rsPrev = new vector<cComplexNumber*>(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		rsPrev->at(i) = cCreateComplexNumber(0.4, 0.9);
		cPower(rsPrev->at(i), i);
	}
	t2 = cCreateComplexNumber(0.0, 0.0);
	t3 = cCreateComplexNumber(0.0, 0.0);

	success = false;

	eigenVectorsReference->matrices = new vector<ComplexMatrix*>(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		eigenVectorsReference->matrices->at(i) = CreateComplexMatrix(n, 1.0);
	}

	for (i = 0.0; i < maxIterations && !success; i = i + 1.0) {
		for (j = 0.0; j < n; j = j + 1.0) {
			xn1 = rsPrev->at(j);

			t1 = pEvaluateComplex(p, xn1);
			cAssignComplexByValues(t2, 1.0, 0.0);
			for (k = 0.0; k < n; k = k + 1.0) {
				if (k < j) {
					cAssignComplex(t3, xn1);
					cSub(t3, rs->at(k));
					cMul(t2, t3);
				}
				if (k > j) {
					cAssignComplex(t3, xn1);
					cSub(t3, rsPrev->at(k));
					cMul(t2, t3);
				}
			}
			cDiv(t1, t2);
			cAssignComplex(rs->at(j), xn1);
			cSub(rs->at(j), t1);

			delete t1;
		}
		withinPrecision = 0.0;
		for (j = 0.0; j < n; j = j + 1.0) {
			eigenValue = rs->at(j);

			/* Calculate the eigenvector corresponding to the eigenvalue. */
			inverseSuccess = InversePowerMethodComplex(M, eigenValue, i + 1.0, evecReference);
			if (inverseSuccess) {
				for (k = 0.0; k < n; k = k + 1.0) {
					eigenVectorsReference->matrices->at(j)->r->at(k)->c->at(0) = evecReference->complexNumbers->at(k);
				}

				/* Check eigenpair agains precision. */
				eigenVector = eigenVectorsReference->matrices->at(j);

				if (CheckComplexEigenpairPrecision(M, eigenValue, eigenVector, precision)) {
					withinPrecision = withinPrecision + 1.0;
				}
				cAssignComplex(rsPrev->at(j), rs->at(j));
			}
		}
		if (withinPrecision == n) {
			success = true;
		}
	}

	eigenValuesReference->complexNumbers = rs;

	return success;
}
bool CheckComplexEigenpairPrecision(ComplexMatrix* a, cComplexNumber* lambda, ComplexMatrix* e, double precision) {
	ComplexMatrix* vec1, * vec2;
	bool equal;

	vec1 = MultiplyComplexToNew(a, e);
	vec2 = ScalarMultiplyComplexToNew(e, lambda);

	equal = ComplexMatrixEqualsEpsilon(vec1, vec2, precision);

	return equal;
}
void testDeterminant(NumberReference* failures) {
	Matrix* a;
	double d;

	a = createExample4x4_1();

	d = Determinant(a);

	AssertTrue(mathEpsilonCompare(d, 88.0, 0.00000000001), failures);
}
void testInverse(NumberReference* failures) {
	Matrix* a, * b, * f;

	a = createExample4x4_1();

	b = CreateSquareMatrix(4.0);
	Inverse(a, b);

	f = createAnswerInverse4x4_1();

	AssertTrue(MatrixEqualsEpsilon(b, f, 0.00000000001), failures);
}
Matrix* createAnswerInverse4x4_1() {
	Matrix* a;
	MatrixReference* aref;
	vector<wchar_t>* s;

	s = new vector<wchar_t>(0.0);
	s = AppendString(s, toVector(L"-0.13636363636363635,  0.8636363636363636,   -0.6818181818181818,   -0.4090909090909091 \n"));
	s = AppendString(s, toVector(L"-0.6363636363636364,   2.3636363636363638,   -0.9318181818181818,   -0.6590909090909091 \n"));
	s = AppendString(s, toVector(L" 0.045454545454545456, 0.045454545454545456, -0.022727272727272728, -0.11363636363636363 \n"));
	s = AppendString(s, toVector(L" 0.045454545454545456, 0.045454545454545456,  0.22727272727272727,   0.13636363636363635 \n"));

	aref = new MatrixReference();
	ParseMatrixFromString(aref, s, new StringReference());
	a = aref->matrix;

	delete aref;

	return a;
}
Matrix* createExample4x4_1() {
	Matrix* a;
	MatrixReference* aref;
	vector<wchar_t>* s;

	s = new vector<wchar_t>(0.0);
	s = AppendString(s, toVector(L" 5, -2,  2, 7 \n"));
	s = AppendString(s, toVector(L" 1,  0,  0, 3 \n"));
	s = AppendString(s, toVector(L"-3,  1,  5, 0 \n"));
	s = AppendString(s, toVector(L" 3, -1, -9, 4 \n"));

	aref = new MatrixReference();
	ParseMatrixFromString(aref, s, new StringReference());
	a = aref->matrix;

	delete aref;

	return a;
}
void testInverseAlt2(NumberReference* failures) {
	Matrix* a, * b, * f;

	a = createExample4x4_1();

	b = InverseUsingCharacteristicPolynomial(a);

	f = createAnswerInverse4x4_1();

	AssertTrue(MatrixEqualsEpsilon(b, f, 0.00000000001), failures);
}
void TestCholesky(NumberReference* failures) {
	Matrix* a, * b, * f;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L"4, 1,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"1, 5,  3"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L"1, 3, 15"));

	b = CreateSquareMatrix(3.0);
	Cholesky(a, b);

	f = CreateSquareMatrix(3.0);

	f->r->at(0)->c = nStringToNumberArray(toVector(L"2,    0,       0"));
	f->r->at(1)->c = nStringToNumberArray(toVector(L"0.5,  2.17945, 0"));
	f->r->at(2)->c = nStringToNumberArray(toVector(L"0.5,  1.26179, 3.62738"));

	AssertTrue(MatrixEqualsEpsilon(b, f, 0.001), failures);
}
void testCharacteristicPolynomial(NumberReference* failures) {
	Matrix* a;
	vector<double>* coeffs, * answerCoeffs;
	double m;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 0, -4, -6"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"-1,  0, -3"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 1,  2,  5"));

	coeffs = CharacteristicPolynomial(a);

	answerCoeffs = nStringToNumberArray(toVector(L"-4, 8, -5, 1"));

	for (m = 0.0; m < 4.0; m = m + 1.0) {
		AssertTrue(mathEpsilonCompare(answerCoeffs->at(m), coeffs->at(m), 0.000001), failures);
	}
}
void testCharacteristicPolynomial2(NumberReference* failures) {
	Matrix* a;
	vector<double>* coeffs, * answerCoeffs;
	double m;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 3,  1,  5"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 3,  3,  1"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 4,  6,  4"));

	coeffs = CharacteristicPolynomial(a);

	answerCoeffs = new vector<double>(NumberOfRows(a) + 1.0);
	answerCoeffs->at(0) = -40.0;
	answerCoeffs->at(1) = 4.0;
	answerCoeffs->at(2) = -10.0;
	answerCoeffs->at(3) = 1.0;

	for (m = 0.0; m < 4.0; m = m + 1.0) {
		AssertTrue(mathEpsilonCompare(answerCoeffs->at(m), coeffs->at(m), 0.000001), failures);
	}
}
void testCharacteristicPolynomial3(NumberReference* failures) {
	Matrix* a;
	vector<double>* coeffs, * answerCoeffs;
	double m;

	failures = CreateNumberReference(0.0);

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 1, -3,  3"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 3, -5,  3"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 6, -6,  4"));

	coeffs = CharacteristicPolynomial(a);

	answerCoeffs = nStringToNumberArray(toVector(L"-16, -12, 0, 1"));

	for (m = 0.0; m < 4.0; m = m + 1.0) {
		AssertTrue(mathEpsilonCompare(answerCoeffs->at(m), coeffs->at(m), 0.000001), failures);
	}
}
void testEigenValues(NumberReference* failures) {
	Matrix* a;
	NumberArrayReference* evsReference;
	vector<double>* evs, * answerEvs;
	double i, j, m;
	bool found, success;

	a = CreateSquareMatrix(2.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 0,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"-2, -3"));

	evsReference = new NumberArrayReference();
	success = Eigenvalues(a, evsReference);
	AssertTrue(success, failures);

	if (success) {
		evs = evsReference->numberArray;

		answerEvs = new vector<double>(NumberOfRows(a));
		answerEvs->at(0) = -1.0;
		answerEvs->at(1) = -2.0;

		for (i = 0.0; i < 2.0; i = i + 1.0) {
			found = false;
			for (j = 0.0; j < 2.0; j = j + 1.0) {
				if (mathEpsilonCompare(evs->at(i), answerEvs->at(j), 0.00001)) {
					found = true;
				}
			}
			AssertTrue(found, failures);
		}
	}
}
void testEigenValues2(NumberReference* failures) {
	Matrix* a;
	NumberArrayReference* evsReference;
	vector<double>* evs, * answerEvs;
	bool success;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 1, -3,  3"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 3, -5,  3"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 6, -6,  4"));

	evsReference = new NumberArrayReference();
	success = Eigenvalues(a, evsReference);

	AssertTrue(success, failures);

	if (success) {
		evs = evsReference->numberArray;

		answerEvs = new vector<double>(NumberOfRows(a));
		answerEvs->at(0) = -2.0;
		answerEvs->at(1) = -2.0;
		answerEvs->at(2) = 4.0;

		SetsEqualEpsilon(evs, answerEvs, 0.001, failures);
	}
}
void testGaussianEliminationInto(NumberReference* failures) {
	Matrix* a, * f;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 1,  3,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 1,  1, -1"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 3, 11,  5"));

	GaussianElimination(a);

	f = CreateSquareMatrix(3.0);

	f->r->at(0)->c = nStringToNumberArray(toVector(L" 3,  11,  5"));
	f->r->at(1)->c = nStringToNumberArray(toVector(L" 0,  -2.6666666666666, -2.6666666666666"));
	f->r->at(2)->c = nStringToNumberArray(toVector(L" 0, 0,  0"));

	AssertTrue(MatrixEqualsEpsilon(a, f, 0.001), failures);
}
void testEigenpairs(NumberReference* failures) {
	Matrix* a;
	double i, eval;
	NumberArrayReference* eigenvaluesReference;
	MatrixArrayReference* eigenvectorsReference, * eigenvectorsReference2;
	Matrix* evec;
	bool success, equal;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 1,  2,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 6, -1,  0"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L"-1, -2, -1"));

	eigenvectorsReference = new MatrixArrayReference();
	eigenvectorsReference2 = new MatrixArrayReference();
	eigenvaluesReference = new NumberArrayReference();
	success = Eigenpairs(a, eigenvaluesReference, eigenvectorsReference);

	AssertTrue(success, failures);

	if (success) {
		/* Check */
		for (i = 0.0; i < 3.0; i = i + 1.0) {
			evec = eigenvectorsReference->matrices->at(i);
			eval = eigenvaluesReference->numberArray->at(i);

			equal = CheckEigenpairPrecision(a, eval, evec, 0.001);

			AssertTrue(equal, failures);
		}
	}

	success = Eigenvectors(a, eigenvectorsReference2);

	/*System.out.println(new String(MatrixArrayToString(eigenvectorsReference2.matrices, 3d))); */
	AssertTrue(success, failures);

	AssertTrue(MatrixEqualsEpsilon(eigenvectorsReference->matrices->at(0), eigenvectorsReference2->matrices->at(0), 0.001), failures);
	AssertTrue(MatrixEqualsEpsilon(eigenvectorsReference->matrices->at(1), eigenvectorsReference2->matrices->at(1), 0.001), failures);
	AssertTrue(MatrixEqualsEpsilon(eigenvectorsReference->matrices->at(2), eigenvectorsReference2->matrices->at(2), 0.001), failures);
}
void testQRDecomposition(NumberReference* failures) {
	Matrix* a, * fq, * fr, * r, * q;

	r = CreateSquareMatrix(3.0);
	q = CreateSquareMatrix(3.0);

	/* Create example: */
	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L"12, -51,   4"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 6, 167, -68"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L"-4,  24, -41"));

	QRDecomposition(a, q, r);

	fq = CreateSquareMatrix(3.0);

	fq->r->at(0)->c->at(0) = 6.0 / 7.0;
	fq->r->at(0)->c->at(1) = -69.0 / 175.0;
	fq->r->at(0)->c->at(2) = -58.0 / 175.0;

	fq->r->at(1)->c->at(0) = 3.0 / 7.0;
	fq->r->at(1)->c->at(1) = 158.0 / 175.0;
	fq->r->at(1)->c->at(2) = 6.0 / 175.0;

	fq->r->at(2)->c->at(0) = -2.0 / 7.0;
	fq->r->at(2)->c->at(1) = 6.0 / 35.0;
	fq->r->at(2)->c->at(2) = -33.0 / 35.0;

	fr = CreateSquareMatrix(3.0);

	fr->r->at(0)->c = nStringToNumberArray(toVector(L"14,  21, -14"));
	fr->r->at(1)->c = nStringToNumberArray(toVector(L" 0, 175, -70"));
	fr->r->at(2)->c = nStringToNumberArray(toVector(L" 0,   0,  35"));

	AssertTrue(MatrixEqualsEpsilon(r, fr, 0.000000000001), failures);
	AssertTrue(MatrixEqualsEpsilon(q, fq, 0.000000000001), failures);
}
void testQRDecomposition2(NumberReference* failures) {
	Matrix* a, * fq, * fr, * r, * q;

	r = CreateSquareMatrix(2.0);
	q = CreateSquareMatrix(2.0);

	/* Create example: */
	a = CreateSquareMatrix(2.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 0,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"-2, -3"));

	QRDecomposition(a, q, r);

	fq = CreateSquareMatrix(2.0);

	fq->r->at(0)->c->at(0) = 0.0;
	fq->r->at(0)->c->at(1) = 1.0;

	fq->r->at(1)->c->at(0) = -1.0;
	fq->r->at(1)->c->at(1) = 0.0;

	fr = CreateSquareMatrix(2.0);

	fr->r->at(0)->c = nStringToNumberArray(toVector(L"2, 3"));
	fr->r->at(1)->c = nStringToNumberArray(toVector(L"0, 1"));

	AssertTrue(MatrixEqualsEpsilon(r, fr, 0.0000000000001), failures);
	AssertTrue(MatrixEqualsEpsilon(q, fq, 0.0000000000001), failures);
}
void testQRDecomposition3(NumberReference* failures) {
	Matrix* a, * fq, * fr, * r, * q;

	r = CreateSquareMatrix(3.0);
	q = CreateMatrix(4.0, 3.0);

	/* Create example: */
	a = CreateMatrix(4.0, 3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" -1, -1, 1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"  1,  3, 3"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" -1, -1, 5"));
	a->r->at(3)->c = nStringToNumberArray(toVector(L"  1,  3, 7"));

	QRDecomposition(a, q, r);

	fq = CreateMatrix(4.0, 3.0);

	fq->r->at(0)->c = nStringToNumberArray(toVector(L"-0.5, 0.5, -0.5"));
	fq->r->at(1)->c = nStringToNumberArray(toVector(L" 0.5, 0.5, -0.5"));
	fq->r->at(2)->c = nStringToNumberArray(toVector(L"-0.5, 0.5,  0.5"));
	fq->r->at(3)->c = nStringToNumberArray(toVector(L" 0.5, 0.5,  0.5"));

	fr = CreateSquareMatrix(3.0);

	fr->r->at(0)->c = nStringToNumberArray(toVector(L" 2, 4, 2"));
	fr->r->at(1)->c = nStringToNumberArray(toVector(L" 0, 2, 8"));
	fr->r->at(2)->c = nStringToNumberArray(toVector(L" 0, 0, 4"));

	AssertTrue(MatrixEqualsEpsilon(r, fr, 0.0000000000001), failures);
	AssertTrue(MatrixEqualsEpsilon(q, fq, 0.0000000000001), failures);
}
void testEigenValues3QR(NumberReference* failures) {
	Matrix* a, * x, * q, * r;
	vector<double>* evs, * answerEvs;
	bool success;

	a = CreateSquareMatrix(2.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 0,  1"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L"-2, -3"));

	x = CreateMatrix(2.0, 2.0);
	q = CreateMatrix(2.0, 2.0);
	r = CreateMatrix(2.0, 2.0);
	success = QRAlgorithm(a, r, x, q, 0.00001, 100.0);
	AssertTrue(success, failures);

	if (success) {
		evs = new vector<double>(2.0);
		evs->at(0) = x->r->at(0)->c->at(0);
		evs->at(1) = x->r->at(1)->c->at(1);

		answerEvs = new vector<double>(NumberOfRows(a));
		answerEvs->at(0) = -1.0;
		answerEvs->at(1) = -2.0;

		SetsEqualEpsilon(evs, answerEvs, 0.00001, failures);
	}
}
void SetsEqualEpsilon(vector<double>* a, vector<double>* b, double epsilon, NumberReference* failures) {
	double i, j;
	bool found;

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		found = false;
		for (j = 0.0; j < b->size(); j = j + 1.0) {
			if (mathEpsilonCompare(a->at(i), b->at(j), epsilon)) {
				found = true;
			}
		}
		AssertTrue(found, failures);
	}
}
void testLUDecomposition(NumberReference* failures) {
	Matrix* a, * l, * u, * fl, * fu, * ui, * li, * i, * fi;
	bool success;

	a = CreateSquareMatrix(2.0);
	l = CreateSquareMatrix(2.0);
	u = CreateSquareMatrix(2.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 4,  3"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 6,  3"));

	success = LUDecomposition(a, l, u);
	AssertTrue(success, failures);

	fl = CreateSquareMatrix(2.0);
	fl->r->at(0)->c = nStringToNumberArray(toVector(L" 1,  0"));
	fl->r->at(1)->c = nStringToNumberArray(toVector(L" 1.5,  1"));

	fu = CreateSquareMatrix(2.0);
	fu->r->at(0)->c = nStringToNumberArray(toVector(L" 4,  3"));
	fu->r->at(1)->c = nStringToNumberArray(toVector(L" 0,  -1.5"));

	AssertTrue(MatrixEqualsEpsilon(l, fl, 0.001), failures);
	AssertTrue(MatrixEqualsEpsilon(l, fl, 0.001), failures);

	li = CreateSquareMatrix(2.0);
	ui = CreateSquareMatrix(2.0);

	InvertLowerTriangularMatrix(l, li);
	InvertUpperTriangularMatrix(u, ui);

	i = CreateSquareMatrix(2.0);
	fi = CreateIdentityMatrix(2.0);

	Multiply(i, l, li);
	AssertTrue(MatrixEqualsEpsilon(i, fi, 0.001), failures);
	Multiply(i, u, ui);
	AssertTrue(MatrixEqualsEpsilon(i, fi, 0.001), failures);
}
void TestSingularValueDecomposition1(NumberReference* failures) {
	Matrix* A;
	MatrixReference* URef, * SRef, * VRef;
	bool success;

	A = CreateSquareMatrix(3.0);

	A->r->at(0)->c = nStringToNumberArray(toVector(L"-1.000, 0.000, 0.000"));
	A->r->at(1)->c = nStringToNumberArray(toVector(L"-0.830, -0.557, 0.000"));
	A->r->at(2)->c = nStringToNumberArray(toVector(L"-0.691, -0.685, 0.232"));

	URef = new MatrixReference();
	SRef = new MatrixReference();
	VRef = new MatrixReference();

	success = SingularValueDecomposition(A, URef, SRef, VRef);

	CheckSingularValueDecomposition(A, URef, SRef, VRef, success, failures);
}
void CheckSingularValueDecomposition(Matrix* a, MatrixReference* URef, MatrixReference* SRef, MatrixReference* VRef, bool success, NumberReference* failures) {
	Matrix* U, * S, * V, * Ac, * Vt, * I, * Ut, * Ic;
	double i;

	AssertTrue(success, failures);

	if (success) {
		U = URef->matrix;
		S = SRef->matrix;
		V = VRef->matrix;

		Ac = MultiplyToNew(U, S);
		Vt = TransposeToNew(V);
		Ac = MultiplyToNew(Ac, Vt);

		AssertTrue(MatrixEqualsEpsilon(a, Ac, 0.0000001), failures);

		I = CreateSquareMatrix(NumberOfRows(U));
		for (i = 0.0; i < fmin(NumberOfRows(U), NumberOfColumns(U)); i = i + 1.0) {
			I->r->at(i)->c->at(i) = 1.0;
		}
		Ut = TransposeToNew(U);
		Ic = MultiplyToNew(U, Ut);

		AssertTrue(MatrixEqualsEpsilon(I, Ic, 0.0000001), failures);

		I = CreateIdentityMatrix(NumberOfRows(V));
		Ic = MultiplyToNew(V, Vt);

		AssertTrue(MatrixEqualsEpsilon(I, Ic, 0.0000001), failures);
	}
}
void TestSingularValueDecomposition2(NumberReference* failures) {
	Matrix* A;
	MatrixReference* URef, * SRef, * VRef;
	bool success;

	A = CreateMatrix(2.0, 3.0);

	A->r->at(0)->c = nStringToNumberArray(toVector(L"3, 2, 2"));
	A->r->at(1)->c = nStringToNumberArray(toVector(L"2, 3, -2"));

	URef = new MatrixReference();
	SRef = new MatrixReference();
	VRef = new MatrixReference();

	success = SingularValueDecomposition(A, URef, SRef, VRef);

	CheckSingularValueDecomposition(A, URef, SRef, VRef, success, failures);
}
void TestSingularValueDecomposition3(NumberReference* failures) {
	Matrix* A;
	MatrixReference* URef, * SRef, * VRef;
	bool success;

	A = CreateSquareMatrix(3.0);

	A->r->at(0)->c = nStringToNumberArray(toVector(L"2, 0, 0"));
	A->r->at(1)->c = nStringToNumberArray(toVector(L"2, 1, 0"));
	A->r->at(2)->c = nStringToNumberArray(toVector(L"0, -2, 0"));

	URef = new MatrixReference();
	SRef = new MatrixReference();
	VRef = new MatrixReference();

	success = SingularValueDecomposition(A, URef, SRef, VRef);

	CheckSingularValueDecomposition(A, URef, SRef, VRef, success, failures);
}
void TestSingularValueDecomposition4(NumberReference* failures) {
	Matrix* A;
	MatrixReference* URef, * SRef, * VRef;
	bool success;

	A = CreateMatrix(2.0, 3.0);

	A->r->at(0)->c = nStringToNumberArray(toVector(L"1, 1, 0"));
	A->r->at(1)->c = nStringToNumberArray(toVector(L"1, 0, 1"));

	URef = new MatrixReference();
	SRef = new MatrixReference();
	VRef = new MatrixReference();

	success = SingularValueDecomposition(A, URef, SRef, VRef);

	CheckSingularValueDecomposition(A, URef, SRef, VRef, success, failures);
}
void TestSingularValueDecomposition5(NumberReference* failures) {
	Matrix* A;
	MatrixReference* URef, * SRef, * VRef;
	bool success;

	A = CreateMatrix(3.0, 2.0);

	A->r->at(0)->c = nStringToNumberArray(toVector(L"1, 1"));
	A->r->at(1)->c = nStringToNumberArray(toVector(L"1, 0"));
	A->r->at(2)->c = nStringToNumberArray(toVector(L"0, 1"));

	URef = new MatrixReference();
	SRef = new MatrixReference();
	VRef = new MatrixReference();

	success = SingularValueDecomposition(A, URef, SRef, VRef);

	CheckSingularValueDecomposition(A, URef, SRef, VRef, success, failures);
}
Matrix* CreateMatrixWithComplexEigenvalues() {
	Matrix* a;

	a = CreateSquareMatrix(3.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 0, 1, 0"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 0, 0, 1"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 1, 0, 0"));

	return a;
}
Matrix* CreateMatrixWithRepeatedEigenvalues() {
	Matrix* a;

	a = CreateSquareMatrix(4.0);

	a->r->at(0)->c = nStringToNumberArray(toVector(L" 2, 0, 0, 0"));
	a->r->at(1)->c = nStringToNumberArray(toVector(L" 1, 2, 0, 0"));
	a->r->at(2)->c = nStringToNumberArray(toVector(L" 0, 1, 3, 0"));
	a->r->at(3)->c = nStringToNumberArray(toVector(L" 0, 0, 1, 3"));

	return a;
}
double test() {
	NumberReference* failures;
	Matrix* m;
	vector<double>* diag;

	failures = CreateNumberReference(0.0);

	testDeterminant(failures);
	testInverse(failures);
	testInverseAlt2(failures);
	TestCholesky(failures);
	testCharacteristicPolynomial(failures);
	testCharacteristicPolynomial2(failures);
	testCharacteristicPolynomial3(failures);
	testComplexCharacteristicPolynomial(failures);
	testComplexCharacteristicPolynomial2(failures);
	testEigenValues(failures);
	testEigenValues2(failures);
	testEigenValues3QR(failures);
	testComplexEigenValues(failures);
	testComplexEigenValues2(failures);
	testGaussianEliminationInto(failures);
	testEigenpairs(failures);
	testComplexMatrices(failures);

	testQRDecomposition(failures);
	testQRDecomposition2(failures);
	testQRDecomposition3(failures);
	testComplexInversePowerMethod(failures);
	testComplexEigenpairs(failures);
	testLUDecomposition(failures);

	TestSingularValueDecomposition1(failures);
	TestSingularValueDecomposition2(failures);
	TestSingularValueDecomposition3(failures);
	TestSingularValueDecomposition4(failures);
	TestSingularValueDecomposition5(failures);

	m = CreateSquareMatrix(2.0);
	m->r->at(0)->c->at(0) = 1.0;
	m->r->at(0)->c->at(1) = 2.0;
	m->r->at(1)->c->at(0) = -3.0;
	m->r->at(1)->c->at(1) = 2.0;

	Transpose(m);
	ElementWisePower(m, 2.0);
	m = CreateDiagonalMatrixFromArray(m->r->at(0)->c);
	diag = new vector<double>(2.0);
	ExtractDiagonal(m, diag);

	return failures->numberValue;
}
void testComplexMatrices(NumberReference* failures) {
	testComplexDeterminant(failures);
	testComplexInverse(failures);
}
void testComplexDeterminant(NumberReference* failures) {
	Matrix* a;
	ComplexMatrix* ca;
	cComplexNumber* d;

	a = createExample4x4_1();
	ca = CreateComplexMatrixFromMatrix(a);

	d = DeterminantComplex(ca);

	AssertTrue(mathEpsilonCompare(d->re, 88.0, 0.00000000001), failures);
	AssertTrue(mathEpsilonCompare(d->im, 0.0, 0.00000000001), failures);

	ca = CreateExampleComplexMatrix2x2();

	d = DeterminantComplex(ca);

	AssertTrue(mathEpsilonCompare(d->re, 0.0, 0.00000000001), failures);
	AssertTrue(mathEpsilonCompare(d->im, -1.0, 0.00000000001), failures);
}
ComplexMatrix* CreateExampleComplexMatrix2x2() {
	ComplexMatrix* ca;
	ca = CreateComplexMatrix(2.0, 2.0);
	IndexComplex(ca, 0.0, 0.0)->re = 0.0;
	IndexComplex(ca, 0.0, 0.0)->im = 1.0;
	IndexComplex(ca, 0.0, 1.0)->re = 1.0;
	IndexComplex(ca, 0.0, 1.0)->im = 0.0;
	IndexComplex(ca, 1.0, 0.0)->re = 0.0;
	IndexComplex(ca, 1.0, 0.0)->im = 0.0;
	IndexComplex(ca, 1.0, 1.0)->re = -1.0;
	IndexComplex(ca, 1.0, 1.0)->im = 0.0;
	return ca;
}
ComplexMatrix* CreateExampleComplexMatrix2x2Inverse() {
	ComplexMatrix* ca;
	ca = CreateComplexMatrix(2.0, 2.0);
	IndexComplex(ca, 0.0, 0.0)->re = 0.0;
	IndexComplex(ca, 0.0, 0.0)->im = -1.0;
	IndexComplex(ca, 0.0, 1.0)->re = 0.0;
	IndexComplex(ca, 0.0, 1.0)->im = -1.0;
	IndexComplex(ca, 1.0, 0.0)->re = 0.0;
	IndexComplex(ca, 1.0, 0.0)->im = 0.0;
	IndexComplex(ca, 1.0, 1.0)->re = -1.0;
	IndexComplex(ca, 1.0, 1.0)->im = 0.0;
	return ca;
}
void testComplexInverse(NumberReference* failures) {
	Matrix* a, * f;
	ComplexMatrix* ca, * b, * cf;

	a = createExample4x4_1();
	ca = CreateComplexMatrixFromMatrix(a);

	b = CreateSquareComplexMatrix(4.0);
	InverseComplex(ca, b);

	f = createAnswerInverse4x4_1();
	cf = CreateComplexMatrixFromMatrix(f);

	AssertTrue(ComplexMatrixEqualsEpsilon(b, cf, 0.00000000001), failures);

	ca = CreateExampleComplexMatrix2x2();

	b = CreateSquareComplexMatrix(2.0);
	InverseComplex(ca, b);

	cf = CreateExampleComplexMatrix2x2Inverse();

	AssertTrue(ComplexMatrixEqualsEpsilon(b, cf, 0.00000000001), failures);
}
void testComplexCharacteristicPolynomial(NumberReference* failures) {
	ComplexMatrix* a;
	pComplexPolynomial* coeffs, * answerCoeffs;
	double m;
	bool success;

	a = CreateComplexMatrix(3.0, 3.0);
	IndexComplex(a, 0.0, 0.0)->re = 0.0;
	IndexComplex(a, 0.0, 0.0)->im = 0.0;
	IndexComplex(a, 0.0, 1.0)->re = -4.0;
	IndexComplex(a, 0.0, 1.0)->im = 0.0;
	IndexComplex(a, 0.0, 2.0)->re = -6.0;
	IndexComplex(a, 0.0, 2.0)->im = 0.0;
	IndexComplex(a, 1.0, 0.0)->re = -1.0;
	IndexComplex(a, 1.0, 0.0)->im = 0.0;
	IndexComplex(a, 1.0, 1.0)->re = 0.0;
	IndexComplex(a, 1.0, 1.0)->im = 0.0;
	IndexComplex(a, 1.0, 2.0)->re = -3.0;
	IndexComplex(a, 1.0, 2.0)->im = 0.0;
	IndexComplex(a, 2.0, 0.0)->re = 1.0;
	IndexComplex(a, 2.0, 0.0)->im = 0.0;
	IndexComplex(a, 2.0, 1.0)->re = 2.0;
	IndexComplex(a, 2.0, 1.0)->im = 0.0;
	IndexComplex(a, 2.0, 2.0)->re = 5.0;
	IndexComplex(a, 2.0, 2.0)->im = 0.0;

	coeffs = new pComplexPolynomial();
	ComplexCharacteristicPolynomial(a, coeffs);

	answerCoeffs = pCreateComplexPolynomial(3.0);
	answerCoeffs->cs->at(0) = cCreateComplexNumber(-4.0, 0.0);
	answerCoeffs->cs->at(1) = cCreateComplexNumber(8.0, 0.0);
	answerCoeffs->cs->at(2) = cCreateComplexNumber(-5.0, 0.0);
	answerCoeffs->cs->at(3) = cCreateComplexNumber(1.0, 0.0);

	for (m = 0.0; m < 4.0; m = m + 1.0) {
		AssertTrue(cEpsilonCompareComplex(answerCoeffs->cs->at(m), coeffs->cs->at(m), 0.000001), failures);
	}
}
void testComplexCharacteristicPolynomial2(NumberReference* failures) {
	ComplexMatrix* a;
	pComplexPolynomial* p, * ap;
	double m;

	a = CreateComplexMatrix(2.0, 2.0);
	IndexComplex(a, 0.0, 0.0)->re = 0.0;
	IndexComplex(a, 0.0, 0.0)->im = 1.0;
	IndexComplex(a, 0.0, 1.0)->re = 0.0;
	IndexComplex(a, 0.0, 1.0)->im = 0.0;
	IndexComplex(a, 1.0, 0.0)->re = 0.0;
	IndexComplex(a, 1.0, 0.0)->im = 1.0;
	IndexComplex(a, 1.0, 1.0)->re = 1.0;
	IndexComplex(a, 1.0, 1.0)->im = 0.0;

	p = new pComplexPolynomial();
	ComplexCharacteristicPolynomial(a, p);

	ap = pCreateComplexPolynomial(2.0);
	ap->cs->at(0) = cCreateComplexNumber(0.0, 1.0);
	ap->cs->at(1) = cCreateComplexNumber(-1.0, -1.0);
	ap->cs->at(2) = cCreateComplexNumber(1.0, 0.0);

	for (m = 0.0; m < 3.0; m = m + 1.0) {
		AssertTrue(cEpsilonCompareComplex(ap->cs->at(m), p->cs->at(m), 0.000001), failures);
	}
}
void testComplexEigenValues(NumberReference* failures) {
	ComplexMatrix* a;
	cComplexNumberArrayReference* evsReference;
	vector<cComplexNumber*>* evs, * answerEvs;
	double i, j;
	bool found, success;

	a = CreateComplexMatrix(2.0, 2.0);
	IndexComplex(a, 0.0, 0.0)->re = 0.0;
	IndexComplex(a, 0.0, 0.0)->im = 0.0;
	IndexComplex(a, 0.0, 1.0)->re = 1.0;
	IndexComplex(a, 0.0, 1.0)->im = 0.0;
	IndexComplex(a, 1.0, 0.0)->re = -2.0;
	IndexComplex(a, 1.0, 0.0)->im = 0.0;
	IndexComplex(a, 1.0, 1.0)->re = -3.0;
	IndexComplex(a, 1.0, 1.0)->im = 0.0;

	evsReference = new cComplexNumberArrayReference();
	success = EigenvaluesComplex(a, evsReference);
	AssertTrue(success, failures);

	if (success) {
		evs = evsReference->complexNumbers;

		answerEvs = new vector<cComplexNumber*>(NumberOfRowsComplex(a));
		answerEvs->at(0) = cCreateComplexNumber(-1.0, 0.0);
		answerEvs->at(1) = cCreateComplexNumber(-2.0, 0.0);

		for (i = 0.0; i < 2.0; i = i + 1.0) {
			found = false;
			for (j = 0.0; j < 2.0; j = j + 1.0) {
				if (cEpsilonCompareComplex(evs->at(i), answerEvs->at(j), 0.00001)) {
					found = true;
				}
			}
			AssertTrue(found, failures);
		}
	}
}
void testComplexEigenValues2(NumberReference* failures) {
	ComplexMatrix* a;
	cComplexNumberArrayReference* evsReference;
	vector<cComplexNumber*>* evs, * answerEvs;
	double i, j;
	bool found, success;

	a = CreateComplex3x3Matrix();

	evsReference = new cComplexNumberArrayReference();
	success = EigenvaluesComplex(a, evsReference);
	AssertTrue(success, failures);

	if (success) {
		evs = evsReference->complexNumbers;

		answerEvs = new vector<cComplexNumber*>(NumberOfRowsComplex(a));
		answerEvs->at(0) = cCreateComplexNumber(1.0, 0.0);
		answerEvs->at(1) = cCreateComplexNumber(-0.5, 0.8660254037844386467637231);
		answerEvs->at(2) = cCreateComplexNumber(-0.5, -0.8660254037844386467637231);

		for (i = 0.0; i < 3.0; i = i + 1.0) {
			found = false;
			for (j = 0.0; j < 3.0; j = j + 1.0) {
				if (cEpsilonCompareComplex(evs->at(i), answerEvs->at(j), 0.00001)) {
					found = true;
				}
			}
			AssertTrue(found, failures);
		}
	}
}
ComplexMatrix* CreateComplex3x3Matrix() {
	ComplexMatrix* a;
	a = CreateComplexMatrix(3.0, 3.0);
	IndexComplex(a, 0.0, 0.0)->re = 0.0;
	IndexComplex(a, 0.0, 0.0)->im = 0.0;
	IndexComplex(a, 0.0, 1.0)->re = 1.0;
	IndexComplex(a, 0.0, 1.0)->im = 0.0;
	IndexComplex(a, 0.0, 2.0)->re = 0.0;
	IndexComplex(a, 0.0, 2.0)->im = 0.0;
	IndexComplex(a, 1.0, 0.0)->re = 0.0;
	IndexComplex(a, 1.0, 0.0)->im = 0.0;
	IndexComplex(a, 1.0, 1.0)->re = 0.0;
	IndexComplex(a, 1.0, 1.0)->im = 0.0;
	IndexComplex(a, 1.0, 2.0)->re = 1.0;
	IndexComplex(a, 1.0, 2.0)->im = 0.0;
	IndexComplex(a, 2.0, 0.0)->re = 1.0;
	IndexComplex(a, 2.0, 0.0)->im = 0.0;
	IndexComplex(a, 2.0, 1.0)->re = 0.0;
	IndexComplex(a, 2.0, 1.0)->im = 0.0;
	IndexComplex(a, 2.0, 2.0)->re = 0.0;
	IndexComplex(a, 2.0, 2.0)->im = 0.0;
	return a;
}
void testComplexInversePowerMethod(NumberReference* failures) {
	ComplexMatrix* a, * evecMatrix;
	cComplexNumberArrayReference* evsReference;
	cComplexNumber* eval, * evalApprox, * answerEvs;
	double i, N;
	bool found, success, equals;
	cComplexNumberArrayReference* eigenvector;

	a = CreateComplex3x3Matrix();

	eval = cCreateComplexNumber(-0.5, 0.8660254037844386467637231);
	evalApprox = cCreateComplexNumber(-0.5 * 1.1, 0.8660254037844386467637231 * 1.1);

	eigenvector = new cComplexNumberArrayReference();
	success = InversePowerMethodComplex(a, evalApprox, 100.0, eigenvector);

	AssertTrue(success, failures);

	N = NumberOfRowsComplex(a);
	evecMatrix = CreateComplexMatrix(N, 1.0);
	for (i = 0.0; i < N; i = i + 1.0) {
		cAssignComplex(evecMatrix->r->at(i)->c->at(0), eigenvector->complexNumbers->at(i));
	}

	equals = CheckComplexEigenpairPrecision(a, eval, evecMatrix, 0.00001);

	AssertTrue(equals, failures);
}
void testComplexEigenpairs(NumberReference* failures) {
	ComplexMatrix* a;
	double i;
	cComplexNumber* eval;
	cComplexNumberArrayReference* eigenvaluesReference;
	ComplexMatrixArrayReference* eigenvectorsReference, * eigenvectorsReference2;
	ComplexMatrix* evec;
	bool success, equal;

	a = CreateSquareComplexMatrix(3.0);

	IndexComplex(a, 0.0, 0.0)->re = 1.0;
	IndexComplex(a, 0.0, 0.0)->im = 0.0;
	IndexComplex(a, 0.0, 1.0)->re = 2.0;
	IndexComplex(a, 0.0, 1.0)->im = 0.0;
	IndexComplex(a, 0.0, 2.0)->re = 1.0;
	IndexComplex(a, 0.0, 2.0)->im = 0.0;
	IndexComplex(a, 1.0, 0.0)->re = 6.0;
	IndexComplex(a, 1.0, 0.0)->im = 0.0;
	IndexComplex(a, 1.0, 1.0)->re = -1.0;
	IndexComplex(a, 1.0, 1.0)->im = 0.0;
	IndexComplex(a, 1.0, 2.0)->re = 0.0;
	IndexComplex(a, 1.0, 2.0)->im = 0.0;
	IndexComplex(a, 2.0, 0.0)->re = -1.0;
	IndexComplex(a, 2.0, 0.0)->im = 0.0;
	IndexComplex(a, 2.0, 1.0)->re = -2.0;
	IndexComplex(a, 2.0, 1.0)->im = 0.0;
	IndexComplex(a, 2.0, 2.0)->re = -1.0;
	IndexComplex(a, 2.0, 2.0)->im = 0.0;

	eigenvectorsReference = new ComplexMatrixArrayReference();
	eigenvaluesReference = new cComplexNumberArrayReference();
	eigenvectorsReference2 = new ComplexMatrixArrayReference();
	success = EigenpairsComplex(a, eigenvaluesReference, eigenvectorsReference);

	AssertTrue(success, failures);

	if (success) {
		/* Check */
		for (i = 0.0; i < 3.0; i = i + 1.0) {
			evec = eigenvectorsReference->matrices->at(i);
			eval = eigenvaluesReference->complexNumbers->at(i);

			equal = CheckComplexEigenpairPrecision(a, eval, evec, 0.001);

			AssertTrue(equal, failures);
		}
	}

	success = EigenvectorsComplex(a, eigenvectorsReference2);

	AssertTrue(success, failures);

	AssertTrue(ComplexMatrixEqualsEpsilon(eigenvectorsReference->matrices->at(0), eigenvectorsReference2->matrices->at(0), 0.001), failures);
	AssertTrue(ComplexMatrixEqualsEpsilon(eigenvectorsReference->matrices->at(1), eigenvectorsReference2->matrices->at(1), 0.001), failures);
	AssertTrue(ComplexMatrixEqualsEpsilon(eigenvectorsReference->matrices->at(2), eigenvectorsReference2->matrices->at(2), 0.001), failures);
}
bool FindRoots(vector<double>* p, NumberArrayReference* rootsReference) {
	return DurandKernerMethod(p, 0.000001, 100.0, rootsReference);
}
bool LaguerresMethodWithRepeatedDivision(vector<double>* p, double maxIterations, double precision, double guess, NumberArrayReference* rootsReference) {
	double n, nr, xk;
	vector<double>* x;
	vector<double>* q, * r, * d;
	bool success;
	NumberReference* xkReference;

	n = pDegree(p);

	x = new vector<double>(n);

	q = pCreatePolynomial(n);
	r = pCreatePolynomial(n);
	d = pCreatePolynomial(n);

	success = true;
	xkReference = CreateNumberReference(0.0);

	for (nr = 0.0; nr < n && success; nr = nr + 1.0) {
		success = LaguerresMethod(p, guess, maxIterations, precision, xkReference);

		if (success) {
			xk = xkReference->numberValue;
			x->at(nr) = xk;

			pFill(d, 0.0);
			d->at(0) = -xk;
			d->at(1) = 1.0;
			pDivide(q, r, p, d);
			pAssign(p, q);
		}
	}

	delete q;
	delete r;
	delete d;

	rootsReference->numberArray = x;

	return success;
}
bool LaguerresMethod(vector<double>* p, double guess, double maxIterations, double precision, NumberReference* rootReference) {
	double k, a, G, H, xk, denom1, denom2, denom, t1, n;
	bool success;

	n = pDegree(p);
	success = true;

	xk = guess;

	for (k = 0.0; (k < maxIterations) && (abs(pEvaluate(p, xk)) >= precision) && success; k = k + 1.0) {
		G = pEvaluateDerivative(p, xk, 1.0) / pEvaluate(p, xk);
		H = pow(G, 2.0) - pEvaluateDerivative(p, xk, 2.0) / pEvaluate(p, xk);
		t1 = (n - 1.0) * (n * H - pow(G, 2.0));
		if (t1 >= 0.0) {
			denom = sqrt(t1);
			denom1 = G + denom;
			denom2 = G - denom;
			if (abs(denom1) >= abs(denom2)) {
				denom = denom1;
			}
			else {
				denom = denom2;
			}
			a = n / denom;

			xk = xk - a;
		}
		else {
			success = false;
		}
	}

	if (k == maxIterations) {
		success = false;
	}

	if (abs(pEvaluate(p, xk)) >= precision) {
		success = false;
	}

	rootReference->numberValue = xk;

	return success;
}
bool DurandKernerMethod(vector<double>* p, double precision, double maxIterations, NumberArrayReference* rootsReference) {
	bool success;
	double n, i, j, k, t1, t2, xn1, withinPrecision;
	vector<double>* rs, * rsPrev;

	n = p->size() - 1.0;
	rs = new vector<double>(n);
	rsPrev = new vector<double>(n);

	for (i = 0.0; i < n; i = i + 1.0) {
		rsPrev->at(i) = i;
	}

	success = false;

	for (i = 0.0; i < maxIterations && !success; i = i + 1.0) {
		for (j = 0.0; j < n; j = j + 1.0) {
			xn1 = rsPrev->at(j);

			t1 = pEvaluate(p, xn1);
			t2 = 1.0;
			for (k = 0.0; k < n; k = k + 1.0) {
				if (k < j) {
					t2 = t2 * (xn1 - rs->at(k));
				}
				if (k > j) {
					t2 = t2 * (xn1 - rsPrev->at(k));
				}
			}
			t1 = t1 / t2;
			rs->at(j) = xn1 - t1;
		}
		withinPrecision = 0.0;
		for (j = 0.0; j < n; j = j + 1.0) {
			if (mathEpsilonCompare(rsPrev->at(j), rs->at(j), precision)) {
				withinPrecision = withinPrecision + 1.0;
			}
			rsPrev->at(j) = rs->at(j);
		}
		if (withinPrecision == n) {
			success = true;
		}
	}

	rootsReference->numberArray = rs;

	return success;
}
bool FindRootsComplex(pComplexPolynomial* p, cComplexNumberArrayReference* rootsReference) {
	return DurandKernerMethodComplex(p, 0.000001, 100.0, rootsReference);
}
bool DurandKernerMethodComplex(pComplexPolynomial* p, double precision, double maxIterations, cComplexNumberArrayReference* rootsReference) {
	bool success;
	double n, i, j, k, withinPrecision;
	cComplexNumber* t1, * t2, * t3, * xn1;
	vector<cComplexNumber*>* rs, * rsPrev;

	n = p->cs->size() - 1.0;
	rs = new vector<cComplexNumber*>(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		rs->at(i) = cCreateComplexNumber(0.0, 0.0);
	}
	rsPrev = new vector<cComplexNumber*>(n);
	for (i = 0.0; i < n; i = i + 1.0) {
		rsPrev->at(i) = cCreateComplexNumber(0.4, 0.9);
		cPower(rsPrev->at(i), i);
	}
	t2 = cCreateComplexNumber(0.0, 0.0);
	t3 = cCreateComplexNumber(0.0, 0.0);

	success = false;

	for (i = 0.0; i < maxIterations && !success; i = i + 1.0) {
		for (j = 0.0; j < n; j = j + 1.0) {
			xn1 = rsPrev->at(j);

			t1 = pEvaluateComplex(p, xn1);
			cAssignComplexByValues(t2, 1.0, 0.0);
			for (k = 0.0; k < n; k = k + 1.0) {
				if (k < j) {
					cAssignComplex(t3, xn1);
					cSub(t3, rs->at(k));
					cMul(t2, t3);
				}
				if (k > j) {
					cAssignComplex(t3, xn1);
					cSub(t3, rsPrev->at(k));
					cMul(t2, t3);
				}
			}
			cDiv(t1, t2);
			cAssignComplex(rs->at(j), xn1);
			cSub(rs->at(j), t1);
		}
		withinPrecision = 0.0;
		for (j = 0.0; j < n; j = j + 1.0) {
			if (cEpsilonCompareComplex(rsPrev->at(j), rs->at(j), precision)) {
				withinPrecision = withinPrecision + 1.0;
			}
			cAssignComplex(rsPrev->at(j), rs->at(j));
		}
		if (withinPrecision == n) {
			success = true;
		}
	}

	rootsReference->complexNumbers = rs;

	return success;
}
BooleanReference* CreateBooleanReference(bool value) {
	BooleanReference* ref;

	ref = new BooleanReference();
	ref->booleanValue = value;

	return ref;
}
BooleanArrayReference* CreateBooleanArrayReference(vector<bool>* value) {
	BooleanArrayReference* ref;

	ref = new BooleanArrayReference();
	ref->booleanArray = value;

	return ref;
}
BooleanArrayReference* CreateBooleanArrayReferenceLengthValue(double length, bool value) {
	BooleanArrayReference* ref;
	double i;

	ref = new BooleanArrayReference();
	ref->booleanArray = new vector<bool>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		ref->booleanArray->at(i) = value;
	}

	return ref;
}
void FreeBooleanArrayReference(BooleanArrayReference* booleanArrayReference) {
	delete booleanArrayReference->booleanArray;
	delete booleanArrayReference;
}
CharacterReference* CreateCharacterReference(wchar_t value) {
	CharacterReference* ref;

	ref = new CharacterReference();
	ref->characterValue = value;

	return ref;
}
NumberReference* CreateNumberReference(double value) {
	NumberReference* ref;

	ref = new NumberReference();
	ref->numberValue = value;

	return ref;
}
NumberArrayReference* CreateNumberArrayReference(vector<double>* value) {
	NumberArrayReference* ref;

	ref = new NumberArrayReference();
	ref->numberArray = value;

	return ref;
}
NumberArrayReference* CreateNumberArrayReferenceLengthValue(double length, double value) {
	NumberArrayReference* ref;
	double i;

	ref = new NumberArrayReference();
	ref->numberArray = new vector<double>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		ref->numberArray->at(i) = value;
	}

	return ref;
}
void FreeNumberArrayReference(NumberArrayReference* numberArrayReference) {
	delete numberArrayReference->numberArray;
	delete numberArrayReference;
}
StringReference* CreateStringReference(vector<wchar_t>* value) {
	StringReference* ref;

	ref = new StringReference();
	ref->string = value;

	return ref;
}
StringReference* CreateStringReferenceLengthValue(double length, wchar_t value) {
	StringReference* ref;
	double i;

	ref = new StringReference();
	ref->string = new vector<wchar_t>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		ref->string->at(i) = value;
	}

	return ref;
}
void FreeStringReference(StringReference* stringReference) {
	delete stringReference->string;
	delete stringReference;
}
StringArrayReference* CreateStringArrayReference(vector<StringReference*>* strings) {
	StringArrayReference* ref;

	ref = new StringArrayReference();
	ref->stringArray = strings;

	return ref;
}
StringArrayReference* CreateStringArrayReferenceLengthValue(double length, vector<wchar_t>* value) {
	StringArrayReference* ref;
	double i;

	ref = new StringArrayReference();
	ref->stringArray = new vector<StringReference*>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		ref->stringArray->at(i) = CreateStringReference(value);
	}

	return ref;
}
void FreeStringArrayReference(StringArrayReference* stringArrayReference) {
	double i;

	for (i = 0.0; i < stringArrayReference->stringArray->size(); i = i + 1.0) {
		delete stringArrayReference->stringArray->at(i);
	}
	delete stringArrayReference->stringArray;
	delete stringArrayReference;
}
double mathNegate(double x) {
	return  -x;
}
double mathPositive(double x) {
	return  +x;
}
double mathFactorial(double x) {
	double i, f;

	f = 1.0;

	for (i = 2.0; i <= x; i = i + 1.0) {
		f = f * i;
	}

	return f;
}
double mathRound(double x) {
	return floor(x + 0.5);
}
double mathBankersRound(double x) {
	double r;

	if (mathAbsolute(x - mathTruncate(x)) == 0.5) {
		if (!mathDivisibleBy(mathRound(x), 2.0)) {
			r = mathRound(x) - 1.0;
		}
		else {
			r = mathRound(x);
		}
	}
	else {
		r = mathRound(x);
	}

	return r;
}
double mathCeil(double x) {
	return ceil(x);
}
double mathFloor(double x) {
	return floor(x);
}
double mathTruncate(double x) {
	double t;

	if (x >= 0.0) {
		t = floor(x);
	}
	else {
		t = ceil(x);
	}

	return t;
}
double mathAbsolute(double x) {
	return abs(x);
}
double mathLogarithm(double x) {
	return log10(x);
}
double mathNaturalLogarithm(double x) {
	return log(x);
}
double mathSin(double x) {
	return sin(x);
}
double mathCos(double x) {
	return cos(x);
}
double mathTan(double x) {
	return tan(x);
}
double mathAsin(double x) {
	return asin(x);
}
double mathAcos(double x) {
	return acos(x);
}
double mathAtan(double x) {
	return atan(x);
}
double mathAtan2(double y, double x) {
	double a;

	/* Atan2 is an invalid operation when x = 0 and y = 0, but this method does not return errors. */
	a = 0.0;

	if (x > 0.0) {
		a = mathAtan(y / x);
	}
	else if (x < 0.0 && y >= 0.0) {
		a = mathAtan(y / x) + 3.141596;
	}
	else if (x < 0.0 && y < 0.0) {
		a = mathAtan(y / x) - 3.141596;
	}
	else if (x == 0.0 && y > 0.0) {
		a = 3.141596 / 2.0;
	}
	else if (x == 0.0 && y < 0.0) {
		a = -3.141596 / 2.0;
	}

	return a;
}
double mathSquareroot(double x) {
	return sqrt(x);
}
double mathExp(double x) {
	return exp(x);
}
bool mathDivisibleBy(double a, double b) {
	return ((fmod(a, b)) == 0.0);
}
double mathCombinations(double n, double k) {
	double i, j, c;

	c = 1.0;
	j = 1.0;
	i = n - k + 1.0;

	for (; i <= n; ) {
		c = c * i;
		c = c / j;

		i = i + 1.0;
		j = j + 1.0;
	}

	return c;
}
double mathPermutations(double n, double k) {
	double i, c;

	c = 1.0;

	for (i = n - k + 1.0; i <= n; i = i + 1.0) {
		c = c * i;
	}

	return c;
}
bool mathEpsilonCompare(double a, double b, double epsilon) {
	return abs(a - b) < epsilon;
}
double mathGreatestCommonDivisor(double a, double b) {
	double t;

	for (; b != 0.0; ) {
		t = b;
		b = fmod(a, b);
		a = t;
	}

	return a;
}
double mathGCDWithSubtraction(double a, double b) {
	double g;

	if (a == 0.0) {
		g = b;
	}
	else {
		for (; b != 0.0; ) {
			if (a > b) {
				a = a - b;
			}
			else {
				b = b - a;
			}
		}

		g = a;
	}

	return g;
}
bool mathIsInteger(double a) {
	return (a - floor(a)) == 0.0;
}
bool mathGreatestCommonDivisorWithCheck(double a, double b, NumberReference* gcdReference) {
	bool success;
	double gcd;

	if (mathIsInteger(a) && mathIsInteger(b)) {
		gcd = mathGreatestCommonDivisor(a, b);
		gcdReference->numberValue = gcd;
		success = true;
	}
	else {
		success = false;
	}

	return success;
}
double mathLeastCommonMultiple(double a, double b) {
	double lcm;

	if (a > 0.0 && b > 0.0) {
		lcm = abs(a * b) / mathGreatestCommonDivisor(a, b);
	}
	else {
		lcm = 0.0;
	}

	return lcm;
}
double mathSign(double a) {
	double s;

	if (a > 0.0) {
		s = 1.0;
	}
	else if (a < 0.0) {
		s = -1.0;
	}
	else {
		s = 0.0;
	}

	return s;
}
double mathMax(double a, double b) {
	return fmax(a, b);
}
double mathMin(double a, double b) {
	return fmin(a, b);
}
double mathPower(double a, double b) {
	return pow(a, b);
}
double mathGamma(double x) {
	return mathLanczosApproximation(x);
}
double mathLogGamma(double x) {
	return log(mathGamma(x));
}
double mathLanczosApproximation(double z) {
	vector<double>* p;
	double i, y, t, x;

	p = new vector<double>(8.0);
	p->at(0) = 676.5203681218851;
	p->at(1) = -1259.1392167224028;
	p->at(2) = 771.32342877765313;
	p->at(3) = -176.61502916214059;
	p->at(4) = 12.507343278686905;
	p->at(5) = -0.13857109526572012;
	p->at(6) = 9.9843695780195716e-6;
	p->at(7) = 1.5056327351493116e-7;

	if (z < 0.5) {
		y = 3.141596 / (sin(3.141596 * z) * mathLanczosApproximation(1.0 - z));
	}
	else {
		z = z - 1.0;
		x = 0.99999999999980993;
		for (i = 0.0; i < p->size(); i = i + 1.0) {
			x = x + p->at(i) / (z + i + 1.0);
		}
		t = z + p->size() - 0.5;
		y = sqrt(2.0 * 3.141596) * pow(t, z + 0.5) * exp(-t) * x;
	}

	return y;
}
double mathBeta(double x, double y) {
	return mathGamma(x) * mathGamma(y) / mathGamma(x + y);
}
double mathSinh(double x) {
	return (exp(x) - exp(-x)) / 2.0;
}
double mathCosh(double x) {
	return (exp(x) + exp(-x)) / 2.0;
}
double mathTanh(double x) {
	return mathSinh(x) / mathCosh(x);
}
double mathCot(double x) {
	return 1.0 / tan(x);
}
double mathSec(double x) {
	return 1.0 / cos(x);
}
double mathCsc(double x) {
	return 1.0 / sin(x);
}
double mathCoth(double x) {
	return mathCosh(x) / mathSinh(x);
}
double mathSech(double x) {
	return 1.0 / mathCosh(x);
}
double mathCsch(double x) {
	return 1.0 / mathSinh(x);
}
double mathError(double x) {
	double y, t, tau, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

	if (x == 0.0) {
		y = 0.0;
	}
	else if (x < 0.0) {
		y = -mathError(-x);
	}
	else {
		c1 = -1.26551223;
		c2 = +1.00002368;
		c3 = +0.37409196;
		c4 = +0.09678418;
		c5 = -0.18628806;
		c6 = +0.27886807;
		c7 = -1.13520398;
		c8 = +1.48851587;
		c9 = -0.82215223;
		c10 = +0.17087277;

		t = 1.0 / (1.0 + 0.5 * abs(x));

		tau = t * exp(-pow(x, 2.0) + c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * (c6 + t * (c7 + t * (c8 + t * (c9 + t * c10)))))))));

		y = 1.0 - tau;
	}

	return y;
}
double mathErrorInverse(double x) {
	double y, a, t;

	a = (8.0 * (3.141596 - 3.0)) / (3.0 * 3.141596 * (4.0 - 3.141596));

	t = 2.0 / (3.141596 * a) + log(1.0 - pow(x, 2.0)) / 2.0;
	y = mathSign(x) * sqrt(sqrt(pow(t, 2.0) - log(1.0 - pow(x, 2.0)) / a) - t);

	return y;
}
double mathFallingFactorial(double x, double n) {
	double k, y;

	y = 1.0;

	for (k = 0.0; k <= n - 1.0; k = k + 1.0) {
		y = y * (x - k);
	}

	return y;
}
double mathRisingFactorial(double x, double n) {
	double k, y;

	y = 1.0;

	for (k = 0.0; k <= n - 1.0; k = k + 1.0) {
		y = y * (x + k);
	}

	return y;
}
double mathHypergeometric(double a, double b, double c, double z, double maxIterations, double precision) {
	double y;

	if (abs(z) >= 0.5) {
		y = pow(1.0 - z, -a) * mathHypergeometricDirect(a, c - b, c, z / (z - 1.0), maxIterations, precision);
	}
	else {
		y = mathHypergeometricDirect(a, b, c, z, maxIterations, precision);
	}

	return y;
}
double mathHypergeometricDirect(double a, double b, double c, double z, double maxIterations, double precision) {
	double y, yp, n;
	bool done;

	y = 0.0;
	done = false;

	for (n = 0.0; n < maxIterations && !done; n = n + 1.0) {
		yp = mathRisingFactorial(a, n) * mathRisingFactorial(b, n) / mathRisingFactorial(c, n) * pow(z, n) / mathFactorial(n);
		if (abs(yp) < precision) {
			done = true;
		}
		y = y + yp;
	}

	return y;
}
double mathBernouilliNumber(double n) {
	return mathAkiyamaTanigawaAlgorithm(n);
}
double mathAkiyamaTanigawaAlgorithm(double n) {
	double m, j, B;
	vector<double>* A;

	A = new vector<double>(n + 1.0);

	for (m = 0.0; m <= n; m = m + 1.0) {
		A->at(m) = 1.0 / (m + 1.0);
		for (j = m; j >= 1.0; j = j - 1.0) {
			A->at(j - 1.0) = j * (A->at(j - 1.0) - A->at(j));
		}
	}

	B = A->at(0);

	delete A;

	return B;
}
cComplexNumber* cCreateComplexNumber(double re, double im) {
	cComplexNumber* z;

	z = new cComplexNumber();
	z->re = re;
	z->im = im;

	return z;
}
cPolarComplexNumber* cCreatePolarComplexNumber(double r, double phi) {
	cPolarComplexNumber* p;

	p = new cPolarComplexNumber();
	p->r = r;
	p->phi = phi;

	return p;
}
void cAdd(cComplexNumber* z1, cComplexNumber* z2) {
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	z1->re = a + c;
	z1->im = b + d;
}
cComplexNumber* cAddToNew(cComplexNumber* z1, cComplexNumber* z2) {
	cComplexNumber* x;
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	x = new cComplexNumber();

	x->re = a + c;
	x->im = b + d;

	return x;
}
void cSub(cComplexNumber* z1, cComplexNumber* z2) {
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	z1->re = a - c;
	z1->im = b - d;
}
cComplexNumber* cSubToNew(cComplexNumber* z1, cComplexNumber* z2) {
	cComplexNumber* x;
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	x = new cComplexNumber();

	x->re = a - c;
	x->im = b - d;

	return x;
}
void cMul(cComplexNumber* z1, cComplexNumber* z2) {
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	z1->re = a * c - b * d;
	z1->im = b * c + a * d;
}
cComplexNumber* cMulToNew(cComplexNumber* z1, cComplexNumber* z2) {
	cComplexNumber* x;
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	x = new cComplexNumber();

	x->re = a * c - b * d;
	x->im = b * c + a * d;

	return x;
}
void cDiv(cComplexNumber* z1, cComplexNumber* z2) {
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	z1->re = (a * c + b * d) / (pow(c, 2.0) + pow(d, 2.0));
	z1->im = (b * c - a * d) / (pow(c, 2.0) + pow(d, 2.0));
}
cComplexNumber* cDivToNew(cComplexNumber* z1, cComplexNumber* z2) {
	cComplexNumber* x;
	double a, b, c, d;

	a = z1->re;
	b = z1->im;
	c = z2->re;
	d = z2->im;

	x = new cComplexNumber();

	x->re = (a * c + b * d) / (pow(c, 2.0) + pow(d, 2.0));
	x->im = (b * c - a * d) / (pow(c, 2.0) + pow(d, 2.0));

	return x;
}
void cConjugate(cComplexNumber* z) {
	z->im = -z->im;
}
cComplexNumber* cConjugateToNew(cComplexNumber* z) {
	cComplexNumber* x;

	x = new cComplexNumber();

	x->re = z->re;
	x->im = -z->im;

	return x;
}
double cAbs(cComplexNumber* z) {
	double x;

	x = sqrt(pow(z->re, 2.0) + pow(z->im, 2.0));

	return x;
}
double cArg(cComplexNumber* z) {
	double x;

	x = mathAtan2(z->im, z->re);

	return x;
}
cPolarComplexNumber* cCreatePolarFromComplexNumber(cComplexNumber* z) {
	cPolarComplexNumber* x;

	x = new cPolarComplexNumber();

	x->r = cAbs(z);
	x->phi = cArg(z);

	return x;
}
cComplexNumber* cCreateComplexFromPolar(cPolarComplexNumber* p) {
	cComplexNumber* z;

	z = new cComplexNumber();

	z->re = p->r * cos(p->phi);
	z->im = p->r * sin(p->phi);

	return z;
}
double cRe(cComplexNumber* z) {
	return z->re;
}
double cIm(cComplexNumber* z) {
	return z->im;
}
void cAddPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	cComplexNumber* z1, * z2;

	z1 = cCreateComplexFromPolar(p1);
	z2 = cCreateComplexFromPolar(p2);

	cAdd(z1, z2);

	x = cCreatePolarFromComplexNumber(z1);

	p1->r = x->r;
	p1->phi = x->phi;

	delete z1;
	delete z2;
	delete x;
}
cPolarComplexNumber* cAddPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	cComplexNumber* z1, * z2;

	z1 = cCreateComplexFromPolar(p1);
	z2 = cCreateComplexFromPolar(p2);

	cAdd(z1, z2);

	x = cCreatePolarFromComplexNumber(z1);

	delete z1;
	delete z2;

	return x;
}
void cSubPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	cComplexNumber* z1, * z2;

	z1 = cCreateComplexFromPolar(p1);
	z2 = cCreateComplexFromPolar(p2);

	cSub(z1, z2);

	x = cCreatePolarFromComplexNumber(z1);

	p1->r = x->r;
	p1->phi = x->phi;

	delete z1;
	delete z2;
	delete x;
}
cPolarComplexNumber* cSubPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	cComplexNumber* z1, * z2;

	z1 = cCreateComplexFromPolar(p1);
	z2 = cCreateComplexFromPolar(p2);

	cSub(z1, z2);

	x = cCreatePolarFromComplexNumber(z1);

	delete z1;
	delete z2;

	return x;
}
void cMulPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	double r1, r2, phi1, phi2;

	r1 = p1->r;
	r2 = p2->r;
	phi1 = p1->phi;
	phi2 = p2->phi;

	p1->r = r1 * r2;
	p1->phi = phi1 + phi2;
}
cPolarComplexNumber* cMulPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	double r1, r2, phi1, phi2;

	r1 = p1->r;
	r2 = p2->r;
	phi1 = p1->phi;
	phi2 = p2->phi;

	x = new cPolarComplexNumber();

	x->r = r1 * r2;
	x->phi = phi1 + phi2;

	return x;
}
void cDivPolar(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	double r1, r2, phi1, phi2;

	r1 = p1->r;
	r2 = p2->r;
	phi1 = p1->phi;
	phi2 = p2->phi;

	p1->r = r1 / r2;
	p1->phi = phi1 - phi2;
}
cPolarComplexNumber* cDivPolarToNew(cPolarComplexNumber* p1, cPolarComplexNumber* p2) {
	cPolarComplexNumber* x;
	double r1, r2, phi1, phi2;

	r1 = p1->r;
	r2 = p2->r;
	phi1 = p1->phi;
	phi2 = p2->phi;

	x = new cPolarComplexNumber();

	x->r = r1 / r2;
	x->phi = phi1 - phi2;

	return x;
}
void cSquareRoot(cComplexNumber* z) {
	double a, b, m;

	a = z->re;
	b = z->im;

	m = sqrt(pow(a, 2.0) + pow(b, 2.0));

	z->re = sqrt((m + a) / 2.0);
	z->im = mathSign(b) * sqrt((m - a) / 2.0);
}
void cPowerPolar(cPolarComplexNumber* p, double n) {
	p->r = pow(p->r, n);
	p->phi = p->phi * n;
}
cComplexNumber* cPowerToNew(cComplexNumber* z, double n) {
	cPolarComplexNumber* p;
	cComplexNumber* zp;

	p = cCreatePolarFromComplexNumber(z);
	cPowerPolar(p, n);
	zp = cCreateComplexFromPolar(p);

	delete p;

	return zp;
}
void cPower(cComplexNumber* z, double n) {
	cComplexNumber* zp;

	zp = cPowerToNew(z, n);
	z->re = zp->re;
	z->im = zp->im;

	delete zp;
}
void cNegate(cComplexNumber* z) {
	z->re = mathNegate(z->re);
	z->im = mathNegate(z->im);
}
void cAssignComplexByValues(cComplexNumber* s, double re, double im) {
	s->re = re;
	s->im = im;
}
void cAssignComplex(cComplexNumber* a, cComplexNumber* b) {
	a->re = b->re;
	a->im = b->im;
}
bool cEpsilonCompareComplex(cComplexNumber* a, cComplexNumber* b, double epsilon) {
	return mathEpsilonCompare(a->re, b->re, epsilon) && mathEpsilonCompare(a->im, b->im, epsilon);
}
void cExpComplex(cComplexNumber* x) {
	double re, im;

	re = exp(x->re) * cos(x->im);
	im = exp(x->re) * sin(x->im);
	x->re = re;
	x->im = im;
}
void cSineComplex(cComplexNumber* x) {
	double re, im;

	re = sin(x->re) * mathCosh(x->im);
	im = cos(x->re) * mathSinh(x->im);
	x->re = re;
	x->im = im;
}
void cCosineComplex(cComplexNumber* x) {
	double re, im;

	re = cos(x->re) * mathCosh(x->im);
	im = sin(x->re) * mathSinh(x->im);
	x->re = re;
	x->im = im;
}
vector<wchar_t>* cComplexToString(cComplexNumber* a) {
	vector<wchar_t>* str, * number;
	LinkedListCharacters* ll;
	double i;

	ll = CreateLinkedListCharacter();

	number = nCreateStringDecimalFromNumber(a->re);

	for (i = 0.0; i < number->size(); i = i + 1.0) {
		LinkedListAddCharacter(ll, number->at(i));
	}

	delete number;

	if (a->im < 0.0) {
		LinkedListAddCharacter(ll, '-');
		number = nCreateStringDecimalFromNumber(-a->im);
	}
	else {
		LinkedListAddCharacter(ll, '+');
		number = nCreateStringDecimalFromNumber(a->im);
	}

	for (i = 0.0; i < number->size(); i = i + 1.0) {
		LinkedListAddCharacter(ll, number->at(i));
	}

	delete number;

	LinkedListAddCharacter(ll, 'i');

	str = LinkedListCharactersToArray(ll);
	FreeLinkedListCharacter(ll);

	return str;
}
vector<wchar_t>* pPolynomialToTextDirect(vector<double>* p, vector<wchar_t>* x) {
	LinkedListCharacters* ll;
	vector<wchar_t>* str;
	double i, c, j;
	StringReference* buffer;

	buffer = new StringReference();

	ll = CreateLinkedListCharacter();

	if (p->size() == 0.0) {
		LinkedListAddCharacter(ll, '0');
	}
	else {
		for (i = 0.0; i < p->size(); i = i + 1.0) {
			c = p->at(i);
			if (c < 0.0) {
				LinkedListAddCharacter(ll, '-');
			}
			else {
				LinkedListAddCharacter(ll, '+');
			}

			nCreateStringFromNumberWithCheck(abs(c), 10.0, buffer);
			pLinkedListCharactersAddString(ll, buffer->string);
			delete buffer->string;

			pLinkedListCharactersAddString(ll, x);
			LinkedListAddCharacter(ll, '^');

			nCreateStringFromNumberWithCheck(i, 10.0, buffer);
			pLinkedListCharactersAddString(ll, buffer->string);
			delete buffer->string;
		}
	}

	str = LinkedListCharactersToArray(ll);
	FreeLinkedListCharacter(ll);

	return str;
}
void pLinkedListCharactersAddString(LinkedListCharacters* ll, vector<wchar_t>* str) {
	double j;

	for (j = 0.0; j < str->size(); j = j + 1.0) {
		LinkedListAddCharacter(ll, str->at(j));
	}
}
void pGenerateCommonRenderSpecification(vector<double>* p, BooleanArrayReference* showCoefficient, StringReference* sign, NumberArrayReference* coefficient, BooleanArrayReference* showPower, BooleanArrayReference* showX) {
	double i, zeros, c;
	bool setZero;

	setZero = false;

	if (p->size() == 0.0) {
		setZero = true;
	}
	else {
		zeros = 0.0;
		for (i = 0.0; i < p->size(); i = i + 1.0) {
			if (p->at(i) == 0.0) {
				zeros = zeros + 1.0;
			}
		}

		if (zeros == p->size()) {
			setZero = true;
		}
		else {
			showCoefficient->booleanArray = new vector<bool>(p->size());
			sign->string = new vector<wchar_t>(p->size());
			coefficient->numberArray = new vector<double>(p->size());
			showPower->booleanArray = new vector<bool>(p->size());
			showX->booleanArray = new vector<bool>(p->size());

			for (i = 0.0; i < p->size(); i = i + 1.0) {
				c = p->at(i);

				if (c < 0.0) {
					sign->string->at(i) = '-';
				}
				else {
					sign->string->at(i) = '+';
				}
				coefficient->numberArray->at(i) = abs(p->at(i));
				if (c == 0.0) {
					showCoefficient->booleanArray->at(i) = false;
				}
				else {
					if (abs(c) == 1.0 && i > 0.0) {
						showCoefficient->booleanArray->at(i) = false;
					}
					else {
						showCoefficient->booleanArray->at(i) = true;
					}

					if (i == 0.0) {
						showX->booleanArray->at(i) = false;
						showPower->booleanArray->at(i) = false;
					}
					else {
						showX->booleanArray->at(i) = true;
						if (i == 1.0) {
							showPower->booleanArray->at(i) = false;
						}
						else {
							showPower->booleanArray->at(i) = true;
						}
					}
				}
			}
		}
	}

	if (setZero) {
		showCoefficient->booleanArray = new vector<bool>(1.0);
		sign->string = new vector<wchar_t>(1.0);
		coefficient->numberArray = new vector<double>(1.0);
		showPower->booleanArray = new vector<bool>(1.0);
		showX->booleanArray = new vector<bool>(1.0);

		showCoefficient->booleanArray->at(0) = true;
		sign->string->at(0) = '+';
		coefficient->numberArray->at(0) = 0.0;
		showPower->booleanArray->at(0) = true;
		showX->booleanArray->at(0) = false;
	}
}
vector<wchar_t>* pPolynomialToText(vector<double>* p, vector<wchar_t>* x) {
	LinkedListCharacters* ll;
	vector<wchar_t>* str;
	double i, c;
	StringReference* buffer;
	BooleanArrayReference* showCoefficient, * showPower, * showX;
	StringReference* sign;
	NumberArrayReference* coefficient;
	bool hasPrinted;

	showCoefficient = CreateBooleanArrayReferenceLengthValue(0.0, false);
	showPower = CreateBooleanArrayReferenceLengthValue(0.0, false);
	showX = CreateBooleanArrayReferenceLengthValue(0.0, false);
	sign = CreateStringReferenceLengthValue(0.0, ' ');
	coefficient = CreateNumberArrayReferenceLengthValue(0.0, 0.0);

	pGenerateCommonRenderSpecification(p, showCoefficient, sign, coefficient, showPower, showX);

	buffer = CreateStringReferenceLengthValue(0.0, ' ');

	ll = CreateLinkedListCharacter();

	hasPrinted = false;
	for (i = 0.0; i < showCoefficient->booleanArray->size(); i = i + 1.0) {
		c = p->at(i);

		if (showCoefficient->booleanArray->at(i) || showX->booleanArray->at(i)) {
			if (!hasPrinted && c >= 0.0) {
			}
			else {
				LinkedListAddCharacter(ll, sign->string->at(i));
			}

			if (showCoefficient->booleanArray->at(i)) {
				nCreateStringFromNumberWithCheck(coefficient->numberArray->at(i), 10.0, buffer);
				pLinkedListCharactersAddString(ll, buffer->string);
				delete buffer->string;
				hasPrinted = true;
			}

			if (showX->booleanArray->at(i)) {
				pLinkedListCharactersAddString(ll, x);
				hasPrinted = true;
				if (showPower->booleanArray->at(i)) {
					LinkedListAddCharacter(ll, '^');
					nCreateStringFromNumberWithCheck(i, 10.0, buffer);
					pLinkedListCharactersAddString(ll, buffer->string);
					delete buffer->string;
				}
			}
		}
	}

	str = LinkedListCharactersToArray(ll);
	FreeLinkedListCharacter(ll);

	return str;
}
vector<wchar_t>* pComplexPolynomialToTextDirect(pComplexPolynomial* p, vector<wchar_t>* x) {
	LinkedListCharacters* ll;
	vector<wchar_t>* str, * number;
	double i;
	cComplexNumber* c;
	StringReference* buffer;

	buffer = new StringReference();

	ll = CreateLinkedListCharacter();

	if (p->cs->size() == 0.0) {
		LinkedListAddCharacter(ll, '0');
	}
	else {
		for (i = 0.0; i < p->cs->size(); i = i + 1.0) {
			c = p->cs->at(i);
			if (i > 0.0) {
				LinkedListAddCharacter(ll, '+');
			}

			number = cComplexToString(c);
			LinkedListAddCharacter(ll, '(');
			pLinkedListCharactersAddString(ll, number);
			LinkedListAddCharacter(ll, ')');
			delete number;

			pLinkedListCharactersAddString(ll, x);
			LinkedListAddCharacter(ll, '^');

			nCreateStringFromNumberWithCheck(i, 10.0, buffer);
			pLinkedListCharactersAddString(ll, buffer->string);
			delete buffer->string;
		}
	}

	str = LinkedListCharactersToArray(ll);
	FreeLinkedListCharacter(ll);

	return str;
}
void pAdd(vector<double>* a, vector<double>* b) {
	double i, nr;

	nr = b->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		a->at(i) = a->at(i) + b->at(i);
	}
}
void pSubtract(vector<double>* a, vector<double>* b) {
	double i, nr;

	nr = b->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		a->at(i) = a->at(i) - b->at(i);
	}
}
void pMultiply(vector<double>* c, vector<double>* a, vector<double>* b) {
	double k, n, m, i;
	double av, bv;

	n = pDegree(a);
	m = pDegree(b);

	pFill(c, 0.0);

	for (i = 0.0; i <= n + m; i = i + 1.0) {
		c->at(i) = 0.0;
		for (k = 0.0; k <= i && k < a->size() && i - k < b->size(); k = k + 1.0) {
			av = a->at(k);
			bv = b->at(i - k);
			c->at(i) = c->at(i) + av * bv;
		}
	}
}
void pDivide(vector<double>* q, vector<double>* r, vector<double>* n, vector<double>* d) {
	vector<double>* t, * t1;
	double tcoff, tdegree, i;
	double deg;
	double rd, dd;

	pFill(q, 0.0);
	pAssign(r, n);
	deg = pDegree(n);
	t = pCreatePolynomial(deg);
	t1 = pCreatePolynomial(deg);

	rd = pDegree(r);
	dd = pDegree(d);
	for (i = 0.0; i < deg + 1.0 && !pIsZero(r) && rd - i >= dd; i = i + 1.0) {
		pFill(t, 0.0);
		tdegree = rd - i - dd;
		tcoff = r->at(rd - i) / d->at(dd);
		t->at(tdegree) = tcoff;
		pAdd(q, t);
		pFill(t1, 0.0);
		pMultiply(t1, t, d);
		pSubtract(r, t1);
	}

	delete t;
	delete t1;
}
bool pIsZero(vector<double>* a) {
	double i;
	bool itIsZero;

	itIsZero = true;

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		if (a->at(i) != 0.0) {
			itIsZero = false;
		}
	}

	return itIsZero;
}
void pAssign(vector<double>* a, vector<double>* b) {
	double i, nr;

	nr = b->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		a->at(i) = b->at(i);
	}
}
vector<double>* pCreatePolynomial(double deg) {
	vector<double>* p;

	p = new vector<double>(deg + 1.0);

	pFill(p, 0.0);

	return p;
}
void pFill(vector<double>* p, double value) {
	double i;

	for (i = 0.0; i < p->size(); i = i + 1.0) {
		p->at(i) = value;
	}
}
double pDegree(vector<double>* A) {
	double i;
	double deg;
	bool done;

	done = false;
	deg = 0.0;
	for (i = A->size() - 1.0; i >= 0.0 && !done; i = i - 1.0) {
		if (A->at(i) != 0.0) {
			deg = i;
			done = true;
		}
	}

	return deg;
}
double pLead(vector<double>* A) {
	double deg;

	deg = pDegree(A);

	return A->at(deg);
}
double pEvaluate(vector<double>* A, double x) {
	return pEvaluateWithHornersMethod(A, x);
}
double pEvaluateWithHornersMethod(vector<double>* A, double x) {
	double r, i;

	r = 0.0;

	for (i = A->size() - 1.0; i >= 0.0; i = i - 1.0) {
		r = r * x;
		r = A->at(i) + r;
	}

	return r;
}
double pEvaluateWithPowers(vector<double>* A, double x) {
	double r, i;

	r = 0.0;

	for (i = 0.0; i < A->size(); i = i + 1.0) {
		r = r + A->at(i) * pow(x, i);
	}

	return r;
}
double pEvaluateDerivative(vector<double>* A, double x, double n) {
	double r, i, v;

	r = 0.0;

	for (i = 0.0; i < A->size(); i = i + 1.0) {
		if (i - n >= 0.0) {
			/*v = A[(int) i] * mFactorial(n) * mCombinations(i, n) * pow(x, i - n); */
			v = A->at(i) * mathPermutations(i, n) * pow(x, i - n);
			r = r + v;
		}
	}

	return r;
}
void pDerivative(vector<double>* A) {
	double i, degree;

	degree = 0.0;
	for (i = 1.0; i < A->size(); i = i + 1.0) {
		degree = degree + 1.0;
		A->at(i - 1.0) = degree * A->at(i);
	}

	A->at(A->size() - 1.0) = 0.0;
}
void pAddComplex(pComplexPolynomial* a, pComplexPolynomial* b) {
	double i, nr;

	nr = a->cs->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		cAdd(a->cs->at(i), b->cs->at(i));
	}
}
void pSubtractComplex(pComplexPolynomial* a, pComplexPolynomial* b) {
	double i, nr;

	nr = a->cs->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		cSub(a->cs->at(i), b->cs->at(i));
	}
}
bool pIsZeroComplex(pComplexPolynomial* a) {
	double i;
	bool itIsZero;

	itIsZero = true;

	for (i = 0.0; i < a->cs->size(); i = i + 1.0) {
		if (a->cs->at(i)->re != 0.0 && a->cs->at(i)->im != 0.0) {
			itIsZero = false;
		}
	}

	return itIsZero;
}
void pAssignComplex(pComplexPolynomial* a, pComplexPolynomial* b) {
	double i, nr;

	nr = b->cs->size();

	for (i = 0.0; i < nr; i = i + 1.0) {
		cAssignComplex(a->cs->at(i), b->cs->at(i));
	}
}
pComplexPolynomial* pCreateComplexPolynomial(double deg) {
	pComplexPolynomial* p;
	double i;

	p = new pComplexPolynomial();
	p->cs = new vector<cComplexNumber*>(deg + 1.0);

	for (i = 0.0; i < deg + 1.0; i = i + 1.0) {
		p->cs->at(i) = new cComplexNumber();
	}

	pFillComplex(p, 0.0, 0.0);

	return p;
}
void pFillComplex(pComplexPolynomial* p, double re, double im) {
	double i;
	cComplexNumber* c;

	c = cCreateComplexNumber(re, im);

	for (i = 0.0; i < p->cs->size(); i = i + 1.0) {
		cAssignComplex(p->cs->at(i), c);
	}

	delete c;
}
double pDegreeComplex(pComplexPolynomial* A) {
	double i;
	double deg;
	bool done;

	done = false;
	deg = 0.0;
	for (i = A->cs->size() - 1.0; i >= 0.0 && !done; i = i - 1.0) {
		if (A->cs->at(i)->re != 0.0 && A->cs->at(i)->im != 0.0) {
			deg = i;
			done = true;
		}
	}

	return deg;
}
cComplexNumber* pLeadComplex(pComplexPolynomial* A) {
	double deg;

	deg = pDegreeComplex(A);

	return A->cs->at(deg);
}
cComplexNumber* pEvaluateComplex(pComplexPolynomial* A, cComplexNumber* x) {
	double i;
	cComplexNumber* r, * t;

	r = cCreateComplexNumber(0.0, 0.0);
	t = cCreateComplexNumber(0.0, 0.0);

	for (i = 0.0; i < A->cs->size(); i = i + 1.0) {
		cAssignComplex(t, x);
		cPower(t, i);
		cMul(t, A->cs->at(i));
		cAdd(r, t);
	}

	return r;
}
double pTotalNumberOfRoots(vector<double>* p) {
	return pDegree(p);
}
void WriteStringToStingStream(vector<wchar_t>* stream, NumberReference* index, vector<wchar_t>* src) {
	double i;

	for (i = 0.0; i < src->size(); i = i + 1.0) {
		stream->at(index->numberValue + i) = src->at(i);
	}
	index->numberValue = index->numberValue + src->size();
}
void WriteCharacterToStingStream(vector<wchar_t>* stream, NumberReference* index, wchar_t src) {
	stream->at(index->numberValue) = src;
	index->numberValue = index->numberValue + 1.0;
}
void WriteBooleanToStingStream(vector<wchar_t>* stream, NumberReference* index, bool src) {
	if (src) {
		WriteStringToStingStream(stream, index, toVector(L"true"));
	}
	else {
		WriteStringToStingStream(stream, index, toVector(L"false"));
	}
}
bool SubstringWithCheck(vector<wchar_t>* string, double from, double to, StringReference* stringReference) {
	bool success;

	if (from >= 0.0 && from <= string->size() && to >= 0.0 && to <= string->size() && from <= to) {
		stringReference->string = Substring(string, from, to);
		success = true;
	}
	else {
		success = false;
	}

	return success;
}
vector<wchar_t>* Substring(vector<wchar_t>* string, double from, double to) {
	vector<wchar_t>* n;
	double i, length;

	length = to - from;

	n = new vector<wchar_t>(length);

	for (i = from; i < to; i = i + 1.0) {
		n->at(i - from) = string->at(i);
	}

	return n;
}
vector<wchar_t>* AppendString(vector<wchar_t>* s1, vector<wchar_t>* s2) {
	vector<wchar_t>* newString;

	newString = ConcatenateString(s1, s2);

	delete s1;

	return newString;
}
vector<wchar_t>* ConcatenateString(vector<wchar_t>* s1, vector<wchar_t>* s2) {
	vector<wchar_t>* newString;
	double i;

	newString = new vector<wchar_t>(s1->size() + s2->size());

	for (i = 0.0; i < s1->size(); i = i + 1.0) {
		newString->at(i) = s1->at(i);
	}

	for (i = 0.0; i < s2->size(); i = i + 1.0) {
		newString->at(s1->size() + i) = s2->at(i);
	}

	return newString;
}
vector<wchar_t>* AppendCharacter(vector<wchar_t>* string, wchar_t c) {
	vector<wchar_t>* newString;

	newString = ConcatenateCharacter(string, c);

	delete string;

	return newString;
}
vector<wchar_t>* ConcatenateCharacter(vector<wchar_t>* string, wchar_t c) {
	vector<wchar_t>* newString;
	double i;
	newString = new vector<wchar_t>(string->size() + 1.0);

	for (i = 0.0; i < string->size(); i = i + 1.0) {
		newString->at(i) = string->at(i);
	}

	newString->at(string->size()) = c;

	return newString;
}
vector<StringReference*>* SplitByCharacter(vector<wchar_t>* toSplit, wchar_t splitBy) {
	vector<StringReference*>* split;
	vector<wchar_t>* stringToSplitBy;

	stringToSplitBy = new vector<wchar_t>(1.0);
	stringToSplitBy->at(0) = splitBy;

	split = SplitByString(toSplit, stringToSplitBy);

	delete stringToSplitBy;

	return split;
}
bool IndexOfCharacter(vector<wchar_t>* string, wchar_t character, NumberReference* indexReference) {
	double i;
	bool found;

	found = false;
	for (i = 0.0; i < string->size() && !found; i = i + 1.0) {
		if (string->at(i) == character) {
			found = true;
			indexReference->numberValue = i;
		}
	}

	return found;
}
bool SubstringEqualsWithCheck(vector<wchar_t>* string, double from, vector<wchar_t>* substring, BooleanReference* equalsReference) {
	bool success;

	if (from < string->size()) {
		success = true;
		equalsReference->booleanValue = SubstringEquals(string, from, substring);
	}
	else {
		success = false;
	}

	return success;
}
bool SubstringEquals(vector<wchar_t>* string, double from, vector<wchar_t>* substring) {
	double i;
	bool equal;

	equal = true;
	for (i = 0.0; i < substring->size() && equal; i = i + 1.0) {
		if (string->at(from + i) != substring->at(i)) {
			equal = false;
		}
	}

	return equal;
}
bool IndexOfString(vector<wchar_t>* string, vector<wchar_t>* substring, NumberReference* indexReference) {
	double i;
	bool found;

	found = false;
	for (i = 0.0; i < string->size() - substring->size() + 1.0 && !found; i = i + 1.0) {
		if (SubstringEquals(string, i, substring)) {
			found = true;
			indexReference->numberValue = i;
		}
	}

	return found;
}
bool ContainsCharacter(vector<wchar_t>* string, wchar_t character) {
	double i;
	bool found;

	found = false;
	for (i = 0.0; i < string->size() && !found; i = i + 1.0) {
		if (string->at(i) == character) {
			found = true;
		}
	}

	return found;
}
bool ContainsString(vector<wchar_t>* string, vector<wchar_t>* substring) {
	return IndexOfString(string, substring, new NumberReference());
}
void ToUpperCase(vector<wchar_t>* string) {
	double i;

	for (i = 0.0; i < string->size(); i = i + 1.0) {
		string->at(i) = charToUpperCase(string->at(i));
	}
}
void ToLowerCase(vector<wchar_t>* string) {
	double i;

	for (i = 0.0; i < string->size(); i = i + 1.0) {
		string->at(i) = charToLowerCase(string->at(i));
	}
}
bool EqualsIgnoreCase(vector<wchar_t>* a, vector<wchar_t>* b) {
	bool equal;
	double i;

	if (a->size() == b->size()) {
		equal = true;
		for (i = 0.0; i < a->size() && equal; i = i + 1.0) {
			if (charToLowerCase(a->at(i)) != charToLowerCase(b->at(i))) {
				equal = false;
			}
		}
	}
	else {
		equal = false;
	}

	return equal;
}
vector<wchar_t>* ReplaceString(vector<wchar_t>* string, vector<wchar_t>* toReplace, vector<wchar_t>* replaceWith) {
	vector<wchar_t>* result;
	double i;
	BooleanReference* equalsReference;
	bool success;

	equalsReference = new BooleanReference();
	result = new vector<wchar_t>(0.0);

	for (i = 0.0; i < string->size(); ) {
		success = SubstringEqualsWithCheck(string, i, toReplace, equalsReference);
		if (success) {
			success = equalsReference->booleanValue;
		}

		if (success && toReplace->size() > 0.0) {
			result = ConcatenateString(result, replaceWith);
			i = i + toReplace->size();
		}
		else {
			result = ConcatenateCharacter(result, string->at(i));
			i = i + 1.0;
		}
	}

	return result;
}
vector<wchar_t>* ReplaceCharacter(vector<wchar_t>* string, wchar_t toReplace, wchar_t replaceWith) {
	vector<wchar_t>* result;
	double i;

	result = new vector<wchar_t>(0.0);

	for (i = 0.0; i < string->size(); i = i + 1.0) {
		if (string->at(i) == toReplace) {
			result = ConcatenateCharacter(result, replaceWith);
		}
		else {
			result = ConcatenateCharacter(result, string->at(i));
		}
	}

	return result;
}
vector<wchar_t>* Trim(vector<wchar_t>* string) {
	vector<wchar_t>* result;
	double i, lastWhitespaceLocationStart, lastWhitespaceLocationEnd;
	bool firstNonWhitespaceFound;

	/* Find whitepaces at the start. */
	lastWhitespaceLocationStart = -1.0;
	firstNonWhitespaceFound = false;
	for (i = 0.0; i < string->size() && !firstNonWhitespaceFound; i = i + 1.0) {
		if (charIsWhiteSpace(string->at(i))) {
			lastWhitespaceLocationStart = i;
		}
		else {
			firstNonWhitespaceFound = true;
		}
	}

	/* Find whitepaces at the end. */
	lastWhitespaceLocationEnd = string->size();
	firstNonWhitespaceFound = false;
	for (i = string->size() - 1.0; i >= 0.0 && !firstNonWhitespaceFound; i = i - 1.0) {
		if (charIsWhiteSpace(string->at(i))) {
			lastWhitespaceLocationEnd = i;
		}
		else {
			firstNonWhitespaceFound = true;
		}
	}

	if (lastWhitespaceLocationStart < lastWhitespaceLocationEnd) {
		result = Substring(string, lastWhitespaceLocationStart + 1.0, lastWhitespaceLocationEnd);
	}
	else {
		result = new vector<wchar_t>(0.0);
	}

	return result;
}
bool StartsWith(vector<wchar_t>* string, vector<wchar_t>* start) {
	bool startsWithString;

	startsWithString = false;
	if (string->size() >= start->size()) {
		startsWithString = SubstringEquals(string, 0.0, start);
	}

	return startsWithString;
}
bool EndsWith(vector<wchar_t>* string, vector<wchar_t>* end) {
	bool endsWithString;

	endsWithString = false;
	if (string->size() >= end->size()) {
		endsWithString = SubstringEquals(string, string->size() - end->size(), end);
	}

	return endsWithString;
}
vector<StringReference*>* SplitByString(vector<wchar_t>* toSplit, vector<wchar_t>* splitBy) {
	vector<StringReference*>* split;
	vector<wchar_t>* next;
	double i;
	wchar_t c;
	StringReference* n;

	split = new vector<StringReference*>(0.0);

	next = new vector<wchar_t>(0.0);
	for (i = 0.0; i < toSplit->size(); ) {
		c = toSplit->at(i);

		if (SubstringEquals(toSplit, i, splitBy)) {
			if (split->size() != 0.0 || i != 0.0) {
				n = new StringReference();
				n->string = next;
				split = AddString(split, n);
				next = new vector<wchar_t>(0.0);
				i = i + splitBy->size();
			}
		}
		else {
			next = AppendCharacter(next, c);
			i = i + 1.0;
		}
	}

	if (next->size() > 0.0) {
		n = new StringReference();
		n->string = next;
		split = AddString(split, n);
	}

	return split;
}
bool StringIsBefore(vector<wchar_t>* a, vector<wchar_t>* b) {
	bool before, equal, done;
	double i;

	before = false;
	equal = true;
	done = false;

	if (a->size() == 0.0 && b->size() > 0.0) {
		before = true;
	}
	else {
		for (i = 0.0; i < a->size() && i < b->size() && !done; i = i + 1.0) {
			if (a->at(i) != b->at(i)) {
				equal = false;
			}
			if (charCharacterIsBefore(a->at(i), b->at(i))) {
				before = true;
			}
			if (charCharacterIsBefore(b->at(i), a->at(i))) {
				done = true;
			}
		}

		if (equal) {
			if (a->size() < b->size()) {
				before = true;
			}
		}
	}

	return before;
}
vector<wchar_t>* nCreateStringScientificNotationDecimalFromNumber(double decimal) {
	StringReference* mantissaReference, * exponentReference;
	double multiplier, inc;
	double exponent;
	bool done, isPositive;
	vector<wchar_t>* result;

	mantissaReference = new StringReference();
	exponentReference = new StringReference();
	result = new vector<wchar_t>(0.0);
	done = false;
	exponent = 0.0;

	if (decimal < 0.0) {
		isPositive = false;
		decimal = -decimal;
	}
	else {
		isPositive = true;
	}

	if (decimal == 0.0) {
		done = true;
	}

	if (!done) {
		multiplier = 0.0;
		inc = 0.0;

		if (decimal < 1.0) {
			multiplier = 10.0;
			inc = -1.0;
		}
		else if (decimal >= 10.0) {
			multiplier = 0.1;
			inc = 1.0;
		}
		else {
			done = true;
		}

		if (!done) {
			for (; decimal >= 10.0 || decimal < 1.0; ) {
				decimal = decimal * multiplier;
				exponent = exponent + inc;
			}
		}
	}

	nCreateStringFromNumberWithCheck(decimal, 10.0, mantissaReference);

	nCreateStringFromNumberWithCheck(exponent, 10.0, exponentReference);

	if (!isPositive) {
		result = AppendString(result, toVector(L"-"));
	}

	result = AppendString(result, mantissaReference->string);
	result = AppendString(result, toVector(L"e"));
	result = AppendString(result, exponentReference->string);

	return result;
}
vector<wchar_t>* nCreateStringDecimalFromNumber(double decimal) {
	StringReference* stringReference;

	stringReference = new StringReference();

	/* This will succeed because base = 10. */
	nCreateStringFromNumberWithCheck(decimal, 10.0, stringReference);

	return stringReference->string;
}
bool nCreateStringFromNumberWithCheck(double decimal, double base, StringReference* stringReference) {
	vector<wchar_t>* string;
	double maximumDigits;
	double digitPosition;
	bool hasPrintedPoint, isPositive;
	double i, d;
	bool success;
	CharacterReference* characterReference;
	wchar_t c;

	isPositive = true;

	if (decimal < 0.0) {
		isPositive = false;
		decimal = -decimal;
	}

	if (decimal == 0.0) {
		stringReference->string = toVector(L"0");
		success = true;
	}
	else {
		characterReference = new CharacterReference();

		if (mathIsInteger(base)) {
			success = true;

			string = new vector<wchar_t>(0.0);

			maximumDigits = nGetMaximumDigitsForBase(base);

			digitPosition = nGetFirstDigitPosition(decimal, base);

			decimal = round(decimal * pow(base, maximumDigits - digitPosition - 1.0));

			hasPrintedPoint = false;

			if (!isPositive) {
				string = AppendCharacter(string, '-');
			}

			/* Print leading zeros. */
			if (digitPosition < 0.0) {
				string = AppendCharacter(string, '0');
				string = AppendCharacter(string, '.');
				hasPrintedPoint = true;
				for (i = 0.0; i < -digitPosition - 1.0; i = i + 1.0) {
					string = AppendCharacter(string, '0');
				}
			}

			/* Print number. */
			for (i = 0.0; i < maximumDigits && success; i = i + 1.0) {
				d = floor(decimal / pow(base, maximumDigits - i - 1.0));

				if (d >= base) {
					d = base - 1.0;
				}

				if (!hasPrintedPoint && digitPosition - i + 1.0 == 0.0) {
					if (decimal != 0.0) {
						string = AppendCharacter(string, '.');
					}
					hasPrintedPoint = true;
				}

				if (decimal == 0.0 && hasPrintedPoint) {
				}
				else {
					success = nGetSingleDigitCharacterFromNumberWithCheck(d, base, characterReference);
					if (success) {
						c = characterReference->characterValue;
						string = AppendCharacter(string, c);
					}
				}

				if (success) {
					decimal = decimal - d * pow(base, maximumDigits - i - 1.0);
				}
			}

			if (success) {
				/* Print trailing zeros. */
				for (i = 0.0; i < digitPosition - maximumDigits + 1.0; i = i + 1.0) {
					string = AppendCharacter(string, '0');
				}

				stringReference->string = string;
			}
		}
		else {
			success = false;
		}
	}

	/* Done */
	return success;
}
double nGetMaximumDigitsForBase(double base) {
	double t;

	t = pow(10.0, 15.0);
	return floor(log10(t) / log10(base));
}
double nGetFirstDigitPosition(double decimal, double base) {
	double power;
	double t;

	power = ceil(log10(decimal) / log10(base));

	t = decimal * pow(base, -power);
	if (t < base && t >= 1.0) {
	}
	else if (t >= base) {
		power = power + 1.0;
	}
	else if (t < 1.0) {
		power = power - 1.0;
	}

	return power;
}
bool nGetSingleDigitCharacterFromNumberWithCheck(double c, double base, CharacterReference* characterReference) {
	vector<wchar_t>* numberTable;
	bool success;

	numberTable = nGetDigitCharacterTable();

	if (c < base || c < numberTable->size()) {
		success = true;
		characterReference->characterValue = numberTable->at(c);
	}
	else {
		success = false;
	}

	return success;
}
vector<wchar_t>* nGetDigitCharacterTable() {
	vector<wchar_t>* numberTable;

	numberTable = toVector(L"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");

	return numberTable;
}
bool nCreateNumberFromDecimalStringWithCheck(vector<wchar_t>* string, NumberReference* decimalReference, StringReference* errorMessage) {
	return nCreateNumberFromStringWithCheck(string, 10.0, decimalReference, errorMessage);
}
double nCreateNumberFromDecimalString(vector<wchar_t>* string) {
	NumberReference* doubleReference;
	StringReference* stringReference;
	double number;

	doubleReference = CreateNumberReference(0.0);
	stringReference = CreateStringReference(toVector(L""));
	nCreateNumberFromStringWithCheck(string, 10.0, doubleReference, stringReference);
	number = doubleReference->numberValue;

	delete doubleReference;
	delete stringReference;

	return number;
}
bool nCreateNumberFromStringWithCheck(vector<wchar_t>* string, double base, NumberReference* numberReference, StringReference* errorMessage) {
	bool success;
	BooleanReference* numberIsPositive, * exponentIsPositive;
	NumberArrayReference* beforePoint, * afterPoint, * exponent;

	numberIsPositive = CreateBooleanReference(true);
	exponentIsPositive = CreateBooleanReference(true);
	beforePoint = new NumberArrayReference();
	afterPoint = new NumberArrayReference();
	exponent = new NumberArrayReference();

	if (base >= 2.0 && base <= 36.0) {
		success = nExtractPartsFromNumberString(string, base, numberIsPositive, beforePoint, afterPoint, exponentIsPositive, exponent, errorMessage);

		if (success) {
			numberReference->numberValue = nCreateNumberFromParts(base, numberIsPositive->booleanValue, beforePoint->numberArray, afterPoint->numberArray, exponentIsPositive->booleanValue, exponent->numberArray);
		}
	}
	else {
		success = false;
		errorMessage->string = toVector(L"Base must be from 2 to 36.");
	}

	return success;
}
double nCreateNumberFromParts(double base, bool numberIsPositive, vector<double>* beforePoint, vector<double>* afterPoint, bool exponentIsPositive, vector<double>* exponent) {
	double n, i, p, e;

	n = 0.0;

	for (i = 0.0; i < beforePoint->size(); i = i + 1.0) {
		p = beforePoint->at(beforePoint->size() - i - 1.0);

		n = n + p * pow(base, i);
	}

	for (i = 0.0; i < afterPoint->size(); i = i + 1.0) {
		p = afterPoint->at(i);

		n = n + p * pow(base, -(i + 1.0));
	}

	if (exponent->size() > 0.0) {
		e = 0.0;
		for (i = 0.0; i < exponent->size(); i = i + 1.0) {
			p = exponent->at(exponent->size() - i - 1.0);

			e = e + p * pow(base, i);
		}

		if (!exponentIsPositive) {
			e = -e;
		}

		n = n * pow(base, e);
	}

	if (!numberIsPositive) {
		n = -n;
	}

	return n;
}
bool nExtractPartsFromNumberString(vector<wchar_t>* n, double base, BooleanReference* numberIsPositive, NumberArrayReference* beforePoint, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages) {
	double i;
	bool success;

	i = 0.0;

	if (i < n->size()) {
		if (n->at(i) == '-') {
			numberIsPositive->booleanValue = false;
			i = i + 1.0;
		}
		else if (n->at(i) == '+') {
			numberIsPositive->booleanValue = true;
			i = i + 1.0;
		}

		success = nExtractPartsFromNumberStringFromSign(n, base, i, beforePoint, afterPoint, exponentIsPositive, exponent, errorMessages);
	}
	else {
		success = false;
		errorMessages->string = toVector(L"Number cannot have length zero.");
	}

	return success;
}
bool nExtractPartsFromNumberStringFromSign(vector<wchar_t>* n, double base, double i, NumberArrayReference* beforePoint, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages) {
	bool success, done;
	double count, j;

	done = false;
	count = 0.0;
	for (; i + count < n->size() && !done; ) {
		if (nCharacterIsNumberCharacterInBase(n->at(i + count), base)) {
			count = count + 1.0;
		}
		else {
			done = true;
		}
	}

	if (count >= 1.0) {
		beforePoint->numberArray = new vector<double>(count);

		for (j = 0.0; j < count; j = j + 1.0) {
			beforePoint->numberArray->at(j) = nGetNumberFromNumberCharacterForBase(n->at(i + j), base);
		}

		i = i + count;

		if (i < n->size()) {
			success = nExtractPartsFromNumberStringFromPointOrExponent(n, base, i, afterPoint, exponentIsPositive, exponent, errorMessages);
		}
		else {
			afterPoint->numberArray = new vector<double>(0.0);
			exponent->numberArray = new vector<double>(0.0);
			success = true;
		}
	}
	else {
		success = false;
		errorMessages->string = toVector(L"Number must have at least one number after the optional sign.");
	}

	return success;
}
bool nExtractPartsFromNumberStringFromPointOrExponent(vector<wchar_t>* n, double base, double i, NumberArrayReference* afterPoint, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages) {
	bool success, done;
	double count, j;

	if (n->at(i) == '.') {
		i = i + 1.0;

		if (i < n->size()) {
			done = false;
			count = 0.0;
			for (; i + count < n->size() && !done; ) {
				if (nCharacterIsNumberCharacterInBase(n->at(i + count), base)) {
					count = count + 1.0;
				}
				else {
					done = true;
				}
			}

			if (count >= 1.0) {
				afterPoint->numberArray = new vector<double>(count);

				for (j = 0.0; j < count; j = j + 1.0) {
					afterPoint->numberArray->at(j) = nGetNumberFromNumberCharacterForBase(n->at(i + j), base);
				}

				i = i + count;

				if (i < n->size()) {
					success = nExtractPartsFromNumberStringFromExponent(n, base, i, exponentIsPositive, exponent, errorMessages);
				}
				else {
					exponent->numberArray = new vector<double>(0.0);
					success = true;
				}
			}
			else {
				success = false;
				errorMessages->string = toVector(L"There must be at least one digit after the decimal point.");
			}
		}
		else {
			success = false;
			errorMessages->string = toVector(L"There must be at least one digit after the decimal point.");
		}
	}
	else if (base <= 14.0 && (n->at(i) == 'e' || n->at(i) == 'E')) {
		if (i < n->size()) {
			success = nExtractPartsFromNumberStringFromExponent(n, base, i, exponentIsPositive, exponent, errorMessages);
			afterPoint->numberArray = new vector<double>(0.0);
		}
		else {
			success = false;
			errorMessages->string = toVector(L"There must be at least one digit after the exponent.");
		}
	}
	else {
		success = false;
		errorMessages->string = toVector(L"Expected decimal point or exponent symbol.");
	}

	return success;
}
bool nExtractPartsFromNumberStringFromExponent(vector<wchar_t>* n, double base, double i, BooleanReference* exponentIsPositive, NumberArrayReference* exponent, StringReference* errorMessages) {
	bool success, done;
	double count, j;

	if (base <= 14.0 && (n->at(i) == 'e' || n->at(i) == 'E')) {
		i = i + 1.0;

		if (i < n->size()) {
			if (n->at(i) == '-') {
				exponentIsPositive->booleanValue = false;
				i = i + 1.0;
			}
			else if (n->at(i) == '+') {
				exponentIsPositive->booleanValue = true;
				i = i + 1.0;
			}

			if (i < n->size()) {
				done = false;
				count = 0.0;
				for (; i + count < n->size() && !done; ) {
					if (nCharacterIsNumberCharacterInBase(n->at(i + count), base)) {
						count = count + 1.0;
					}
					else {
						done = true;
					}
				}

				if (count >= 1.0) {
					exponent->numberArray = new vector<double>(count);

					for (j = 0.0; j < count; j = j + 1.0) {
						exponent->numberArray->at(j) = nGetNumberFromNumberCharacterForBase(n->at(i + j), base);
					}

					i = i + count;

					if (i == n->size()) {
						success = true;
					}
					else {
						success = false;
						errorMessages->string = toVector(L"There cannot be any characters past the exponent of the number.");
					}
				}
				else {
					success = false;
					errorMessages->string = toVector(L"There must be at least one digit after the decimal point.");
				}
			}
			else {
				success = false;
				errorMessages->string = toVector(L"There must be at least one digit after the exponent symbol.");
			}
		}
		else {
			success = false;
			errorMessages->string = toVector(L"There must be at least one digit after the exponent symbol.");
		}
	}
	else {
		success = false;
		errorMessages->string = toVector(L"Expected exponent symbol.");
	}

	return success;
}
double nGetNumberFromNumberCharacterForBase(wchar_t c, double base) {
	vector<wchar_t>* numberTable;
	double i;
	double position;

	numberTable = nGetDigitCharacterTable();
	position = 0.0;

	for (i = 0.0; i < base; i = i + 1.0) {
		if (numberTable->at(i) == c) {
			position = i;
		}
	}

	return position;
}
bool nCharacterIsNumberCharacterInBase(wchar_t c, double base) {
	vector<wchar_t>* numberTable;
	double i;
	bool found;

	numberTable = nGetDigitCharacterTable();
	found = false;

	for (i = 0.0; i < base; i = i + 1.0) {
		if (numberTable->at(i) == c) {
			found = true;
		}
	}

	return found;
}
vector<double>* nStringToNumberArray(vector<wchar_t>* str) {
	NumberArrayReference* numberArrayReference;
	StringReference* stringReference;
	vector<double>* numbers;

	numberArrayReference = new NumberArrayReference();
	stringReference = new StringReference();

	nStringToNumberArrayWithCheck(str, numberArrayReference, stringReference);

	numbers = numberArrayReference->numberArray;

	delete numberArrayReference;
	delete stringReference;

	return numbers;
}
bool nStringToNumberArrayWithCheck(vector<wchar_t>* str, NumberArrayReference* numberArrayReference, StringReference* errorMessage) {
	vector<StringReference*>* numberStrings;
	vector<double>* numbers;
	double i;
	vector<wchar_t>* numberString, * trimmedNumberString;
	bool success;
	NumberReference* numberReference;

	numberStrings = SplitByString(str, toVector(L","));

	numbers = new vector<double>(numberStrings->size());
	success = true;
	numberReference = new NumberReference();

	for (i = 0.0; i < numberStrings->size(); i = i + 1.0) {
		numberString = numberStrings->at(i)->string;
		trimmedNumberString = Trim(numberString);
		success = nCreateNumberFromDecimalStringWithCheck(trimmedNumberString, numberReference, errorMessage);
		numbers->at(i) = numberReference->numberValue;

		FreeStringReference(numberStrings->at(i));
		delete trimmedNumberString;
	}

	delete numberStrings;
	delete numberReference;

	numberArrayReference->numberArray = numbers;

	return success;
}
void AssertFalse(bool b, NumberReference* failures) {
	if (b) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertTrue(bool b, NumberReference* failures) {
	if (!b) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertEquals(double a, double b, NumberReference* failures) {
	if (a != b) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertBooleansEqual(bool a, bool b, NumberReference* failures) {
	if (a != b) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertCharactersEqual(wchar_t a, wchar_t b, NumberReference* failures) {
	if (a != b) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertStringEquals(vector<wchar_t>* a, vector<wchar_t>* b, NumberReference* failures) {
	if (!StringsEqual(a, b)) {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertNumberArraysEqual(vector<double>* a, vector<double>* b, NumberReference* failures) {
	double i;

	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size(); i = i + 1.0) {
			AssertEquals(a->at(i), b->at(i), failures);
		}
	}
	else {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertBooleanArraysEqual(vector<bool>* a, vector<bool>* b, NumberReference* failures) {
	double i;

	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size(); i = i + 1.0) {
			AssertBooleansEqual(a->at(i), b->at(i), failures);
		}
	}
	else {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
void AssertStringArraysEqual(vector<StringReference*>* a, vector<StringReference*>* b, NumberReference* failures) {
	double i;

	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size(); i = i + 1.0) {
			AssertStringEquals(a->at(i)->string, b->at(i)->string, failures);
		}
	}
	else {
		failures->numberValue = failures->numberValue + 1.0;
	}
}
vector<double>* AddNumber(vector<double>* list, double a) {
	vector<double>* newlist;
	double i;

	newlist = new vector<double>(list->size() + 1.0);
	for (i = 0.0; i < list->size(); i = i + 1.0) {
		newlist->at(i) = list->at(i);
	}
	newlist->at(list->size()) = a;

	delete list;

	return newlist;
}
void AddNumberRef(NumberArrayReference* list, double i) {
	list->numberArray = AddNumber(list->numberArray, i);
}
vector<double>* RemoveNumber(vector<double>* list, double n) {
	vector<double>* newlist;
	double i;

	newlist = new vector<double>(list->size() - 1.0);

	if (n >= 0.0 && n < list->size()) {
		for (i = 0.0; i < list->size(); i = i + 1.0) {
			if (i < n) {
				newlist->at(i) = list->at(i);
			}
			if (i > n) {
				newlist->at(i - 1.0) = list->at(i);
			}
		}

		delete list;
	}
	else {
		delete newlist;
	}

	return newlist;
}
double GetNumberRef(NumberArrayReference* list, double i) {
	return list->numberArray->at(i);
}
void RemoveNumberRef(NumberArrayReference* list, double i) {
	list->numberArray = RemoveNumber(list->numberArray, i);
}
vector<StringReference*>* AddString(vector<StringReference*>* list, StringReference* a) {
	vector<StringReference*>* newlist;
	double i;

	newlist = new vector<StringReference*>(list->size() + 1.0);

	for (i = 0.0; i < list->size(); i = i + 1.0) {
		newlist->at(i) = list->at(i);
	}
	newlist->at(list->size()) = a;

	delete list;

	return newlist;
}
void AddStringRef(StringArrayReference* list, StringReference* i) {
	list->stringArray = AddString(list->stringArray, i);
}
vector<StringReference*>* RemoveString(vector<StringReference*>* list, double n) {
	vector<StringReference*>* newlist;
	double i;

	newlist = new vector<StringReference*>(list->size() - 1.0);

	if (n >= 0.0 && n < list->size()) {
		for (i = 0.0; i < list->size(); i = i + 1.0) {
			if (i < n) {
				newlist->at(i) = list->at(i);
			}
			if (i > n) {
				newlist->at(i - 1.0) = list->at(i);
			}
		}

		delete list;
	}
	else {
		delete newlist;
	}

	return newlist;
}
StringReference* GetStringRef(StringArrayReference* list, double i) {
	return list->stringArray->at(i);
}
void RemoveStringRef(StringArrayReference* list, double i) {
	list->stringArray = RemoveString(list->stringArray, i);
}
vector<bool>* AddBoolean(vector<bool>* list, bool a) {
	vector<bool>* newlist;
	double i;

	newlist = new vector<bool>(list->size() + 1.0);
	for (i = 0.0; i < list->size(); i = i + 1.0) {
		newlist->at(i) = list->at(i);
	}
	newlist->at(list->size()) = a;

	delete list;

	return newlist;
}
void AddBooleanRef(BooleanArrayReference* list, bool i) {
	list->booleanArray = AddBoolean(list->booleanArray, i);
}
vector<bool>* RemoveBoolean(vector<bool>* list, double n) {
	vector<bool>* newlist;
	double i;

	newlist = new vector<bool>(list->size() - 1.0);

	if (n >= 0.0 && n < list->size()) {
		for (i = 0.0; i < list->size(); i = i + 1.0) {
			if (i < n) {
				newlist->at(i) = list->at(i);
			}
			if (i > n) {
				newlist->at(i - 1.0) = list->at(i);
			}
		}

		delete list;
	}
	else {
		delete newlist;
	}

	return newlist;
}
bool GetBooleanRef(BooleanArrayReference* list, double i) {
	return list->booleanArray->at(i);
}
void RemoveDecimalRef(BooleanArrayReference* list, double i) {
	list->booleanArray = RemoveBoolean(list->booleanArray, i);
}
LinkedListStrings* CreateLinkedListString() {
	LinkedListStrings* ll;

	ll = new LinkedListStrings();
	ll->first = new LinkedListNodeStrings();
	ll->last = ll->first;
	ll->last->end = true;

	return ll;
}
void LinkedListAddString(LinkedListStrings* ll, vector<wchar_t>* value) {
	ll->last->end = false;
	ll->last->value = value;
	ll->last->next = new LinkedListNodeStrings();
	ll->last->next->end = true;
	ll->last = ll->last->next;
}
vector<StringReference*>* LinkedListStringsToArray(LinkedListStrings* ll) {
	vector<StringReference*>* array;
	double length, i;
	LinkedListNodeStrings* node;

	node = ll->first;

	length = LinkedListStringsLength(ll);

	array = new vector<StringReference*>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		array->at(i) = new StringReference();
		array->at(i)->string = node->value;
		node = node->next;
	}

	return array;
}
double LinkedListStringsLength(LinkedListStrings* ll) {
	double l;
	LinkedListNodeStrings* node;

	l = 0.0;
	node = ll->first;
	for (; !node->end; ) {
		node = node->next;
		l = l + 1.0;
	}

	return l;
}
void FreeLinkedListString(LinkedListStrings* ll) {
	LinkedListNodeStrings* node, * prev;

	node = ll->first;

	for (; !node->end; ) {
		prev = node;
		node = node->next;
		delete prev;
	}

	delete node;
}
LinkedListNumbers* CreateLinkedListNumbers() {
	LinkedListNumbers* ll;

	ll = new LinkedListNumbers();
	ll->first = new LinkedListNodeNumbers();
	ll->last = ll->first;
	ll->last->end = true;

	return ll;
}
vector<LinkedListNumbers*>* CreateLinkedListNumbersArray(double length) {
	vector<LinkedListNumbers*>* lls;
	double i;

	lls = new vector<LinkedListNumbers*>(length);
	for (i = 0.0; i < lls->size(); i = i + 1.0) {
		lls->at(i) = CreateLinkedListNumbers();
	}

	return lls;
}
void LinkedListAddNumber(LinkedListNumbers* ll, double value) {
	ll->last->end = false;
	ll->last->value = value;
	ll->last->next = new LinkedListNodeNumbers();
	ll->last->next->end = true;
	ll->last = ll->last->next;
}
double LinkedListNumbersLength(LinkedListNumbers* ll) {
	double l;
	LinkedListNodeNumbers* node;

	l = 0.0;
	node = ll->first;
	for (; !node->end; ) {
		node = node->next;
		l = l + 1.0;
	}

	return l;
}
double LinkedListNumbersIndex(LinkedListNumbers* ll, double index) {
	double i;
	LinkedListNodeNumbers* node;

	node = ll->first;
	for (i = 0.0; i < index; i = i + 1.0) {
		node = node->next;
	}

	return node->value;
}
void LinkedListInsertNumber(LinkedListNumbers* ll, double index, double value) {
	double i;
	LinkedListNodeNumbers* node, * tmp;

	if (index == 0.0) {
		tmp = ll->first;
		ll->first = new LinkedListNodeNumbers();
		ll->first->next = tmp;
		ll->first->value = value;
		ll->first->end = false;
	}
	else {
		node = ll->first;
		for (i = 0.0; i < index - 1.0; i = i + 1.0) {
			node = node->next;
		}

		tmp = node->next;
		node->next = new LinkedListNodeNumbers();
		node->next->next = tmp;
		node->next->value = value;
		node->next->end = false;
	}
}
void LinkedListSet(LinkedListNumbers* ll, double index, double value) {
	double i;
	LinkedListNodeNumbers* node;

	node = ll->first;
	for (i = 0.0; i < index; i = i + 1.0) {
		node = node->next;
	}

	node->next->value = value;
}
void LinkedListRemoveNumber(LinkedListNumbers* ll, double index) {
	double i;
	LinkedListNodeNumbers* node, * prev;

	node = ll->first;
	prev = ll->first;

	for (i = 0.0; i < index; i = i + 1.0) {
		prev = node;
		node = node->next;
	}

	if (index == 0.0) {
		ll->first = prev->next;
	}
	if (!prev->next->end) {
		prev->next = prev->next->next;
	}
}
void FreeLinkedListNumbers(LinkedListNumbers* ll) {
	LinkedListNodeNumbers* node, * prev;

	node = ll->first;

	for (; !node->end; ) {
		prev = node;
		node = node->next;
		delete prev;
	}

	delete node;
}
void FreeLinkedListNumbersArray(vector<LinkedListNumbers*>* lls) {
	double i;

	for (i = 0.0; i < lls->size(); i = i + 1.0) {
		FreeLinkedListNumbers(lls->at(i));
	}
	delete lls;
}
vector<double>* LinkedListNumbersToArray(LinkedListNumbers* ll) {
	vector<double>* array;
	double length, i;
	LinkedListNodeNumbers* node;

	node = ll->first;

	length = LinkedListNumbersLength(ll);

	array = new vector<double>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		array->at(i) = node->value;
		node = node->next;
	}

	return array;
}
LinkedListNumbers* ArrayToLinkedListNumbers(vector<double>* array) {
	LinkedListNumbers* ll;
	double i;

	ll = CreateLinkedListNumbers();

	for (i = 0.0; i < array->size(); i = i + 1.0) {
		LinkedListAddNumber(ll, array->at(i));
	}

	return ll;
}
bool LinkedListNumbersEqual(LinkedListNumbers* a, LinkedListNumbers* b) {
	bool equal, done;
	LinkedListNodeNumbers* an, * bn;

	an = a->first;
	bn = b->first;

	equal = true;
	done = false;
	for (; equal && !done; ) {
		if (an->end == bn->end) {
			if (an->end) {
				done = true;
			}
			else if (an->value == bn->value) {
				an = an->next;
				bn = bn->next;
			}
			else {
				equal = false;
			}
		}
		else {
			equal = false;
		}
	}

	return equal;
}
LinkedListCharacters* CreateLinkedListCharacter() {
	LinkedListCharacters* ll;

	ll = new LinkedListCharacters();
	ll->first = new LinkedListNodeCharacters();
	ll->last = ll->first;
	ll->last->end = true;

	return ll;
}
void LinkedListAddCharacter(LinkedListCharacters* ll, wchar_t value) {
	ll->last->end = false;
	ll->last->value = value;
	ll->last->next = new LinkedListNodeCharacters();
	ll->last->next->end = true;
	ll->last = ll->last->next;
}
vector<wchar_t>* LinkedListCharactersToArray(LinkedListCharacters* ll) {
	vector<wchar_t>* array;
	double length, i;
	LinkedListNodeCharacters* node;

	node = ll->first;

	length = LinkedListCharactersLength(ll);

	array = new vector<wchar_t>(length);

	for (i = 0.0; i < length; i = i + 1.0) {
		array->at(i) = node->value;
		node = node->next;
	}

	return array;
}
double LinkedListCharactersLength(LinkedListCharacters* ll) {
	double l;
	LinkedListNodeCharacters* node;

	l = 0.0;
	node = ll->first;
	for (; !node->end; ) {
		node = node->next;
		l = l + 1.0;
	}

	return l;
}
void FreeLinkedListCharacter(LinkedListCharacters* ll) {
	LinkedListNodeCharacters* node, * prev;

	node = ll->first;

	for (; !node->end; ) {
		prev = node;
		node = node->next;
		delete prev;
	}

	delete node;
}
DynamicArrayNumbers* CreateDynamicArrayNumbers() {
	DynamicArrayNumbers* da;

	da = new DynamicArrayNumbers();
	da->array = new vector<double>(10.0);
	da->length = 0.0;

	return da;
}
DynamicArrayNumbers* CreateDynamicArrayNumbersWithInitialCapacity(double capacity) {
	DynamicArrayNumbers* da;

	da = new DynamicArrayNumbers();
	da->array = new vector<double>(capacity);
	da->length = 0.0;

	return da;
}
void DynamicArrayAddNumber(DynamicArrayNumbers* da, double value) {
	if (da->length == da->array->size()) {
		DynamicArrayNumbersIncreaseSize(da);
	}

	da->array->at(da->length) = value;
	da->length = da->length + 1.0;
}
void DynamicArrayNumbersIncreaseSize(DynamicArrayNumbers* da) {
	double newLength, i;
	vector<double>* newArray;

	newLength = round(da->array->size() * 3.0 / 2.0);
	newArray = new vector<double>(newLength);

	for (i = 0.0; i < da->array->size(); i = i + 1.0) {
		newArray->at(i) = da->array->at(i);
	}

	delete da->array;

	da->array = newArray;
}
bool DynamicArrayNumbersDecreaseSizeNecessary(DynamicArrayNumbers* da) {
	bool needsDecrease;

	needsDecrease = false;

	if (da->length > 10.0) {
		needsDecrease = da->length <= round(da->array->size() * 2.0 / 3.0);
	}

	return needsDecrease;
}
void DynamicArrayNumbersDecreaseSize(DynamicArrayNumbers* da) {
	double newLength, i;
	vector<double>* newArray;

	newLength = round(da->array->size() * 2.0 / 3.0);
	newArray = new vector<double>(newLength);

	for (i = 0.0; i < da->array->size(); i = i + 1.0) {
		newArray->at(i) = da->array->at(i);
	}

	delete da->array;

	da->array = newArray;
}
double DynamicArrayNumbersIndex(DynamicArrayNumbers* da, double index) {
	return da->array->at(index);
}
double DynamicArrayNumbersLength(DynamicArrayNumbers* da) {
	return da->length;
}
void DynamicArrayInsertNumber(DynamicArrayNumbers* da, double index, double value) {
	double i;

	if (da->length == da->array->size()) {
		DynamicArrayNumbersIncreaseSize(da);
	}

	for (i = da->length; i > index; i = i - 1.0) {
		da->array->at(i) = da->array->at(i - 1.0);
	}

	da->array->at(index) = value;

	da->length = da->length + 1.0;
}
void DynamicArraySet(DynamicArrayNumbers* da, double index, double value) {
	da->array->at(index) = value;
}
void DynamicArrayRemoveNumber(DynamicArrayNumbers* da, double index) {
	double i;

	for (i = index; i < da->length - 1.0; i = i + 1.0) {
		da->array->at(i) = da->array->at(i + 1.0);
	}

	da->length = da->length - 1.0;

	if (DynamicArrayNumbersDecreaseSizeNecessary(da)) {
		DynamicArrayNumbersDecreaseSize(da);
	}
}
void FreeDynamicArrayNumbers(DynamicArrayNumbers* da) {
	delete da->array;
	delete da;
}
vector<double>* DynamicArrayNumbersToArray(DynamicArrayNumbers* da) {
	vector<double>* array;
	double i;

	array = new vector<double>(da->length);

	for (i = 0.0; i < da->length; i = i + 1.0) {
		array->at(i) = da->array->at(i);
	}

	return array;
}
DynamicArrayNumbers* ArrayToDynamicArrayNumbersWithOptimalSize(vector<double>* array) {
	DynamicArrayNumbers* da;
	double i;
	double c, n, newCapacity;

	/*
		   c = 10*(3/2)^n
		   log(c) = log(10*(3/2)^n)
		   log(c) = log(10) + log((3/2)^n)
		   log(c) = 1 + log((3/2)^n)
		   log(c) - 1 = log((3/2)^n)
		   log(c) - 1 = n*log(3/2)
		   n = (log(c) - 1)/log(3/2)
		   */
	c = array->size();
	n = (log(c) - 1.0) / log(3.0 / 2.0);
	newCapacity = floor(n) + 1.0;

	da = CreateDynamicArrayNumbersWithInitialCapacity(newCapacity);

	for (i = 0.0; i < array->size(); i = i + 1.0) {
		da->array->at(i) = array->at(i);
	}

	return da;
}
DynamicArrayNumbers* ArrayToDynamicArrayNumbers(vector<double>* array) {
	DynamicArrayNumbers* da;

	da = new DynamicArrayNumbers();
	da->array = CopyNumberArray(array);
	da->length = array->size();

	return da;
}
bool DynamicArrayNumbersEqual(DynamicArrayNumbers* a, DynamicArrayNumbers* b) {
	bool equal;
	double i;

	equal = true;
	if (a->length == b->length) {
		for (i = 0.0; i < a->length && equal; i = i + 1.0) {
			if (a->array->at(i) != b->array->at(i)) {
				equal = false;
			}
		}
	}
	else {
		equal = false;
	}

	return equal;
}
LinkedListNumbers* DynamicArrayNumbersToLinkedList(DynamicArrayNumbers* da) {
	LinkedListNumbers* ll;
	double i;

	ll = CreateLinkedListNumbers();

	for (i = 0.0; i < da->length; i = i + 1.0) {
		LinkedListAddNumber(ll, da->array->at(i));
	}

	return ll;
}
DynamicArrayNumbers* LinkedListToDynamicArrayNumbers(LinkedListNumbers* ll) {
	DynamicArrayNumbers* da;
	double i;
	LinkedListNodeNumbers* node;

	node = ll->first;

	da = new DynamicArrayNumbers();
	da->length = LinkedListNumbersLength(ll);

	da->array = new vector<double>(da->length);

	for (i = 0.0; i < da->length; i = i + 1.0) {
		da->array->at(i) = node->value;
		node = node->next;
	}

	return da;
}
vector<wchar_t>* AddCharacter(vector<wchar_t>* list, wchar_t a) {
	vector<wchar_t>* newlist;
	double i;

	newlist = new vector<wchar_t>(list->size() + 1.0);
	for (i = 0.0; i < list->size(); i = i + 1.0) {
		newlist->at(i) = list->at(i);
	}
	newlist->at(list->size()) = a;

	delete list;

	return newlist;
}
void AddCharacterRef(StringReference* list, wchar_t i) {
	list->string = AddCharacter(list->string, i);
}
vector<wchar_t>* RemoveCharacter(vector<wchar_t>* list, double n) {
	vector<wchar_t>* newlist;
	double i;

	newlist = new vector<wchar_t>(list->size() - 1.0);

	if (n >= 0.0 && n < list->size()) {
		for (i = 0.0; i < list->size(); i = i + 1.0) {
			if (i < n) {
				newlist->at(i) = list->at(i);
			}
			if (i > n) {
				newlist->at(i - 1.0) = list->at(i);
			}
		}

		delete list;
	}
	else {
		delete newlist;
	}

	return newlist;
}
wchar_t GetCharacterRef(StringReference* list, double i) {
	return list->string->at(i);
}
void RemoveCharacterRef(StringReference* list, double i) {
	list->string = RemoveCharacter(list->string, i);
}
wchar_t charToLowerCase(wchar_t character) {
	wchar_t toReturn;

	toReturn = character;
	if (character == 'A') {
		toReturn = 'a';
	}
	else if (character == 'B') {
		toReturn = 'b';
	}
	else if (character == 'C') {
		toReturn = 'c';
	}
	else if (character == 'D') {
		toReturn = 'd';
	}
	else if (character == 'E') {
		toReturn = 'e';
	}
	else if (character == 'F') {
		toReturn = 'f';
	}
	else if (character == 'G') {
		toReturn = 'g';
	}
	else if (character == 'H') {
		toReturn = 'h';
	}
	else if (character == 'I') {
		toReturn = 'i';
	}
	else if (character == 'J') {
		toReturn = 'j';
	}
	else if (character == 'K') {
		toReturn = 'k';
	}
	else if (character == 'L') {
		toReturn = 'l';
	}
	else if (character == 'M') {
		toReturn = 'm';
	}
	else if (character == 'N') {
		toReturn = 'n';
	}
	else if (character == 'O') {
		toReturn = 'o';
	}
	else if (character == 'P') {
		toReturn = 'p';
	}
	else if (character == 'Q') {
		toReturn = 'q';
	}
	else if (character == 'R') {
		toReturn = 'r';
	}
	else if (character == 'S') {
		toReturn = 's';
	}
	else if (character == 'T') {
		toReturn = 't';
	}
	else if (character == 'U') {
		toReturn = 'u';
	}
	else if (character == 'V') {
		toReturn = 'v';
	}
	else if (character == 'W') {
		toReturn = 'w';
	}
	else if (character == 'X') {
		toReturn = 'x';
	}
	else if (character == 'Y') {
		toReturn = 'y';
	}
	else if (character == 'Z') {
		toReturn = 'z';
	}

	return toReturn;
}
wchar_t charToUpperCase(wchar_t character) {
	wchar_t toReturn;

	toReturn = character;
	if (character == 'a') {
		toReturn = 'A';
	}
	else if (character == 'b') {
		toReturn = 'B';
	}
	else if (character == 'c') {
		toReturn = 'C';
	}
	else if (character == 'd') {
		toReturn = 'D';
	}
	else if (character == 'e') {
		toReturn = 'E';
	}
	else if (character == 'f') {
		toReturn = 'F';
	}
	else if (character == 'g') {
		toReturn = 'G';
	}
	else if (character == 'h') {
		toReturn = 'H';
	}
	else if (character == 'i') {
		toReturn = 'I';
	}
	else if (character == 'j') {
		toReturn = 'J';
	}
	else if (character == 'k') {
		toReturn = 'K';
	}
	else if (character == 'l') {
		toReturn = 'L';
	}
	else if (character == 'm') {
		toReturn = 'M';
	}
	else if (character == 'n') {
		toReturn = 'N';
	}
	else if (character == 'o') {
		toReturn = 'O';
	}
	else if (character == 'p') {
		toReturn = 'P';
	}
	else if (character == 'q') {
		toReturn = 'Q';
	}
	else if (character == 'r') {
		toReturn = 'R';
	}
	else if (character == 's') {
		toReturn = 'S';
	}
	else if (character == 't') {
		toReturn = 'T';
	}
	else if (character == 'u') {
		toReturn = 'U';
	}
	else if (character == 'v') {
		toReturn = 'V';
	}
	else if (character == 'w') {
		toReturn = 'W';
	}
	else if (character == 'x') {
		toReturn = 'X';
	}
	else if (character == 'y') {
		toReturn = 'Y';
	}
	else if (character == 'z') {
		toReturn = 'Z';
	}

	return toReturn;
}
bool charIsUpperCase(wchar_t character) {
	bool isUpper;

	isUpper = false;
	if (character == 'A') {
		isUpper = true;
	}
	else if (character == 'B') {
		isUpper = true;
	}
	else if (character == 'C') {
		isUpper = true;
	}
	else if (character == 'D') {
		isUpper = true;
	}
	else if (character == 'E') {
		isUpper = true;
	}
	else if (character == 'F') {
		isUpper = true;
	}
	else if (character == 'G') {
		isUpper = true;
	}
	else if (character == 'H') {
		isUpper = true;
	}
	else if (character == 'I') {
		isUpper = true;
	}
	else if (character == 'J') {
		isUpper = true;
	}
	else if (character == 'K') {
		isUpper = true;
	}
	else if (character == 'L') {
		isUpper = true;
	}
	else if (character == 'M') {
		isUpper = true;
	}
	else if (character == 'N') {
		isUpper = true;
	}
	else if (character == 'O') {
		isUpper = true;
	}
	else if (character == 'P') {
		isUpper = true;
	}
	else if (character == 'Q') {
		isUpper = true;
	}
	else if (character == 'R') {
		isUpper = true;
	}
	else if (character == 'S') {
		isUpper = true;
	}
	else if (character == 'T') {
		isUpper = true;
	}
	else if (character == 'U') {
		isUpper = true;
	}
	else if (character == 'V') {
		isUpper = true;
	}
	else if (character == 'W') {
		isUpper = true;
	}
	else if (character == 'X') {
		isUpper = true;
	}
	else if (character == 'Y') {
		isUpper = true;
	}
	else if (character == 'Z') {
		isUpper = true;
	}

	return isUpper;
}
bool charIsLowerCase(wchar_t character) {
	bool isLower;

	isLower = false;
	if (character == 'a') {
		isLower = true;
	}
	else if (character == 'b') {
		isLower = true;
	}
	else if (character == 'c') {
		isLower = true;
	}
	else if (character == 'd') {
		isLower = true;
	}
	else if (character == 'e') {
		isLower = true;
	}
	else if (character == 'f') {
		isLower = true;
	}
	else if (character == 'g') {
		isLower = true;
	}
	else if (character == 'h') {
		isLower = true;
	}
	else if (character == 'i') {
		isLower = true;
	}
	else if (character == 'j') {
		isLower = true;
	}
	else if (character == 'k') {
		isLower = true;
	}
	else if (character == 'l') {
		isLower = true;
	}
	else if (character == 'm') {
		isLower = true;
	}
	else if (character == 'n') {
		isLower = true;
	}
	else if (character == 'o') {
		isLower = true;
	}
	else if (character == 'p') {
		isLower = true;
	}
	else if (character == 'q') {
		isLower = true;
	}
	else if (character == 'r') {
		isLower = true;
	}
	else if (character == 's') {
		isLower = true;
	}
	else if (character == 't') {
		isLower = true;
	}
	else if (character == 'u') {
		isLower = true;
	}
	else if (character == 'v') {
		isLower = true;
	}
	else if (character == 'w') {
		isLower = true;
	}
	else if (character == 'x') {
		isLower = true;
	}
	else if (character == 'y') {
		isLower = true;
	}
	else if (character == 'z') {
		isLower = true;
	}

	return isLower;
}
bool charIsLetter(wchar_t character) {
	return charIsUpperCase(character) || charIsLowerCase(character);
}
bool charIsNumber(wchar_t character) {
	bool isNumberx;

	isNumberx = false;
	if (character == '0') {
		isNumberx = true;
	}
	else if (character == '1') {
		isNumberx = true;
	}
	else if (character == '2') {
		isNumberx = true;
	}
	else if (character == '3') {
		isNumberx = true;
	}
	else if (character == '4') {
		isNumberx = true;
	}
	else if (character == '5') {
		isNumberx = true;
	}
	else if (character == '6') {
		isNumberx = true;
	}
	else if (character == '7') {
		isNumberx = true;
	}
	else if (character == '8') {
		isNumberx = true;
	}
	else if (character == '9') {
		isNumberx = true;
	}

	return isNumberx;
}
bool charIsWhiteSpace(wchar_t character) {
	bool isWhiteSpacex;

	isWhiteSpacex = false;
	if (character == ' ') {
		isWhiteSpacex = true;
	}
	else if (character == '\t') {
		isWhiteSpacex = true;
	}
	else if (character == '\n') {
		isWhiteSpacex = true;
	}
	else if (character == '\r') {
		isWhiteSpacex = true;
	}

	return isWhiteSpacex;
}
bool charIsSymbol(wchar_t character) {
	bool isSymbolx;

	isSymbolx = false;
	if (character == '!') {
		isSymbolx = true;
	}
	else if (character == '\"') {
		isSymbolx = true;
	}
	else if (character == '#') {
		isSymbolx = true;
	}
	else if (character == '$') {
		isSymbolx = true;
	}
	else if (character == '%') {
		isSymbolx = true;
	}
	else if (character == '&') {
		isSymbolx = true;
	}
	else if (character == '\'') {
		isSymbolx = true;
	}
	else if (character == '(') {
		isSymbolx = true;
	}
	else if (character == ')') {
		isSymbolx = true;
	}
	else if (character == '*') {
		isSymbolx = true;
	}
	else if (character == '+') {
		isSymbolx = true;
	}
	else if (character == ',') {
		isSymbolx = true;
	}
	else if (character == '-') {
		isSymbolx = true;
	}
	else if (character == '.') {
		isSymbolx = true;
	}
	else if (character == '/') {
		isSymbolx = true;
	}
	else if (character == ':') {
		isSymbolx = true;
	}
	else if (character == ';') {
		isSymbolx = true;
	}
	else if (character == '<') {
		isSymbolx = true;
	}
	else if (character == '=') {
		isSymbolx = true;
	}
	else if (character == '>') {
		isSymbolx = true;
	}
	else if (character == '?') {
		isSymbolx = true;
	}
	else if (character == '@') {
		isSymbolx = true;
	}
	else if (character == '[') {
		isSymbolx = true;
	}
	else if (character == '\\') {
		isSymbolx = true;
	}
	else if (character == ']') {
		isSymbolx = true;
	}
	else if (character == '^') {
		isSymbolx = true;
	}
	else if (character == '_') {
		isSymbolx = true;
	}
	else if (character == '`') {
		isSymbolx = true;
	}
	else if (character == '{') {
		isSymbolx = true;
	}
	else if (character == '|') {
		isSymbolx = true;
	}
	else if (character == '}') {
		isSymbolx = true;
	}
	else if (character == '~') {
		isSymbolx = true;
	}

	return isSymbolx;
}
bool charCharacterIsBefore(wchar_t a, wchar_t b) {
	double ad, bd;

	ad = a;
	bd = b;

	return ad < bd;
}
vector<double>* StringToNumberArray(vector<wchar_t>* string) {
	double i;
	vector<double>* array;

	array = new vector<double>(string->size());

	for (i = 0.0; i < string->size(); i = i + 1.0) {
		array->at(i) = string->at(i);
	}
	return array;
}
vector<wchar_t>* NumberArrayToString(vector<double>* array) {
	double i;
	vector<wchar_t>* string;

	string = new vector<wchar_t>(array->size());

	for (i = 0.0; i < array->size(); i = i + 1.0) {
		string->at(i) = array->at(i);
	}
	return string;
}
bool NumberArraysEqual(vector<double>* a, vector<double>* b) {
	bool equal;
	double i;

	equal = true;
	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size() && equal; i = i + 1.0) {
			if (a->at(i) != b->at(i)) {
				equal = false;
			}
		}
	}
	else {
		equal = false;
	}

	return equal;
}
bool BooleanArraysEqual(vector<bool>* a, vector<bool>* b) {
	bool equal;
	double i;

	equal = true;
	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size() && equal; i = i + 1.0) {
			if (a->at(i) != b->at(i)) {
				equal = false;
			}
		}
	}
	else {
		equal = false;
	}

	return equal;
}
bool StringsEqual(vector<wchar_t>* a, vector<wchar_t>* b) {
	bool equal;
	double i;

	equal = true;
	if (a->size() == b->size()) {
		for (i = 0.0; i < a->size() && equal; i = i + 1.0) {
			if (a->at(i) != b->at(i)) {
				equal = false;
			}
		}
	}
	else {
		equal = false;
	}

	return equal;
}
void FillNumberArray(vector<double>* a, double value) {
	double i;

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		a->at(i) = value;
	}
}
void FillString(vector<wchar_t>* a, wchar_t value) {
	double i;

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		a->at(i) = value;
	}
}
void FillBooleanArray(vector<bool>* a, bool value) {
	double i;

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		a->at(i) = value;
	}
}
bool FillNumberArrayRange(vector<double>* a, double value, double from, double to) {
	double i, length;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		for (i = 0.0; i < length; i = i + 1.0) {
			a->at(from + i) = value;
		}

		success = true;
	}
	else {
		success = false;
	}

	return success;
}
bool FillBooleanArrayRange(vector<bool>* a, bool value, double from, double to) {
	double i, length;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		for (i = 0.0; i < length; i = i + 1.0) {
			a->at(from + i) = value;
		}

		success = true;
	}
	else {
		success = false;
	}

	return success;
}
bool FillStringRange(vector<wchar_t>* a, wchar_t value, double from, double to) {
	double i, length;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		for (i = 0.0; i < length; i = i + 1.0) {
			a->at(from + i) = value;
		}

		success = true;
	}
	else {
		success = false;
	}

	return success;
}
vector<double>* CopyNumberArray(vector<double>* a) {
	double i;
	vector<double>* n;

	n = new vector<double>(a->size());

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		n->at(i) = a->at(i);
	}

	return n;
}
vector<bool>* CopyBooleanArray(vector<bool>* a) {
	double i;
	vector<bool>* n;

	n = new vector<bool>(a->size());

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		n->at(i) = a->at(i);
	}

	return n;
}
vector<wchar_t>* CopyString(vector<wchar_t>* a) {
	double i;
	vector<wchar_t>* n;

	n = new vector<wchar_t>(a->size());

	for (i = 0.0; i < a->size(); i = i + 1.0) {
		n->at(i) = a->at(i);
	}

	return n;
}
bool CopyNumberArrayRange(vector<double>* a, double from, double to, NumberArrayReference* copyReference) {
	double i, length;
	vector<double>* n;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		n = new vector<double>(length);

		for (i = 0.0; i < length; i = i + 1.0) {
			n->at(i) = a->at(from + i);
		}

		copyReference->numberArray = n;
		success = true;
	}
	else {
		success = false;
	}

	return success;
}
bool CopyBooleanArrayRange(vector<bool>* a, double from, double to, BooleanArrayReference* copyReference) {
	double i, length;
	vector<bool>* n;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		n = new vector<bool>(length);

		for (i = 0.0; i < length; i = i + 1.0) {
			n->at(i) = a->at(from + i);
		}

		copyReference->booleanArray = n;
		success = true;
	}
	else {
		success = false;
	}

	return success;
}
bool CopyStringRange(vector<wchar_t>* a, double from, double to, StringReference* copyReference) {
	double i, length;
	vector<wchar_t>* n;
	bool success;

	if (from >= 0.0 && from <= a->size() && to >= 0.0 && to <= a->size() && from <= to) {
		length = to - from;
		n = new vector<wchar_t>(length);

		for (i = 0.0; i < length; i = i + 1.0) {
			n->at(i) = a->at(from + i);
		}

		copyReference->string = n;
		success = true;
	}
	else {
		success = false;
	}

	return success;
}
bool IsLastElement(double length, double index) {
	return index + 1.0 == length;
}
vector<double>* CreateNumberArray(double length, double value) {
	vector<double>* array;

	array = new vector<double>(length);
	FillNumberArray(array, value);

	return array;
}
vector<bool>* CreateBooleanArray(double length, bool value) {
	vector<bool>* array;

	array = new vector<bool>(length);
	FillBooleanArray(array, value);

	return array;
}
vector<wchar_t>* CreateString(double length, wchar_t value) {
	vector<wchar_t>* array;

	array = new vector<wchar_t>(length);
	FillString(array, value);

	return array;
}
void SwapElementsOfNumberArray(vector<double>* A, double ai, double bi) {
	double tmp;

	tmp = A->at(ai);
	A->at(ai) = A->at(bi);
	A->at(bi) = tmp;
}
void SwapElementsOfStringArray(StringArrayReference* A, double ai, double bi) {
	StringReference* tmp;

	tmp = A->stringArray->at(ai);
	A->stringArray->at(ai) = A->stringArray->at(bi);
	A->stringArray->at(bi) = tmp;
}
void ReverseNumberArray(vector<double>* array) {
	double i;

	for (i = 0.0; i < array->size() / 2.0; i = i + 1.0) {
		SwapElementsOfNumberArray(array, i, array->size() - i - 1.0);
	}
}

