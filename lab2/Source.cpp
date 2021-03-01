#include <iostream>
#include <windows.h>
#include <immintrin.h>
#include <functional>
#include <sysinfoapi.h>
#include <math.h>

using namespace std;

#define l_VAL 8
#define m_VAL 8
#define n_VAL 8

int block_n;
int block_m;
int block_l;

struct Interval
{
	int start;
	int end;
};

struct CustomSize
{
	Interval width;
	Interval height;
};

class Matrix;

class SmallMatrix
{
	friend class Matrix;
public:
	double** data;
	int row;
	int column;
public:
	SmallMatrix(int row, int column)
	{
		this->data = new double* [row];
		for (int i = 0; i < row; i++)
		{
			this->data[i] = new double[column];
			for (int j = 0; j < column; j++)
				this->data[i][j] = 0.0;
		}
		this->row = row;
		this->column = column;
	}
	~SmallMatrix()
	{
		for (int i = 0; i < row; i++)
			delete[]data[i];
		delete[] data;
	}
	void InitSmallMatrix(int row, int column)
	{
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < column; j++)
			{
				double first = (rand() % 2000) - 1000;
				double second = rand() % 1000;
				this->data[i][j] = first + second / 1000;
			}
		}
	}
	SmallMatrix(SmallMatrix& obj)
	{
		this->column = obj.column;
		this->row = obj.row;
		this->data = new double* [this->row];
		for (int i = 0; i < this->row; i++)
		{
			this->data[i] = new double[this->column];
			for (int j = 0; j < this->column; j++)
				this->data[i][j] = obj.data[i][j];
		}
	}
	void setValueInMatrix(int row, int column, double value)
	{
		this->data[row][column] = value;
	}
	double getValueInMatrix(int row, int column)
	{
		return this->data[row][column];
	}
	friend void MulMatricesWithIntrin(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC);
	friend void SumArrayOfMatricesWithIntrin(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result);
	friend void MulMatrices(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC);
	friend void SumArrayOfMatrices(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result);
	friend void MulMatrix(function<void(SmallMatrix&, SmallMatrix&, SmallMatrix&)>Mulfunction, function<void(SmallMatrix**, long, SmallMatrix&)>SumArrayfunction, Matrix& objA, Matrix& objB, Matrix& objC);

};

class Matrix
{
	friend class SmallMatrix;
public:
	SmallMatrix*** matrix;
	int row_matrix;
	int column_matrix;
	int row_smallmatrix;
	int column_smallmatrix;

public:
	Matrix(int row_m, int column_m, int row_s, int column_s)
	{
		this->matrix = new SmallMatrix * *[row_m];
		for (int i = 0; i < row_m; i++)
		{
			this->matrix[i] = new SmallMatrix * [column_m];
			for (int j = 0; j < column_m; j++)
				this->matrix[i][j] = new SmallMatrix(row_s, column_s);
		}
		this->row_matrix = row_m;
		this->column_matrix = column_m;
		this->row_smallmatrix = row_s;
		this->column_smallmatrix = column_s;
	}
	void InitMatrix()
	{
		for (int i = 0; i < row_matrix; i++)
		{
			for (int j = 0; j < column_matrix; j++)
				matrix[i][j]->InitSmallMatrix(row_smallmatrix, column_smallmatrix);
		}
	}
	Matrix(Matrix& obj)
	{
		this->matrix = new SmallMatrix * *[obj.row_matrix];
		for (int i = 0; i < obj.row_matrix; i++)
		{
			this->matrix[i] = new SmallMatrix * [obj.column_matrix];
			for (int j = 0; j < obj.column_matrix; j++)
			{
				this->matrix[i][j] = new SmallMatrix(obj.row_smallmatrix, obj.column_smallmatrix);
				for (int x = 0; x < obj.row_smallmatrix; x++)
				{
					for (int y = 0; y < obj.column_smallmatrix; y++)
						this->matrix[i][j]->data[x][y] = obj.matrix[i][j]->data[x][y];
				}
			}
		}
		this->row_matrix = obj.row_matrix;
		this->column_matrix = obj.column_matrix;
		this->row_smallmatrix = obj.row_smallmatrix;
		this->column_smallmatrix = obj.column_smallmatrix;
	}
	bool operator ==(Matrix& obj)
	{
		if (obj.row_matrix == this->row_matrix && this->column_matrix == obj.column_matrix && this->row_smallmatrix == obj.row_smallmatrix && this->column_smallmatrix == obj.column_smallmatrix)
		{
			for (int i = 0; i < obj.row_matrix; i++)
			{
				for (int j = 0; j < obj.column_matrix; j++)
				{
					for (int x = 0; x < obj.row_smallmatrix; x++)
					{
						for (int y = 0; y < obj.column_smallmatrix; y++)
						{
							long compareable1 = (long)(this->matrix[i][j]->data[x][y] * 100.0), compareable2 = ((long)(obj.matrix[i][j]->data[x][y] * 100.0));
							if (compareable1 - compareable2 < -1 && compareable1 - compareable2>1)
							{
								return false;
							}
						}
					}
				}
			}
			return true;
		}
		else
			return false;
	}
	SmallMatrix& getSmallMatrix(int row, int column) { return *(this->matrix[row][column]); }
	void setSmallMatrix(SmallMatrix& obj, int row, int column)
	{
		for (int i = 0; i < obj.row; i++)
		{
			for (int j = 0; j < obj.column; j++)
				this->matrix[row][column]->data[i][j] = obj.data[i][j];
		}
	}

	~Matrix()
	{
		for (int i = 0; i < row_matrix; i++)
		{
			for (int j = 0; j < column_matrix; j++)
			{
				for (int x = 0; x < row_smallmatrix; x++)
					delete[] matrix[i][j]->data[x];
				delete[] matrix[i][j]->data;
			}
			delete[] matrix[i];
		}
		delete[] matrix;
	}
	friend void MulMatricesWithIntrin(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC);
	friend void SumArrayOfMatricesWithIntrin(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result);
	friend void MulMatrices(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC);
	friend void SumArrayOfMatrices(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result);
	friend void MulMatrix(function<void(SmallMatrix&, SmallMatrix&, SmallMatrix&)>Mulfunction, function<void(SmallMatrix**, long, SmallMatrix&)>SumArrayfunction, Matrix& objA, Matrix& objB, Matrix& objC);
};





void MulMatricesWithIntrin(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC)
{
	for (int l = 0; l < l_VAL; l++)
	{
		double bufferWithZero = 0.0;
		for (int n = 0; n < n_VAL; n += 4)
		{
			__m256d c = _mm256_broadcast_sd(&bufferWithZero);

			for (int m = 0; m < m_VAL; m++)
			{
				__m256d b = _mm256_load_pd(&B.data[m][n]);
				__m256d a = _mm256_broadcast_sd(&A.data[l][m]);
				c = _mm256_fmadd_pd(a, b, c);
			}

			_mm256_store_pd(&resultC.data[l][n], c);
		}
	}
}

void SumMatricesWithIntrin(SmallMatrix& summand, SmallMatrix& result)
{
	for (int l = 0; l < l_VAL; l++)
	{
		for (int n = 0; n < n_VAL; n += 4)
		{
			__m256d c = _mm256_load_pd(&result.data[l][n]);
			__m256d a = _mm256_load_pd(&(summand.data[l][n]));//
			c = _mm256_add_pd(a, c);

			_mm256_store_pd(&result.data[l][n], c);
		}
	}

}

void SumArrayOfMatricesWithIntrin(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result)
{
	for (int l = 0; l < l_VAL; l++)
	{
		double bufferWithZero = 0.0;
		for (int n = 0; n < n_VAL; n += 4)
		{
			__m256d c = _mm256_broadcast_sd(&bufferWithZero);
			for (int i = 0; i < numberOfSummands; i++)
			{
				__m256d a = _mm256_load_pd(&(arrayOfSummands[i]->data[l][n]));//
				c = _mm256_add_pd(a, c);
			}

			_mm256_store_pd(&result.data[l][n], c);
		}
	}
}

void MulMatrices(SmallMatrix& A, SmallMatrix& B, SmallMatrix& resultC)
{

	for (int l = 0; l < l_VAL; l++)
	{
		double bufferWithZero = 0.0;
		for (int i = 0; i < m_VAL; i++)
		{
			resultC.data[l][i] = bufferWithZero;
		}

		for (int m = 0; m < m_VAL; m++)
		{
			double a = A.data[l][m];
			for (int i = 0; i < n_VAL; i++)
			{
				resultC.data[l][i] += a * B.data[m][i];
			}
		}

	}

}

void SumMatrices(SmallMatrix& summand, SmallMatrix& result)
{
	for (int l = 0; l < l_VAL; l++)
	{
		double** dataPointer = summand.data;
		//#pragma loop(no_vector)
		for (int n = 0; n < n_VAL; n++)
		{
			result.data[l][n] += dataPointer[l][n];
		}
	}
}

void SumArrayOfMatrices(SmallMatrix** arrayOfSummands, long numberOfSummands, SmallMatrix& result)
{
	for (int l = 0; l < l_VAL; l++)
	{
		for (int n = 0; n < n_VAL; n++)
		{
			result.data[l][n] = 0.0;
		}

		for (int i = 0; i < numberOfSummands; i++)
		{
			double** dataPointer = arrayOfSummands[i]->data;
			//#pragma loop(no_vector)
			for (int n = 0; n < n_VAL; n++)
			{
				result.data[l][n] += dataPointer[l][n];
			}
		}
	}
}

void MulCustomMatrix(Matrix& A, Matrix& B, Matrix& C, Interval lSize, Interval nSize, Interval mSize)
{
	SmallMatrix* bufferMatrix = new SmallMatrix(l_VAL, n_VAL);
	for (int l = lSize.start; l <= lSize.end; l++)
	{
		for (int n = nSize.start; n <= nSize.end; n++)
		{
			for (int m = mSize.start; m <= mSize.end; m++)
			{
				MulMatricesWithIntrin(*A.matrix[l][m], *B.matrix[m][n], *bufferMatrix);
				SumMatricesWithIntrin(*bufferMatrix, *C.matrix[l][n]);
			}
		}
	}
}

void MulMatrixWithBlocks(Matrix& objA, Matrix& objB, Matrix& objC)
{
	for (int l = 0; l < objA.row_matrix; l += block_l)
	{
		for (int n = 0; n < objC.column_matrix; n += block_n)
		{
			for (int m = 0; m < objB.row_matrix; m += block_m)
			{
				Interval l_interv{ l, l + block_l - 1 }, m_interv{ m, m + block_m - 1 }, n_interv{ n, n + block_n - 1 };

				MulCustomMatrix(objA, objB, objC, l_interv, n_interv, m_interv);
			}
		}
	}
}

void MulMatrix(function<void(SmallMatrix&, SmallMatrix&, SmallMatrix&)>Mulfunction, function<void(SmallMatrix**, long, SmallMatrix&)>SumArrayfunction, Matrix& objA, Matrix& objB, Matrix& objC)
{
	SmallMatrix result(l_VAL, n_VAL);
	SmallMatrix** buffers = new SmallMatrix * [objB.row_matrix];

	for (int M = 0; M < objB.row_matrix; M++)
		buffers[M] = new SmallMatrix(objC.row_smallmatrix, objC.column_smallmatrix);

	for (int L = 0; L < objC.row_matrix; L++)
	{
		for (int N = 0; N < objC.column_matrix; N++)
		{
			for (int M = 0; M < objB.row_matrix; M++)
				Mulfunction(objA.getSmallMatrix(L, M), objB.getSmallMatrix(M, N), *buffers[M]);

			SumArrayfunction(buffers, objB.row_matrix, result);
			objC.setSmallMatrix(result, L, N);
		}
	}

	for (int i = 0; i < objB.row_matrix; i++)
	{
		for (int y = 0; y < objC.row_smallmatrix; y++)
			delete[] buffers[i]->data[y];
		delete[] buffers[i]->data;
	}

}

int main()
{
	srand(NULL);
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	int l, m, n;
	cout << "Введите параметр L для матрицы матриц" << endl;
	cin >> l;
	cout << "Введите параметр M для матрицы матриц" << endl;
	cin >> m;
	cout << "Введите параметр N для матрицы матриц" << endl;
	cin >> n;

	cout << "Введите параметр l для блока" << endl;
	cin >> block_l;
	cout << "Введите параметр m для блока" << endl;
	cin >> block_m;
	cout << "Введите параметр n для блока" << endl;
	cin >> block_n;

	Matrix objA(l, m, l_VAL, m_VAL);
	Matrix objB(m, n, m_VAL, n_VAL);
	Matrix objCManual(l, n, l_VAL, n_VAL);
	Matrix objCAuto(l, n, l_VAL, n_VAL);
	Matrix objCCache(l, n, l_VAL, n_VAL);
	objA.InitMatrix();
	objB.InitMatrix();

	//MulMatrixWithBlocks(objA, objB, objCCache);

	ULONGLONG startTime, endTime, timeOfExecutionAuto, timeOfExecutionManual, timeOfExecutionCacheFriendly;
	startTime = GetTickCount64();

	//умножение без ручной векторизации

	MulMatrix(MulMatrices, SumArrayOfMatrices, objA, objB, objCAuto);


	endTime = GetTickCount64();
	timeOfExecutionAuto = endTime - startTime;
	cout << "Время работы без использованием ручной векторизации: " << timeOfExecutionAuto << endl;


	startTime = GetTickCount64();

	//умножение с ручной векторизацией

	MulMatrix(MulMatricesWithIntrin, SumArrayOfMatricesWithIntrin, objA, objB, objCManual);

	endTime = GetTickCount64();
	timeOfExecutionManual = endTime - startTime;
	cout << "Время работы с использованием ручной векторизации: " << timeOfExecutionManual << endl;

	startTime = GetTickCount64();

	//умножение с ручной векторизацией

	//MulMatrix(MulMatricesWithIntrin, SumArrayOfMatricesWithIntrin, objA, objB, objCManual);
	MulMatrixWithBlocks(objA, objB, objCCache);

	endTime = GetTickCount64();
	timeOfExecutionCacheFriendly = endTime - startTime;
	cout << "Время работы с при рационализации хранения в кэше: " << timeOfExecutionCacheFriendly << endl;

	if (objCAuto == objCManual && objCAuto == objCCache)
	{
		cout << "Матрицы равны " << endl;
	}
	else
	{
		cout << "Матрицы не равны " << endl;
	}

	if (timeOfExecutionManual == 0)
		return 0;

	cout << "Коэффициент Auto/Manual:" << (double)timeOfExecutionAuto / timeOfExecutionManual << endl;


	return 0;
}