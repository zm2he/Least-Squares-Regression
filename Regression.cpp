//#include "stdafx.h"
#include <iostream>

/*
	This program performs polynomial least-squares regressions for data sets of ordered x and y pairs where y is dependent upon x.

	A coefficient matrix c for the regressed polynomial of best fit is calculated using the formula c = ((Vt)(V))^-1(Vt)(y) 
	where V is the vandermonde matrix calculated from the column vector of x values, Vt is the transpose of the vandermonde matrix, 
	and y is the columrn vector of the y values. 

	The first entry of the coefficient matrix will be the coefficient of the highest order term, the second, that of the next highest term
	all the way until the last entry, which will be that of the constant term.
*/


/*
	Matrices are represented by a Matrix object which consists of the row and column dimensions 
	and a 1-dimensional array tat stores matrix terms in row major order
*/

// Matrix Class
class Matrix
{
public:
	Matrix(int r, int c, double *a);
	~Matrix();
	int getRow();
	int getCol();
	double* getArray();

private:
	int row;
	int col;
	double* array;
};

// Constructor of the Matrix object
Matrix::Matrix(int r, int c, double *a)
{
	row = r;
	col = c;
	array = a;
}

Matrix::~Matrix() {}

// Returns row dimension of Matrix object
int Matrix::getRow() 
{ 
	return row; 
}

// Returns column dimension of Matrix object
int Matrix::getCol() 
{ 
	return col; 
}

// Returns array storing matrix
double* Matrix::getArray() 
{ 
	return array; 
}

// Function accepts an array storing a matrix and row and column dimensions and prints out the array 
void print(double *a, int row, int col)
{
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
		{
			std::cout << " " << a[r * col + c] << " ";
		}
		std::cout << "\n";
	}
}

// Function accepts an an array storing a column matrix and row and column dimensions and returns a pointer to an array
// storing the vandermonde matrix for that array
double* vandermonde(double *x, int size, int order)
{
	double *vandermonde = new double[size * (order + 1)];

	for (int r = 0; r < size; r++)
	{
		double entry = 1;

		// start from the end of each row in the vandermonde matrix 
		// and move to the beginning, increasing power of the x value by 1
		// each time
		for (int c = order; c >= 0; c--)
		{
			vandermonde[r * (order + 1) + c] = entry;
			entry = entry * x[r];
		}
	}

	return vandermonde;
}

// Function accepts an array storing a matrix and row and column dimensions and returns a pointer to an array
// storing the transpose of that array
double* transpose(double *a, int row, int col)
{
	double *transpose = new double[col * row];
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
		{
			transpose[c * row + r] = a[r * col + c];
		}
	}
	return transpose;
}

// Function accepts two matrices, a and b, represented by arrays and dimensions (array of a, row a, column a, array of b, row b, column b)
// and returns a pointer to an array storing the the product of a and b (that is (a)(b), where a is on the left)
double* multiply(double *a, int ra, int ca, double *b, int rb, int cb)
{
	double *product = new double[ra * cb];
	for (int r = 0; r < ra; r++)
	{
		for (int c = 0; c < cb; c++)
		{
			product[r*cb + c] = 0;
			for (int k = 0; k < ca; k++)
			{
				product[r * cb + c] += (a[r * ca + k]) * (b[c + cb * k]);
			}
		}
	}
	return product;
}


// Function accepts row and column dimensions and returns a pointer to an array storing a matrix
// where all entries are zero, exept for entries on the principal diagonal, which are one
double *identity(int row, int col)
{
	int PrincipalDiagonal = 0; // position of principal diagonal
	double *identity = new double[row*col];
	for (int r = 0; r < row; r++)
	{
		for (int c = 0; c < col; c++)
		{
			if (c == PrincipalDiagonal)
			{
				identity[r * col + c] = 1;
			}
			else
			{
				identity[r * col + c] = 0;
			}
		}
		PrincipalDiagonal++;
	}
	return identity;
}


// Function accepts an array stroring a square matrix and its dimension and returns the determinant of 
// the matrix through recursive calculation
double determinant(double *a, int n)
{
	double det = 0;
	
	// base case for 1x1
	if (n == 1)
	{
		det = 1;
	}

	// base case for 2x2
	if (n == 2)
	{
		det = (a[0] * a[3]) - (a[1] * a[2]);
	}

	// recursive case
	else
	{
		double *submatrix = new double[(n - 1) * (n - 1)];
		
		// perform the Laplacian expansion along the top row
		for (int x = 0; x < n; x++)
		{
			// -1 multiplier
			double term = 1;
			if (x % 2 == 1)
			{
				term = term * (-1);
			}

			// top row term
			term = term * a[x];

			// copy out submatrix
			int i = 0;
			for (int r = 1; r < n; r++)
			{
				for (int c = 0; c < n; c++)
				{
					if (c != x)
					{
						submatrix[i] = a[r * n + c];
						i++;
					}
				}
			}

			// take determinant of submatrix
			term = term * determinant(submatrix, n - 1);
			//std::cout << term;
			//std::cout << "\n";
			det = det + term;
		}
		delete[] submatrix;
	}
	return det;
}

// Function accepts an array storing a square matrix and its dimension and returns a pointer to an array storing
// its matrix of cofactors
double* matrix_of_cofactors(double *a, int n)
{
	double *cofactor_matrix = new double[n*n];

	//case for dimension 1x1
	if (n == 1)
	{
		cofactor_matrix[0] = a[0];
	}

	//case for dimension 2x2
	else if (n == 2)
	{
		cofactor_matrix[0] = a[3];
		cofactor_matrix[1] = -1 * a[2];
		cofactor_matrix[2] = -1 * a[1];
		cofactor_matrix[3] = a[0];
	}

	// larger dimensions
	else
	{
		double *submatrix = new double[(n - 1)*(n - 1)];
		
		// loop through each entry
		for (int r = 0; r < n; r++)
		{
			for (int c = 0; c < n; c++)
			{
				int i = 0;

				// copy out submatrix
				for (int sr = 0; sr < n; sr++)
				{
					for (int sc = 0; sc < n; sc++)
					{
						if (sr != r && sc != c)
						{
							submatrix[i] = a[sr * n + sc];
							i++;
						}
					}
				}

				// calculate determinant
				cofactor_matrix[r * n + c] = determinant(submatrix, n - 1);
				
				// add -1 exponent to find cofactor
				if ((r + c + 2) % 2 == 1)
				{
					cofactor_matrix[r * n + c] *= -1;
				}
			}
		}
		delete[] submatrix;
	}
	return cofactor_matrix;
}

// Function accepts a square matrix stored in an array and its dimension and returns its inverse
double* inverse(double *a, int n)
{
	double det = determinant(a, n);
	// return null if the determinant is 0 and the matrix is not invertible
	if (det == 0)
	{
		std::cout << "The Matrix not invertible. \n";
		return nullptr;
	}
	else
	{
		// use inverse = 1/determinant*transpose(matrix of cofactors)
		double *cofactor_matrix = matrix_of_cofactors(a, n);
		double *inv = transpose(cofactor_matrix, n, n);
		delete[] cofactor_matrix;
		for (int x = 0; x < n * n; x++)
		{
			inv[x] = inv[x] / det;
		}

		return inv;
	}
}


int main()
{
	// get number of data points and order of polynomial
	int numberPoints, order;
	std::cout << "Number of data points:\n";
	std::cin >> numberPoints;
	std::cout << "Order of regressed polynomial:\n";
	std::cin >> order;
	
	// read in data
	double *xValues = new double [numberPoints];
	double *yValues = new double[numberPoints];

	
	std::cout << "Enter all x values:\n";
	for (int count = 0; count < numberPoints; count++)
	{
		double a;
		std::cin >> a;
		xValues[count] = a;
	}

	std::cout << "Enter all y values:\n";
	for (int count = 0; count < numberPoints; count++)
	{
		double a;
		std::cin >> a;
		yValues[count] = a;
	}

	// we use formula c = ((Vt)V)^-1(Vt)(y)

	// create matrix objects for column matrices containing x and y values
	Matrix x_values = Matrix(numberPoints, 1, xValues);
	Matrix y_values = Matrix(numberPoints, 1, yValues);

	// create matrix object for vandermonde matrix and its transpose
	Matrix V = Matrix(numberPoints, order + 1, vandermonde(x_values.getArray(), numberPoints, order));
	Matrix Vt = Matrix(V.getCol(), V.getRow(), transpose(V.getArray(), V.getRow(), V.getCol()));
	
	// create matrix objects for (Vt)V and its inverse
	Matrix VtV = Matrix(Vt.getRow(), V.getCol(), multiply(Vt.getArray(), Vt.getRow(), Vt.getCol(), V.getArray(), V.getRow(), V.getCol()));
	Matrix VtVinverse = Matrix(VtV.getRow(), VtV.getCol(), inverse(VtV.getArray(), VtV.getRow()));

	// if (Vt)V is not invertible, calculation cannot proceed; display error messave
	if (VtVinverse.getArray() == nullptr)
	{
		std::cout << "The polynomial of best fit cannot be calculated. \n";
	}
	// if (Vt)V is invertible
	else
	{
		// create matrix object for (Vt)(y_values_
		Matrix Vty = Matrix(Vt.getRow(), y_values.getCol(), multiply(Vt.getArray(), Vt.getRow(), Vt.getCol(), y_values.getArray(), y_values.getRow(), y_values.getCol()));

		// create matrix object for final column vectors
		Matrix c = Matrix( VtVinverse.getRow(), Vty.getCol(), multiply(VtVinverse.getArray(), VtVinverse.getRow(), VtVinverse.getCol(), Vty.getArray(), Vty.getRow(), Vty.getCol()));
		
		// delete pointer
		delete[] Vty.getArray();

		// display coefficient column matrix
		print(c.getArray(), c.getRow(), c.getCol());
		// delete pointer
		delete[] c.getArray();
		
		int x = 0;
		std::cin >> x;
	}

	// delete pointers
	delete[] x_values.getArray();
	delete[] y_values.getArray();
	delete[] V.getArray();
	delete[] Vt.getArray();
	delete[] VtV.getArray();
	delete[] VtVinverse.getArray();

	int x = 0;
	std::cin >> x;
	return 0;
}