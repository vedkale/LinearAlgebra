// LinearAlgebra.cpp : Defines the entry point for the application.
//

#include "LinearAlgebra.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

const bool CONSOLE_OUTPUT = true;
vector<double**> matrices;
vector<pair<int, int>> sizes;

void defineMatrix();
void displayMatrix(double** mat, int rows, int columns);
void addMatrix(double** mat1, int rows1, int columns1, double** mat2, int rows2, int columns2);
void subtractMatrix(double** mat1, int rows1, int columns1, double** mat2, int rows2, int columns2);
int determinantOfMatrix(double** mat, int n);


int main()
{
	//matrices = new vector<double**>();

	if (CONSOLE_OUTPUT)
	{
		bool exit = false;
		while (!exit)
		{
			int choice = -1;
			cout << "" //"\n** You choose from the following **\n\n"
				<< "0 - Define Matrix\n"
				<< "1 - Display Matrix\n"
				<< "2 - Addition\n"
				<< "3 - Substraction\n"
				<< "4 - Determinant\n"
				<< "5 - Quit\n"
				<< "Enter your choice: ";

			while (!(cin >> choice)) {
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\nERROR: Please enter valid input!\n"
					<< "Enter your choice again: ";
			}

			if (choice == 0)
			{
				cout << "Define Matrix " << matrices.size() << " " << endl;
				defineMatrix();
			}
			else if (choice == 1)
			{
				int ans;
				cout << "Matrix to display: ";
				cin >> ans;
				displayMatrix(matrices.at(ans), sizes.at(ans).first, sizes.at(ans).second);
			}
			else if (choice == 2)
			{
				int a1, a2;
				cout << "Input two matrix numbers to add: " << endl;
				cout << "Matrix A: ";
				cin >> a1;
				cout << "Matrix B: ";
				cin >> a2;
				cout << "Addition of Matrix A and Matrix B is:\n";
				addMatrix(matrices.at(a1), sizes.at(a1).first, sizes.at(a1).second, matrices.at(a2), sizes.at(a2).first, sizes.at(a2).second);
			}
			else if (choice == 3)
			{
				int a1, a2;
				cout << "Input two matrix numbers to subtract (A - B) : " << endl;
				cout << "Matrix A: ";
				cin >> a1;
				cout << "Matrix B: ";
				cin >> a2;
				cout << "Subtraction of Matrix A and Matrix B is:\n";
				subtractMatrix(matrices.at(a1), sizes.at(a1).first, sizes.at(a1).second, matrices.at(a2), sizes.at(a2).first, sizes.at(a2).second);
			}
			else if (choice == 4)
			{
				int a1, a2;
				cout << "Input matrix number to find determinant of: " << endl;
				cout << "Matrix: ";
				cin >> a1;
				if (sizes.at(a1).first != sizes.at(a1).second)
					cout << "ERROR! Dimensions are not equal" << endl;
				else
				{
					cout << "Determinant of Matrix is:\n";
					cout << determinantOfMatrix(matrices.at(a1), sizes.at(a1).first) << endl;
				}
			}
			else if (choice == 5)
			{
				cout << "Exiting...\n";
				exit = true;
			}
			else
				cout << "That's not a choice.\n";

			cout << "Number of Matrices: " << matrices.size() << endl;
		}
	}
	return 0;
}

void defineMatrix()
{
	int r, c;

	cout << "Enter number of rows: ";
	cin >> r;

	cout << "Enter number of columns: ";
	cin >> c;

	double** matrix;
	matrix = new double* [r];

	for (int i = 0; i < r; i++) // loop to create 2d matrix thats dynamic
	{
		matrix[i] = new double[c];
	}

	for (int i = 0; i < r; i++) // input stuff in
	{
		for (int j = 0; j < c; j++)
		{
			cin >> matrix[i][j];
		}
	}

	matrices.push_back(matrix);
	sizes.push_back(make_pair(r, c));

}

void displayMatrix(double** mat, int rows, int columns)
{
	double** matrix = mat;
	for (int i = 0; i < rows; i++)
	{
		//cout << "\t";
		for (int j = 0; j < columns; j++)
		{
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void addMatrix(double** mat1, int rows1, int columns1, double** mat2, int rows2, int columns2)
{
	if (rows1 != rows2 || columns1 != columns2)
	{
		cout << "ERROR: Cannot add matrices. Dimensions mismatch." << endl;
	}

	double** matrix;
	matrix = new double* [rows1];

	for (int i = 0; i < rows1; i++) // loop to create 2d matrix thats dynamic
	{
		matrix[i] = new double[columns1];
	}

	for (int i = 0; i < rows1; i++)
	{
		for (int j = 0; j < columns1; j++)
		{
			matrix[i][j] = mat1[i][j] + mat2[i][j];
		}
		cout << endl;
	}

	displayMatrix(matrix, rows1, columns1);

}

void subtractMatrix(double** mat1, int rows1, int columns1, double** mat2, int rows2, int columns2)
{
	if (rows1 != rows2 || columns1 != columns2)
	{
		cout << "ERROR: Cannot add matrices. Dimensions mismatch." << endl;
	}

	double** matrix;
	matrix = new double* [rows1];

	for (int i = 0; i < rows1; i++) // loop to create 2d matrix thats dynamic
	{
		matrix[i] = new double[columns1];
	}

	for (int i = 0; i < rows1; i++)
	{
		for (int j = 0; j < columns1; j++)
		{
			matrix[i][j] = mat1[i][j] - mat2[i][j];
		}
		cout << endl;
	}

	displayMatrix(matrix, rows1, columns1);

}

void getCofactor(double** mat, double** temp, int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp[i][j++] = mat[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of mat[][]. */
int determinantOfMatrix(double** mat, int n)
{
	int D = 0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1)
		return mat[0][0];

	double** temp;
	temp = new double* [n];
	for (int i = 0; i < n; i++) // loop to create 2d matrix thats dynamic
	{
		temp[i] = new double[n];
	}
	// To store cofactors 

	int sign = 1;  // To store sign multiplier 

	 // Iterate for each element of first row 
	for (int i = 0; i < n; i++)
	{
		// Getting Cofactor of mat[0][f] 
		getCofactor(mat, temp, 0, i, n);
		D += sign * mat[0][i] * determinantOfMatrix(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}