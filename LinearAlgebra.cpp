////////////////////////////////////////////////////////////////////////////////
// File Name:      LinearAlgebra.cpp
//
// Author:         Ved Kale & Amairani Zepeda
// CS email:       vpkale@wisc.edu
// 				   zepeda2@wisc.edu
//
// Description:    A matrix calculator that computes math operations (addition,
//				   subtraction, and multiplication), matrix functions
//                 (determinant and transpose of matrices), and provides various
//                 display options.
//
////////////////////////////////////////////////////////////////////////////////

#include "LinearAlgebra.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

using namespace std;

const bool CONSOLE_OUTPUT = true;
vector<double**> matrices;
vector<pair<int, int>> sizes;

/**
 *  Checks if the given matrix exists. If it does not exist, it prompts the user
 *  for a valid matrix and displays all valid matrices.
 *
 *  @param matrix is the user given matrix
 *  @return a valid matrix
 */
int checkMatrix(int matrix) {
	bool valid = false;
	while (!valid) {
		int matSize = matrices.size();
		for (int mat = 0; mat < matSize; mat++) {
			if (mat == matrix) {
				valid = true;
			}
		}
		if (!valid) {
			cout << "\nMatrix " << matrix << " does not exist.\n"
					<< "Valid matrices are: ";
			for (unsigned i = 0; i < matrices.size(); i++) {
				cout << i << " ";
			}
			cout << "\nPlease enter a valid matrix: ";
			cin >> matrix;
		}
	}
	return matrix;
}

/**
 *  Displays the matrix that is chosen by the user.
 *
 *  @param mat is the chosen matrix
 *  @param rows is the number of rows of the chosen matrix
 *  @param columns is the number of columns of the chosen matrix
 */
void displayMatrix(double** mat, int rows, int columns) {
	double** matrix = mat;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

/**
 *  Defines a new matrix by prompting the user to enter the matrix dimensions.
 *  It prompts the user for individual matrix element input then displays the
 *  matrix.
 */
void defineMatrix() {
	int r, c;

	cout << "Enter number of rows: ";
	cin >> r;

	cout << "Enter number of columns: ";
	cin >> c;

	double** matrix;
	matrix = new double*[r];

	for (int i = 0; i < r; i++) {	// loop to create 2d matrix thats dynamic
		matrix[i] = new double[c];
	}

	cout << "Enter " << r * c << " numbers:\n";

	for (int i = 0; i < r; i++) {		// input stuff in
		for (int j = 0; j < c; j++) {
			cout << "  Position a" << i + 1 << j + 1 << ": ";
			cin >> matrix[i][j];
		}
	}

	matrices.push_back(matrix);
	sizes.push_back(make_pair(r, c));

	cout << "Matrix " << matrices.size() << " created.\n";
	displayMatrix(matrix, r, c);
}

/**
 *  Adds two given matrices. If the given matrices are of different dimensions,
 *  it displays an error message.
 *
 *  @param mat1 is the first matrix
 *  @param rows1 is the number of rows of the first matrix
 *  @param columns1 is the number of columns of the first matrix
 *  @param mat2 is the second matrix
 *  @param rows1 is the number of rows of the second matrix
 *  @param columns2 is the number of columns of the second matrix
 */
void addMatrix(double** mat1, int rows1, int columns1, double** mat2, int rows2,
		int columns2) {
	if (rows1 != rows2 || columns1 != columns2) {
		cout << "ERROR: Cannot add matrices. Dimensions mismatch." << endl;
	} else {
		double** matrix;
		matrix = new double*[rows1];

		for (int i = 0; i < rows1; i++) // loop to create 2d matrix thats dynamic
				{
			matrix[i] = new double[columns1];
		}

		for (int i = 0; i < rows1; i++) {
			for (int j = 0; j < columns1; j++) {
				matrix[i][j] = mat1[i][j] + mat2[i][j];
			}
			cout << endl;
		}
		displayMatrix(matrix, rows1, columns1);
	}
}

/**
 *  Subtracts two given matrices. If the given matrices are of different
 *  dimensions, it displays an error message.
 *
 *  @param mat1 is the first matrix
 *  @param rows1 is the number of rows of the first matrix
 *  @param columns1 is the number of columns of the first matrix
 *  @param mat2 is the second matrix
 *  @param rows1 is the number of rows of the second matrix
 *  @param columns2 is the number of columns of the second matrix
 */
void subtractMatrix(double** mat1, int rows1, int columns1, double** mat2,
		int rows2, int columns2) {
	if (rows1 != rows2 || columns1 != columns2) {
		cout << "ERROR: Cannot subtract matrices. Dimensions mismatch." << endl;
	} else {
		double** matrix;
		matrix = new double*[rows1];

		for (int i = 0; i < rows1; i++) { // loop to create 2d matrix thats dynamic
			matrix[i] = new double[columns1];
		}

		for (int i = 0; i < rows1; i++) {
			for (int j = 0; j < columns1; j++) {
				matrix[i][j] = mat1[i][j] - mat2[i][j];
			}
			cout << endl;
		}

		displayMatrix(matrix, rows1, columns1);
	}
}

/**
 *  Multiplies two given matrices. If the given matrices are of invalid
 *  dimensions, it displays an error message.
 *
 *  @param mat1 is the first matrix
 *  @param rows1 is the number of rows of the first matrix
 *  @param columns1 is the number of columns of the first matrix
 *  @param mat2 is the second matrix
 *  @param rows1 is the number of rows of the second matrix
 *  @param columns2 is the number of columns of the second matrix
 */
void multiplyMatrix(double** mat1, int rows1, int columns1, double** mat2,
		int rows2, int columns2) {
	if (columns1 != rows2) {
		cout
				<< "ERROR: Cannot multiply matrices. Dimensions must be (m x k) and (k x n)."
				<< endl;
	} else {
		double** matrix;
		matrix = new double*[rows1];

		for (int i = 0; i < rows1; i++) { // loop to create 2d matrix thats dynamic
			matrix[i] = new double[columns2];
		}

		for (int i = 0; i < rows1; i++) {
			for (int j = 0; j < columns2; j++) {
				for (int k = 0; k < rows2; k++) {
					matrix[i][j] = matrix[i][j] + mat1[i][k] * mat2[k][j];
				}
			}
			cout << endl;
		}
		displayMatrix(matrix, rows1, columns2);
	}
}

/**
 *  A helper method for determining the determinant of a square matrix. It gets
 *  the cofactor of a submatrix.
 *
 *  @param mat is the given matrix
 *  @param temp is a 2D matrix that is dynamic
 *  @param p is zero
 *  @param q is the column element of the first row
 *  @param n is the dimension of the square matrix
 */
void getCofactor(double** mat, double** temp, int p, int q, int n) {
	int i = 0, j = 0;

	// looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				temp[i][j++] = mat[row][col];

				// row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

/**
 *  A recursive function for finding the determinant of a square matrix.
 *
 *  @param mat is the given square matrix
 *  @param n is the dimension of the square matrix.
 */
int determinantOfMatrix(double** mat, int n) {
	int D = 0; // Initialize result

	// base case : if matrix contains single element
	if (n == 1)
		return mat[0][0];

	double** temp;
	temp = new double*[n];
	for (int i = 0; i < n; i++) { // loop to create 2d matrix thats dynamic
		temp[i] = new double[n];
	}
	// to store cofactors

	int sign = 1;  // to store sign multiplier

	// iteration for each element of first row
	for (int i = 0; i < n; i++) {
		// getts cofactor of mat[0][f]
		getCofactor(mat, temp, 0, i, n);
		D += sign * mat[0][i] * determinantOfMatrix(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

/**
 *  Finds the transpose of a given matrix.
 *
 *  @param mat is the given matrix
 *  @param rows is the number of rows of the matrix
 *  @param columns is the number of columns of the matrix
 */
void transposeMatrix(double** mat, int rows, int columns) {
	double** matrix = mat;

	cout << "\nOrignial Matrix: " << endl;
	displayMatrix(mat, rows, columns);

	cout << "\nTransposed Matrix: " << endl;
	for (int i = 0; i < columns; i++) {
		for (int j = 0; j < rows; j++) {
			cout << matrix[j][i] << "\t";
		}
		cout << endl;
	}
}
/**
 *  Determines the singularity of a matrix.
 *
 *  @param mat is the given matrix
 *  @param rows is the number of rows of the given matrix
 *  @param columns is the number of columns of the given matrix
 *  @param matrix is the name of the matrix
 */
void matrixSingularity(double** mat, int rows, int columns, int matrix) {
	int determinant = determinantOfMatrix(mat, rows);

	if (rows != columns) {
		cout << "Matrix " << matrix
				<< " is not singular because dimensions are not equal.\n";
	} else if (determinant != 0) {
		cout << "Matrix " << matrix << " is singular.\n";
	} else {
		cout << "Matrix " << matrix << " is not singular.\n";
	}
}

/**
 *  A helper method for the sub menu Display Matrix. There are 3 display matrix
 *  options: number of matrices defined, a single matrix, and multiple matrices.
 *
 *  @param choice is the users choice
 */
void displayMatrixOptions(int choice) {
	int ans;

	if (choice == 0) {

		cout << "Number of Matrices: " << matrices.size() << "\n" << endl;

	} else if (choice == 1) {

		cout << "Matrix to display: ";
		cin >> ans;
		ans = checkMatrix(ans);
		displayMatrix(matrices.at(ans), sizes.at(ans).first,
				sizes.at(ans).second);

	} else if (choice == 2) {

		cout << "\nMatrices to display (enter -1 to stop)\n"
				<< "Matrix to display: ";
		cin >> ans;

		while (ans != -1) {
			if (ans != -1) {
				ans = checkMatrix(ans);
				displayMatrix(matrices.at(ans), sizes.at(ans).first,
						sizes.at(ans).second);
				cout << "\nMatrix to display: ";
				cin >> ans;
			}
		}
	}
}

/**
 *  A helper method for the sub menu Math Operations. There are 3 math operation
 *  options: addition, subtraction, and multiplication. It checks that the given
 *  matrices are valid and executes the matrix operation if input is valid.
 *
 *  @param choice is the users choice
 */
void mathOperations(int choice) {
	int a1, a2;
	if (choice == 0) {			// addition

		cout << "Input two matrix numbers to add: " << endl;
		cout << "Matrix A: ";
		cin >> a1;
		a1 = checkMatrix(a1);

		cout << "Matrix B: ";
		cin >> a2;
		a2 = checkMatrix(a2);

		cout << "Addition of Matrix A and Matrix B is:\n";
		addMatrix(matrices.at(a1), sizes.at(a1).first, sizes.at(a1).second,
				matrices.at(a2), sizes.at(a2).first, sizes.at(a2).second);

	} else if (choice == 1) {	// subtraction

		cout << "Input two matrix numbers to subtract (A - B) : " << endl;
		cout << "Matrix A: ";
		cin >> a1;
		a1 = checkMatrix(a1);

		cout << "Matrix B: ";
		cin >> a2;
		a2 = checkMatrix(a2);

		cout << "Subtraction of Matrix A and Matrix B is:\n";
		subtractMatrix(matrices.at(a1), sizes.at(a1).first, sizes.at(a1).second,
				matrices.at(a2), sizes.at(a2).first, sizes.at(a2).second);

	} else if (choice == 2) {	// multiplication

		cout << "Input two matrix numbers to multiply (A * B) : " << endl;
		cout << "Matrix A: ";
		cin >> a1;
		a1 = checkMatrix(a1);

		cout << "Matrix B: ";
		cin >> a2;
		a2 = checkMatrix(a2);

		cout << "Multiplication of Matrix A and Matrix B is:\n";
		multiplyMatrix(matrices.at(a1), sizes.at(a1).first, sizes.at(a1).second,
				matrices.at(a2), sizes.at(a2).first, sizes.at(a2).second);

	} else
		cout << "That's not a choice.\n";
}

/**
 *  A helper method for the sub menu Matrix Functions. There are 2 matrix
 *  function options: determinant and transpose. It checks that the given
 *  matrix is valid and executes the matrix operation if input is valid.
 *
 *  @param choice is the users choice
 */
void matrixFunctions(int choice) {
	int ans = -1;

	if (choice == 0) {				// matrix determinant

		cout << "Input matrix number to find determinant of: " << endl;
		cout << "Matrix: ";
		cin >> ans;
		ans = checkMatrix(ans);

		if (sizes.at(ans).first != sizes.at(ans).second)
			cout << "ERROR! Dimensions are not equal" << endl;
		else {
			cout << "Determinant of Matrix is:\n";
			cout << determinantOfMatrix(matrices.at(ans), sizes.at(ans).first)
					<< endl;
		}

	} else if (choice == 1) {		// transpose matrix

		cout << "Matrix to transpose: ";
		cin >> ans;
		ans = checkMatrix(ans);

		transposeMatrix(matrices.at(ans), sizes.at(ans).first,
				sizes.at(ans).second);

	} else if (choice == 2) {

		cout << "Matrix to find singularity: ";
		cin >> ans;
		ans = checkMatrix(ans);

		matrixSingularity(matrices.at(ans), sizes.at(ans).first,
				sizes.at(ans).second, ans);

	} else
		cout << "That's not a choice.\n";
}

/**
 *  Executes the main program. Displays a menu with various options to choose
 *  from.
 */
int main() {

	if (CONSOLE_OUTPUT) {
		bool exit = false;

		// displays title
		cout << "*************************\n" << "**  MATRIX CALCULATOR  **\n"
				<< "*************************" << endl;

		while (!exit) {
			int choice = -1;

			// menu options
			cout << "\nMain Menu\n" << "0 - Quit\n" << "1 - Define Matrix\n"
					<< "2 - Display Matrix\n" << "3 - Math Operations\n"
					<< "4 - Matrix Functions\n" << "Enter your choice: ";

			while (!(cin >> choice)) {

				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "\nERROR: Please enter valid input!\n"
						<< "Enter your choice again: ";
			}

			if (choice == 0) {		// exits program

				cout << "Exiting Matrix Calculator...\n" << "Goodbye!";
				exit = true;

			} else if (choice == 1) {				// defines a matrix

				cout << "Define Matrix " << matrices.size() << " " << endl;
				defineMatrix();

			} else if (matrices.empty()) {	// if no matrix has been defined

				cout << "Please define a matrix first." << endl;

			} else if (choice == 2) {		// displays a matrix

				cout << "\nDisplay Options\n"
						<< "0 - Total number of matrices\n"
						<< "1 - Single matrix\n" << "2 - Multiple matrices\n"
						<< "Enter your choice: ";

				while (!(cin >> choice)) {

					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "\nERROR: Please enter valid input!\n"
							<< "Enter your choice again: ";
				}

				displayMatrixOptions(choice);

			} else if (choice == 3) {		// basic operations

				cout << "\nMath Operations\n" << "0 - Addition\n"
						<< "1 - Subtraction\n" << "2 - Multiplication\n"
						<< "Enter your choice: ";

				while (!(cin >> choice)) {

					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "\nERROR: Please enter valid input!\n"
							<< "Enter your choice again: ";
				}
				mathOperations(choice);

			} else if (choice == 4) {		// other matrix options

				cout << "\nMatrix Functions\n" << "0 - Determinant\n"
						<< "1 - Transpose\n" << "2 - Singularity\n"
						<< "Enter your choice: ";

				while (!(cin >> choice)) {

					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "\nERROR: Please enter valid input!\n"
							<< "Enter your choice again: ";
				}

				matrixFunctions(choice);

			} else
				cout << "That's not a choice.\n";

		}
	}
	return 0;
}
