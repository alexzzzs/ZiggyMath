using System;
using ZiggyMath.Core;

namespace ZiggyMath.LinearAlgebra
{
    /// <summary>
    /// Gaussian elimination with partial pivoting
    /// </summary>
    public static class GaussianElimination
    {
        /// <summary>
        /// Solve linear system Ax = b
        /// </summary>
        public static void Solve(Matrix A, Vector b, out Vector x)
        {
            if (!A.IsSquare)
                throw new ArgumentException("Matrix A must be square");
            if (A.Rows != b.Length)
                throw new ArgumentException("Matrix rows must equal vector length");

            int n = A.Rows;
            x = new Vector(n);

            // Create augmented matrix [A|b]
            Matrix augmented = CreateAugmentedMatrix(A, b);

            // Forward elimination with partial pivoting
            ForwardElimination(augmented);

            // Back substitution
            BackSubstitution(augmented, x);
        }

        /// <summary>
        /// Solve linear system AX = B
        /// </summary>
        public static void Solve(Matrix A, Matrix B, out Matrix X)
        {
            if (!A.IsSquare)
                throw new ArgumentException("Matrix A must be square");
            if (A.Rows != B.Rows)
                throw new ArgumentException("Matrix A rows must equal matrix B rows");

            int n = A.Rows;
            int m = B.Columns;
            X = new Matrix(n, m);

            // Solve each column of B
            for (int j = 0; j < m; j++)
            {
                Vector b = new Vector(n);
                for (int i = 0; i < n; i++)
                {
                    b[i] = B[i, j];
                }

                Vector x_col = new Vector(n);
                Solve(A, b, out x_col);

                for (int i = 0; i < n; i++)
                {
                    X[i, j] = x_col[i];
                }
            }
        }

        /// <summary>
        /// Invert matrix using Gaussian elimination
        /// </summary>
        public static Matrix Invert(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for inversion");

            int n = matrix.Rows;
            Matrix identity = Matrix.Identity(n);
            Matrix inverse = new Matrix(n);

            Solve(matrix, identity, out inverse);
            return inverse;
        }

        /// <summary>
        /// Calculate matrix determinant
        /// </summary>
        public static double Determinant(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for determinant");

            int n = matrix.Rows;
            if (n == 1)
                return matrix[0, 0];
            if (n == 2)
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];

            // Create copy for elimination
            Matrix work = new Matrix(matrix.Rows, matrix.Columns);
            matrix.CopyTo(work);

            // Forward elimination
            ForwardElimination(work);

            // Calculate determinant from diagonal
            double det = 1.0;
            for (int i = 0; i < n; i++)
            {
                det *= work[i, i];
            }

            return det;
        }

        /// <summary>
        /// Create augmented matrix [A|b]
        /// </summary>
        private static Matrix CreateAugmentedMatrix(Matrix A, Vector b)
        {
            int n = A.Rows;
            int m = A.Columns + 1;
            Matrix augmented = new Matrix(n, m);

            // Copy A
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < A.Columns; j++)
                {
                    augmented[i, j] = A[i, j];
                }
            }

            // Copy b
            for (int i = 0; i < n; i++)
            {
                augmented[i, A.Columns] = b[i];
            }

            return augmented;
        }

        /// <summary>
        /// Forward elimination with partial pivoting - the heart of Gaussian elimination.
        /// This is where we transform the matrix into upper triangular form.
        /// Think of it as solving a system of equations by eliminating variables one by one.
        /// </summary>
        /// <param name="augmented">The augmented matrix [A|b] to be transformed</param>
        private static void ForwardElimination(Matrix augmented)
        {
            int n = augmented.Rows;

            // Process each column from left to right
            for (int i = 0; i < n; i++)
            {
                // PARTIAL PIVOTING: Find the row with the largest absolute value in column i
                // This prevents division by small numbers and improves numerical stability
                // Because floating-point arithmetic can be surprisingly fragile!
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmented[k, i]) > Math.Abs(augmented[maxRow, i]))
                    {
                        maxRow = k;  // Found a better pivot element
                    }
                }

                // Swap rows if we found a better pivot
                // This is like choosing the best player for your basketball team!
                if (maxRow != i)
                {
                    SwapRows(augmented, i, maxRow);
                }

                // Check for singular matrix (no unique solution)
                // If the pivot is too small, the matrix might be singular
                if (Math.Abs(augmented[i, i]) < 1e-10)
                {
                    throw new InvalidOperationException("Matrix is singular - no unique solution exists. " +
                        "This is like trying to divide by zero, but for linear systems!");
                }

                // ELIMINATION STEP: Make all elements below the pivot zero
                // For each row below the current row...
                for (int k = i + 1; k < n; k++)
                {
                    // Calculate the elimination factor
                    double factor = augmented[k, i] / augmented[i, i];

                    // Subtract factor * pivot_row from current row
                    // This eliminates the element at position (k, i)
                    for (int j = i; j < augmented.Columns; j++)
                    {
                        augmented[k, j] -= factor * augmented[i, j];
                    }
                }
            }

            // After this process, the matrix is in upper triangular form
            // Ready for back substitution to find the solution!
        }

        /// <summary>
        /// Back substitution
        /// </summary>
        private static void BackSubstitution(Matrix augmented, Vector x)
        {
            int n = augmented.Rows;

            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = augmented[i, n];

                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= augmented[i, j] * x[j];
                }

                x[i] /= augmented[i, i];
            }
        }

        /// <summary>
        /// Swap two rows in matrix
        /// </summary>
        private static void SwapRows(Matrix matrix, int row1, int row2)
        {
            int columns = matrix.Columns;
            for (int j = 0; j < columns; j++)
            {
                double temp = matrix[row1, j];
                matrix[row1, j] = matrix[row2, j];
                matrix[row2, j] = temp;
            }
        }
    }
}