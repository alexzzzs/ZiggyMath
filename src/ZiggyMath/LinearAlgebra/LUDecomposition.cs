using System;
using ZiggyMath.Core;

namespace ZiggyMath.LinearAlgebra
{
    /// <summary>
    /// LU Decomposition with partial pivoting
    /// </summary>
    public static class LUDecomposition
    {
        /// <summary>
        /// Decompose matrix into L, U, and permutation matrices
        /// </summary>
        public static (Matrix L, Matrix U, int[] permutations) Decompose(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for LU decomposition");

            int n = matrix.Rows;
            Matrix L = new Matrix(n);
            Matrix U = new Matrix(matrix.Rows, matrix.Columns);
            int[] permutations = new int[n];

            // Initialize permutations
            for (int i = 0; i < n; i++)
            {
                permutations[i] = i;
            }

            // Copy matrix to U
            matrix.CopyTo(U);

            // LU decomposition with partial pivoting
            for (int i = 0; i < n; i++)
            {
                // Partial pivoting
                int maxRow = i;
                for (int k = i; k < n; k++)
                {
                    if (Math.Abs(U[k, i]) > Math.Abs(U[maxRow, i]))
                    {
                        maxRow = k;
                    }
                }

                // Swap rows in U and permutations
                if (maxRow != i)
                {
                    SwapRows(U, i, maxRow);
                    (permutations[i], permutations[maxRow]) = (permutations[maxRow], permutations[i]);
                }

                // Check for singular matrix
                if (Math.Abs(U[i, i]) < 1e-10)
                {
                    throw new InvalidOperationException("Matrix is singular");
                }

                // Gaussian elimination
                for (int k = i + 1; k < n; k++)
                {
                    double factor = U[k, i] / U[i, i];
                    L[k, i] = factor;

                    for (int j = i; j < n; j++)
                    {
                        U[k, j] -= factor * U[i, j];
                    }
                }
            }

            // Extract L (lower triangular with unit diagonal)
            for (int i = 0; i < n; i++)
            {
                L[i, i] = 1.0;
                for (int j = 0; j < i; j++)
                {
                    L[i, j] = U[i, j];
                    U[i, j] = 0.0;
                }
            }

            return (L, U, permutations);
        }

        /// <summary>
        /// Solve system using LU decomposition
        /// </summary>
        public static void Solve(Matrix L, Matrix U, Vector b, int[] permutations, out Vector x)
        {
            int n = L.Rows;
            x = new Vector(n);

            // Apply permutations to b
            Vector pb = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                pb[i] = b[permutations[i]];
            }

            // Forward substitution: L * y = P * b
            Vector y = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                {
                    sum += L[i, j] * y[j];
                }
                y[i] = (pb[i] - sum) / L[i, i];
            }

            // Back substitution: U * x = y
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                {
                    sum += U[i, j] * x[j];
                }
                x[i] = (y[i] - sum) / U[i, i];
            }
        }

        /// <summary>
        /// Solve system using pre-computed LU decomposition
        /// </summary>
        public static Vector Solve(Matrix L, Matrix U, int[] permutations, Vector b)
        {
            Vector x;
            Solve(L, U, b, permutations, out x);
            return x;
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