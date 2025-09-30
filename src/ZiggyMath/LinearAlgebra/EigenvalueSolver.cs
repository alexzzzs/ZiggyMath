using System;
using ZiggyMath.Core;

namespace ZiggyMath.LinearAlgebra
{
    /// <summary>
    /// Eigenvalue and eigenvector computation
    /// </summary>
    public static class EigenvalueSolver
    {
        /// <summary>
        /// Compute eigenvalues and eigenvectors using QR algorithm
        /// </summary>
        public static (Vector eigenvalues, Matrix eigenvectors) Compute(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for eigenvalue computation");

            int n = matrix.Rows;
            Matrix A = new Matrix(matrix.Rows, matrix.Columns);
            matrix.CopyTo(A);

            Matrix Q_total = Matrix.Identity(n);

            const int maxIterations = 1000;
            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // QR decomposition
                var (Q, R) = QRDecomposition.Decompose(A);

                // Update A = R * Q
                A = Matrix.Multiply(R, Q);

                // Accumulate Q
                Q_total = Matrix.Multiply(Q_total, Q);

                // Check convergence
                bool converged = true;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        if (Math.Abs(A[i, j]) > tolerance)
                        {
                            converged = false;
                            break;
                        }
                    }
                    if (!converged)
                        break;
                }

                if (converged)
                    break;
            }

            // Extract eigenvalues (diagonal elements)
            Vector eigenvalues = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i, i];
            }

            return (eigenvalues, Q_total);
        }

        /// <summary>
        /// Power method for finding dominant eigenvalue and eigenvector
        /// </summary>
        public static Vector PowerMethod(Matrix matrix, Vector initialGuess, int maxIterations = 1000)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for power method");
            if (matrix.Rows != initialGuess.Length)
                throw new ArgumentException("Matrix size must match vector length");

            Vector x = VectorMath.Normalize(initialGuess);
            Vector x_prev = new Vector(x.Length);

            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                x_prev.CopyTo(x);

                // Multiply matrix by vector
                x = Matrix.Multiply(matrix, x);

                // Normalize
                x = VectorMath.Normalize(x);

                // Check convergence
                double diff = VectorMath.Max(VectorMath.Abs(x - x_prev));
                if (diff < tolerance)
                    break;
            }

            return x;
        }

        /// <summary>
        /// Compute eigenvalues using simple iterative method
        /// </summary>
        public static Vector ComputeEigenvalues(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square for eigenvalue computation");

            int n = matrix.Rows;
            Matrix A = new Matrix(matrix.Rows, matrix.Columns);
            matrix.CopyTo(A);

            const int maxIterations = 100;
            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // QR decomposition
                var (Q, R) = QRDecomposition.Decompose(A);

                // Update A = R * Q
                A = Matrix.Multiply(R, Q);

                // Check convergence
                bool converged = true;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        if (Math.Abs(A[i, j]) > tolerance)
                        {
                            converged = false;
                            break;
                        }
                    }
                    if (!converged)
                        break;
                }

                if (converged)
                    break;
            }

            // Extract eigenvalues (diagonal elements)
            Vector eigenvalues = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i, i];
            }

            return eigenvalues;
        }

        /// <summary>
        /// Compute spectral radius (largest absolute eigenvalue)
        /// </summary>
        public static double SpectralRadius(Matrix matrix)
        {
            Vector eigenvalues = ComputeEigenvalues(matrix);
            double maxMagnitude = 0.0;

            for (int i = 0; i < eigenvalues.Length; i++)
            {
                double magnitude = Math.Abs(eigenvalues[i]);
                if (magnitude > maxMagnitude)
                    maxMagnitude = magnitude;
            }

            return maxMagnitude;
        }

        /// <summary>
        /// Check if matrix is positive definite
        /// </summary>
        public static bool IsPositiveDefinite(Matrix matrix)
        {
            if (!matrix.IsSquare)
                throw new ArgumentException("Matrix must be square");

            try
            {
                Vector eigenvalues = ComputeEigenvalues(matrix);
                for (int i = 0; i < eigenvalues.Length; i++)
                {
                    if (eigenvalues[i] <= 0)
                        return false;
                }
                return true;
            }
            catch
            {
                return false;
            }
        }
    }
}