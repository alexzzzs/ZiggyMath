using System;
using ZiggyMath.Core;

namespace ZiggyMath.LinearAlgebra
{
    /// <summary>
    /// QR Decomposition using Householder reflections
    /// </summary>
    public static class QRDecomposition
    {
        /// <summary>
        /// Decompose matrix into Q and R matrices
        /// </summary>
        public static (Matrix Q, Matrix R) Decompose(Matrix matrix)
        {
            int m = matrix.Rows;
            int n = matrix.Columns;
            int k = Math.Min(m, n);

            // Create working matrices
            Matrix R = new Matrix(matrix.Rows, matrix.Columns);
            matrix.CopyTo(R);

            Matrix Q = Matrix.Identity(m);

            // Householder reflections
            for (int j = 0; j < k; j++)
            {
                // Compute Householder vector
                Vector x = new Vector(m - j);
                for (int i = 0; i < m - j; i++)
                {
                    x[i] = R[i + j, j];
                }

                double norm = x.Magnitude();
                if (norm == 0)
                    continue;

                // Create Householder vector
                Vector v = new Vector(m - j);
                v[0] = x[0] >= 0 ? x[0] + norm : x[0] - norm;

                for (int i = 1; i < m - j; i++)
                {
                    v[i] = x[i];
                }

                double vNorm = v.Magnitude();
                if (vNorm == 0)
                    continue;

                // Normalize v
                v = Vector.Multiply(v, 1.0 / vNorm);

                // Apply reflection to R
                for (int p = j; p < n; p++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < m - j; i++)
                    {
                        dot += v[i] * R[i + j, p];
                    }
                    dot *= 2.0;

                    for (int i = 0; i < m - j; i++)
                    {
                        R[i + j, p] -= dot * v[i];
                    }
                }

                // Apply reflection to Q
                Matrix QTemp = new Matrix(m);
                Q.CopyTo(QTemp);

                for (int p = 0; p < m; p++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < m - j; i++)
                    {
                        dot += v[i] * QTemp[i + j, p];
                    }
                    dot *= 2.0;

                    for (int i = 0; i < m - j; i++)
                    {
                        Q[i + j, p] = QTemp[i + j, p] - dot * v[i];
                    }
                }
            }

            // Extract upper triangular R
            Matrix R_upper = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    R_upper[i, j] = R[i, j];
                }
            }

            return (Q, R_upper);
        }

        /// <summary>
        /// Solve system using QR decomposition
        /// </summary>
        public static void Solve(Matrix Q, Matrix R, Vector b, out Vector x)
        {
            int n = R.Rows;

            // Solve R * x = Q^T * b
            Vector Qt_b = new Vector(n);

            // Compute Q^T * b
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < Q.Rows; j++)
                {
                    sum += Q[j, i] * b[j];
                }
                Qt_b[i] = sum;
            }

            // Solve R * x = Q^T * b using back substitution
            x = new Vector(n);
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = Qt_b[i];

                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= R[i, j] * x[j];
                }

                if (Math.Abs(R[i, i]) < 1e-10)
                    throw new InvalidOperationException("Matrix is singular");

                x[i] /= R[i, i];
            }
        }

        /// <summary>
        /// Solve system using pre-computed QR decomposition
        /// </summary>
        public static Vector Solve(Matrix Q, Matrix R, Vector b)
        {
            Vector x;
            Solve(Q, R, b, out x);
            return x;
        }
    }
}