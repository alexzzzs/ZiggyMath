using System;
using ZiggyMath.Core;
using ZiggyAlloc;

namespace ZiggyMath.LinearAlgebra
{
    /// <summary>
    /// High-performance linear algebra operations using ZiggyAlloc optimization
    /// </summary>
    public static class OptimizedLinearAlgebra
    {
        /// <summary>
        /// Gaussian elimination with partial pivoting using ZiggyAlloc optimization
        /// </summary>
        public static Vector SolveLinearSystem(Matrix A, Vector b, MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.LinearAlgebra);

            // Use optimized matrix for row operations
            using var optimizedA = MatrixLayoutOptimizer.CreateOptimizedMatrix(
                A.Rows, A.Columns, MatrixOperation.RowWise, ctx);
            A.CopyTo(optimizedA);

            using var x = ctx.MemoryPool.Rent<double>(A.Rows);
            using var bCopy = ctx.MemoryPool.Rent<double>(b.Length);

            // Copy b to modifiable buffer
            b.AsSpan().CopyTo(bCopy.AsSpan());

            // Forward elimination with partial pivoting
            for (int i = 0; i < A.Rows - 1; i++)
            {
                // Find pivot
                int pivot = i;
                for (int j = i + 1; j < A.Rows; j++)
                {
                    if (Math.Abs(optimizedA[j, i]) > Math.Abs(optimizedA[pivot, i]))
                        pivot = j;
                }

                // Swap rows if needed
                if (pivot != i)
                    SwapRows(optimizedA, bCopy, i, pivot, ctx.Allocator);

                // Eliminate
                for (int j = i + 1; j < A.Rows; j++)
                {
                    double factor = optimizedA[j, i] / optimizedA[i, i];
                    for (int k = i; k < A.Columns; k++)
                    {
                        optimizedA[j, k] -= factor * optimizedA[i, k];
                    }
                    bCopy[j] -= factor * bCopy[i];
                }
            }

            // Back substitution
            for (int i = A.Rows - 1; i >= 0; i--)
            {
                x[i] = bCopy[i];
                for (int j = i + 1; j < A.Rows; j++)
                {
                    x[i] -= optimizedA[i, j] * x[j];
                }
                x[i] /= optimizedA[i, i];
            }

            return new Vector(x.AsSpan(), ctx.Allocator);
        }

        private static void SwapRows(Matrix matrix, UnmanagedBuffer<double> vector, int row1, int row2, IUnmanagedMemoryAllocator allocator)
        {
            // Efficient row swapping with ZiggyAlloc optimization
            using var tempRow = allocator.Allocate<double>(matrix.Columns);
            using var tempVec = allocator.Allocate<double>(1);

            // Swap matrix rows using GetRow/SetRow pattern
            var row1Data = new double[matrix.Columns];
            var row2Data = new double[matrix.Columns];
            for (int c = 0; c < matrix.Columns; c++)
            {
                row1Data[c] = matrix[row1, c];
                row2Data[c] = matrix[row2, c];
            }

            for (int c = 0; c < matrix.Columns; c++)
            {
                matrix[row1, c] = row2Data[c];
                matrix[row2, c] = row1Data[c];
            }

            // Swap vector elements
            tempVec[0] = vector[row1];
            vector[row1] = vector[row2];
            vector[row2] = tempVec[0];
        }

        /// <summary>
        /// LU Decomposition with partial pivoting using ZiggyAlloc optimization
        /// </summary>
        public static (Matrix L, Matrix U, int[] permutations) DecomposeLU(Matrix A, MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.LinearAlgebra);

            using var optimizedA = MatrixLayoutOptimizer.CreateOptimizedMatrix(
                A.Rows, A.Columns, MatrixOperation.RowWise, ctx);
            A.CopyTo(optimizedA);

            int n = A.Rows;
            var permutations = new int[n];
            for (int i = 0; i < n; i++) permutations[i] = i;

            // LU decomposition with partial pivoting
            for (int i = 0; i < n - 1; i++)
            {
                // Find pivot
                int pivot = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(optimizedA[permutations[j], i]) > Math.Abs(optimizedA[permutations[pivot], i]))
                        pivot = j;
                }

                // Swap permutations
                (permutations[i], permutations[pivot]) = (permutations[pivot], permutations[i]);

                // Elimination
                for (int j = i + 1; j < n; j++)
                {
                    optimizedA[permutations[j], i] /= optimizedA[permutations[i], i];
                    for (int k = i + 1; k < n; k++)
                    {
                        optimizedA[permutations[j], k] -=
                            optimizedA[permutations[j], i] * optimizedA[permutations[i], k];
                    }
                }
            }

            // Extract L and U matrices
            using var L = MatrixLayoutOptimizer.CreateOptimizedMatrix(n, n, MatrixOperation.ElementWise, ctx);
            using var U = MatrixLayoutOptimizer.CreateOptimizedMatrix(n, n, MatrixOperation.ElementWise, ctx);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i > j)
                    {
                        L[i, j] = optimizedA[permutations[i], j];
                        U[i, j] = 0;
                    }
                    else if (i == j)
                    {
                        L[i, j] = 1;
                        U[i, j] = optimizedA[permutations[i], j];
                    }
                    else
                    {
                        L[i, j] = 0;
                        U[i, j] = optimizedA[permutations[i], j];
                    }
                }
            }

            return (new Matrix(L.AsSpan(), n, n, ctx.Allocator),
                    new Matrix(U.AsSpan(), n, n, ctx.Allocator), permutations);
        }
    }
}