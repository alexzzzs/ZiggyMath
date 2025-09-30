using System;
using ZiggyAlloc;

namespace ZiggyMath.Core
{
    /// <summary>
    /// Optimizes matrix memory layout for specific operations
    /// </summary>
    public static class MatrixLayoutOptimizer
    {
        /// <summary>
        /// Creates a matrix optimized for the specified operation type
        /// </summary>
        public static Matrix CreateOptimizedMatrix(
            int rows,
            int columns,
            MatrixOperation operation,
            MathComputationContext context = null)
        {
            var allocator = operation switch
            {
                MatrixOperation.RowWise =>
                    new SlabAllocator(context?.Allocator ?? Z.DefaultAllocator, slabSize: 512 * 1024),
                MatrixOperation.ColumnWise =>
                    new SystemMemoryAllocator(),
                MatrixOperation.BlockWise =>
                    new LargeBlockAllocator(context?.Allocator ?? Z.DefaultAllocator),
                MatrixOperation.ElementWise =>
                    new UnmanagedMemoryPool(context?.Allocator ?? Z.DefaultAllocator),
                _ => context?.Allocator ?? Z.DefaultAllocator
            };

            return new Matrix(rows, columns, allocator);
        }

        /// <summary>
        /// Optimizes matrix multiplication based on dimensions
        /// </summary>
        public static Matrix MultiplyOptimized(Matrix a, Matrix b, MathComputationContext context = null)
        {
            if (a.Columns != b.Rows)
                throw new ArgumentException("Matrix dimensions don't match for multiplication");

            // Choose optimization strategy based on matrix characteristics
            var strategy = ChooseMultiplicationStrategy(a, b);

            return strategy switch
            {
                MultiplicationStrategy.Classic => MultiplyClassic(a, b, context),
                MultiplicationStrategy.Blocked => MultiplyBlocked(a, b, context),
                MultiplicationStrategy.Strassen => MultiplyStrassen(a, b, context),
                _ => MultiplyClassic(a, b, context)
            };
        }

        private static Matrix MultiplyClassic(Matrix a, Matrix b, MathComputationContext context)
        {
            using var result = CreateOptimizedMatrix(a.Rows, b.Columns, MatrixOperation.RowWise, context);

            // Classic O(n^3) multiplication with optimized memory access
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < b.Columns; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < a.Columns; k++)
                    {
                        sum += a[i, k] * b[k, j];
                    }
                    result[i, j] = sum;
                }
            }

            return result;
        }

        private static Matrix MultiplyBlocked(Matrix a, Matrix b, MathComputationContext context)
        {
            using var result = CreateOptimizedMatrix(a.Rows, b.Columns, MatrixOperation.BlockWise, context);

            const int blockSize = 64; // Tune based on cache size

            for (int i = 0; i < a.Rows; i += blockSize)
            {
                for (int j = 0; j < b.Columns; j += blockSize)
                {
                    for (int k = 0; k < a.Columns; k += blockSize)
                    {
                        MultiplyBlock(a, b, result, i, j, k, blockSize);
                    }
                }
            }

            return result;
        }

        private static void MultiplyBlock(Matrix a, Matrix b, Matrix result,
                                        int i, int j, int k, int blockSize)
        {
            int iEnd = Math.Min(i + blockSize, a.Rows);
            int jEnd = Math.Min(j + blockSize, b.Columns);
            int kEnd = Math.Min(k + blockSize, a.Columns);

            for (int ii = i; ii < iEnd; ii++)
            {
                for (int jj = j; jj < jEnd; jj++)
                {
                    double sum = result[ii, jj];
                    for (int kk = k; kk < kEnd; kk++)
                    {
                        sum += a[ii, kk] * b[kk, jj];
                    }
                    result[ii, jj] = sum;
                }
            }
        }

        private static Matrix MultiplyStrassen(Matrix a, Matrix b, MathComputationContext context)
        {
            // For now, fall back to classic multiplication
            // Strassen algorithm would be implemented here for very large matrices
            return MultiplyClassic(a, b, context);
        }

        private enum MultiplicationStrategy { Classic, Blocked, Strassen }

        private static MultiplicationStrategy ChooseMultiplicationStrategy(Matrix a, Matrix b)
        {
            int n = Math.Max(Math.Max(a.Rows, a.Columns), Math.Max(b.Rows, b.Columns));

            return n switch
            {
                < 64 => MultiplicationStrategy.Classic,
                < 512 => MultiplicationStrategy.Blocked,
                _ => MultiplicationStrategy.Strassen
            };
        }
    }
}