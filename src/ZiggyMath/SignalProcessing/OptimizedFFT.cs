using System;
using ZiggyMath.Core;
using ZiggyAlloc;

namespace ZiggyMath.SignalProcessing
{
    /// <summary>
    /// High-performance FFT implementation using ZiggyAlloc optimization
    /// </summary>
    public static class OptimizedFFT
    {
        /// <summary>
        /// Fast Fourier Transform using Cooley-Tukey algorithm with ZiggyAlloc optimization
        /// </summary>
        public static ComplexVector Transform(ComplexVector input, MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.SignalProcessing);
            int n = input.Length;

            if (n <= 1) return input;

            // Check if power of 2
            if ((n & (n - 1)) != 0)
                throw new ArgumentException("Input size must be a power of 2");

            // Use optimized memory layout for FFT
            using var optimizedInput = CreateOptimizedComplexVector(n, ctx);
            for (int i = 0; i < n; i++)
            {
                optimizedInput[i] = input[i];
            }

            // Recursive FFT implementation
            return RecursiveFFT(optimizedInput, ctx);
        }

        private static ComplexVector RecursiveFFT(ComplexVector input, MathComputationContext context)
        {
            int n = input.Length;
            if (n <= 1) return input;

            // Split into even and odd components - avoid using statements to prevent premature disposal
            var even = ExtractEvenElements(input, context);
            var odd = ExtractOddElements(input, context);

            // Recursive calls
            var fftEven = RecursiveFFT(even, context);
            var fftOdd = RecursiveFFT(odd, context);

            // Combine results with optimized memory usage
            var result = new ComplexVector(n, context.Allocator);

            for (int k = 0; k < n / 2; k++)
            {
                double angle = -2 * Math.PI * k / n;
                Complex twiddleFactor = new Complex(Math.Cos(angle), Math.Sin(angle));

                var t = twiddleFactor * fftOdd[k];
                result[k] = fftEven[k] + t;
                result[k + n / 2] = fftEven[k] - t;
            }

            // Dispose the intermediate vectors - even and odd are managed by caller
            fftEven.Dispose();
            fftOdd.Dispose();

            return result;
        }

        private static ComplexVector ExtractEvenElements(ComplexVector input, MathComputationContext context)
        {
            int n = input.Length / 2;
            var result = new ComplexVector(n, context.Allocator);

            for (int i = 0; i < n; i++)
            {
                result[i] = input[2 * i];
            }

            return result;
        }

        private static ComplexVector ExtractOddElements(ComplexVector input, MathComputationContext context)
        {
            int n = input.Length / 2;
            var result = new ComplexVector(n, context.Allocator);

            for (int i = 0; i < n; i++)
            {
                result[i] = input[2 * i + 1];
            }

            return result;
        }

        private static ComplexVector CreateOptimizedComplexVector(int size, MathComputationContext context)
        {
            // Use the context's optimized allocator for FFT operations
            return new ComplexVector(size, context.Allocator);
        }
    }
}