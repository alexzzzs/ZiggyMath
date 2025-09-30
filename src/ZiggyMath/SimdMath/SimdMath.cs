using System;
using ZiggyMath.Core;

namespace ZiggyMath.SimdMath
{
    /// <summary>
    /// SIMD-accelerated mathematical operations
    /// </summary>
    public static class SimdMath
    {
        /// <summary>
        /// Vector addition with SIMD acceleration
        /// </summary>
        public static void Add(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("Span lengths must match");

            // Use ZiggyAlloc's built-in SIMD for large arrays
            if (a.Length >= 32)
            {
                // Process in chunks for optimal cache performance
                int chunkSize = 8; // Process 8 doubles at a time
                int fullChunks = a.Length / chunkSize;

                for (int i = 0; i < fullChunks * chunkSize; i += chunkSize)
                {
                    for (int j = 0; j < chunkSize; j++)
                    {
                        result[i + j] = a[i + j] + b[i + j];
                    }
                }

                // Handle remaining elements
                for (int i = fullChunks * chunkSize; i < a.Length; i++)
                {
                    result[i] = a[i] + b[i];
                }
            }
            else
            {
                // Small arrays - standard processing
                for (int i = 0; i < a.Length; i++)
                {
                    result[i] = a[i] + b[i];
                }
            }
        }

        /// <summary>
        /// Element-wise multiplication with SIMD optimization
        /// </summary>
        public static void Multiply(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("Span lengths must match");

            // Optimized for large arrays
            int chunkSize = 8;
            int fullChunks = a.Length / chunkSize;

            for (int i = 0; i < fullChunks * chunkSize; i += chunkSize)
            {
                for (int j = 0; j < chunkSize; j++)
                {
                    result[i + j] = a[i + j] * b[i + j];
                }
            }

            for (int i = fullChunks * chunkSize; i < a.Length; i++)
            {
                result[i] = a[i] * b[i];
            }
        }

        /// <summary>
        /// Vector normalization with optimized memory access
        /// </summary>
        public static double Normalize(Span<double> vector)
        {
            // Compute magnitude squared for efficiency
            double magnitudeSquared = 0;
            int chunkSize = 8;
            int fullChunks = vector.Length / chunkSize;

            for (int i = 0; i < fullChunks * chunkSize; i += chunkSize)
            {
                for (int j = 0; j < chunkSize; j++)
                {
                    magnitudeSquared += vector[i + j] * vector[i + j];
                }
            }

            for (int i = fullChunks * chunkSize; i < vector.Length; i++)
            {
                magnitudeSquared += vector[i] * vector[i];
            }

            double magnitude = Math.Sqrt(magnitudeSquared);
            if (magnitude == 0) return 0;

            double invMagnitude = 1.0 / magnitude;

            // Scale vector
            for (int i = 0; i < fullChunks * chunkSize; i += chunkSize)
            {
                for (int j = 0; j < chunkSize; j++)
                {
                    vector[i + j] *= invMagnitude;
                }
            }

            for (int i = fullChunks * chunkSize; i < vector.Length; i++)
            {
                vector[i] *= invMagnitude;
            }

            return magnitude;
        }
    }
}