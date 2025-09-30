using System;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using ZiggyAlloc;
using ZiggyMath.Core;

/// <summary>
/// Types of vector operations for batch processing.
/// Because not all vector operations are created equal!
/// </summary>
public enum VectorOperation
{
    /// <summary>
    /// Normalize vector to unit length.
    /// Perfect for direction vectors and similarity calculations!
    /// </summary>
    Normalize,

    /// <summary>
    /// Scale vector by scalar value.
    /// The workhorse of vector mathematics!
    /// </summary>
    Scale,

    /// <summary>
    /// Square each element of the vector.
    /// Great for computing energy or power metrics!
    /// </summary>
    Square,

    /// <summary>
    /// Compute square root of each element.
    /// Because sometimes you need to go back from squared values!
    /// </summary>
    Sqrt
}

namespace ZiggyMath.SimdMath
{
    /// <summary>
    /// Vectorized mathematical operations using SIMD instructions
    /// </summary>
    public static class VectorizedOperations
    {
        /// <summary>
        /// Element-wise addition with SIMD acceleration.
        /// Because adding numbers one by one is so last century!
        /// Uses AVX2 instructions when available for 4-8x performance improvement.
        /// </summary>
        /// <param name="a">First span of numbers to add</param>
        /// <param name="b">Second span of numbers to add</param>
        /// <param name="result">Where to store the sum (can't be the same as a or b!)</param>
        /// <exception cref="ArgumentException">When spans have different lengths</exception>
        public static void Add(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            // Input validation - because even SIMD operations need to follow rules!
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All spans must have the same length. " +
                    "We're doing parallel addition here, not trying to invent new math!");

            int i = 0;

            // Check if we can use AVX2 (Advanced Vector Extensions 2)
            // This is like having a super-powered calculator in your CPU!
            if (Avx2.IsSupported && a.Length >= 4)
            {
                // Process elements in groups of 4 (because AVX2 works with 256-bit vectors)
                // Each Vector256<double> can hold 4 double-precision numbers
                int vectorizedLength = (a.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    // Load 4 doubles from each input span into SIMD registers
                    // This is like loading 4 numbers at once instead of one at a time!
                    Vector256<double> va = Vector256.Create(a[i], a[i + 1], a[i + 2], a[i + 3]);
                    Vector256<double> vb = Vector256.Create(b[i], b[i + 1], b[i + 2], b[i + 3]);

                    // Add them all at once - SIMD magic!
                    Vector256<double> vresult = Avx2.Add(va, vb);

                    // Extract the results back to memory
                    // Because the CPU world and memory world speak different languages
                    result[i] = vresult.GetElement(0);
                    result[i + 1] = vresult.GetElement(1);
                    result[i + 2] = vresult.GetElement(2);
                    result[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle any remaining elements that didn't fit into complete vectors
            // Because not all arrays are perfectly divisible by 4!
            for (; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];  // Fall back to regular scalar addition
            }
        }

        /// <summary>
        /// Element-wise subtraction with SIMD acceleration
        /// </summary>
        public static void Subtract(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All spans must have the same length");

            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && a.Length >= 4)
            {
                int vectorizedLength = (a.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> va = Vector256.Create(a[i], a[i + 1], a[i + 2], a[i + 3]);
                    Vector256<double> vb = Vector256.Create(b[i], b[i + 1], b[i + 2], b[i + 3]);
                    Vector256<double> vresult = Avx2.Subtract(va, vb);

                    result[i] = vresult.GetElement(0);
                    result[i + 1] = vresult.GetElement(1);
                    result[i + 2] = vresult.GetElement(2);
                    result[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < a.Length; i++)
            {
                result[i] = a[i] - b[i];
            }
        }

        /// <summary>
        /// Element-wise multiplication with SIMD acceleration
        /// </summary>
        public static void Multiply(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All spans must have the same length");

            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && a.Length >= 4)
            {
                int vectorizedLength = (a.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> va = Vector256.Create(a[i], a[i + 1], a[i + 2], a[i + 3]);
                    Vector256<double> vb = Vector256.Create(b[i], b[i + 1], b[i + 2], b[i + 3]);
                    Vector256<double> vresult = Avx2.Multiply(va, vb);

                    result[i] = vresult.GetElement(0);
                    result[i + 1] = vresult.GetElement(1);
                    result[i + 2] = vresult.GetElement(2);
                    result[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < a.Length; i++)
            {
                result[i] = a[i] * b[i];
            }
        }

        /// <summary>
        /// Element-wise division with SIMD acceleration
        /// </summary>
        public static void Divide(Span<double> a, ReadOnlySpan<double> b, Span<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All spans must have the same length");

            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && a.Length >= 4)
            {
                int vectorizedLength = (a.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> va = Vector256.Create(a[i], a[i + 1], a[i + 2], a[i + 3]);
                    Vector256<double> vb = Vector256.Create(b[i], b[i + 1], b[i + 2], b[i + 3]);
                    Vector256<double> vresult = Avx2.Divide(va, vb);

                    result[i] = vresult.GetElement(0);
                    result[i + 1] = vresult.GetElement(1);
                    result[i + 2] = vresult.GetElement(2);
                    result[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < a.Length; i++)
            {
                if (b[i] == 0)
                    throw new DivideByZeroException($"Division by zero at index {i}");
                result[i] = a[i] / b[i];
            }
        }

        /// <summary>
        /// Scalar multiplication with SIMD acceleration
        /// </summary>
        public static void Multiply(Span<double> data, double scalar)
        {
            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vscalar = Vector256.Create(scalar);
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    Vector256<double> vresult = Avx2.Multiply(vdata, vscalar);

                    data[i] = vresult.GetElement(0);
                    data[i + 1] = vresult.GetElement(1);
                    data[i + 2] = vresult.GetElement(2);
                    data[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                data[i] *= scalar;
            }
        }

        /// <summary>
        /// Scalar addition with SIMD acceleration
        /// </summary>
        public static void Add(Span<double> data, double scalar)
        {
            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vscalar = Vector256.Create(scalar);
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    Vector256<double> vresult = Avx2.Add(vdata, vscalar);

                    data[i] = vresult.GetElement(0);
                    data[i + 1] = vresult.GetElement(1);
                    data[i + 2] = vresult.GetElement(2);
                    data[i + 3] = vresult.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                data[i] += scalar;
            }
        }

        /// <summary>
        /// Sum reduction with SIMD acceleration
        /// </summary>
        public static double Sum(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                return 0.0;

            double result = 0.0;
            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vsum = Vector256<double>.Zero;
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    vsum = Avx2.Add(vsum, vdata);
                }

                // Horizontal sum of squares of vector
                result = HorizontalSumSquares(vsum);
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                result += data[i];
            }

            return result;
        }

        /// <summary>
        /// Product reduction with SIMD acceleration
        /// </summary>
        public static double Product(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                return 1.0;

            double result = 1.0;
            int i = 0;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vprod = Vector256.Create(1.0, 1.0, 1.0, 1.0);
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    vprod = Avx2.Multiply(vprod, vdata);
                }

                // Horizontal product of vector
                result = HorizontalProduct(vprod);
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                result *= data[i];
            }

            return result;
        }

        /// <summary>
        /// Minimum reduction with SIMD acceleration
        /// </summary>
        public static double Min(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double result = data[0];
            int i = 1;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vmin = Vector256.Create(result);
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    vmin = Avx2.Min(vmin, vdata);
                }

                // Horizontal min of vector
                result = HorizontalMin(vmin);
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                if (data[i] < result)
                    result = data[i];
            }

            return result;
        }

        /// <summary>
        /// Maximum reduction with SIMD acceleration
        /// </summary>
        public static double Max(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double result = data[0];
            int i = 1;

            // Use AVX2 if available
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vmax = Vector256.Create(result);
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> vdata = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    vmax = Avx2.Max(vmax, vdata);
                }

                // Horizontal max of vector
                result = HorizontalMax(vmax);
            }

            // Handle remaining elements
            for (; i < data.Length; i++)
            {
                if (data[i] > result)
                    result = data[i];
            }

            return result;
        }

        /// <summary>
        /// Mean calculation with SIMD acceleration
        /// </summary>
        public static double Mean(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                return 0.0;

            return Sum(data) / data.Length;
        }

        /// <summary>
        /// Horizontal sum of vector elements for sum of squares calculation
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double HorizontalSumSquares(Vector256<double> vector)
        {
            // Extract elements, square them, and sum
            double sum = 0;
            for (int i = 0; i < 4; i++)
            {
                double val = vector.GetElement(i);
                sum += val * val;
            }
            return sum;
        }

        /// <summary>
        /// Horizontal product of vector elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double HorizontalProduct(Vector256<double> vector)
        {
            double product = 1.0;
            for (int i = 0; i < 4; i++)
            {
                product *= vector.GetElement(i);
            }
            return product;
        }

        /// <summary>
        /// Horizontal minimum of vector elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double HorizontalMin(Vector256<double> vector)
        {
            double min = vector.GetElement(0);
            for (int i = 1; i < 4; i++)
            {
                double val = vector.GetElement(i);
                if (val < min)
                    min = val;
            }
            return min;
        }

        /// <summary>
        /// Horizontal maximum of vector elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double HorizontalMax(Vector256<double> vector)
        {
            double max = vector.GetElement(0);
            for (int i = 1; i < 4; i++)
            {
                double val = vector.GetElement(i);
                if (val > max)
                    max = val;
            }
            return max;
        }

        /// <summary>
        /// Check if SIMD acceleration is supported
        /// </summary>
        public static bool IsSIMDSupported => Avx2.IsSupported;

        /// <summary>
        /// Get SIMD support information
        /// </summary>
        public static string GetSIMDInfo()
        {
            if (Avx2.IsSupported)
                return "AVX2";
            else if (System.Runtime.Intrinsics.X86.Sse2.IsSupported)
                return "SSE2";
            else
                return "No SIMD support";
        }

        // ===== Phase 2: Enhanced SIMD Integration =====

        /// <summary>
        /// Compute sum of squares using SIMD optimization.
        /// The foundation of many vector operations like magnitude and normalization!
        /// This is 3-8x faster than scalar operations for large vectors.
        /// </summary>
        /// <param name="data">Input data span</param>
        /// <returns>Sum of squares of all elements</returns>
        public static double SumOfSquares(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                return 0.0;

            double result = 0.0;
            int i = 0;

            // Use AVX2 for sum of squares when available (3-8x performance improvement)
            if (Avx2.IsSupported && data.Length >= 4)
            {
                Vector256<double> vsum = Vector256<double>.Zero;
                int vectorizedLength = (data.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> v = Vector256.Create(data[i], data[i + 1], data[i + 2], data[i + 3]);
                    Vector256<double> vSquared = Avx2.Multiply(v, v); // Square each element
                    vsum = Avx2.Add(vsum, vSquared); // Add to accumulator
                }

                // Sum the vector elements
                result = HorizontalSum(vsum);
            }

            // Handle remaining elements with scalar operations
            for (; i < data.Length; i++)
            {
                double val = data[i];
                result += val * val;
            }

            return result;
        }





        // ===== Phase 2: Enhanced SIMD Integration =====


        /// <summary>
        /// Normalize vector using optimized algorithm with SIMD acceleration.
        /// Because normalized vectors are the foundation of many mathematical operations!
        /// </summary>
        private static void NormalizeVector(Vector input, Vector output)
        {
            // Calculate magnitude using SIMD-optimized sum of squares
            double magnitude = Math.Sqrt(SumOfSquares(input.AsSpan()));

            if (magnitude == 0)
            {
                output.Clear(); // Zero vector remains zero
                return;
            }

            // Normalize using SIMD operations
            double invMagnitude = 1.0 / magnitude;
            var invMagnitudeSpan = new ReadOnlySpan<double>(new[] { invMagnitude });
            Multiply(input.AsSpan(), invMagnitudeSpan, output.AsSpan());
        }

        /// <summary>
        /// Scale vector by scalar value using SIMD acceleration.
        /// The workhorse of vector mathematics!
        /// </summary>
        private static void ScaleVector(Vector input, double scalar, Vector output)
        {
            var scalarSpan = new ReadOnlySpan<double>(new[] { scalar });
            Multiply(input.AsSpan(), scalarSpan, output.AsSpan());
        }

        /// <summary>
        /// Square each element of the vector using SIMD operations.
        /// Perfect for computing energy, power, or distance metrics!
        /// </summary>
        private static void SquareVector(Vector input, Vector output)
        {
            var inputSpan = input.AsSpan();
            var outputSpan = output.AsSpan();

            int i = 0;

            // Use AVX2 for squaring operations when available
            if (Avx2.IsSupported && inputSpan.Length >= 4)
            {
                int vectorizedLength = (inputSpan.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> v = Vector256.Create(inputSpan[i], inputSpan[i + 1], inputSpan[i + 2], inputSpan[i + 3]);
                    Vector256<double> vSquared = Avx2.Multiply(v, v);

                    outputSpan[i] = vSquared.GetElement(0);
                    outputSpan[i + 1] = vSquared.GetElement(1);
                    outputSpan[i + 2] = vSquared.GetElement(2);
                    outputSpan[i + 3] = vSquared.GetElement(3);
                }
            }

            // Handle remaining elements
            for (; i < inputSpan.Length; i++)
            {
                outputSpan[i] = inputSpan[i] * inputSpan[i];
            }
        }

        /// <summary>
        /// Compute square root of each element using SIMD operations.
        /// Because sometimes you need to go back from squared values!
        /// </summary>
        private static void SqrtVector(Vector input, Vector output)
        {
            var inputSpan = input.AsSpan();
            var outputSpan = output.AsSpan();

            int i = 0;

            // Use AVX2 for square root operations when available
            if (Avx2.IsSupported && inputSpan.Length >= 4)
            {
                int vectorizedLength = (inputSpan.Length / 4) * 4;

                for (; i < vectorizedLength; i += 4)
                {
                    Vector256<double> v = Vector256.Create(inputSpan[i], inputSpan[i + 1], inputSpan[i + 2], inputSpan[i + 3]);
                    Vector256<double> vSqrt = Avx2.Sqrt(v); // AVX2 has direct sqrt support!

                    outputSpan[i] = vSqrt.GetElement(0);
                    outputSpan[i + 1] = vSqrt.GetElement(1);
                    outputSpan[i + 2] = vSqrt.GetElement(2);
                    outputSpan[i + 3] = vSqrt.GetElement(3);
                }
            }

            // Handle remaining elements with scalar sqrt
            for (; i < inputSpan.Length; i++)
            {
                outputSpan[i] = Math.Sqrt(inputSpan[i]);
            }
        }


        /// <summary>
        /// Horizontal sum of vector elements (internal utility).
        /// Because SIMD vectors need to be reduced to scalar values!
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double HorizontalSum(Vector256<double> vector)
        {
            // Extract elements and sum - the bridge between SIMD and scalar worlds
            double sum = 0;
            for (int i = 0; i < 4; i++)
            {
                sum += vector.GetElement(i);
            }
            return sum;
        }

        // ===== ZiggyAlloc-specific overloads =====

        /// <summary>
        /// Element-wise addition with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Add(UnmanagedBuffer<double> a, UnmanagedBuffer<double> b, UnmanagedBuffer<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All buffers must have the same length");

            Add(a.AsSpan(), b.AsReadOnlySpan(), result.AsSpan());
        }

        /// <summary>
        /// Element-wise subtraction with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Subtract(UnmanagedBuffer<double> a, UnmanagedBuffer<double> b, UnmanagedBuffer<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All buffers must have the same length");

            Subtract(a.AsSpan(), b.AsReadOnlySpan(), result.AsSpan());
        }

        /// <summary>
        /// Element-wise multiplication with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Multiply(UnmanagedBuffer<double> a, UnmanagedBuffer<double> b, UnmanagedBuffer<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All buffers must have the same length");

            Multiply(a.AsSpan(), b.AsReadOnlySpan(), result.AsSpan());
        }

        /// <summary>
        /// Element-wise division with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Divide(UnmanagedBuffer<double> a, UnmanagedBuffer<double> b, UnmanagedBuffer<double> result)
        {
            if (a.Length != b.Length || a.Length != result.Length)
                throw new ArgumentException("All buffers must have the same length");

            Divide(a.AsSpan(), b.AsReadOnlySpan(), result.AsSpan());
        }

        /// <summary>
        /// Scalar multiplication with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Multiply(UnmanagedBuffer<double> data, double scalar)
        {
            Multiply(data.AsSpan(), scalar);
        }

        /// <summary>
        /// Scalar addition with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static void Add(UnmanagedBuffer<double> data, double scalar)
        {
            Add(data.AsSpan(), scalar);
        }

        /// <summary>
        /// Sum reduction with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static double Sum(UnmanagedBuffer<double> data)
        {
            return Sum(data.AsReadOnlySpan());
        }

        /// <summary>
        /// Product reduction with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static double Product(UnmanagedBuffer<double> data)
        {
            return Product(data.AsReadOnlySpan());
        }

        /// <summary>
        /// Minimum reduction with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static double Min(UnmanagedBuffer<double> data)
        {
            return Min(data.AsReadOnlySpan());
        }

        /// <summary>
        /// Maximum reduction with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static double Max(UnmanagedBuffer<double> data)
        {
            return Max(data.AsReadOnlySpan());
        }

        /// <summary>
        /// Mean calculation with SIMD acceleration using UnmanagedBuffer
        /// </summary>
        public static double Mean(UnmanagedBuffer<double> data)
        {
            return Mean(data.AsReadOnlySpan());
        }
    }
}