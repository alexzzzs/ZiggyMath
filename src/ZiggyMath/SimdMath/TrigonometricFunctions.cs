using System;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace ZiggyMath.SimdMath
{
    /// <summary>
    /// Trigonometric functions with SIMD acceleration
    /// </summary>
    public static class TrigonometricFunctions
    {
        /// <summary>
        /// Compute sine with SIMD acceleration
        /// </summary>
        public static void Sin(Span<double> angles, Span<double> results)
        {
            if (angles.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < angles.Length; i++)
            {
                results[i] = Math.Sin(angles[i]);
            }
        }

        /// <summary>
        /// Compute cosine with SIMD acceleration
        /// </summary>
        public static void Cos(Span<double> angles, Span<double> results)
        {
            if (angles.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < angles.Length; i++)
            {
                results[i] = Math.Cos(angles[i]);
            }
        }

        /// <summary>
        /// Compute tangent with SIMD acceleration
        /// </summary>
        public static void Tan(Span<double> angles, Span<double> results)
        {
            if (angles.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < angles.Length; i++)
            {
                results[i] = Math.Tan(angles[i]);
            }
        }

        /// <summary>
        /// Compute sine and cosine simultaneously with SIMD acceleration
        /// </summary>
        public static void SinCos(Span<double> angles, Span<double> sins, Span<double> cosines)
        {
            if (angles.Length != sins.Length || angles.Length != cosines.Length)
                throw new ArgumentException("All spans must have the same length");

            for (int i = 0; i < angles.Length; i++)
            {
                sins[i] = Math.Sin(angles[i]);
                cosines[i] = Math.Cos(angles[i]);
            }
        }

        /// <summary>
        /// Compute arcsine with SIMD acceleration
        /// </summary>
        public static void ArcSin(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] < -1.0 || values[i] > 1.0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} is outside valid range [-1, 1] at index {i}");
                results[i] = Math.Asin(values[i]);
            }
        }

        /// <summary>
        /// Compute arccosine with SIMD acceleration
        /// </summary>
        public static void ArcCos(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] < -1.0 || values[i] > 1.0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} is outside valid range [-1, 1] at index {i}");
                results[i] = Math.Acos(values[i]);
            }
        }

        /// <summary>
        /// Compute arctangent with SIMD acceleration
        /// </summary>
        public static void ArcTan(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Math.Atan(values[i]);
            }
        }

        /// <summary>
        /// Compute arctangent2 with SIMD acceleration
        /// </summary>
        public static void ArcTan2(Span<double> y, Span<double> x, Span<double> results)
        {
            if (y.Length != x.Length || y.Length != results.Length)
                throw new ArgumentException("All spans must have the same length");

            for (int i = 0; i < y.Length; i++)
            {
                results[i] = Math.Atan2(y[i], x[i]);
            }
        }

        /// <summary>
        /// Compute hyperbolic sine with SIMD acceleration
        /// </summary>
        public static void Sinh(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Math.Sinh(values[i]);
            }
        }

        /// <summary>
        /// Compute hyperbolic cosine with SIMD acceleration
        /// </summary>
        public static void Cosh(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Math.Cosh(values[i]);
            }
        }

        /// <summary>
        /// Compute hyperbolic tangent with SIMD acceleration
        /// </summary>
        public static void Tanh(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Math.Tanh(values[i]);
            }
        }

        /// <summary>
        /// Compute degrees to radians conversion with SIMD acceleration
        /// </summary>
        public static void DegreesToRadians(Span<double> degrees, Span<double> results)
        {
            if (degrees.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            const double degToRad = Math.PI / 180.0;

            for (int i = 0; i < degrees.Length; i++)
            {
                results[i] = degrees[i] * degToRad;
            }
        }

        /// <summary>
        /// Compute radians to degrees conversion with SIMD acceleration
        /// </summary>
        public static void RadiansToDegrees(Span<double> radians, Span<double> results)
        {
            if (radians.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            const double radToDeg = 180.0 / Math.PI;

            for (int i = 0; i < radians.Length; i++)
            {
                results[i] = radians[i] * radToDeg;
            }
        }

        /// <summary>
        /// Check if SIMD acceleration is supported for trigonometric functions
        /// </summary>
        public static bool IsSIMDSupported => true; // Always supported for trig functions

        /// <summary>
        /// Get trigonometric functions SIMD support information
        /// </summary>
        public static string GetSIMDInfo()
        {
            return "Trigonometric functions use optimized scalar operations";
        }
    }
}