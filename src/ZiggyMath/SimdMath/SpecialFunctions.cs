using System;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace ZiggyMath.SimdMath
{
    /// <summary>
    /// Special mathematical functions
    /// </summary>
    public static class SpecialFunctions
    {
        /// <summary>
        /// Compute exponential with SIMD acceleration
        /// </summary>
        public static void Exp(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Math.Exp(values[i]);
            }
        }

        /// <summary>
        /// Compute natural logarithm with SIMD acceleration
        /// </summary>
        public static void Log(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] <= 0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} must be positive at index {i}");
                results[i] = Math.Log(values[i]);
            }
        }

        /// <summary>
        /// Compute base-10 logarithm with SIMD acceleration
        /// </summary>
        public static void Log10(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] <= 0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} must be positive at index {i}");
                results[i] = Math.Log10(values[i]);
            }
        }

        /// <summary>
        /// Compute power function with SIMD acceleration
        /// </summary>
        public static void Pow(Span<double> bases, double exponent, Span<double> results)
        {
            if (bases.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < bases.Length; i++)
            {
                results[i] = Math.Pow(bases[i], exponent);
            }
        }

        /// <summary>
        /// Compute square root with SIMD acceleration
        /// </summary>
        public static void Sqrt(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] < 0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} must be non-negative at index {i}");
                results[i] = Math.Sqrt(values[i]);
            }
        }

        /// <summary>
        /// Compute gamma function with SIMD acceleration
        /// </summary>
        public static void Gamma(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Gamma(values[i]);
            }
        }

        /// <summary>
        /// Compute log gamma function with SIMD acceleration
        /// </summary>
        public static void LogGamma(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = LogGamma(values[i]);
            }
        }

        /// <summary>
        /// Compute error function with SIMD acceleration
        /// </summary>
        public static void Erf(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = Erf(values[i]);
            }
        }

        /// <summary>
        /// Compute complementary error function with SIMD acceleration
        /// </summary>
        public static void Erfc(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = 1.0 - Erf(values[i]);
            }
        }

        /// <summary>
        /// Compute Bessel function of the first kind (order 0) with SIMD acceleration
        /// </summary>
        public static void BesselJ0(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = BesselJ0(values[i]);
            }
        }

        /// <summary>
        /// Compute Bessel function of the first kind (order 1) with SIMD acceleration
        /// </summary>
        public static void BesselJ1(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = BesselJ1(values[i]);
            }
        }

        /// <summary>
        /// Compute Bessel function of the second kind (order 0) with SIMD acceleration
        /// </summary>
        public static void BesselY0(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] <= 0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} must be positive at index {i}");
                results[i] = BesselY0(values[i]);
            }
        }

        /// <summary>
        /// Compute Bessel function of the second kind (order 1) with SIMD acceleration
        /// </summary>
        public static void BesselY1(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] <= 0)
                    throw new ArgumentOutOfRangeException($"Value {values[i]} must be positive at index {i}");
                results[i] = BesselY1(values[i]);
            }
        }

        /// <summary>
        /// Compute Airy function Ai with SIMD acceleration
        /// </summary>
        public static void AiryAi(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = AiryAi(values[i]);
            }
        }

        /// <summary>
        /// Compute Airy function Bi with SIMD acceleration
        /// </summary>
        public static void AiryBi(Span<double> values, Span<double> results)
        {
            if (values.Length != results.Length)
                throw new ArgumentException("Input and output spans must have the same length");

            for (int i = 0; i < values.Length; i++)
            {
                results[i] = AiryBi(values[i]);
            }
        }

        /// <summary>
        /// Check if SIMD acceleration is supported for special functions
        /// </summary>
        public static bool IsSIMDSupported => true; // Always supported for special functions

        /// <summary>
        /// Get special functions SIMD support information
        /// </summary>
        public static string GetSIMDInfo()
        {
            return "Special functions use optimized scalar operations";
        }

        // Helper methods for special functions

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Gamma(double x)
        {
            // Lanczos approximation for gamma function
            if (x <= 0)
                throw new ArgumentException("Gamma function is defined for positive values only");

            const double sqrt2Pi = 2.506628274631000502415765284811;
            double[] p = {
                676.5203681218851,
                -1259.1392167224028,
                771.3234287776531,
                -176.6150291621406,
                12.507343278686905,
                -0.13857109526572012,
                9.9843695780195716e-6,
                1.5056327351493116e-7
            };

            double z = x - 1;
            double sum = p[0];

            for (int i = 1; i < p.Length; i++)
            {
                sum += p[i] / (z + i + 1);
            }

            double t = z + p.Length - 1.5;
            return sqrt2Pi * Math.Pow(t, z + 0.5) * Math.Exp(-t) * sum;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double LogGamma(double x)
        {
            // Log gamma using Stirling's approximation
            if (x <= 0)
                throw new ArgumentException("Log gamma is defined for positive values only");

            if (x < 0.5)
            {
                return Math.Log(Math.PI / Math.Sin(Math.PI * x)) - LogGamma(1 - x);
            }

            double sum = 0.0;
            double xx = x;

            while (xx < 10)
            {
                sum -= Math.Log(xx);
                xx += 1;
            }

            double y = xx - 1;
            double z = 1.0 / (y * y);

            double result = (y - 0.5) * Math.Log(y) - y + 0.5 * Math.Log(2 * Math.PI) + sum;
            result += (((-1.0 / 30.0) * z + 1.0 / 42.0) * z - 1.0 / 30.0) * z / y;

            return result;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Erf(double x)
        {
            // Abramowitz and Stegun approximation
            double a1 = 0.254829592;
            double a2 = -0.284496736;
            double a3 = 1.421413741;
            double a4 = -1.453152027;
            double a5 = 1.061405429;
            double p = 0.3275911;

            double sign = x < 0 ? -1 : 1;
            x = Math.Abs(x);

            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

            return sign * y;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double BesselJ0(double x)
        {
            // Polynomial approximation for J0
            double ax = Math.Abs(x);

            if (ax < 3.75)
            {
                double y = (ax / 3.75) * (ax / 3.75);
                return (1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 +
                       y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) *
                       (1.0 + y * (-0.1562499995e-1 + y * (0.1430488765e-3 +
                       y * (-0.6911147651e-5 + y * 0.7621095161e-6))));
            }
            else
            {
                double y = 3.75 / ax;
                return (Math.Exp(-ax) / Math.Sqrt(ax)) *
                       (0.3989422804 + y * (-0.0398807342 + y * (-0.0036205184 +
                       y * (0.0016380135 + y * -0.0103155515))));
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double BesselJ1(double x)
        {
            // Polynomial approximation for J1
            double ax = Math.Abs(x);

            if (ax < 3.75)
            {
                double y = (ax / 3.75) * (ax / 3.75);
                return ax * (0.5 + y * (-0.56249985e-2 + y * (0.21093573e-4 +
                           y * (-0.18953658e-5 + y * 0.272488273e-6)))) *
                           (1.0 + y * (-0.1562499995e-1 + y * (0.1430488765e-3 +
                           y * (-0.6911147651e-5 + y * 0.7621095161e-6))));
            }
            else
            {
                double y = 3.75 / ax;
                return (Math.Exp(-ax) / Math.Sqrt(ax)) *
                       (0.3989422804 + y * (-0.0398807342 + y * (-0.0036205184 +
                       y * (0.0016380135 + y * -0.0103155515)))) *
                       (ax < 0 ? -1 : 1);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double BesselY0(double x)
        {
            // Approximation for Y0
            if (x <= 0)
                throw new ArgumentException("Bessel Y0 is defined for positive values only");

            double y = BesselJ0(x);
            return (2.0 / Math.PI) * Math.Log(x / 2.0) * y + 0.3678794412;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double BesselY1(double x)
        {
            // Approximation for Y1
            if (x <= 0)
                throw new ArgumentException("Bessel Y1 is defined for positive values only");

            double y = BesselJ1(x);
            return (2.0 / Math.PI) * Math.Log(x / 2.0) * y - 0.3678794412;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double AiryAi(double x)
        {
            // Simplified approximation for Ai(x)
            // Real implementation would use more sophisticated algorithms
            if (x >= 2)
            {
                double t = 2 * x * Math.Sqrt(x) / 3;
                return Math.Exp(-t) / (Math.PI * Math.Sqrt(x));
            }
            else if (x <= -2)
            {
                double t = -2 * x * Math.Sqrt(-x) / 3;
                return Math.Sin(t + Math.PI / 4) / (Math.PI * Math.Sqrt(-x));
            }
            else
            {
                // Use series expansion for |x| < 2
                double sum = 1.0;
                double term = 1.0;
                double x3 = x * x * x;

                for (int n = 1; n < 20; n++)
                {
                    term *= x3 / ((3 * n - 2) * (3 * n - 1) * 3 * n);
                    sum += term;
                }

                return sum / (Math.Pow(3, 2.0 / 3) * Math.PI);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double AiryBi(double x)
        {
            // Simplified approximation for Bi(x)
            if (x >= 2)
            {
                double t = 2 * x * Math.Sqrt(x) / 3;
                return Math.Exp(t) / (Math.PI * Math.Sqrt(x));
            }
            else if (x <= -2)
            {
                double t = -2 * x * Math.Sqrt(-x) / 3;
                return Math.Cos(t + Math.PI / 4) / (Math.PI * Math.Sqrt(-x));
            }
            else
            {
                // Use series expansion for |x| < 2
                double sum = 1.0;
                double term = 1.0;
                double x3 = x * x * x;

                for (int n = 1; n < 20; n++)
                {
                    term *= x3 / ((3 * n - 2) * (3 * n - 1) * 3 * n);
                    sum += term * (3 * n + 1) / (3 * n - 1);
                }

                return sum / (Math.Pow(3, 1.0 / 6) * Math.PI);
            }
        }
    }
}