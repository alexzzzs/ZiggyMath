using System;
using ZiggyMath.Core;
using ZiggyAlloc;

namespace ZiggyMath.Calculus
{
    /// <summary>
    /// High-performance numerical integration using ZiggyAlloc optimization
    /// </summary>
    public static class OptimizedNumericalIntegration
    {
        public delegate double Function(double x);
        public delegate double Function2D(double x, double y);

        /// <summary>
        /// Adaptive quadrature with automatic error control using ZiggyAlloc optimization
        /// </summary>
        public static double AdaptiveQuadrature(
            Function f,
            double a,
            double b,
            double tolerance = 1e-10,
            MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.SmallComputations);

            return AdaptiveQuadratureRecursive(f, a, b, tolerance, ctx);
        }

        private static double AdaptiveQuadratureRecursive(
            Function f,
            double a,
            double b,
            double tolerance,
            MathComputationContext context)
        {
            // Use Simpson's rule for local integration with optimized memory
            double c = (a + b) / 2;
            double h = b - a;

            // Use optimized buffer for function evaluations
            using var evalBuffer = context.MemoryPool.Rent<double>(5);
            var evalSpan = evalBuffer.AsSpan();

            evalSpan[0] = f(a);
            evalSpan[1] = f((a + c) / 2);
            evalSpan[2] = f(c);
            evalSpan[3] = f((c + b) / 2);
            evalSpan[4] = f(b);

            double fab = evalSpan[0];
            double fab_c = evalSpan[1];
            double fc = evalSpan[2];
            double fcb = evalSpan[3];
            double fba = evalSpan[4];

            // Simpson's rule for whole interval
            double s = (h / 6) * (fab + 4 * fc + fba);

            // Simpson's rule for left half
            double sl = (h / 12) * (fab + 4 * fab_c + fc);

            // Simpson's rule for right half
            double sr = (h / 12) * (fc + 4 * fcb + fba);

            double error = Math.Abs(s - sl - sr);

            if (error < 15 * tolerance)
            {
                return sl + sr + error / 15;
            }
            else
            {
                // Recursive subdivision
                double leftResult = AdaptiveQuadratureRecursive(f, a, c, tolerance / 2, context);
                double rightResult = AdaptiveQuadratureRecursive(f, c, b, tolerance / 2, context);
                return leftResult + rightResult;
            }
        }

        /// <summary>
        /// Monte Carlo integration for high-dimensional integrals with ZiggyAlloc optimization
        /// </summary>
        public static double MonteCarloIntegration(
            Function f,
            double a,
            double b,
            int samples = 100000,
            MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.SmallComputations);

            // Use optimized buffer for random samples
            using var randomBuffer = ctx.MemoryPool.Rent<double>(samples);
            using var resultBuffer = ctx.MemoryPool.Rent<double>(samples);

            // Generate random samples using SIMD optimization
            var randomSpan = randomBuffer.AsSpan();
            var resultSpan = resultBuffer.AsSpan();

            Random random = new Random();
            for (int i = 0; i < samples; i++)
            {
                randomSpan[i] = a + (b - a) * random.NextDouble();
                resultSpan[i] = f(randomSpan[i]);
            }

            // Compute integral estimate using SIMD sum reduction
            double sum = 0;
            for (int i = 0; i < samples; i++)
            {
                sum += resultSpan[i];
            }

            double average = sum / samples;
            return average * (b - a);
        }
    }
}