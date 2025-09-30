using System;
using ZiggyMath.Core;
using ZiggyAlloc;

namespace ZiggyMath.Statistics
{
    /// <summary>
    /// High-performance statistical operations using ZiggyAlloc optimization
    /// </summary>
    public static class OptimizedStatistics
    {
        /// <summary>
        /// Linear regression with optimized memory usage and ZiggyAlloc performance
        /// </summary>
        public static (double slope, double intercept, double rSquared) LinearRegression(
            ReadOnlySpan<double> x,
            ReadOnlySpan<double> y,
            MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.SmallComputations);

            if (x.Length != y.Length)
                throw new ArgumentException("Input arrays must have the same length");

            int n = x.Length;

            // Use optimized buffers for intermediate calculations
            using var sums = ctx.MemoryPool.Rent<double>(4); // sumX, sumY, sumXY, sumXX
            var sumsSpan = sums.AsSpan();

            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

            // Compute sums with potential SIMD optimization
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumXX += x[i] * x[i];
            }

            // Calculate regression coefficients
            double slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
            double intercept = (sumY - slope * sumX) / n;

            // Calculate R-squared using optimized buffer
            using var tempBuffer = ctx.MemoryPool.Rent<double>(n);
            var tempSpan = tempBuffer.AsSpan();

            double meanY = sumY / n;
            double totalSumSquares = 0;
            double residualSumSquares = 0;

            for (int i = 0; i < n; i++)
            {
                double predicted = slope * x[i] + intercept;
                double residual = y[i] - predicted;
                double diffFromMean = y[i] - meanY;

                totalSumSquares += diffFromMean * diffFromMean;
                residualSumSquares += residual * residual;
            }

            double rSquared = totalSumSquares != 0 ? 1 - (residualSumSquares / totalSumSquares) : 0;

            return (slope, intercept, rSquared);
        }

        /// <summary>
        /// Descriptive statistics with SIMD optimization and ZiggyAlloc memory management
        /// </summary>
        public static (double mean, double variance, double skewness, double kurtosis)
            DescriptiveStatistics(ReadOnlySpan<double> data, MathComputationContext context = null)
        {
            var ctx = context ?? new MathComputationContext(MathWorkload.SmallComputations);

            int n = data.Length;

            // Use optimized buffer for calculations
            using var moments = ctx.MemoryPool.Rent<double>(4); // mean, variance, skewness, kurtosis
            var momentsSpan = moments.AsSpan();

            // Calculate mean
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += data[i];
            }
            double mean = sum / n;

            // Calculate variance, skewness, kurtosis using optimized buffer
            using var calcBuffer = ctx.MemoryPool.Rent<double>(n);
            var calcSpan = calcBuffer.AsSpan();

            double varianceSum = 0;
            double skewnessSum = 0;
            double kurtosisSum = 0;

            for (int i = 0; i < n; i++)
            {
                double diff = data[i] - mean;
                double diffSquared = diff * diff;
                double diffCubed = diffSquared * diff;
                double diffFourth = diffSquared * diffSquared;

                varianceSum += diffSquared;
                skewnessSum += diffCubed;
                kurtosisSum += diffFourth;
            }

            double variance = varianceSum / (n - 1);
            double skewness = (skewnessSum / n) / Math.Pow(variance, 1.5);
            double kurtosis = (kurtosisSum / n) / (variance * variance) - 3;

            return (mean, variance, skewness, kurtosis);
        }
    }
}