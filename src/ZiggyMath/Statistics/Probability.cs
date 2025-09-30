using System;
using ZiggyMath.Core;

namespace ZiggyMath.Statistics
{
    /// <summary>
    /// Probability distributions and functions
    /// </summary>
    public static class Probability
    {
        /// <summary>
        /// Normal probability density function
        /// </summary>
        public static double NormalPDF(double x, double mean, double stdDev)
        {
            if (stdDev <= 0)
                throw new ArgumentException("Standard deviation must be positive");

            double diff = x - mean;
            double variance = stdDev * stdDev;
            double exponent = -diff * diff / (2 * variance);

            return Math.Exp(exponent) / Math.Sqrt(2 * Math.PI * variance);
        }

        /// <summary>
        /// Normal cumulative distribution function (approximation)
        /// </summary>
        public static double NormalCDF(double x, double mean, double stdDev)
        {
            if (stdDev <= 0)
                throw new ArgumentException("Standard deviation must be positive");

            double z = (x - mean) / stdDev;
            return 0.5 * (1 + Erf(z / Math.Sqrt(2)));
        }

        /// <summary>
        /// Inverse normal CDF (approximation)
        /// </summary>
        public static double InverseNormalCDF(double p, double mean, double stdDev)
        {
            if (p <= 0 || p >= 1)
                throw new ArgumentException("Probability must be between 0 and 1");
            if (stdDev <= 0)
                throw new ArgumentException("Standard deviation must be positive");

            // Approximation using rational function
            double z;
            if (p < 0.5)
            {
                z = -InverseErf(2 * p);
            }
            else
            {
                z = InverseErf(2 * (1 - p));
            }

            return mean + stdDev * z * Math.Sqrt(2);
        }

        /// <summary>
        /// Chi-square test for goodness of fit
        /// </summary>
        public static double ChiSquareTest(double[] observed, double[] expected)
        {
            if (observed == null || expected == null)
                throw new ArgumentException("Observed and expected arrays cannot be null");
            if (observed.Length != expected.Length)
                throw new ArgumentException("Arrays must have the same length");

            double chiSquare = 0.0;
            int degreesOfFreedom = 0;

            for (int i = 0; i < observed.Length; i++)
            {
                if (expected[i] > 0)
                {
                    chiSquare += (observed[i] - expected[i]) * (observed[i] - expected[i]) / expected[i];
                    degreesOfFreedom++;
                }
            }

            return chiSquare;
        }

        /// <summary>
        /// Student's t-test for two samples
        /// </summary>
        public static double TTest(double[] sample1, double[] sample2)
        {
            if (sample1 == null || sample2 == null)
                throw new ArgumentException("Samples cannot be null");
            if (sample1.Length < 2 || sample2.Length < 2)
                throw new ArgumentException("Each sample must have at least 2 observations");

            double mean1 = DescriptiveStats.Mean(sample1);
            double mean2 = DescriptiveStats.Mean(sample2);
            double var1 = DescriptiveStats.Variance(sample1.AsSpan());
            double var2 = DescriptiveStats.Variance(sample2.AsSpan());

            double pooledVar = ((sample1.Length - 1) * var1 + (sample2.Length - 1) * var2) /
                              (sample1.Length + sample2.Length - 2);

            double tStatistic = (mean1 - mean2) / Math.Sqrt(pooledVar * (1.0 / sample1.Length + 1.0 / sample2.Length));

            return tStatistic;
        }

        /// <summary>
        /// F-test for equality of variances
        /// </summary>
        public static double FTest(double[] sample1, double[] sample2)
        {
            if (sample1 == null || sample2 == null)
                throw new ArgumentException("Samples cannot be null");
            if (sample1.Length < 2 || sample2.Length < 2)
                throw new ArgumentException("Each sample must have at least 2 observations");

            double var1 = DescriptiveStats.Variance(sample1.AsSpan());
            double var2 = DescriptiveStats.Variance(sample2.AsSpan());

            if (var1 == 0 && var2 == 0)
                return 1.0;

            double fStatistic = var1 > var2 ? var1 / var2 : var2 / var1;

            return fStatistic;
        }

        /// <summary>
        /// Generate random samples from normal distribution
        /// </summary>
        public static double[] RandomNormal(int count, double mean, double stdDev)
        {
            if (count <= 0)
                throw new ArgumentException("Count must be positive");
            if (stdDev <= 0)
                throw new ArgumentException("Standard deviation must be positive");

            Random random = new Random();
            double[] samples = new double[count];

            // Box-Muller transformation
            for (int i = 0; i < count - 1; i += 2)
            {
                double u1 = random.NextDouble();
                double u2 = random.NextDouble();

                double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                double z1 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);

                samples[i] = mean + stdDev * z0;
                samples[i + 1] = mean + stdDev * z1;
            }

            // Handle odd count
            if (count % 2 == 1)
            {
                double u1 = random.NextDouble();
                double u2 = random.NextDouble();
                double z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                samples[count - 1] = mean + stdDev * z0;
            }

            return samples;
        }

        /// <summary>
        /// Generate random samples from uniform distribution
        /// </summary>
        public static double[] RandomUniform(int count, double min, double max)
        {
            if (count <= 0)
                throw new ArgumentException("Count must be positive");
            if (min >= max)
                throw new ArgumentException("Min must be less than max");

            Random random = new Random();
            double[] samples = new double[count];

            for (int i = 0; i < count; i++)
            {
                samples[i] = min + (max - min) * random.NextDouble();
            }

            return samples;
        }

        /// <summary>
        /// Generate random permutation of integers
        /// </summary>
        public static int[] RandomPermutation(int n)
        {
            if (n <= 0)
                throw new ArgumentException("N must be positive");

            int[] permutation = new int[n];
            for (int i = 0; i < n; i++)
            {
                permutation[i] = i;
            }

            Random random = new Random();
            for (int i = n - 1; i > 0; i--)
            {
                int j = random.Next(i + 1);
                (permutation[i], permutation[j]) = (permutation[j], permutation[i]);
            }

            return permutation;
        }

        /// <summary>
        /// Beta function B(a,b)
        /// </summary>
        public static double Beta(double a, double b)
        {
            if (a <= 0 || b <= 0)
                throw new ArgumentException("Parameters must be positive");

            return Math.Exp(LogGamma(a) + LogGamma(b) - LogGamma(a + b));
        }

        /// <summary>
        /// Incomplete beta function (approximation)
        /// </summary>
        public static double IncompleteBeta(double x, double a, double b)
        {
            if (x < 0 || x > 1)
                throw new ArgumentException("X must be between 0 and 1");
            if (a <= 0 || b <= 0)
                throw new ArgumentException("Parameters must be positive");

            // Continued fraction approximation
            const int maxIterations = 100;
            const double tolerance = 1e-10;

            if (x == 0)
                return 0.0;
            if (x == 1)
                return 1.0;

            double result = Math.Exp(LogGamma(a + b) - LogGamma(a) - LogGamma(b) +
                                   a * Math.Log(x) + b * Math.Log(1 - x));

            // This is a simplified implementation
            // Real implementation would use more sophisticated algorithms
            return result;
        }

        /// <summary>
        /// Error function
        /// </summary>
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

        /// <summary>
        /// Inverse error function (approximation)
        /// </summary>
        private static double InverseErf(double x)
        {
            // Approximation for inverse error function
            double a = 0.886226899;
            double b = -1.645349621;
            double c = 0.914624893;
            double d = -0.140543331;

            double y = Math.Sign(x) * Math.Sqrt(-Math.Log((1.0 - Math.Abs(x)) / 2.0));

            return Math.Sign(x) * Math.Sqrt(-Math.Log((1.0 - Math.Abs(x)) / 2.0)) *
                   (a + b * y + c * y * y + d * y * y * y);
        }

        /// <summary>
        /// Log gamma function (approximation)
        /// </summary>
        private static double LogGamma(double x)
        {
            if (x <= 0)
                throw new ArgumentException("X must be positive");

            // Stirling's approximation for large x
            if (x > 10)
            {
                return (x - 0.5) * Math.Log(x) - x + 0.5 * Math.Log(2 * Math.PI) + 1.0 / (12 * x);
            }

            // For smaller x, use recursive relation
            double result = 0.0;
            double y = x;

            while (y < 10)
            {
                result -= Math.Log(y);
                y += 1;
            }

            return result + (y - 0.5) * Math.Log(y) - y + 0.5 * Math.Log(2 * Math.PI) + 1.0 / (12 * y);
        }

        /// <summary>
        /// Generate random samples from exponential distribution
        /// </summary>
        public static double[] RandomExponential(int count, double rate)
        {
            if (count <= 0)
                throw new ArgumentException("Count must be positive");
            if (rate <= 0)
                throw new ArgumentException("Rate must be positive");

            Random random = new Random();
            double[] samples = new double[count];

            for (int i = 0; i < count; i++)
            {
                samples[i] = -Math.Log(random.NextDouble()) / rate;
            }

            return samples;
        }

        /// <summary>
        /// Generate random samples from gamma distribution
        /// </summary>
        public static double[] RandomGamma(int count, double shape, double scale)
        {
            if (count <= 0)
                throw new ArgumentException("Count must be positive");
            if (shape <= 0 || scale <= 0)
                throw new ArgumentException("Shape and scale must be positive");

            Random random = new Random();
            double[] samples = new double[count];

            if (shape >= 1)
            {
                // Marsaglia and Tsang method
                for (int i = 0; i < count; i++)
                {
                    double d = shape - 1.0 / 3.0;
                    double c = 1.0 / Math.Sqrt(9.0 * d);

                    double x, v;
                    do
                    {
                        do
                        {
                            double u1 = random.NextDouble();
                            double u2 = random.NextDouble();
                            x = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                            v = 1.0 + c * x;
                        } while (v <= 0);

                        v = v * v * v;
                        double u = random.NextDouble();

                        double x2 = x * x;
                        double d_inv = 1.0 / d;

                        if (u < 1.0 - 0.0331 * x2 * x2 ||
                            Math.Log(u) < 0.5 * x2 + d * (1.0 - v + Math.Log(v)))
                        {
                            samples[i] = d * v * scale;
                        }
                    } while (samples[i] == 0);
                }
            }
            else
            {
                // Use exponential distribution for shape < 1
                double[] exponentialSamples = RandomExponential(count, 1.0);
                for (int i = 0; i < count; i++)
                {
                    samples[i] = Math.Pow(random.NextDouble(), 1.0 / shape) * exponentialSamples[i] * scale;
                }
            }

            return samples;
        }

        /// <summary>
        /// Kolmogorov-Smirnov test statistic
        /// </summary>
        public static double KolmogorovSmirnovTest(double[] sample, Func<double, double> cdf)
        {
            if (sample == null || sample.Length == 0)
                throw new ArgumentException("Sample cannot be null or empty");

            double[] sortedSample = new double[sample.Length];
            Array.Copy(sample, sortedSample, sample.Length);
            Array.Sort(sortedSample);

            double maxDiff = 0.0;

            for (int i = 0; i < sortedSample.Length; i++)
            {
                double empiricalCDF = (i + 1.0) / sortedSample.Length;
                double theoreticalCDF = cdf(sortedSample[i]);
                double diff = Math.Abs(empiricalCDF - theoreticalCDF);

                if (diff > maxDiff)
                    maxDiff = diff;
            }

            return maxDiff;
        }

        /// <summary>
        /// Anderson-Darling test statistic
        /// </summary>
        public static double AndersonDarlingTest(double[] sample, Func<double, double> cdf)
        {
            if (sample == null || sample.Length == 0)
                throw new ArgumentException("Sample cannot be null or empty");

            double[] sortedSample = new double[sample.Length];
            Array.Copy(sample, sortedSample, sample.Length);
            Array.Sort(sortedSample);

            double adStatistic = 0.0;
            int n = sortedSample.Length;

            for (int i = 0; i < n; i++)
            {
                double empiricalCDF = (i + 1.0) / (n + 1); // Adjusted for AD test
                double theoreticalCDF = cdf(sortedSample[i]);

                if (theoreticalCDF > 0 && theoreticalCDF < 1)
                {
                    adStatistic += (2 * (i + 1) - 1) * Math.Log(theoreticalCDF) +
                                  (2 * (n - i) - 1) * Math.Log(1 - theoreticalCDF);
                }
            }

            adStatistic = -n - adStatistic / n;

            return adStatistic;
        }
    }
}