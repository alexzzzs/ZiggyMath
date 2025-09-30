using System;
using System.Linq;
using ZiggyMath.Core;

namespace ZiggyMath.Statistics
{
    /// <summary>
    /// Descriptive statistics
    /// </summary>
    public static class DescriptiveStats
    {
        /// <summary>
        /// Calculate arithmetic mean
        /// </summary>
        public static double Mean(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double sum = 0.0;
            for (int i = 0; i < data.Length; i++)
            {
                sum += data[i];
            }

            return sum / data.Length;
        }

        /// <summary>
        /// Calculate median
        /// </summary>
        public static double Median(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double[] sortedData = data.ToArray();
            Array.Sort(sortedData);

            int n = sortedData.Length;
            if (n % 2 == 0)
            {
                return (sortedData[n / 2 - 1] + sortedData[n / 2]) / 2.0;
            }
            else
            {
                return sortedData[n / 2];
            }
        }

        /// <summary>
        /// Calculate mode (most frequent value)
        /// </summary>
        public static double Mode(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            // Simple implementation - find most frequent value
            double[] sortedData = data.ToArray();
            Array.Sort(sortedData);

            double currentValue = sortedData[0];
            double mode = currentValue;
            int currentCount = 1;
            int maxCount = 1;

            for (int i = 1; i < sortedData.Length; i++)
            {
                if (Math.Abs(sortedData[i] - currentValue) < 1e-10)
                {
                    currentCount++;
                }
                else
                {
                    if (currentCount > maxCount)
                    {
                        maxCount = currentCount;
                        mode = currentValue;
                    }
                    currentValue = sortedData[i];
                    currentCount = 1;
                }
            }

            if (currentCount > maxCount)
            {
                mode = currentValue;
            }

            return mode;
        }

        /// <summary>
        /// Calculate variance
        /// </summary>
        public static double Variance(ReadOnlySpan<double> data)
        {
            if (data.Length <= 1)
                return 0.0;

            double mean = Mean(data);
            double sumSquaredDiffs = 0.0;

            for (int i = 0; i < data.Length; i++)
            {
                double diff = data[i] - mean;
                sumSquaredDiffs += diff * diff;
            }

            return sumSquaredDiffs / (data.Length - 1); // Sample variance
        }

        /// <summary>
        /// Calculate standard deviation
        /// </summary>
        public static double StandardDeviation(ReadOnlySpan<double> data)
        {
            return Math.Sqrt(Variance(data));
        }

        /// <summary>
        /// Calculate skewness
        /// </summary>
        public static double Skewness(ReadOnlySpan<double> data)
        {
            if (data.Length <= 2)
                return 0.0;

            double mean = Mean(data);
            double stdDev = StandardDeviation(data);

            if (stdDev == 0)
                return 0.0;

            double sumCubedDiffs = 0.0;

            for (int i = 0; i < data.Length; i++)
            {
                double diff = (data[i] - mean) / stdDev;
                sumCubedDiffs += diff * diff * diff;
            }

            return sumCubedDiffs / data.Length;
        }

        /// <summary>
        /// Calculate kurtosis
        /// </summary>
        public static double Kurtosis(ReadOnlySpan<double> data)
        {
            if (data.Length <= 3)
                return 0.0;

            double mean = Mean(data);
            double stdDev = StandardDeviation(data);

            if (stdDev == 0)
                return 0.0;

            double sumFourthPowerDiffs = 0.0;

            for (int i = 0; i < data.Length; i++)
            {
                double diff = (data[i] - mean) / stdDev;
                sumFourthPowerDiffs += diff * diff * diff * diff;
            }

            return (sumFourthPowerDiffs / data.Length) - 3.0; // Excess kurtosis
        }

        /// <summary>
        /// Find minimum value
        /// </summary>
        public static double Min(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double min = data[0];
            for (int i = 1; i < data.Length; i++)
            {
                if (data[i] < min)
                    min = data[i];
            }

            return min;
        }

        /// <summary>
        /// Find maximum value
        /// </summary>
        public static double Max(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double max = data[0];
            for (int i = 1; i < data.Length; i++)
            {
                if (data[i] > max)
                    max = data[i];
            }

            return max;
        }

        /// <summary>
        /// Calculate range (max - min)
        /// </summary>
        public static double Range(ReadOnlySpan<double> data)
        {
            return Max(data) - Min(data);
        }

        /// <summary>
        /// Calculate quartiles (Q1 and Q3)
        /// </summary>
        public static (double q1, double q3) Quartiles(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double[] sortedData = data.ToArray();
            Array.Sort(sortedData);

            int n = sortedData.Length;
            int q1Index = n / 4;
            int q3Index = 3 * n / 4;

            double q1 = sortedData[q1Index];
            double q3 = sortedData[q3Index];

            // Linear interpolation for more accurate quartiles
            if (n % 4 == 0 && q1Index > 0)
            {
                q1 = (sortedData[q1Index - 1] + sortedData[q1Index]) / 2.0;
            }

            if (3 * n % 4 == 0 && q3Index < n - 1)
            {
                q3 = (sortedData[q3Index] + sortedData[q3Index + 1]) / 2.0;
            }

            return (q1, q3);
        }

        /// <summary>
        /// Calculate interquartile range
        /// </summary>
        public static double InterquartileRange(ReadOnlySpan<double> data)
        {
            var (q1, q3) = Quartiles(data);
            return q3 - q1;
        }

        /// <summary>
        /// Calculate median absolute deviation
        /// </summary>
        public static double MedianAbsoluteDeviation(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double median = Median(data);
            double[] deviations = new double[data.Length];

            for (int i = 0; i < data.Length; i++)
            {
                deviations[i] = Math.Abs(data[i] - median);
            }

            return Median(deviations);
        }

        /// <summary>
        /// Calculate percentiles
        /// </summary>
        public static double Percentile(ReadOnlySpan<double> data, double percentile)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");
            if (percentile < 0 || percentile > 100)
                throw new ArgumentException("Percentile must be between 0 and 100");

            double[] sortedData = data.ToArray();
            Array.Sort(sortedData);

            double index = (percentile / 100.0) * (sortedData.Length - 1);
            int lowerIndex = (int)Math.Floor(index);
            int upperIndex = (int)Math.Ceiling(index);

            if (lowerIndex == upperIndex)
            {
                return sortedData[lowerIndex];
            }
            else
            {
                double weight = index - lowerIndex;
                return sortedData[lowerIndex] * (1 - weight) + sortedData[upperIndex] * weight;
            }
        }

        /// <summary>
        /// Calculate trimmed mean (remove outliers)
        /// </summary>
        public static double TrimmedMean(ReadOnlySpan<double> data, double trimPercentage)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");
            if (trimPercentage < 0 || trimPercentage > 50)
                throw new ArgumentException("Trim percentage must be between 0 and 50");

            double[] sortedData = data.ToArray();
            Array.Sort(sortedData);

            int trimCount = (int)(sortedData.Length * trimPercentage / 100.0);
            int startIndex = trimCount;
            int endIndex = sortedData.Length - trimCount;

            double sum = 0.0;
            int count = 0;

            for (int i = startIndex; i < endIndex; i++)
            {
                sum += sortedData[i];
                count++;
            }

            return count > 0 ? sum / count : 0.0;
        }

        /// <summary>
        /// Calculate geometric mean
        /// </summary>
        public static double GeometricMean(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double product = 1.0;

            for (int i = 0; i < data.Length; i++)
            {
                if (data[i] <= 0)
                    throw new ArgumentException("All values must be positive for geometric mean");
                product *= data[i];
            }

            return Math.Pow(product, 1.0 / data.Length);
        }

        /// <summary>
        /// Calculate harmonic mean
        /// </summary>
        public static double HarmonicMean(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double sumReciprocal = 0.0;

            for (int i = 0; i < data.Length; i++)
            {
                if (data[i] == 0)
                    throw new ArgumentException("Cannot compute harmonic mean with zero values");
                sumReciprocal += 1.0 / data[i];
            }

            return data.Length / sumReciprocal;
        }

        /// <summary>
        /// Calculate root mean square
        /// </summary>
        public static double RootMeanSquare(ReadOnlySpan<double> data)
        {
            if (data.Length == 0)
                throw new ArgumentException("Data cannot be empty");

            double sumSquares = 0.0;

            for (int i = 0; i < data.Length; i++)
            {
                sumSquares += data[i] * data[i];
            }

            return Math.Sqrt(sumSquares / data.Length);
        }

        /// <summary>
        /// Calculate coefficient of variation
        /// </summary>
        public static double CoefficientOfVariation(ReadOnlySpan<double> data)
        {
            double mean = Mean(data);
            double stdDev = StandardDeviation(data);

            if (mean == 0)
                return 0.0;

            return stdDev / Math.Abs(mean);
        }
    }
}