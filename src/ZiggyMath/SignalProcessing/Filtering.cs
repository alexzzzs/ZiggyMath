using System;
using ZiggyMath.Core;

namespace ZiggyMath.SignalProcessing
{
    /// <summary>
    /// Digital filter design and application
    /// </summary>
    public static class Filtering
    {
        /// <summary>
        /// Apply FIR filter to input signal
        /// </summary>
        public static double[] ApplyFIRFilter(double[] input, double[] coefficients)
        {
            if (input == null || input.Length == 0)
                throw new ArgumentException("Input signal cannot be null or empty");
            if (coefficients == null || coefficients.Length == 0)
                throw new ArgumentException("Filter coefficients cannot be null or empty");

            return Convolution.Convolve(input, coefficients);
        }

        /// <summary>
        /// Apply IIR filter to input signal
        /// </summary>
        public static double[] ApplyIIRFilter(double[] input, double[] aCoefficients, double[] bCoefficients)
        {
            if (input == null || input.Length == 0)
                throw new ArgumentException("Input signal cannot be null or empty");
            if (aCoefficients == null || aCoefficients.Length == 0)
                throw new ArgumentException("A coefficients cannot be null or empty");
            if (bCoefficients == null || bCoefficients.Length == 0)
                throw new ArgumentException("B coefficients cannot be null or empty");

            int order = Math.Max(aCoefficients.Length, bCoefficients.Length) - 1;
            double[] output = new double[input.Length + order];

            // Initialize delay lines
            double[] xDelay = new double[order];
            double[] yDelay = new double[order];

            for (int n = 0; n < input.Length; n++)
            {
                // Shift input delay line
                for (int i = order - 1; i > 0; i--)
                {
                    xDelay[i] = xDelay[i - 1];
                }
                xDelay[0] = input[n];

                // Compute output using IIR difference equation
                double y = 0.0;

                // Feed-forward part (zeros)
                for (int i = 0; i < Math.Min(bCoefficients.Length, order + 1); i++)
                {
                    if (n - i >= 0)
                    {
                        y += bCoefficients[i] * (i == 0 ? input[n] : xDelay[i - 1]);
                    }
                }

                // Feedback part (poles)
                for (int i = 1; i < Math.Min(aCoefficients.Length, order + 1); i++)
                {
                    if (n - i >= 0)
                    {
                        y -= aCoefficients[i] * (i == 0 ? output[n] : yDelay[i - 1]);
                    }
                }

                // Normalize by a[0]
                if (aCoefficients.Length > 0 && aCoefficients[0] != 0)
                {
                    y /= aCoefficients[0];
                }

                output[n] = y;

                // Shift output delay line
                for (int i = order - 1; i > 0; i--)
                {
                    yDelay[i] = yDelay[i - 1];
                }
                if (order > 0)
                {
                    yDelay[0] = y;
                }
            }

            // Return only the valid output samples
            double[] result = new double[input.Length];
            Array.Copy(output, result, input.Length);

            return result;
        }

        /// <summary>
        /// Design low-pass FIR filter
        /// </summary>
        public static double[] DesignLowPassFIR(int order, double cutoffFrequency)
        {
            if (order <= 0)
                throw new ArgumentException("Filter order must be positive");
            if (cutoffFrequency <= 0 || cutoffFrequency >= 0.5)
                throw new ArgumentException("Cutoff frequency must be between 0 and 0.5");

            double[] coefficients = new double[order + 1];
            double normalizedCutoff = 2.0 * cutoffFrequency;

            // Simple rectangular window design
            for (int n = 0; n <= order; n++)
            {
                if (n == order / 2)
                {
                    coefficients[n] = normalizedCutoff;
                }
                else
                {
                    double sincArg = normalizedCutoff * (n - order / 2.0);
                    coefficients[n] = Math.Sin(sincArg) / sincArg;
                }
            }

            return coefficients;
        }

        /// <summary>
        /// Design high-pass FIR filter
        /// </summary>
        public static double[] DesignHighPassFIR(int order, double cutoffFrequency)
        {
            double[] lowPass = DesignLowPassFIR(order, cutoffFrequency);
            double[] highPass = new double[lowPass.Length];

            // High-pass = Low-pass inverted with center tap adjustment
            for (int i = 0; i < lowPass.Length; i++)
            {
                highPass[i] = -lowPass[i];
            }

            highPass[order / 2] += 1.0;

            return highPass;
        }

        /// <summary>
        /// Design band-pass FIR filter
        /// </summary>
        public static double[] DesignBandPassFIR(int order, double lowCutoff, double highCutoff)
        {
            if (lowCutoff >= highCutoff)
                throw new ArgumentException("Low cutoff must be less than high cutoff");

            double[] lowPass = DesignLowPassFIR(order, highCutoff);
            double[] highPass = DesignHighPassFIR(order, lowCutoff);
            double[] bandPass = new double[lowPass.Length];

            // Band-pass = High-pass + Low-pass
            for (int i = 0; i < lowPass.Length; i++)
            {
                bandPass[i] = highPass[i] + lowPass[i];
            }

            return bandPass;
        }

        /// <summary>
        /// Hamming window function
        /// </summary>
        public static double[] HammingWindow(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Window length must be positive");

            double[] window = new double[length];

            for (int n = 0; n < length; n++)
            {
                window[n] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * n / (length - 1));
            }

            return window;
        }

        /// <summary>
        /// Hann window function
        /// </summary>
        public static double[] HannWindow(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Window length must be positive");

            double[] window = new double[length];

            for (int n = 0; n < length; n++)
            {
                window[n] = 0.5 * (1 - Math.Cos(2 * Math.PI * n / (length - 1)));
            }

            return window;
        }

        /// <summary>
        /// Blackman window function
        /// </summary>
        public static double[] BlackmanWindow(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Window length must be positive");

            double[] window = new double[length];

            for (int n = 0; n < length; n++)
            {
                double a0 = 0.42;
                double a1 = 0.5;
                double a2 = 0.08;

                window[n] = a0 - a1 * Math.Cos(2 * Math.PI * n / (length - 1)) +
                           a2 * Math.Cos(4 * Math.PI * n / (length - 1));
            }

            return window;
        }

        /// <summary>
        /// Kaiser window function
        /// </summary>
        public static double[] KaiserWindow(int length, double beta)
        {
            if (length <= 0)
                throw new ArgumentException("Window length must be positive");

            double[] window = new double[length];
            double denominator = BesselI0(beta);

            for (int n = 0; n < length; n++)
            {
                double x = 2.0 * n / (length - 1) - 1.0;
                window[n] = BesselI0(beta * Math.Sqrt(1 - x * x)) / denominator;
            }

            return window;
        }

        /// <summary>
        /// Apply window function to filter coefficients
        /// </summary>
        public static double[] ApplyWindow(double[] coefficients, double[] window)
        {
            if (coefficients == null || coefficients.Length == 0)
                throw new ArgumentException("Coefficients cannot be null or empty");
            if (window == null || window.Length == 0)
                throw new ArgumentException("Window cannot be null or empty");
            if (coefficients.Length != window.Length)
                throw new ArgumentException("Coefficients and window must have the same length");

            double[] result = new double[coefficients.Length];
            for (int i = 0; i < coefficients.Length; i++)
            {
                result[i] = coefficients[i] * window[i];
            }

            return result;
        }

        /// <summary>
        /// Design moving average filter
        /// </summary>
        public static double[] MovingAverageFilter(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Filter length must be positive");

            double[] coefficients = new double[length];
            double value = 1.0 / length;

            for (int i = 0; i < length; i++)
            {
                coefficients[i] = value;
            }

            return coefficients;
        }

        /// <summary>
        /// Design Gaussian filter
        /// </summary>
        public static double[] GaussianFilter(int length, double sigma)
        {
            if (length <= 0)
                throw new ArgumentException("Filter length must be positive");
            if (sigma <= 0)
                throw new ArgumentException("Sigma must be positive");

            double[] coefficients = new double[length];
            int center = length / 2;
            double sum = 0.0;

            for (int i = 0; i < length; i++)
            {
                int x = i - center;
                coefficients[i] = Math.Exp(-(x * x) / (2 * sigma * sigma));
                sum += coefficients[i];
            }

            // Normalize
            for (int i = 0; i < length; i++)
            {
                coefficients[i] /= sum;
            }

            return coefficients;
        }

        /// <summary>
        /// Modified Bessel function of the first kind (order 0)
        /// </summary>
        private static double BesselI0(double x)
        {
            double sum = 1.0;
            double term = 1.0;
            double k = 1.0;

            while (Math.Abs(term) > 1e-12)
            {
                term *= (x * x) / (4 * k * k);
                sum += term;
                k += 1.0;
            }

            return sum;
        }

        /// <summary>
        /// Compute filter frequency response
        /// </summary>
        public static Complex[] FrequencyResponse(double[] coefficients, int numPoints = 512)
        {
            Complex[] response = new Complex[numPoints];

            for (int k = 0; k < numPoints; k++)
            {
                double frequency = 2.0 * Math.PI * k / numPoints;
                Complex sum = new Complex(0, 0);

                for (int n = 0; n < coefficients.Length; n++)
                {
                    double angle = -frequency * n;
                    Complex coeff = new Complex(coefficients[n] * Math.Cos(angle), coefficients[n] * Math.Sin(angle));
                    sum += coeff;
                }

                response[k] = sum;
            }

            return response;
        }

        /// <summary>
        /// Design Butterworth low-pass filter (simplified)
        /// </summary>
        public static (double[] a, double[] b) ButterworthLowPass(int order, double cutoffFrequency)
        {
            if (order <= 0)
                throw new ArgumentException("Filter order must be positive");
            if (cutoffFrequency <= 0 || cutoffFrequency >= 1)
                throw new ArgumentException("Cutoff frequency must be between 0 and 1");

            // Simplified Butterworth design
            double[] a = new double[order + 1];
            double[] b = new double[order + 1];

            // This is a simplified implementation
            // Real Butterworth design would use more sophisticated algorithms
            a[0] = 1.0;
            b[0] = 1.0;

            for (int i = 1; i <= order; i++)
            {
                a[i] = -Math.Pow(cutoffFrequency, i);
                b[i] = Math.Pow(cutoffFrequency, i);
            }

            return (a, b);
        }
    }
}