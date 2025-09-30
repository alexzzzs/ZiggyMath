using System;
using ZiggyMath.Core;

namespace ZiggyMath.SignalProcessing
{
    /// <summary>
    /// Fast Fourier Transform implementation - because sometimes you need to see your signal's frequency domain.
    /// This is the Cooley-Tukey radix-2 algorithm, which runs in O(n log n) time instead of O(n²).
    /// </summary>
    public static class FFT
    {
        /// <summary>
        /// Forward FFT of complex vector - transforms time domain to frequency domain.
        /// Because who doesn't want to know what their signal looks like in frequency space?
        /// </summary>
        /// <param name="input">Complex input signal (must be power of 2 length)</param>
        /// <returns>Frequency domain representation of the input signal</returns>
        /// <exception cref="ArgumentException">When input is invalid or not power of 2</exception>
        public static ComplexVector Forward(ComplexVector input)
        {
            // Input validation - because even FFTs need good inputs!
            if (input == null || input.Length == 0)
                throw new ArgumentException("Input cannot be null or empty. " +
                    "Even the frequency domain needs some data to work with!");

            int n = input.Length;

            // FFT requires power-of-2 lengths for the divide-and-conquer approach
            if (!IsPowerOfTwo(n))
                throw new ArgumentException($"Input length must be a power of 2, but got {n}. " +
                    "This isn't because we're lazy - it's because the algorithm is most efficient this way!");

            // Create result vector - we'll be modifying this in-place
            ComplexVector result = new ComplexVector(n);

            // Copy input data to result (we'll modify result in-place)
            for (int i = 0; i < n; i++)
            {
                result[i] = input[i];
            }

            // Perform the actual FFT using our recursive implementation
            // This is where the magic happens - O(n log n) complexity!
            FFT_Recursive(result, 0, n, 1);

            return result;
        }

        /// <summary>
        /// Inverse FFT of complex vector
        /// </summary>
        public static ComplexVector Inverse(ComplexVector input)
        {
            if (input == null || input.Length == 0)
                throw new ArgumentException("Input cannot be null or empty");

            int n = input.Length;
            if (!IsPowerOfTwo(n))
                throw new ArgumentException("Input length must be a power of 2");

            ComplexVector result = new ComplexVector(n);

            // Copy input to result
            for (int i = 0; i < n; i++)
            {
                result[i] = input[i];
            }

            // Perform inverse FFT
            FFT_Recursive(result, 0, n, -1);

            // Scale by 1/n
            double scale = 1.0 / n;
            for (int i = 0; i < n; i++)
            {
                result[i] = new Complex(result[i].Real * scale, result[i].Imaginary * scale);
            }

            return result;
        }

        /// <summary>
        /// Forward FFT of real-valued signal
        /// </summary>
        public static ComplexVector Forward(double[] realInput)
        {
            if (realInput == null || realInput.Length == 0)
                throw new ArgumentException("Input cannot be null or empty");

            int n = realInput.Length;
            if (!IsPowerOfTwo(n))
                throw new ArgumentException("Input length must be a power of 2");

            ComplexVector complexInput = new ComplexVector(n);
            for (int i = 0; i < n; i++)
            {
                complexInput[i] = new Complex(realInput[i], 0);
            }

            return Forward(complexInput);
        }

        /// <summary>
        /// Inverse FFT returning real and imaginary parts
        /// </summary>
        public static (double[] real, double[] imaginary) InverseToReal(ComplexVector input)
        {
            ComplexVector result = Inverse(input);

            int n = result.Length;
            double[] real = new double[n];
            double[] imaginary = new double[n];

            for (int i = 0; i < n; i++)
            {
                real[i] = result[i].Real;
                imaginary[i] = result[i].Imaginary;
            }

            return (real, imaginary);
        }

        /// <summary>
        /// 2D FFT of complex matrix
        /// </summary>
        public static ComplexMatrix Forward2D(ComplexMatrix input)
        {
            throw new NotImplementedException("2D FFT not yet implemented");
        }

        /// <summary>
        /// 2D inverse FFT of complex matrix
        /// </summary>
        public static ComplexMatrix Inverse2D(ComplexMatrix input)
        {
            throw new NotImplementedException("2D inverse FFT not yet implemented");
        }

        /// <summary>
        /// Get next power of two greater than or equal to n
        /// </summary>
        public static int NextPowerOfTwo(int n)
        {
            if (n <= 0)
                return 1;

            int power = 1;
            while (power < n)
            {
                power *= 2;
            }

            return power;
        }

        /// <summary>
        /// Check if n is a power of two
        /// </summary>
        public static bool IsPowerOfTwo(int n)
        {
            return n > 0 && (n & (n - 1)) == 0;
        }

        /// <summary>
        /// Recursive FFT implementation using the Cooley-Tukey algorithm.
        /// This is the "divide and conquer" approach: break the problem into smaller pieces,
        /// solve recursively, then combine the results. Classic computer science magic!
        /// </summary>
        /// <param name="x">The complex vector to transform (modified in-place)</param>
        /// <param name="start">Starting index in the vector</param>
        /// <param name="length">Length of the sub-array to transform</param>
        /// <param name="sign">1 for forward FFT, -1 for inverse FFT</param>
        private static void FFT_Recursive(ComplexVector x, int start, int length, int sign)
        {
            // Base case: if length is 1 or less, nothing to do
            // This is the "stop dividing when you can't divide anymore" part
            if (length <= 1)
                return;

            // Split the array into two halves of equal size
            int halfLength = length / 2;

            // Bit reversal permutation - this is crucial for the algorithm to work correctly
            // It's like shuffling a deck of cards in a very specific way
            for (int i = 0; i < halfLength; i++)
            {
                int evenIndex = start + i;           // Index in the "even" half
                int oddIndex = start + i + halfLength; // Index in the "odd" half

                // Swap the elements - because the algorithm needs them in a specific order
                Complex temp = x[evenIndex];
                x[evenIndex] = x[oddIndex];
                x[oddIndex] = temp;
            }

            // Recursive calls: transform the two halves separately
            // This is the "divide" part of divide-and-conquer
            FFT_Recursive(x, start, halfLength, sign);           // Transform first half
            FFT_Recursive(x, start + halfLength, halfLength, sign); // Transform second half

            // Combine the results from the two halves
            // This is the "conquer" part - where we use the butterfly operations
            for (int i = 0; i < halfLength; i++)
            {
                int evenIndex = start + i;           // Index in the first half
                int oddIndex = start + i + halfLength; // Index in the second half

                // Get the values from both halves
                Complex even = x[evenIndex];
                Complex odd = x[oddIndex];

                // Calculate the twiddle factor (the "rotation" factor in the unit circle)
                // This is what makes FFT different from regular DFT
                double angle = sign * 2 * Math.PI * i / length;
                Complex w = new Complex(Math.Cos(angle), Math.Sin(angle));

                // The butterfly operation: combine the two halves with the twiddle factor
                // This is where the frequency domain transformation actually happens!
                x[evenIndex] = even + w * odd;  // Top part of the butterfly
                x[oddIndex] = even - w * odd;   // Bottom part of the butterfly
            }
        }

        /// <summary>
        /// Compute power spectrum of signal
        /// </summary>
        public static double[] PowerSpectrum(double[] signal)
        {
            ComplexVector fft = Forward(signal);
            int n = fft.Length;
            double[] powerSpectrum = new double[n / 2 + 1];

            powerSpectrum[0] = fft[0].Magnitude * fft[0].Magnitude;

            for (int i = 1; i < n / 2; i++)
            {
                double magnitude = fft[i].Magnitude;
                powerSpectrum[i] = magnitude * magnitude;
            }

            if (n > 1)
            {
                powerSpectrum[n / 2] = fft[n / 2].Magnitude * fft[n / 2].Magnitude;
            }

            return powerSpectrum;
        }

        /// <summary>
        /// Compute magnitude spectrum of signal
        /// </summary>
        public static double[] MagnitudeSpectrum(double[] signal)
        {
            ComplexVector fft = Forward(signal);
            int n = fft.Length;
            double[] magnitudeSpectrum = new double[n / 2 + 1];

            magnitudeSpectrum[0] = fft[0].Magnitude;

            for (int i = 1; i < n / 2; i++)
            {
                magnitudeSpectrum[i] = fft[i].Magnitude;
            }

            if (n > 1)
            {
                magnitudeSpectrum[n / 2] = fft[n / 2].Magnitude;
            }

            return magnitudeSpectrum;
        }

        /// <summary>
        /// Compute phase spectrum of signal
        /// </summary>
        public static double[] PhaseSpectrum(double[] signal)
        {
            ComplexVector fft = Forward(signal);
            int n = fft.Length;
            double[] phaseSpectrum = new double[n / 2 + 1];

            phaseSpectrum[0] = fft[0].Phase;

            for (int i = 1; i < n / 2; i++)
            {
                phaseSpectrum[i] = fft[i].Phase;
            }

            if (n > 1)
            {
                phaseSpectrum[n / 2] = fft[n / 2].Phase;
            }

            return phaseSpectrum;
        }
    }

    /// <summary>
    /// Complex matrix for 2D operations
    /// </summary>
    public class ComplexMatrix : IDisposable
    {
        private readonly Complex[,] _data;

        public int Rows { get; }
        public int Columns { get; }
        public bool IsValid => _data != null;

        public ComplexMatrix(int rows, int columns)
        {
            Rows = rows;
            Columns = columns;
            _data = new Complex[rows, columns];
        }

        public Complex this[int row, int column]
        {
            get => _data[row, column];
            set => _data[row, column] = value;
        }

        public void Dispose()
        {
            // No unmanaged resources to dispose
        }
    }
}