using System;
using ZiggyMath.Core;

namespace ZiggyMath.SignalProcessing
{
    /// <summary>
    /// Convolution and filtering operations
    /// </summary>
    public static class Convolution
    {
        /// <summary>
        /// 1D convolution of real-valued signals
        /// </summary>
        public static double[] Convolve(double[] signal, double[] kernel)
        {
            if (signal == null || signal.Length == 0)
                throw new ArgumentException("Signal cannot be null or empty");
            if (kernel == null || kernel.Length == 0)
                throw new ArgumentException("Kernel cannot be null or empty");

            int signalLength = signal.Length;
            int kernelLength = kernel.Length;
            int resultLength = signalLength + kernelLength - 1;

            double[] result = new double[resultLength];

            for (int n = 0; n < resultLength; n++)
            {
                double sum = 0.0;
                int start = Math.Max(0, n - kernelLength + 1);
                int end = Math.Min(n, signalLength - 1);

                for (int k = start; k <= end; k++)
                {
                    sum += signal[k] * kernel[n - k];
                }

                result[n] = sum;
            }

            return result;
        }

        /// <summary>
        /// 1D convolution of complex-valued signals
        /// </summary>
        public static ComplexVector Convolve(ComplexVector signal, ComplexVector kernel)
        {
            if (signal == null || signal.Length == 0)
                throw new ArgumentException("Signal cannot be null or empty");
            if (kernel == null || kernel.Length == 0)
                throw new ArgumentException("Kernel cannot be null or empty");

            int signalLength = signal.Length;
            int kernelLength = kernel.Length;
            int resultLength = signalLength + kernelLength - 1;

            ComplexVector result = new ComplexVector(resultLength);

            for (int n = 0; n < resultLength; n++)
            {
                Complex sum = new Complex(0, 0);
                int start = Math.Max(0, n - kernelLength + 1);
                int end = Math.Min(n, signalLength - 1);

                for (int k = start; k <= end; k++)
                {
                    sum += signal[k] * kernel[n - k];
                }

                result[n] = sum;
            }

            return result;
        }

        /// <summary>
        /// 2D convolution of real-valued images/kernels
        /// </summary>
        public static double[,] Convolve2D(double[,] image, double[,] kernel)
        {
            if (image == null || image.GetLength(0) == 0 || image.GetLength(1) == 0)
                throw new ArgumentException("Image cannot be null or empty");
            if (kernel == null || kernel.GetLength(0) == 0 || kernel.GetLength(1) == 0)
                throw new ArgumentException("Kernel cannot be null or empty");

            int imageRows = image.GetLength(0);
            int imageCols = image.GetLength(1);
            int kernelRows = kernel.GetLength(0);
            int kernelCols = kernel.GetLength(1);

            // Check if kernel dimensions are odd
            if (kernelRows % 2 == 0 || kernelCols % 2 == 0)
                throw new ArgumentException("Kernel dimensions must be odd");

            int resultRows = imageRows - kernelRows + 1;
            int resultCols = imageCols - kernelCols + 1;

            if (resultRows <= 0 || resultCols <= 0)
                throw new ArgumentException("Kernel is larger than image");

            double[,] result = new double[resultRows, resultCols];

            int kernelCenterRow = kernelRows / 2;
            int kernelCenterCol = kernelCols / 2;

            for (int i = 0; i < resultRows; i++)
            {
                for (int j = 0; j < resultCols; j++)
                {
                    double sum = 0.0;

                    for (int ki = 0; ki < kernelRows; ki++)
                    {
                        for (int kj = 0; kj < kernelCols; kj++)
                        {
                            int imageRow = i + ki - kernelCenterRow;
                            int imageCol = j + kj - kernelCenterCol;

                            if (imageRow >= 0 && imageRow < imageRows &&
                                imageCol >= 0 && imageCol < imageCols)
                            {
                                sum += image[imageRow, imageCol] * kernel[ki, kj];
                            }
                        }
                    }

                    result[i, j] = sum;
                }
            }

            return result;
        }

        /// <summary>
        /// 2D convolution of complex-valued matrices
        /// </summary>
        public static ComplexMatrix Convolve2D(ComplexMatrix image, ComplexMatrix kernel)
        {
            throw new NotImplementedException("Complex 2D convolution not yet implemented");
        }

        /// <summary>
        /// Separable 2D convolution using row and column kernels
        /// </summary>
        public static double[,] SeparableConvolve2D(double[,] image, double[] rowKernel, double[] colKernel)
        {
            if (image == null || image.GetLength(0) == 0 || image.GetLength(1) == 0)
                throw new ArgumentException("Image cannot be null or empty");
            if (rowKernel == null || rowKernel.Length == 0)
                throw new ArgumentException("Row kernel cannot be null or empty");
            if (colKernel == null || colKernel.Length == 0)
                throw new ArgumentException("Column kernel cannot be null or empty");

            int rows = image.GetLength(0);
            int cols = image.GetLength(1);

            // First convolve with row kernel
            double[,] temp = new double[rows, cols];
            for (int i = 0; i < rows; i++)
            {
                double[] row = new double[cols];
                for (int j = 0; j < cols; j++)
                {
                    row[j] = image[i, j];
                }

                double[] convolvedRow = Convolve(row, rowKernel);
                Array.Copy(convolvedRow, 0, temp, i * cols, Math.Min(cols, convolvedRow.Length));
            }

            // Then convolve with column kernel
            double[,] result = new double[rows, cols];
            for (int j = 0; j < cols; j++)
            {
                double[] col = new double[rows];
                for (int i = 0; i < rows; i++)
                {
                    col[i] = temp[i, j];
                }

                double[] convolvedCol = Convolve(col, colKernel);

                for (int i = 0; i < Math.Min(rows, convolvedCol.Length); i++)
                {
                    result[i, j] = convolvedCol[i];
                }
            }

            return result;
        }

        /// <summary>
        /// Convolution using FFT for efficiency
        /// </summary>
        public static double[] ConvolveFFT(double[] signal, double[] kernel)
        {
            int signalLength = signal.Length;
            int kernelLength = kernel.Length;
            int fftSize = NextPowerOfTwo(signalLength + kernelLength - 1);

            // Pad signals to FFT size
            double[] paddedSignal = new double[fftSize];
            double[] paddedKernel = new double[fftSize];

            Array.Copy(signal, paddedSignal, signalLength);
            Array.Copy(kernel, paddedKernel, kernelLength);

            // Compute FFTs
            ComplexVector signalFFT = FFT.Forward(paddedSignal);
            ComplexVector kernelFFT = FFT.Forward(paddedKernel);

            // Multiply in frequency domain
            ComplexVector productFFT = new ComplexVector(fftSize);
            for (int i = 0; i < fftSize; i++)
            {
                productFFT[i] = signalFFT[i] * kernelFFT[i];
            }

            // Inverse FFT
            var (real, _) = FFT.InverseToReal(productFFT);

            // Extract valid convolution result
            double[] result = new double[signalLength + kernelLength - 1];
            Array.Copy(real, result, result.Length);

            return result;
        }

        /// <summary>
        /// Get next power of two
        /// </summary>
        private static int NextPowerOfTwo(int n)
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
        /// Compute autocorrelation of signal
        /// </summary>
        public static double[] Autocorrelation(double[] signal)
        {
            int n = signal.Length;
            double[] result = new double[n];

            for (int lag = 0; lag < n; lag++)
            {
                double sum = 0.0;
                int count = 0;

                for (int i = 0; i < n - lag; i++)
                {
                    sum += signal[i] * signal[i + lag];
                    count++;
                }

                result[lag] = count > 0 ? sum / count : 0.0;
            }

            return result;
        }

        /// <summary>
        /// Compute cross-correlation between two signals
        /// </summary>
        public static double[] CrossCorrelation(double[] signal1, double[] signal2)
        {
            if (signal1.Length != signal2.Length)
                throw new ArgumentException("Signals must have the same length");

            int n = signal1.Length;
            double[] result = new double[2 * n - 1];

            for (int lag = -n + 1; lag < n; lag++)
            {
                double sum = 0.0;
                int count = 0;

                for (int i = 0; i < n; i++)
                {
                    int j = i - lag;
                    if (j >= 0 && j < n)
                    {
                        sum += signal1[i] * signal2[j];
                        count++;
                    }
                }

                int index = lag + n - 1;
                result[index] = count > 0 ? sum / count : 0.0;
            }

            return result;
        }
    }
}