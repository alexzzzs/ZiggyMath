using System;
using ZiggyMath;
using ZiggyMath.Core;
using ZiggyMath.SimdMath;
using ZiggyMath.LinearAlgebra;
using ZiggyMath.SignalProcessing;
using ZiggyMath.Statistics;
using ZiggyMath.Calculus;

namespace Examples
{
    /// <summary>
    /// Comprehensive ZiggyMath examples organized by functionality
    /// Demonstrates core features, performance optimizations, and advanced usage
    /// </summary>
    public static class ZiggyMathExamples
    {
        /// <summary>
        /// Run all example categories
        /// </summary>
        public static void RunAllExamples()
        {
            Console.WriteLine("🚀 ZiggyMath Comprehensive Examples");
            Console.WriteLine("===================================");

            CoreExamples();
            PerformanceExamples();
            AdvancedExamples();
        }

        /// <summary>
        /// Core mathematical operations and basic usage
        /// </summary>
        public static void CoreExamples()
        {
            Console.WriteLine("\n📐 Core Mathematical Operations");
            Console.WriteLine("-------------------------------");

            VectorOperations();
            MatrixOperations();
            ComplexOperations();
        }

        /// <summary>
        /// Performance optimizations and benchmarks
        /// </summary>
        public static void PerformanceExamples()
        {
            Console.WriteLine("\n⚡ Performance & Optimization");
            Console.WriteLine("----------------------------");

            MemoryOptimization();
            ComputationContexts();
            Benchmarking();
        }

        /// <summary>
        /// Advanced algorithms and specialized operations
        /// </summary>
        public static void AdvancedExamples()
        {
            Console.WriteLine("\n🔬 Advanced Mathematical Features");
            Console.WriteLine("--------------------------------");

            LinearAlgebra();
            SignalProcessing();
            NumericalMethods();
            StatisticalAnalysis();
        }

        #region Core Examples

        private static void VectorOperations()
        {
            Console.WriteLine("\n  Vector Operations:");

            using var context = new MathComputationContext(MathWorkload.SmallComputations);
            using var v1 = context.RentBuffer(5);
            using var v2 = context.RentBuffer(5);
            using var result = context.RentBuffer(5);

            // Initialize vectors
            for (int i = 0; i < 5; i++)
            {
                v1[i] = i + 1;
                v2[i] = i + 2;
            }

            // SIMD-accelerated operations
            SimdMath.Add(v1.AsSpan(), v2.AsSpan(), result.AsSpan());

            Console.WriteLine($"    Vector addition: [{string.Join(", ", result.AsSpan().Slice(0, 5).ToArray())}]");
            Console.WriteLine($"    Vector magnitude squared: {VectorizedOperations.SumOfSquares(v1.AsSpan())}");
        }

        private static void MatrixOperations()
        {
            Console.WriteLine("\n  Matrix Operations:");

            using var A = new Matrix(3, 3);
            using var B = new Matrix(3, 3);

            // Initialize matrices (make A non-singular)
            A[0, 0] = 2; A[0, 1] = 1; A[0, 2] = -1;
            A[1, 0] = -3; A[1, 1] = -1; A[1, 2] = 2;
            A[2, 0] = -2; A[2, 1] = 1; A[2, 2] = 2;

            B[0, 0] = 8; B[0, 1] = 1; B[0, 2] = 3;
            B[1, 0] = -11; B[1, 1] = 4; B[1, 2] = 6;
            B[2, 0] = -3; B[2, 1] = 7; B[2, 2] = 9;

            Console.WriteLine($"    Matrix trace: {A.Trace()}");
            Console.WriteLine($"    Matrix determinant: {GaussianElimination.Determinant(A)}");
        }

        private static void ComplexOperations()
        {
            Console.WriteLine("\n  Complex Number Operations:");

            var c1 = new Complex(3, 4);
            var c2 = new Complex(1, -2);

            Console.WriteLine($"    Complex 1: {c1}");
            Console.WriteLine($"    Complex 2: {c2}");
            Console.WriteLine($"    Sum: {c1 + c2}");
            Console.WriteLine($"    Product: {c1 * c2}");
            Console.WriteLine($"    Magnitude: {c1.Magnitude:F2}");
        }

        #endregion

        #region Performance Examples

        private static void MemoryOptimization()
        {
            Console.WriteLine("\n  Memory Optimization:");

            using var context = new MathComputationContext(MathWorkload.LinearAlgebra, MathProblemSize.Large);
            var monitor = new MathPerformanceMonitor();
            monitor.Start();

            // Demonstrate optimized matrix creation
            using var matrix = MatrixLayoutOptimizer.CreateOptimizedMatrix(
                1000, 1000, MatrixOperation.RowWise, context);

            var report = monitor.Stop();
            Console.WriteLine($"    Optimized matrix created: {matrix.Rows}x{matrix.Columns}");
            Console.WriteLine($"    Memory efficiency: {report.MemoryEfficiency:F2} operations/MB");
        }

        private static void ComputationContexts()
        {
            Console.WriteLine("\n  Computation Contexts:");

            // Demonstrate context switching
            using var context = new MathComputationContext(MathWorkload.GeneralPurpose);

            using (context.EnterWorkload(MathWorkload.SignalProcessing, MathProblemSize.Large))
            {
                Console.WriteLine($"    Switched to: {context.CurrentWorkload} ({context.CurrentSize})");
                using var signalBuffer = context.MemoryPool.Rent<double>(1024);
                Console.WriteLine($"    Allocated signal buffer: {signalBuffer.Length} samples");
            }

            Console.WriteLine($"    Restored context: {context.CurrentWorkload} ({context.CurrentSize})");
        }

        private static void Benchmarking()
        {
            Console.WriteLine("\n  Performance Monitoring:");

            var monitor = new MathPerformanceMonitor();
            monitor.Start();

            // Simulate some work
            using var context = new MathComputationContext(MathWorkload.SmallComputations);
            using var data = context.RentBuffer(1000);

            for (int i = 0; i < 1000; i++)
            {
                data[i] = Math.Sin(i * 0.01);
            }

            var report = monitor.Stop();
            Console.WriteLine($"    Operations monitored: {report.AllocationsMade}");
            Console.WriteLine($"    Total time: {report.TotalTime.TotalMilliseconds:F2}ms");
        }

        #endregion

        #region Advanced Examples

        private static void LinearAlgebra()
        {
            Console.WriteLine("\n  Linear Algebra:");

            using var A = new Matrix(3, 3);
            using var b = new Vector(3);

            // Set up a simple linear system: 2x + y - z = 8, -3x - y + 2z = -11, -2x + y + 2z = -3
            A[0, 0] = 2; A[0, 1] = 1; A[0, 2] = -1; b[0] = 8;
            A[1, 0] = -3; A[1, 1] = -1; A[1, 2] = 2; b[1] = -11;
            A[2, 0] = -2; A[2, 1] = 1; A[2, 2] = 2; b[2] = -3;

            Vector x;
            GaussianElimination.Solve(A, b, out x);

            Console.WriteLine($"    Linear system solved: x ≈ {x[0]:F2}, y ≈ {x[1]:F2}, z ≈ {x[2]:F2}");
        }

        private static void SignalProcessing()
        {
            Console.WriteLine("\n  Signal Processing:");

            // Create a signal with multiple frequency components
            using var signal = new ComplexVector(1024);
            for (int i = 0; i < 1024; i++)
            {
                double t = 2 * Math.PI * i / 1024;
                signal[i] = new Complex(
                    Math.Sin(2 * t) + 0.5 * Math.Sin(10 * t), // Real part
                    0); // Imaginary part
            }

            var fftResult = FFT.Forward(signal);
            Console.WriteLine($"    FFT completed for {signal.Length} samples");
            Console.WriteLine($"    First few frequency components: {fftResult[0].Magnitude:F2}, {fftResult[1].Magnitude:F2}");
        }

        private static void NumericalMethods()
        {
            Console.WriteLine("\n  Numerical Methods:");

            // Integrate f(x) = x² * sin(x) from 0 to π
            double IntegralFunction(double x) => x * x * Math.Sin(x);

            double result = NumericalIntegration.AdaptiveQuadrature(IntegralFunction, 0, Math.PI, 1e-8);
            Console.WriteLine($"    ∫₀^π x²⋅sin(x) dx ≈ {result:F6}");

            // Monte Carlo integration
            double mcResult = NumericalIntegration.MonteCarlo(IntegralFunction, 0, Math.PI, 10000);
            Console.WriteLine($"    Monte Carlo estimate: {mcResult:F6}");
        }

        private static void StatisticalAnalysis()
        {
            Console.WriteLine("\n  Statistical Analysis:");

            // Generate sample data
            var data = new double[1000];
            var random = new Random(42);
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = Math.Sin(i * 0.1) + random.NextDouble() * 0.2;
            }

            double mean = DescriptiveStats.Mean(data);
            double variance = DescriptiveStats.Variance(data);
            double skewness = DescriptiveStats.Skewness(data);
            double kurtosis = DescriptiveStats.Kurtosis(data);

            Console.WriteLine($"    Sample size: {data.Length}");
            Console.WriteLine($"    Mean: {mean:F4}, Variance: {variance:F4}");
            Console.WriteLine($"    Skewness: {skewness:F4}, Kurtosis: {kurtosis:F4}");
        }

        #endregion
    }
}