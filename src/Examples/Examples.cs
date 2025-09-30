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
    /// Examples demonstrating Phase 2: Core Infrastructure Implementation
    /// Shows how to use MathComputationContext, MathPerformanceMonitor,
    /// MatrixLayoutOptimizer, and enhanced SIMD operations
    /// </summary>
    public static class Phase2Examples
    {
        public static void RunAllExamples()
        {
            Console.WriteLine("🚀 ZiggyMath Phase 2 Examples");
            Console.WriteLine("================================");

            BasicUsageExample();
            AdvancedUsageExample();
            PerformanceBenchmarkExample();
        }

        /// <summary>
        /// Basic usage of Phase 2 components
        /// </summary>
        private static void BasicUsageExample()
        {
            Console.WriteLine("\n📋 Basic Usage Example");
            Console.WriteLine("----------------------");

            // Create optimized computational context
            using var context = new MathComputationContext(
                MathWorkload.LinearAlgebra,
                MathProblemSize.Medium);

            // Create optimized matrix
            using var matrix = MatrixLayoutOptimizer.CreateOptimizedMatrix(
                1000, 1000, MatrixOperation.RowWise, context);

            // Monitor performance
            var monitor = new MathPerformanceMonitor();
            monitor.Start();

            using (context.EnterWorkload(MathWorkload.SignalProcessing, MathProblemSize.Large))
            {
                // Signal processing operations with optimized allocator
                using var signalBuffer = context.MemoryPool.Rent<double>(1024);
                SimdMath.Normalize(signalBuffer.AsSpan());
            }

            var report = monitor.Stop();
            Console.WriteLine($"Operations: {report.AllocationsMade}, Memory: {report.TotalBytesAllocated} bytes");
        }

        /// <summary>
        /// Advanced usage with nested contexts and performance monitoring
        /// </summary>
        private static void AdvancedUsageExample()
        {
            Console.WriteLine("\n📊 Advanced Usage Example");
            Console.WriteLine("-------------------------");

            // Performance-optimized matrix operations
            using var context = new MathComputationContext(MathWorkload.LinearAlgebra, MathProblemSize.Large);
            var monitor = new MathPerformanceMonitor();
            monitor.Start();

            using var a = MatrixLayoutOptimizer.CreateOptimizedMatrix(1000, 1000, MatrixOperation.RowWise, context);
            using var b = MatrixLayoutOptimizer.CreateOptimizedMatrix(1000, 1000, MatrixOperation.ColumnWise, context);

            monitor.RecordOperation("MatrixMultiply", () =>
            {
                MatrixLayoutOptimizer.MultiplyOptimized(a, b, context);
            });

            var report = monitor.Stop();
            Console.WriteLine($"Matrix multiplication completed");
            Console.WriteLine($"Memory efficiency: {report.MemoryEfficiency:F2} operations/byte");
        }

        /// <summary>
        /// Performance benchmark comparing different approaches
        /// </summary>
        private static void PerformanceBenchmarkExample()
        {
            Console.WriteLine("\n⚡ Performance Benchmark Example");
            Console.WriteLine("--------------------------------");

            const int size = 10000;

            // Benchmark standard approach
            var standardMonitor = new MathPerformanceMonitor();
            standardMonitor.Start();

            var standardResult = BenchmarkStandardApproach(size);

            var standardReport = standardMonitor.Stop();
            Console.WriteLine($"Standard approach: {standardReport.TotalTime.TotalMilliseconds:F2}ms, " +
                            $"{standardReport.AllocationsMade} allocations");

            // Benchmark optimized approach
            var optimizedMonitor = new MathPerformanceMonitor();
            optimizedMonitor.Start();

            var optimizedResult = BenchmarkOptimizedApproach(size);

            var optimizedReport = optimizedMonitor.Stop();
            Console.WriteLine($"Optimized approach: {optimizedReport.TotalTime.TotalMilliseconds:F2}ms, " +
                            $"{optimizedReport.AllocationsMade} allocations");

            // Calculate improvement
            double improvement = ((standardReport.TotalTime.TotalMilliseconds - optimizedReport.TotalTime.TotalMilliseconds)
                                / standardReport.TotalTime.TotalMilliseconds) * 100;
            Console.WriteLine($"Performance improvement: {improvement:F1}%");
        }

        private static double[] BenchmarkStandardApproach(int size)
        {
            var a = new double[size];
            var b = new double[size];
            var result = new double[size];

            // Initialize with random data
            var random = new Random(42);
            for (int i = 0; i < size; i++)
            {
                a[i] = random.NextDouble();
                b[i] = random.NextDouble();
            }

            // Standard vector addition
            for (int i = 0; i < size; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        private static double[] BenchmarkOptimizedApproach(int size)
        {
            using var context = new MathComputationContext(MathWorkload.SmallComputations, MathProblemSize.Large);

            using var a = context.RentBuffer(size);
            using var b = context.RentBuffer(size);
            using var result = context.RentBuffer(size);

            // Initialize with random data
            var random = new Random(42);
            for (int i = 0; i < size; i++)
            {
                a[i] = random.NextDouble();
                b[i] = random.NextDouble();
            }

            // Optimized SIMD vector addition
            SimdMath.Add(a.AsSpan(), b.AsSpan(), result.AsSpan());

            return result.AsSpan().ToArray();
        }

    }
}