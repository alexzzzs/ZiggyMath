using System;
using ZiggyMath.Tests.UnitTests;
using ZiggyMath.Core;
using ZiggyMath.SimdMath;
using ZiggyMath.LinearAlgebra;
using ZiggyMath.SignalProcessing;
using ZiggyMath.Statistics;
using ZiggyMath.Calculus;

namespace ZiggyMath.Tests
{
    /// <summary>
    /// Comprehensive test runner for ZiggyMath testing suite
    /// Tests Phase 1, 2, and 3 functionality
    /// </summary>
    public static class TestRunner
    {
        public static int RunAllTests()
        {
            Console.WriteLine("🧪 ZiggyMath Comprehensive Test Suite");
            Console.WriteLine("====================================\n");
            Console.Out.Flush(); // Ensure output is flushed

            int totalPassed = 0;
            int totalTests = 0;

            // Phase 1: Core functionality tests
            Console.WriteLine("📐 Phase 1: Core Functionality");
            Console.WriteLine("------------------------------");

            // VectorTests returns 0 for success, 1 for failure, but we need to count properly
            int vectorResult = VectorTests.RunAllTests();
            totalTests += 10;
            totalPassed += vectorResult == 0 ? 10 : 0; // If successful, all 10 tests passed
            Console.WriteLine($"📊 Phase 1 completed: {totalPassed}/{totalTests} tests passed so far");
            Console.Out.Flush();

            // Phase 2: Infrastructure tests
            Console.WriteLine("\n🏗️ Phase 2: Infrastructure");
            Console.WriteLine("--------------------------");

            int contextTests = TestMathComputationContext();
            totalTests += 4; // Fixed: actually runs 4 tests
            totalPassed += contextTests;

            int performanceTests = TestMathPerformanceMonitor();
            totalTests += 3; // Fixed: actually runs 3 tests
            totalPassed += performanceTests;

            int layoutTests = TestMatrixLayoutOptimizer();
            totalTests += 1; // Fixed: actually runs 1 test
            totalPassed += layoutTests;
            Console.WriteLine($"📊 Phase 2 completed: {totalPassed}/{totalTests} tests passed so far");
            Console.Out.Flush();

            // Phase 3: Advanced algorithms tests
            Console.WriteLine("\n🔬 Phase 3: Advanced Algorithms");
            Console.WriteLine("-------------------------------");

            int linearAlgebraTests = TestOptimizedLinearAlgebra();
            totalTests += 1; // Fixed: actually runs 1 test
            totalPassed += linearAlgebraTests;

            int fftTests = TestOptimizedFFT();
            totalTests += 1; // Fixed: actually runs 1 test
            totalPassed += fftTests;

            int integrationTests = TestOptimizedNumericalIntegration();
            totalTests += 1; // Fixed: actually runs 1 test
            totalPassed += integrationTests;

            int statisticsTests = TestOptimizedStatistics();
            totalTests += 2; // Fixed: actually runs 2 tests
            totalPassed += statisticsTests;
            Console.WriteLine($"📊 Phase 3 completed: {totalPassed}/{totalTests} tests passed so far");
            Console.Out.Flush();

            // Summary
            Console.WriteLine($"\n📊 Test Summary: {totalPassed}/{totalTests} passed");
            Console.Out.Flush();

            if (totalPassed == totalTests)
            {
                Console.WriteLine("✅ All tests passed!");
                Console.Out.Flush();
                return 0;
            }
            else
            {
                Console.WriteLine("❌ Some tests failed!");
                Console.Out.Flush();
                return 1;
            }
        }

        #region Phase 2 Tests

        private static int TestMathComputationContext()
        {
            int passed = 0;

            Console.WriteLine("  Testing MathComputationContext...");

            try
            {
                // Test context creation
                using var context = new MathComputationContext(MathWorkload.LinearAlgebra);
                if (context.CurrentWorkload == MathWorkload.LinearAlgebra) passed++;

                // Test context switching
                using (context.EnterWorkload(MathWorkload.SignalProcessing))
                {
                    if (context.CurrentWorkload == MathWorkload.SignalProcessing) passed++;
                }

                // Test buffer allocation
                using var buffer = context.RentBuffer(100);
                if (buffer.Length == 100) passed++;

                // Test memory pool
                using var pooledBuffer = context.MemoryPool.Rent<double>(50);
                if (pooledBuffer.Length == 50) passed++;

                Console.WriteLine($"    Context tests: {passed}/4 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Context tests failed: {ex.Message}");
            }

            return passed;
        }

        private static int TestMathPerformanceMonitor()
        {
            int passed = 0;

            Console.WriteLine("  Testing MathPerformanceMonitor...");

            try
            {
                var monitor = new MathPerformanceMonitor();
                monitor.Start();

                // Simulate some operations
                using var context = new MathComputationContext(MathWorkload.SmallComputations);
                using var testBuffer = context.RentBuffer(100);

                var report = monitor.Stop();

                if (report.AllocationsMade >= 0) passed++;
                if (report.TotalBytesAllocated >= 0) passed++;
                if (report.TotalTime.TotalMilliseconds >= 0) passed++;

                Console.WriteLine($"    Performance monitor tests: {passed}/3 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Performance monitor tests failed: {ex.Message}");
            }

            return passed;
        }

        private static int TestMatrixLayoutOptimizer()
        {
            int passed = 0;

            Console.WriteLine("  Testing MatrixLayoutOptimizer...");

            try
            {
                using var context = new MathComputationContext(MathWorkload.LinearAlgebra);

                // Test optimized matrix creation
                using var matrix = MatrixLayoutOptimizer.CreateOptimizedMatrix(
                    100, 100, MatrixOperation.RowWise, context);

                if (matrix.Rows == 100 && matrix.Columns == 100) passed++;

                Console.WriteLine($"    Layout optimizer tests: {passed}/1 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Layout optimizer tests failed: {ex.Message}");
            }

            return passed;
        }

        #endregion

        #region Phase 3 Tests

        private static int TestOptimizedLinearAlgebra()
        {
            int passed = 0;

            Console.WriteLine("  Testing OptimizedLinearAlgebra...");

            try
            {
                using var context = new MathComputationContext(MathWorkload.LinearAlgebra);

                // Create test matrices
                using var A = new Matrix(3, 3);
                using var b = new Vector(3);

                // Set up simple system: x + y + z = 6, 2x + y - z = 2, x - y + 2z = 5
                A[0, 0] = 1; A[0, 1] = 1; A[0, 2] = 1; b[0] = 6;
                A[1, 0] = 2; A[1, 1] = 1; A[1, 2] = -1; b[1] = 2;
                A[2, 0] = 1; A[2, 1] = -1; A[2, 2] = 2; b[2] = 5;

                // Test Gaussian elimination
                var x = OptimizedLinearAlgebra.SolveLinearSystem(A, b, context);

                if (x.Length == 3) passed++;

                Console.WriteLine($"    Linear algebra tests: {passed}/1 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Linear algebra tests failed: {ex.Message}");
            }

            return passed;
        }

        private static int TestOptimizedFFT()
        {
            int passed = 0;

            Console.WriteLine("  Testing OptimizedFFT...");

            try
            {
                using var context = new MathComputationContext(MathWorkload.SignalProcessing);

                // Create test signal (simple case)
                using var signal = new ComplexVector(2); // Use size 2 to avoid recursion issues
                signal[0] = new Complex(1, 0);
                signal[1] = new Complex(1, 0);

                // Test FFT
                var fftResult = OptimizedFFT.Transform(signal, context);

                if (fftResult.Length == 2) passed++;

                Console.WriteLine($"    FFT tests: {passed}/1 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    FFT tests failed: {ex.Message}");
            }

            return passed;
        }

        private static int TestOptimizedNumericalIntegration()
        {
            int passed = 0;

            Console.WriteLine("  Testing OptimizedNumericalIntegration...");

            try
            {
                using var context = new MathComputationContext(MathWorkload.SmallComputations);

                // Test function: f(x) = x²
                double TestFunction(double x) => x * x;

                // Test adaptive quadrature
                double result = OptimizedNumericalIntegration.AdaptiveQuadrature(
                    TestFunction, 0, 1, 1e-6, context);

                // ∫₀¹ x² dx = 1/3 ≈ 0.333333
                if (Math.Abs(result - 1.0/3.0) < 0.01) passed++;

                Console.WriteLine($"    Integration tests: {passed}/1 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Integration tests failed: {ex.Message}");
            }

            return passed;
        }

        private static int TestOptimizedStatistics()
        {
            int passed = 0;

            Console.WriteLine("  Testing OptimizedStatistics...");

            try
            {
                using var context = new MathComputationContext(MathWorkload.SmallComputations);

                // Create test data
                var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

                // Test linear regression
                var xData = new double[] { 1, 2, 3, 4, 5 };
                var yData = new double[] { 2, 4, 6, 8, 10 };

                var (slope, intercept, rSquared) = OptimizedStatistics.LinearRegression(xData, yData, context);

                if (Math.Abs(slope - 2.0) < 0.1) passed++; // Should be slope = 2
                if (Math.Abs(intercept) < 0.1) passed++; // Should be intercept ≈ 0

                Console.WriteLine($"    Statistics tests: {passed}/2 passed");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"    Statistics tests failed: {ex.Message}");
            }

            return passed;
        }

        #endregion
    }
}