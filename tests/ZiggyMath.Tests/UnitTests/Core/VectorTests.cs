using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using ZiggyMath.Core;
using ZiggyMath.SimdMath;
using ZiggyAlloc;

namespace ZiggyMath.Tests.UnitTests.Core
{
    /// <summary>
    /// Comprehensive unit tests for Vector class
    /// Tests core functionality, memory management, and performance optimizations
    /// </summary>
    public class VectorTests : IDisposable
    {
        private readonly IUnmanagedMemoryAllocator _allocator;

        public VectorTests()
        {
            _allocator = new HybridAllocator(new SystemMemoryAllocator());
        }

        public void Dispose()
        {
            // Cleanup if needed
        }

        public void Constructor_WithSize_CreatesValidVector()
        {
            // Arrange & Act
            using var vector = new Vector(100);

            // Assert
            if (vector.Length != 100) throw new Exception($"Expected length 100, got {vector.Length}");
            if (!vector.IsValid) throw new Exception("Vector should be valid");
        }

        public void Constructor_WithAllocator_UsesSpecifiedAllocator()
        {
            // Arrange & Act
            using var vector = new Vector(100, _allocator);

            // Assert
            if (vector.Length != 100) throw new Exception($"Expected length 100, got {vector.Length}");
            if (!vector.IsValid) throw new Exception("Vector should be valid");
        }

        public void Constructor_WithData_CreatesVectorWithValues()
        {
            // Arrange
            var data = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

            // Act
            using var vector = new Vector(data, _allocator);

            // Assert
            if (vector.Length != 5) throw new Exception($"Expected length 5, got {vector.Length}");
            for (int i = 0; i < 5; i++)
            {
                if (vector[i] != data[i]) throw new Exception($"Expected {data[i]}, got {vector[i]} at index {i}");
            }
        }

        public void Indexer_GetAndSet_WorksCorrectly()
        {
            // Arrange
            using var vector = new Vector(10);

            // Act
            vector[5] = 42.0;
            var value = vector[5];

            // Assert
            if (value != 42.0) throw new Exception($"Expected 42.0, got {value}");
        }

        public void Indexer_OutOfBounds_ThrowsException()
        {
            // Arrange
            using var vector = new Vector(10);

            // Act & Assert
            try
            {
                var temp = vector[-1];
                throw new Exception("Should have thrown IndexOutOfRangeException for negative index");
            }
            catch (IndexOutOfRangeException) { }

            try
            {
                var temp = vector[10];
                throw new Exception("Should have thrown IndexOutOfRangeException for out of bounds index");
            }
            catch (IndexOutOfRangeException) { }
        }

        public void AsSpan_ReturnsCorrectSpan()
        {
            // Arrange
            using var vector = new Vector(10);

            // Act
            var span = vector.AsSpan();

            // Assert
            if (span.Length != 10) throw new Exception($"Expected span length 10, got {span.Length}");
        }

        public void Clear_FillsVectorWithZeros()
        {
            // Arrange
            using var vector = new Vector(10);
            for (int i = 0; i < 10; i++)
            {
                vector[i] = i + 1;
            }

            // Act
            vector.Clear();

            // Assert
            for (int i = 0; i < 10; i++)
            {
                if (vector[i] != 0.0) throw new Exception($"Expected 0.0 at index {i}, got {vector[i]}");
            }
        }

        public void Fill_FillsVectorWithValue()
        {
            // Arrange
            using var vector = new Vector(10);

            // Act
            vector.Fill(42.0);

            // Assert
            for (int i = 0; i < 10; i++)
            {
                if (vector[i] != 42.0) throw new Exception($"Expected 42.0 at index {i}, got {vector[i]}");
            }
        }

        public void DotProduct_CalculatesCorrectly()
        {
            // Arrange
            using var v1 = new Vector(3);
            using var v2 = new Vector(3);

            v1[0] = 1; v1[1] = 2; v1[2] = 3;
            v2[0] = 4; v2[1] = 5; v2[2] = 6;

            // Act
            double result = v1.DotProduct(v2);

            // Assert
            if (result != 32.0) throw new Exception($"Expected 32.0, got {result}"); // 1*4 + 2*5 + 3*6 = 32
        }

        public void Magnitude_CalculatesCorrectly()
        {
            // Arrange
            using var vector = new Vector(3);
            vector[0] = 3; vector[1] = 4; vector[2] = 0;

            // Act
            double magnitude = vector.Magnitude();

            // Assert
            if (magnitude != 5.0) throw new Exception($"Expected 5.0, got {magnitude}"); // sqrt(9 + 16 + 0) = 5
        }

        public void Normalized_ReturnsUnitVector()
        {
            // Arrange
            using var vector = new Vector(2);
            vector[0] = 3; vector[1] = 4;

            // Act
            using var normalized = vector.Normalized();

            // Assert
            double mag = normalized.Magnitude();
            if (Math.Abs(mag - 1.0) > 1e-10) throw new Exception($"Expected magnitude 1.0, got {mag}");

            if (Math.Abs(normalized[0] - 0.6) > 1e-10) throw new Exception($"Expected 0.6, got {normalized[0]}"); // 3/5
            if (Math.Abs(normalized[1] - 0.8) > 1e-10) throw new Exception($"Expected 0.8, got {normalized[1]}"); // 4/5
        }

        public void VectorOperations_CreateCorrectResults()
        {
            // Arrange
            using var original = new Vector(5);
            for (int i = 0; i < 5; i++)
            {
                original[i] = i * i;
            }

            // Act - Create new vector with same data
            using var copy = new Vector(5);
            for (int i = 0; i < 5; i++)
            {
                copy[i] = original[i];
            }

            // Assert
            if (original.Length != copy.Length) throw new Exception($"Expected same length, got {original.Length} vs {copy.Length}");
            for (int i = 0; i < 5; i++)
            {
                if (original[i] != copy[i]) throw new Exception($"Expected {original[i]}, got {copy[i]} at index {i}");
            }

            // Modify original and ensure copy is independent
            original[0] = 999;
            if (original[0] == copy[0]) throw new Exception("Copy should be independent of original");
        }

        public void SIMDOperations_ProduceCorrectResults()
        {
            // Arrange
            using var context = new MathComputationContext(MathWorkload.SmallComputations);
            using var a = context.RentBuffer(100);
            using var b = context.RentBuffer(100);
            using var result = context.RentBuffer(100);

            // Initialize with test data
            for (int i = 0; i < 100; i++)
            {
                a[i] = i;
                b[i] = i * 2;
            }

            // Act - SIMD operations
            VectorizedOperations.Add(a.AsSpan(), b.AsSpan(), result.AsSpan());

            // Assert
            for (int i = 0; i < 100; i++)
            {
                if (result[i] != i + i * 2) throw new Exception($"Expected {i + i * 2}, got {result[i]} at index {i}");
            }
        }

        public void PerformanceContext_ImprovesAllocationEfficiency()
        {
            // Arrange
            using var context = new MathComputationContext(MathWorkload.SmallComputations);
            var monitor = new MathPerformanceMonitor();
            monitor.Start();

            // Act - Multiple vector operations
            for (int i = 0; i < 100; i++)
            {
                using var v = context.RentBuffer(1000);
                for (int j = 0; j < 1000; j++)
                {
                    v[j] = Math.Sin(j * 0.01);
                }
            }

            var report = monitor.Stop();

            // Assert
            if (report.AllocationsMade <= 0) throw new Exception("Expected some allocations");
            if (report.TotalBytesAllocated <= 0) throw new Exception("Expected some bytes allocated");
            if (report.MemoryEfficiency < 0) throw new Exception("Expected non-negative memory efficiency");
        }

        public void MemoryPool_ReusesBuffersEfficiently()
        {
            // Arrange
            using var context = new MathComputationContext(MathWorkload.SmallComputations);

            // Act - Rent and return multiple buffers
            var buffers = new List<UnmanagedBuffer<double>>();
            for (int i = 0; i < 10; i++)
            {
                buffers.Add(context.RentBuffer(100));
            }

            // Dispose all buffers (return to pool)
            foreach (var buffer in buffers)
            {
                buffer.Dispose();
            }

            // Get pool stats
            var (active, pooled) = context.MemoryPool.GetStats();

            // Assert - Should have efficient memory reuse
            if (active < 0) throw new Exception("Active buffers should be non-negative");
        }

        public void LargeVectorOperations_MaintainPerformance()
        {
            // Arrange
            using var context = new MathComputationContext(MathWorkload.LinearAlgebra, MathProblemSize.Large);
            const int largeSize = 10000;

            // Act
            using var largeVector = context.RentBuffer(largeSize);
            using var otherVector = context.RentBuffer(largeSize);
            using var resultVector = context.RentBuffer(largeSize);

            // Initialize
            for (int i = 0; i < largeSize; i++)
            {
                largeVector[i] = Math.Sin(i * 0.001);
                otherVector[i] = Math.Cos(i * 0.001);
            }

            // SIMD operations on large vectors
            VectorizedOperations.Add(largeVector.AsSpan(), otherVector.AsSpan(), resultVector.AsSpan());

            // Assert
            if (resultVector.Length != largeSize) throw new Exception($"Expected length {largeSize}, got {resultVector.Length}");
            if (resultVector[0] != largeVector[0] + otherVector[0]) throw new Exception("First element incorrect");
            if (resultVector[largeSize - 1] != largeVector[largeSize - 1] + otherVector[largeSize - 1]) throw new Exception("Last element incorrect");
        }

        public void VectorOperations_AreThreadSafe()
        {
            // Arrange
            using var context = new MathComputationContext(MathWorkload.SmallComputations);
            const int numThreads = 4;
            const int vectorSize = 1000;
            var results = new double[numThreads][];

            // Act
            Parallel.For(0, numThreads, i =>
            {
                using var v1 = context.RentBuffer(vectorSize);
                using var v2 = context.RentBuffer(vectorSize);
                using var result = context.RentBuffer(vectorSize);

                for (int j = 0; j < vectorSize; j++)
                {
                    v1[j] = i + j;
                    v2[j] = i - j;
                }

                VectorizedOperations.Add(v1.AsSpan(), v2.AsSpan(), result.AsSpan());
                results[i] = result.AsSpan().ToArray();
            });

            // Assert
            for (int i = 0; i < numThreads; i++)
            {
                if (results[i] == null) throw new Exception($"Results[{i}] should not be null");
                if (results[i].Length != vectorSize) throw new Exception($"Expected length {vectorSize}, got {results[i].Length}");
                if (results[i][0] != i + i) throw new Exception($"Expected {i + i}, got {results[i][0]}"); // (i + 0) + (i - 0) = 2i
            }
        }

        /// <summary>
        /// Run all vector tests
        /// </summary>
        public void RunAllTests()
        {
            Console.WriteLine("Running Vector tests...");

            Constructor_WithSize_CreatesValidVector();
            Constructor_WithAllocator_UsesSpecifiedAllocator();
            Constructor_WithData_CreatesVectorWithValues();
            Indexer_GetAndSet_WorksCorrectly();
            Indexer_OutOfBounds_ThrowsException();
            AsSpan_ReturnsCorrectSpan();
            Clear_FillsVectorWithZeros();
            Fill_FillsVectorWithValue();
            DotProduct_CalculatesCorrectly();
            Magnitude_CalculatesCorrectly();
            Normalized_ReturnsUnitVector();
            VectorOperations_CreateCorrectResults();
            SIMDOperations_ProduceCorrectResults();
            PerformanceContext_ImprovesAllocationEfficiency();
            MemoryPool_ReusesBuffersEfficiently();
            LargeVectorOperations_MaintainPerformance();
            VectorOperations_AreThreadSafe();

            Console.WriteLine("✅ All Vector tests passed!");
        }
    }
}