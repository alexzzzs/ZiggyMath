using System;
using ZiggyAlloc;

namespace ZiggyMath
{
    /// <summary>
    /// Main entry point for ZiggyMath library - your gateway to high-performance mathematics!
    /// This class provides optimized allocators and utilities for mathematical computing.
    /// </summary>
    public static class ZiggyMath
    {
        /// <summary>
        /// Default allocator optimized for mathematical computations.
        /// Uses a hybrid strategy that automatically chooses the best allocation method
        /// based on data type and size for optimal performance.
        /// </summary>
        public static readonly IUnmanagedMemoryAllocator DefaultAllocator =
            new HybridAllocator(new SystemMemoryAllocator());

        /// <summary>
        /// Create an allocator optimized for specific mathematical workloads.
        /// Because one size doesn't fit all - different math needs different memory strategies!
        /// </summary>
        /// <param name="workload">The type of mathematical work you'll be doing</param>
        /// <returns>An allocator fine-tuned for your specific use case</returns>
        public static IUnmanagedMemoryAllocator CreateAllocator(MathWorkload workload)
        {
            return CreateAllocator(workload, MathProblemSize.Medium);
        }

        /// <summary>
        /// Create an allocator optimized for specific mathematical workloads and problem sizes.
        /// Because one size doesn't fit all - different math needs different memory strategies!
        /// </summary>
        /// <param name="workload">The type of mathematical work you'll be doing</param>
        /// <param name="size">The expected size of mathematical problems</param>
        /// <returns>An allocator fine-tuned for your specific use case</returns>
        public static IUnmanagedMemoryAllocator CreateAllocator(MathWorkload workload, MathProblemSize size)
        {
            // Start with the system allocator as our base
            var baseAllocator = new SystemMemoryAllocator();

            // Choose the best allocator strategy based on workload type AND problem size
            return (workload, size) switch
            {
                // Linear Algebra optimizations based on size
                (MathWorkload.LinearAlgebra, MathProblemSize.Small) =>
                    new SlabAllocator(baseAllocator, slabSize: 512 * 1024), // Smaller slabs for small problems
                (MathWorkload.LinearAlgebra, MathProblemSize.Medium) =>
                    new SlabAllocator(baseAllocator, slabSize: 1024 * 1024), // Standard slabs for medium problems
                (MathWorkload.LinearAlgebra, MathProblemSize.Large) =>
                    new LargeBlockAllocator(new SlabAllocator(baseAllocator)), // Large block + slab for big problems

                // Signal Processing optimizations
                (MathWorkload.SignalProcessing, MathProblemSize.Small) =>
                    new LargeBlockAllocator(baseAllocator), // Even small signal processing needs large blocks
                (MathWorkload.SignalProcessing, MathProblemSize.Medium) =>
                    new LargeBlockAllocator(baseAllocator),
                (MathWorkload.SignalProcessing, MathProblemSize.Large) =>
                    new LargeBlockAllocator(new SlabAllocator(baseAllocator)), // Slab backing for very large problems

                // Small computations - maximum pooling for best performance
                (MathWorkload.SmallComputations, MathProblemSize.Small) =>
                    new UnmanagedMemoryPool(new SlabAllocator(baseAllocator, slabSize: 256 * 1024)),
                (MathWorkload.SmallComputations, MathProblemSize.Medium) =>
                    new UnmanagedMemoryPool(new SlabAllocator(baseAllocator, slabSize: 512 * 1024)),
                (MathWorkload.SmallComputations, MathProblemSize.Large) =>
                    new UnmanagedMemoryPool(new SlabAllocator(baseAllocator, slabSize: 1024 * 1024)),

                // General purpose: Let the hybrid allocator decide based on size
                (MathWorkload.GeneralPurpose, MathProblemSize.Small) =>
                    new HybridAllocator(baseAllocator),
                (MathWorkload.GeneralPurpose, MathProblemSize.Medium) =>
                    new HybridAllocator(baseAllocator),
                (MathWorkload.GeneralPurpose, MathProblemSize.Large) =>
                    new LargeBlockAllocator(baseAllocator),

                // Fallback to default for unknown combinations
                _ => DefaultAllocator
            };
        }

        /// <summary>
        /// Get library version information
        /// </summary>
        public static string Version => "1.0.0";

        /// <summary>
        /// Check if hardware acceleration is available
        /// </summary>
        public static bool IsHardwareAccelerationSupported =>
            System.Runtime.Intrinsics.X86.Avx2.IsSupported ||
            System.Runtime.Intrinsics.Arm.AdvSimd.IsSupported;
    }

    /// <summary>
    /// Mathematical workload types for allocator optimization
    /// </summary>
    public enum MathWorkload
    {
        /// <summary>
        /// General purpose mathematical operations
        /// </summary>
        GeneralPurpose,

        /// <summary>
        /// Linear algebra operations (matrices, vectors)
        /// </summary>
        LinearAlgebra,

        /// <summary>
        /// Signal processing operations (FFT, filtering)
        /// </summary>
        SignalProcessing,

        /// <summary>
        /// Small computational tasks
        /// </summary>
        SmallComputations
    }

    /// <summary>
    /// Mathematical problem sizes for fine-tuned allocator optimization
    /// </summary>
    public enum MathProblemSize
    {
        /// <summary>
        /// Small problems (typically < 1000 elements)
        /// </summary>
        Small,

        /// <summary>
        /// Medium problems (typically 1000-10000 elements)
        /// </summary>
        Medium,

        /// <summary>
        /// Large problems (typically > 10000 elements)
        /// </summary>
        Large
    }

    /// <summary>
    /// Matrix operation types for layout optimization
    /// </summary>
    public enum MatrixOperation
    {
        /// <summary>
        /// Row-wise operations (matrix-vector multiplication, row reductions)
        /// </summary>
        RowWise,

        /// <summary>
        /// Column-wise operations (column reductions, matrix-column operations)
        /// </summary>
        ColumnWise,

        /// <summary>
        /// Block-wise operations (matrix multiplication, block algorithms)
        /// </summary>
        BlockWise,

        /// <summary>
        /// Element-wise operations (Hadamard product, element transformations)
        /// </summary>
        ElementWise
    }
}