using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace ZiggyMath.Core
{
    /// <summary>
    /// Metrics for a single operation execution
    /// </summary>
    public record OperationMetrics
    {
        public string OperationName { get; init; }
        public TimeSpan StartTime { get; init; }
        public TimeSpan EndTime { get; init; }
        public long MemoryUsed { get; init; }
        public int AllocationsMade { get; init; }
        public TimeSpan Duration => EndTime - StartTime;
    }

    /// <summary>
    /// High-precision performance monitor for mathematical computations.
    /// Because you can't optimize what you don't measure - and we measure everything!
    /// </summary>
    public class MathPerformanceMonitor : IDisposable
    {
        /// <summary>
        /// High-precision stopwatch for measuring computation time.
        /// Because milliseconds matter when you're optimizing for performance!
        /// </summary>
        private readonly Stopwatch _totalTime = new();

        /// <summary>
        /// Baseline memory usage when monitoring started.
        /// So we can track exactly how much memory our math operations consume.
        /// </summary>
        private long _baselineMemory;

        /// <summary>
        /// Total number of allocations made during monitoring.
        /// Because allocation count is as important as allocation size!
        /// </summary>
        private int _allocationCount;

        /// <summary>
        /// Total bytes allocated during monitoring period.
        /// The foundation of memory efficiency metrics.
        /// </summary>
        private long _totalBytesAllocated;

        /// <summary>
        /// Whether monitoring is currently active.
        /// Because you shouldn't measure what you're not tracking!
        /// </summary>
        private bool _isMonitoring;

        /// <summary>
        /// Start performance monitoring for mathematical operations.
        /// Let the measurement begin!
        /// </summary>
        public void Start()
        {
            if (_isMonitoring)
                throw new InvalidOperationException("Performance monitoring is already active. Stop it first!");

            _baselineMemory = GC.GetTotalMemory(false);
            _totalTime.Restart();
            _allocationCount = 0;
            _totalBytesAllocated = 0;
            _isMonitoring = true;

            Console.WriteLine("📊 Performance monitoring started");
        }

        /// <summary>
        /// Record a memory allocation for tracking efficiency.
        /// Every allocation tells a story about performance!
        /// </summary>
        /// <param name="byteCount">Number of bytes allocated</param>
        public void RecordAllocation(int byteCount)
        {
            if (!_isMonitoring)
                return; // Silently ignore if not monitoring

            _allocationCount++;
            _totalBytesAllocated += byteCount;
        }

        /// <summary>
        /// Record the execution time and metrics of an operation.
        /// </summary>
        /// <param name="operationName">Name of the operation being recorded</param>
        /// <param name="operation">The operation to execute and measure</param>
        /// <returns>Metrics about the operation execution</returns>
        public OperationMetrics RecordOperation(string operationName, Action operation)
        {
            var startTime = _totalTime.Elapsed;
            var startMemory = GC.GetTotalMemory(false);

            operation();

            var endTime = _totalTime.Elapsed;
            var endMemory = GC.GetTotalMemory(false);

            var metrics = new OperationMetrics
            {
                OperationName = operationName,
                StartTime = startTime,
                EndTime = endTime,
                MemoryUsed = endMemory - startMemory,
                AllocationsMade = _allocationCount
            };

            return metrics;
        }

        /// <summary>
        /// Stop monitoring and generate comprehensive performance report.
        /// Time to see how well our optimizations performed!
        /// </summary>
        /// <returns>Detailed performance metrics and optimization recommendations</returns>
        public MathPerformanceReport Stop()
        {
            if (!_isMonitoring)
                throw new InvalidOperationException("Performance monitoring is not active. Start it first!");

            _totalTime.Stop();
            _isMonitoring = false;

            var managedMemoryUsed = GC.GetTotalMemory(false) - _baselineMemory;

            var report = new MathPerformanceReport
            {
                TotalTime = _totalTime.Elapsed,
                AllocationsMade = _allocationCount,
                TotalBytesAllocated = _totalBytesAllocated,
                ManagedMemoryUsed = managedMemoryUsed,
                MemoryEfficiency = CalculateMemoryEfficiency(),
                AllocationRate = CalculateAllocationRate(),
                Recommendations = GenerateRecommendations()
            };

            Console.WriteLine($"📊 Performance monitoring stopped: {report.AllocationsMade} allocations, {report.TotalBytesAllocated:N0} bytes");
            return report;
        }

        /// <summary>
        /// Calculate memory efficiency metric (bytes per millisecond).
        /// Higher is better - shows how much data you process per unit time.
        /// </summary>
        private double CalculateMemoryEfficiency()
        {
            if (_totalTime.ElapsedMilliseconds == 0)
                return 0;

            return (double)_totalBytesAllocated / _totalTime.ElapsedMilliseconds;
        }

        /// <summary>
        /// Calculate allocation rate (allocations per millisecond).
        /// Lower is better for allocation-heavy workloads.
        /// </summary>
        private double CalculateAllocationRate()
        {
            if (_totalTime.ElapsedMilliseconds == 0)
                return 0;

            return (double)_allocationCount / _totalTime.ElapsedMilliseconds;
        }

        /// <summary>
        /// Generate intelligent recommendations based on performance data.
        /// Because monitoring without insights is just collecting numbers!
        /// </summary>
        private string[] GenerateRecommendations()
        {
            var recommendations = new List<string>();

            // Analyze allocation patterns
            if (_allocationCount > 1000 && _totalTime.ElapsedMilliseconds > 100)
            {
                recommendations.Add("High allocation frequency detected - consider using memory pools");
            }

            // Analyze memory efficiency
            var efficiency = CalculateMemoryEfficiency();
            if (efficiency < 1000) // Less than 1KB/ms
            {
                recommendations.Add("Low memory efficiency - consider larger batch operations");
            }

            // Analyze allocation sizes
            if (_allocationCount > 0)
            {
                double avgAllocationSize = (double)_totalBytesAllocated / _allocationCount;
                if (avgAllocationSize < 100)
                {
                    recommendations.Add("Many small allocations - slab allocator recommended");
                }
                else if (avgAllocationSize > 1000000)
                {
                    recommendations.Add("Large allocations detected - large block allocator recommended");
                }
            }

            return recommendations.ToArray();
        }

        /// <summary>
        /// Clean up monitoring resources.
        /// Because even performance monitors need proper cleanup!
        /// </summary>
        public void Dispose()
        {
            if (_isMonitoring)
            {
                Stop();
            }
        }
    }

    /// <summary>
    /// Comprehensive performance report for mathematical computations.
    /// Because understanding performance requires seeing the big picture!
    /// </summary>
    public class MathPerformanceReport
    {
        /// <summary>
        /// Total time spent in monitored operations.
        /// The foundation of all performance metrics.
        /// </summary>
        public TimeSpan TotalTime { get; init; }

        /// <summary>
        /// Number of memory allocations made during monitoring.
        /// Because allocation count matters as much as allocation size!
        /// </summary>
        public int AllocationsMade { get; init; }

        /// <summary>
        /// Total bytes allocated during monitoring period.
        /// The raw measure of memory consumption.
        /// </summary>
        public long TotalBytesAllocated { get; init; }

        /// <summary>
        /// Managed memory used by the garbage collector.
        /// Shows the GC pressure impact of operations.
        /// </summary>
        public long ManagedMemoryUsed { get; init; }

        /// <summary>
        /// Memory efficiency metric (bytes processed per millisecond).
        /// Higher values indicate better performance.
        /// </summary>
        public double MemoryEfficiency { get; init; }

        /// <summary>
        /// Allocation rate (allocations per millisecond).
        /// Lower values indicate better allocation efficiency.
        /// </summary>
        public double AllocationRate { get; init; }

        /// <summary>
        /// Intelligent recommendations for performance optimization.
        /// Because monitoring without actionable insights is just collecting data!
        /// </summary>
        public string[] Recommendations { get; init; } = Array.Empty<string>();

        /// <summary>
        /// Get a human-readable summary of the performance report.
        /// Because numbers are great, but understanding is better!
        /// </summary>
        public override string ToString()
        {
            return $"Performance Report: {AllocationsMade} allocations, {TotalBytesAllocated:N0} bytes, {TotalTime.TotalMilliseconds:F2}ms";
        }
    }
}