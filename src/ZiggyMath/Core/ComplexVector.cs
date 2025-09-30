using System;
using ZiggyAlloc;

namespace ZiggyMath.Core
{
    /// <summary>
    /// Complex vector for signal processing
    /// </summary>
    public class ComplexVector : IDisposable
    {
        private readonly UnmanagedBuffer<Complex> _buffer;
        private readonly IUnmanagedMemoryAllocator _allocator;

        public int Length { get; }
        public bool IsValid => _buffer.IsValid;

        /// <summary>
        /// Create a new complex vector with specified size
        /// </summary>
        public ComplexVector(int size, IUnmanagedMemoryAllocator allocator = null)
        {
            if (size <= 0)
                throw new ArgumentOutOfRangeException(nameof(size), "Vector size must be positive");

            Length = size;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<Complex>(size);
        }

        /// <summary>
        /// Create a complex vector from existing data
        /// </summary>
        public ComplexVector(ReadOnlySpan<Complex> data, IUnmanagedMemoryAllocator allocator = null)
        {
            Length = data.Length;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<Complex>(Length);
            data.CopyTo(_buffer.AsSpan());
        }

        /// <summary>
        /// Element access with bounds checking
        /// </summary>
        public Complex this[int index]
        {
            get
            {
                if ((uint)index >= (uint)Length)
                    throw new IndexOutOfRangeException($"Index {index} is out of range [0, {Length})");

                return _buffer[index];
            }
            set
            {
                if ((uint)index >= (uint)Length)
                    throw new IndexOutOfRangeException($"Index {index} is out of range [0, {Length})");

                _buffer[index] = value;
            }
        }

        /// <summary>
        /// Get vector as Span<Complex>
        /// </summary>
        public Span<Complex> AsSpan()
        {
            return _buffer.AsSpan();
        }

        /// <summary>
        /// Dispose of unmanaged resources
        /// </summary>
        public void Dispose()
        {
            _buffer.Dispose();
        }

        /// <summary>
        /// Get string representation
        /// </summary>
        public override string ToString()
        {
            return $"ComplexVector(Length={Length}, IsValid={IsValid})";
        }
    }
}