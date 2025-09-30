using System;
using ZiggyAlloc;

namespace ZiggyMath.Core
{
    /// <summary>
    /// High-performance vector with unmanaged memory backing
    /// </summary>
    public readonly struct Vector : IDisposable
    {
        private readonly UnmanagedBuffer<double> _buffer;
        private readonly IUnmanagedMemoryAllocator _allocator;

        public int Length { get; }
        public int SizeInBytes => Length * sizeof(double);
        public bool IsValid => _buffer.IsValid;

        /// <summary>
        /// Create a new vector with specified size
        /// </summary>
        public Vector(int size, IUnmanagedMemoryAllocator allocator = null)
        {
            if (size <= 0)
                throw new ArgumentOutOfRangeException(nameof(size), "Vector size must be positive");

            Length = size;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<double>(size);
        }

        /// <summary>
        /// Create a vector from existing data
        /// </summary>
        public Vector(ReadOnlySpan<double> data, IUnmanagedMemoryAllocator allocator = null)
        {
            Length = data.Length;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<double>(Length);

            // Copy data into the buffer
            data.CopyTo(_buffer.AsSpan());
        }

        /// <summary>
        /// Element access with bounds checking
        /// </summary>
        public double this[int index]
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
        /// Vector addition
        /// </summary>
        public static Vector Add(Vector left, Vector right)
        {
            if (left.Length != right.Length)
                throw new ArgumentException("Vectors must have the same length");

            var result = new Vector(left.Length, left._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < left.Length; i++)
            {
                resultSpan[i] = left[i] + right[i];
            }

            return result;
        }

        /// <summary>
        /// Vector subtraction
        /// </summary>
        public static Vector Subtract(Vector left, Vector right)
        {
            if (left.Length != right.Length)
                throw new ArgumentException("Vectors must have the same length");

            var result = new Vector(left.Length, left._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < left.Length; i++)
            {
                resultSpan[i] = left[i] - right[i];
            }

            return result;
        }

        /// <summary>
        /// Scalar multiplication
        /// </summary>
        public static Vector Multiply(Vector vector, double scalar)
        {
            var result = new Vector(vector.Length, vector._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = vector[i] * scalar;
            }

            return result;
        }

        /// <summary>
        /// Scalar multiplication (commutative)
        /// </summary>
        public static Vector Multiply(double scalar, Vector vector)
        {
            return Multiply(vector, scalar);
        }

        /// <summary>
        /// Dot product
        /// </summary>
        public double DotProduct(Vector other)
        {
            if (Length != other.Length)
                throw new ArgumentException("Vectors must have the same length");

            double result = 0.0;
            for (int i = 0; i < Length; i++)
            {
                result += _buffer[i] * other[i];
            }

            return result;
        }

        /// <summary>
        /// Vector magnitude (L2 norm)
        /// </summary>
        public double Magnitude()
        {
            return Math.Sqrt(MagnitudeSquared());
        }

        /// <summary>
        /// Squared magnitude (more efficient for comparisons)
        /// </summary>
        public double MagnitudeSquared()
        {
            double sum = 0.0;
            for (int i = 0; i < Length; i++)
            {
                sum += _buffer[i] * _buffer[i];
            }
            return sum;
        }

        /// <summary>
        /// Normalize the vector in-place
        /// </summary>
        public void Normalize()
        {
            double magnitude = Magnitude();
            if (magnitude > 0)
            {
                var normalized = Multiply(this, 1.0 / magnitude);
                normalized.CopyTo(this);
            }
        }

        /// <summary>
        /// Get normalized copy of the vector
        /// </summary>
        public Vector Normalized()
        {
            Vector result = this;
            result.Normalize();
            return result;
        }

        /// <summary>
        /// Distance between two vectors
        /// </summary>
        public static double Distance(Vector left, Vector right)
        {
            return (left - right).Magnitude();
        }

        /// <summary>
        /// Angle between two vectors
        /// </summary>
        public static double Angle(Vector left, Vector right)
        {
            double dot = left.DotProduct(right);
            double magnitudeProduct = left.Magnitude() * right.Magnitude();

            if (magnitudeProduct == 0)
                return 0;

            return Math.Acos(dot / magnitudeProduct);
        }

        /// <summary>
        /// Get vector as Span<double>
        /// </summary>
        public Span<double> AsSpan()
        {
            return _buffer.AsSpan();
        }

        /// <summary>
        /// Get vector as ReadOnlySpan<double>
        /// </summary>
        public ReadOnlySpan<double> AsReadOnlySpan()
        {
            return _buffer.AsReadOnlySpan();
        }

        /// <summary>
        /// Copy vector to destination
        /// </summary>
        public void CopyTo(Vector destination)
        {
            if (destination.Length != Length)
                throw new ArgumentException("Destination vector must have the same length");

            _buffer.AsSpan().CopyTo(destination._buffer.AsSpan());
        }

        /// <summary>
        /// Clear all elements to zero
        /// </summary>
        public void Clear()
        {
            _buffer.AsSpan().Clear();
        }

        /// <summary>
        /// Fill vector with specified value
        /// </summary>
        public void Fill(double value)
        {
            _buffer.AsSpan().Fill(value);
        }

        /// <summary>
        /// Implicit operator to subtract vectors
        /// </summary>
        public static Vector operator -(Vector left, Vector right)
        {
            return Subtract(left, right);
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
            return $"Vector(Length={Length}, IsValid={IsValid})";
        }
    }
}