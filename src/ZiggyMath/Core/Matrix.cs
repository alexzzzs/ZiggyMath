using System;
using ZiggyAlloc;

namespace ZiggyMath.Core
{
    /// <summary>
    /// High-performance matrix with unmanaged memory backing
    /// </summary>
    public class Matrix : IDisposable
    {
        private readonly UnmanagedBuffer<double> _buffer;
        private readonly IUnmanagedMemoryAllocator _allocator;

        public int Rows { get; }
        public int Columns { get; }
        public int Size => Rows * Columns;
        public int SizeInBytes => Size * sizeof(double);
        public bool IsValid => _buffer.IsValid;
        public bool IsSquare => Rows == Columns;

        /// <summary>
        /// Create a new matrix with specified dimensions
        /// </summary>
        public Matrix(int rows, int columns, IUnmanagedMemoryAllocator allocator = null)
        {
            if (rows <= 0)
                throw new ArgumentOutOfRangeException(nameof(rows), "Rows must be positive");
            if (columns <= 0)
                throw new ArgumentOutOfRangeException(nameof(columns), "Columns must be positive");

            Rows = rows;
            Columns = columns;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<double>(Size);
        }

        /// <summary>
        /// Create a square matrix
        /// </summary>
        public Matrix(int size, IUnmanagedMemoryAllocator allocator = null) : this(size, size, allocator) { }

        /// <summary>
        /// Create a matrix from 2D array
        /// </summary>
        public Matrix(double[,] data, IUnmanagedMemoryAllocator allocator = null)
        {
            Rows = data.GetLength(0);
            Columns = data.GetLength(1);
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<double>(Size);

            // Copy data in row-major order
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    this[i, j] = data[i, j];
                }
            }
        }

        /// <summary>
        /// Create a matrix from 1D data in row-major order
        /// </summary>
        public Matrix(ReadOnlySpan<double> data, int rows, int columns, IUnmanagedMemoryAllocator allocator = null)
        {
            if (data.Length != rows * columns)
                throw new ArgumentException("Data length must match matrix dimensions");

            Rows = rows;
            Columns = columns;
            _allocator = allocator ?? Z.DefaultAllocator;
            _buffer = _allocator.Allocate<double>(Size);
            data.CopyTo(_buffer.AsSpan());
        }

        /// <summary>
        /// Element access with bounds checking
        /// </summary>
        public double this[int row, int column]
        {
            get
            {
                if ((uint)row >= (uint)Rows)
                    throw new IndexOutOfRangeException($"Row {row} is out of range [0, {Rows})");
                if ((uint)column >= (uint)Columns)
                    throw new IndexOutOfRangeException($"Column {column} is out of range [0, {Columns})");

                return _buffer[row * Columns + column];
            }
            set
            {
                if ((uint)row >= (uint)Rows)
                    throw new IndexOutOfRangeException($"Row {row} is out of range [0, {Rows})");
                if ((uint)column >= (uint)Columns)
                    throw new IndexOutOfRangeException($"Column {column} is out of range [0, {Columns})");

                _buffer[row * Columns + column] = value;
            }
        }

        /// <summary>
        /// Matrix addition
        /// </summary>
        public static Matrix Add(Matrix left, Matrix right)
        {
            if (left.Rows != right.Rows || left.Columns != right.Columns)
                throw new ArgumentException("Matrices must have the same dimensions");

            var result = new Matrix(left.Rows, left.Columns, left._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < left.Size; i++)
            {
                resultSpan[i] = left._buffer[i] + right._buffer[i];
            }

            return result;
        }

        /// <summary>
        /// Matrix subtraction
        /// </summary>
        public static Matrix Subtract(Matrix left, Matrix right)
        {
            if (left.Rows != right.Rows || left.Columns != right.Columns)
                throw new ArgumentException("Matrices must have the same dimensions");

            var result = new Matrix(left.Rows, left.Columns, left._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < left.Size; i++)
            {
                resultSpan[i] = left._buffer[i] - right._buffer[i];
            }

            return result;
        }

        /// <summary>
        /// Matrix multiplication
        /// </summary>
        public static Matrix Multiply(Matrix left, Matrix right)
        {
            if (left.Columns != right.Rows)
                throw new ArgumentException("Left matrix columns must equal right matrix rows");

            var result = new Matrix(left.Rows, right.Columns, left._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < left.Rows; i++)
            {
                for (int j = 0; j < right.Columns; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < left.Columns; k++)
                    {
                        sum += left[i, k] * right[k, j];
                    }
                    resultSpan[i * right.Columns + j] = sum;
                }
            }

            return result;
        }

        /// <summary>
        /// Matrix-vector multiplication
        /// </summary>
        public static Vector Multiply(Matrix matrix, Vector vector)
        {
            if (matrix.Columns != vector.Length)
                throw new ArgumentException("Matrix columns must equal vector length");

            var result = new Vector(matrix.Rows);

            for (int i = 0; i < matrix.Rows; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < matrix.Columns; j++)
                {
                    sum += matrix[i, j] * vector[j];
                }
                result[i] = sum;
            }

            return result;
        }

        /// <summary>
        /// Vector-matrix multiplication
        /// </summary>
        public static Vector Multiply(Vector vector, Matrix matrix)
        {
            if (vector.Length != matrix.Rows)
                throw new ArgumentException("Vector length must equal matrix rows");

            var result = new Vector(matrix.Columns);

            for (int j = 0; j < matrix.Columns; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < matrix.Rows; i++)
                {
                    sum += vector[i] * matrix[i, j];
                }
                result[j] = sum;
            }

            return result;
        }

        /// <summary>
        /// Scalar multiplication
        /// </summary>
        public static Matrix Multiply(Matrix matrix, double scalar)
        {
            var result = new Matrix(matrix.Rows, matrix.Columns, matrix._allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < matrix.Size; i++)
            {
                resultSpan[i] = matrix._buffer[i] * scalar;
            }

            return result;
        }

        /// <summary>
        /// Scalar multiplication (commutative)
        /// </summary>
        public static Matrix Multiply(double scalar, Matrix matrix)
        {
            return Multiply(matrix, scalar);
        }

        /// <summary>
        /// Matrix transpose
        /// </summary>
        public Matrix Transpose()
        {
            var result = new Matrix(Columns, Rows, _allocator);
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    resultSpan[j * Rows + i] = _buffer[i * Columns + j];
                }
            }

            return result;
        }

        /// <summary>
        /// Matrix inverse (placeholder - simplified implementation)
        /// </summary>
        public Matrix Inverse()
        {
            if (!IsSquare)
                throw new InvalidOperationException("Matrix must be square for inversion");

            // Placeholder implementation - would use proper LAPACK/algorithm in real version
            throw new NotImplementedException("Matrix inversion not yet implemented");
        }

        /// <summary>
        /// Matrix determinant (placeholder - simplified implementation)
        /// </summary>
        public double Determinant()
        {
            if (!IsSquare)
                throw new InvalidOperationException("Matrix must be square for determinant");

            // Placeholder implementation - would use proper algorithm in real version
            throw new NotImplementedException("Determinant calculation not yet implemented");
        }

        /// <summary>
        /// Matrix trace
        /// </summary>
        public double Trace()
        {
            if (!IsSquare)
                throw new InvalidOperationException("Matrix must be square for trace");

            double sum = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                sum += _buffer[i * Columns + i];
            }

            return sum;
        }

        /// <summary>
        /// QR decomposition (placeholder)
        /// </summary>
        public (Matrix Q, Matrix R) QRDecomposition()
        {
            throw new NotImplementedException("QR decomposition not yet implemented");
        }

        /// <summary>
        /// LU decomposition (placeholder)
        /// </summary>
        public (Matrix L, Matrix U) LUDecomposition()
        {
            throw new NotImplementedException("LU decomposition not yet implemented");
        }

        /// <summary>
        /// Eigen decomposition (placeholder)
        /// </summary>
        public (Vector eigenvalues, Matrix eigenvectors) EigenDecomposition()
        {
            throw new NotImplementedException("Eigen decomposition not yet implemented");
        }

        /// <summary>
        /// Get matrix as Span<double>
        /// </summary>
        public Span<double> AsSpan()
        {
            return _buffer.AsSpan();
        }

        /// <summary>
        /// Get matrix as ReadOnlySpan<double>
        /// </summary>
        public ReadOnlySpan<double> AsReadOnlySpan()
        {
            return _buffer.AsReadOnlySpan();
        }

        /// <summary>
        /// Convert to 2D array
        /// </summary>
        public double[] ToArray()
        {
            var result = new double[Size];
            _buffer.AsSpan().CopyTo(result);
            return result;
        }

        /// <summary>
        /// Copy matrix to destination
        /// </summary>
        public void CopyTo(Matrix destination)
        {
            if (destination.Rows != Rows || destination.Columns != Columns)
                throw new ArgumentException("Destination matrix must have the same dimensions");

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
        /// Fill matrix with specified value
        /// </summary>
        public void Fill(double value)
        {
            _buffer.AsSpan().Fill(value);
        }

        /// <summary>
        /// Create identity matrix
        /// </summary>
        public static Matrix Identity(int size)
        {
            var result = new Matrix(size);
            for (int i = 0; i < size; i++)
            {
                result[i, i] = 1.0;
            }
            return result;
        }

        /// <summary>
        /// Create zero matrix
        /// </summary>
        public static Matrix Zeros(int rows, int columns)
        {
            return new Matrix(rows, columns);
        }

        /// <summary>
        /// Create matrix filled with ones
        /// </summary>
        public static Matrix Ones(int rows, int columns)
        {
            var result = new Matrix(rows, columns);
            result.Fill(1.0);
            return result;
        }

        /// <summary>
        /// Create random matrix
        /// </summary>
        public static Matrix Random(int rows, int columns, double min = 0, double max = 1)
        {
            var result = new Matrix(rows, columns);
            var random = new Random();

            for (int i = 0; i < result.Size; i++)
            {
                result._buffer[i] = min + (max - min) * random.NextDouble();
            }

            return result;
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
            return $"Matrix({Rows}x{Columns}, IsValid={IsValid})";
        }
    }
}