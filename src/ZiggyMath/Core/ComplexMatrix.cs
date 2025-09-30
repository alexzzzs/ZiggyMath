using System;
using ZiggyAlloc;
using ZiggyMath.SimdMath;

namespace ZiggyMath.Core
{
    /// <summary>
    /// 2D matrix implementation for complex numbers using ZiggyAlloc for high-performance memory management.
    /// Because sometimes you need to do math with imaginary friends... and real performance!
    /// </summary>
    public class ComplexMatrix : IDisposable
    {
        /// <summary>
        /// The underlying unmanaged buffer that stores our complex data.
        /// Thanks to ZiggyAlloc, this bad boy is lightning fast and memory-efficient!
        /// </summary>
        private readonly UnmanagedBuffer<Complex> _buffer;

        /// <summary>
        /// How many rows of complex numbers we're dealing with.
        /// Because every matrix needs to know its dimensions - it's not just a pretty face!
        /// </summary>
        private readonly int _rows;

        /// <summary>
        /// How many columns of complex numbers we have.
        /// Rows and columns: the dynamic duo of linear algebra!
        /// </summary>
        private readonly int _cols;

        /// <summary>
        /// Create a new complex matrix with the specified dimensions.
        /// Perfect for when you need to store complex data in a structured way.
        /// </summary>
        /// <param name="rows">Number of rows (must be positive, or we'll throw a tantrum!)</param>
        /// <param name="cols">Number of columns (also must be positive, no negativity allowed!)</param>
        /// <exception cref="ArgumentException">Thrown when dimensions are zero or negative</exception>
        public ComplexMatrix(int rows, int cols)
        {
            // Input validation - because even mathematicians make mistakes!
            if (rows <= 0 || cols <= 0)
                throw new ArgumentException("Rows and columns must be positive. " +
                    "We're doing math here, not trying to create black holes!");

            _rows = rows;
            _cols = cols;

            // Allocate our high-performance buffer using ZiggyAlloc
            // This is where the magic happens - fast, efficient, and leak-free!
            _buffer = Z.DefaultAllocator.Allocate<Complex>(rows * cols);
        }

        /// <summary>
        /// Create a new complex matrix with the specified dimensions and allocator
        /// </summary>
        public ComplexMatrix(int rows, int cols, IUnmanagedMemoryAllocator allocator)
        {
            if (rows <= 0 || cols <= 0)
                throw new ArgumentException("Rows and columns must be positive");

            _rows = rows;
            _cols = cols;
            _buffer = allocator.Allocate<Complex>(rows * cols);
        }

        /// <summary>
        /// Number of rows in the matrix
        /// </summary>
        public int Rows => _rows;

        /// <summary>
        /// Number of columns in the matrix
        /// </summary>
        public int Columns => _cols;

        /// <summary>
        /// Total number of elements in the matrix
        /// </summary>
        public int Count => _rows * _cols;

        /// <summary>
        /// Get or set a matrix element by row and column index.
        /// Because sometimes you just need to poke a specific complex number!
        /// </summary>
        /// <param name="row">The row index (0-based, like a well-behaved array)</param>
        /// <param name="col">The column index (also 0-based, no favoritism here!)</param>
        /// <returns>A reference to the complex number at the specified position</returns>
        /// <exception cref="IndexOutOfRangeException">When you try to access indices that don't exist</exception>
        public ref Complex this[int row, int col]
        {
            get
            {
                // Bounds checking - because even mathematicians occasionally count from 1 instead of 0!
                if (row < 0 || row >= _rows || col < 0 || col >= _cols)
                    throw new IndexOutOfRangeException($"Index [{row},{col}] is out of range. " +
                        $"Valid range is [0,{_rows-1}] for rows and [0,{_cols-1}] for columns. " +
                        "Remember: arrays are 0-based, not 1-based!");

                // Row-major order: row * cols + col
                // This is like reading a book: left to right, top to bottom!
                return ref _buffer[row * _cols + col];
            }
        }

        /// <summary>
        /// Get or set a matrix element by linear index
        /// </summary>
        public ref Complex this[int index]
        {
            get
            {
                if (index < 0 || index >= Count)
                    throw new IndexOutOfRangeException($"Index {index} is out of range");

                return ref _buffer[index];
            }
        }

        /// <summary>
        /// Get a span of the matrix data for high-performance operations
        /// </summary>
        public Span<Complex> AsSpan() => _buffer.AsSpan();

        /// <summary>
        /// Get a readonly span of the matrix data
        /// </summary>
        public ReadOnlySpan<Complex> AsReadOnlySpan() => _buffer.AsReadOnlySpan();

        /// <summary>
        /// Fill the matrix with a specific complex value
        /// </summary>
        public void Fill(Complex value)
        {
            _buffer.Fill(value);
        }

        /// <summary>
        /// Fill the matrix with zeros
        /// </summary>
        public void Clear()
        {
            _buffer.Clear();
        }

        /// <summary>
        /// Copy data from another matrix
        /// </summary>
        public void CopyFrom(ComplexMatrix source)
        {
            if (source._rows != _rows || source._cols != _cols)
                throw new ArgumentException("Matrix dimensions must match");

            source._buffer.AsSpan().CopyTo(_buffer.AsSpan());
        }

        /// <summary>
        /// Create a copy of this matrix
        /// </summary>
        public ComplexMatrix Clone()
        {
            var clone = new ComplexMatrix(_rows, _cols);
            _buffer.AsSpan().CopyTo(clone._buffer.AsSpan());
            return clone;
        }

        /// <summary>
        /// Create an identity matrix
        /// </summary>
        public static ComplexMatrix CreateIdentity(int size)
        {
            var matrix = new ComplexMatrix(size, size);
            for (int i = 0; i < size; i++)
            {
                matrix[i, i] = new Complex(1.0, 0.0);
            }
            return matrix;
        }

        /// <summary>
        /// Create a matrix filled with a specific value
        /// </summary>
        public static ComplexMatrix CreateFilled(int rows, int cols, Complex value)
        {
            var matrix = new ComplexMatrix(rows, cols);
            matrix.Fill(value);
            return matrix;
        }

        /// <summary>
        /// Matrix addition
        /// </summary>
        public static ComplexMatrix operator +(ComplexMatrix left, ComplexMatrix right)
        {
            if (left._rows != right._rows || left._cols != right._cols)
                throw new ArgumentException("Matrix dimensions must match for addition");

            var result = new ComplexMatrix(left._rows, left._cols);
            var leftSpan = left._buffer.AsSpan();
            var rightSpan = right._buffer.AsSpan();
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < leftSpan.Length; i++)
            {
                resultSpan[i] = leftSpan[i] + rightSpan[i];
            }

            return result;
        }

        /// <summary>
        /// Matrix subtraction
        /// </summary>
        public static ComplexMatrix operator -(ComplexMatrix left, ComplexMatrix right)
        {
            if (left._rows != right._rows || left._cols != right._cols)
                throw new ArgumentException("Matrix dimensions must match for subtraction");

            var result = new ComplexMatrix(left._rows, left._cols);
            var leftSpan = left._buffer.AsSpan();
            var rightSpan = right._buffer.AsSpan();
            var resultSpan = result._buffer.AsSpan();

            for (int i = 0; i < leftSpan.Length; i++)
            {
                resultSpan[i] = leftSpan[i] - rightSpan[i];
            }

            return result;
        }

        /// <summary>
        /// Matrix multiplication - where the magic happens!
        /// This is O(n³) complexity, so be careful with large matrices.
        /// </summary>
        /// <param name="left">Left matrix (the one on the left side of the * operator)</param>
        /// <param name="right">Right matrix (the one on the right side of the * operator)</param>
        /// <returns>A new matrix containing the product</returns>
        /// <exception cref="ArgumentException">When matrix dimensions don't allow multiplication</exception>
        public static ComplexMatrix operator *(ComplexMatrix left, ComplexMatrix right)
        {
            // Dimension compatibility check - because matrix multiplication has standards!
            if (left._cols != right._rows)
                throw new ArgumentException(
                    $"Cannot multiply {left._rows}x{left._cols} matrix with {right._rows}x{right._cols} matrix. " +
                    $"Left matrix columns ({left._cols}) must equal right matrix rows ({right._rows}). " +
                    "It's like trying to multiply apples by oranges - dimensions must match!");

            // Create result matrix with appropriate dimensions
            // Result is (left.Rows x right.Columns) - basic matrix multiplication rule!
            var result = new ComplexMatrix(left._rows, right._cols);

            // The classic triple-nested loop for matrix multiplication
            // This is the "naive" O(n³) algorithm, but it works and is easy to understand
            // TODO: Optimize with SIMD operations for better performance on large matrices
            for (int i = 0; i < left._rows; i++)        // For each row in left matrix
            {
                for (int j = 0; j < right._cols; j++)   // For each column in right matrix
                {
                    // Calculate the dot product of row i from left and column j from right
                    Complex sum = new Complex(0, 0);    // Start with the additive identity

                    for (int k = 0; k < left._cols; k++) // For each element in the row/column
                    {
                        // Multiply corresponding elements and add to sum
                        // This is where the complex arithmetic happens!
                        sum += left[i, k] * right[k, j];
                    }

                    // Store the result - voila, one element of the product matrix!
                    result[i, j] = sum;
                }
            }

            // Return the shiny new product matrix
            // The original matrices are unchanged (immutable operations are good!)
            return result;
        }

        /// <summary>
        /// Scalar multiplication
        /// </summary>
        public static ComplexMatrix operator *(ComplexMatrix matrix, Complex scalar)
        {
            var result = new ComplexMatrix(matrix._rows, matrix._cols);

            // Use vectorized operations for scalar multiplication
            var scalarSpan = new ReadOnlySpan<Complex>(new[] { scalar });
            var resultSpan = result._buffer.AsSpan();
            var matrixSpan = matrix._buffer.AsSpan();

            for (int i = 0; i < matrixSpan.Length; i++)
            {
                resultSpan[i] = matrixSpan[i] * scalar;
            }

            return result;
        }

        /// <summary>
        /// Scalar multiplication (commutative)
        /// </summary>
        public static ComplexMatrix operator *(Complex scalar, ComplexMatrix matrix)
        {
            return matrix * scalar;
        }

        /// <summary>
        /// Get the raw pointer to the underlying buffer (for native interop)
        /// </summary>
        public IntPtr GetRawPointer() => _buffer.RawPointer;

        /// <summary>
        /// Get the conjugate of the matrix
        /// </summary>
        public ComplexMatrix Conjugate()
        {
            var result = new ComplexMatrix(_rows, _cols);
            for (int i = 0; i < Count; i++)
            {
                result._buffer[i] = _buffer[i].Conjugate();
            }
            return result;
        }

        /// <summary>
        /// Get the transpose of the matrix
        /// </summary>
        public ComplexMatrix Transpose()
        {
            var result = new ComplexMatrix(_cols, _rows);
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    result[j, i] = this[i, j];
                }
            }
            return result;
        }

        /// <summary>
        /// Get the Hermitian transpose (conjugate transpose) of the matrix
        /// </summary>
        public ComplexMatrix ConjugateTranspose()
        {
            return Transpose().Conjugate();
        }

        public void Dispose()
        {
            _buffer.Dispose();
        }

        /// <summary>
        /// Get string representation of the matrix
        /// </summary>
        public override string ToString()
        {
            return $"ComplexMatrix({_rows}x{_cols})";
        }
    }
}