using System;
using System.Runtime.CompilerServices;

namespace ZiggyMath.Core
{
    /// <summary>
    /// Advanced vector mathematics operations
    /// </summary>
    public static class VectorMath
    {
        /// <summary>
        /// Cross product of two 3D vectors
        /// </summary>
        public static Vector Cross(Vector left, Vector right)
        {
            if (left.Length != 3 || right.Length != 3)
                throw new ArgumentException("Cross product requires 3D vectors");

            var result = new Vector(3);
            result[0] = left[1] * right[2] - left[2] * right[1];
            result[1] = left[2] * right[0] - left[0] * right[2];
            result[2] = left[0] * right[1] - left[1] * right[0];

            return result;
        }

        /// <summary>
        /// Element-wise multiplication (Hadamard product)
        /// </summary>
        public static Vector HadamardProduct(Vector left, Vector right)
        {
            if (left.Length != right.Length)
                throw new ArgumentException("Vectors must have the same length");

            var result = new Vector(left.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < left.Length; i++)
            {
                resultSpan[i] = left[i] * right[i];
            }

            return result;
        }

        /// <summary>
        /// Element-wise division
        /// </summary>
        public static Vector ElementWiseDivide(Vector left, Vector right)
        {
            if (left.Length != right.Length)
                throw new ArgumentException("Vectors must have the same length");

            var result = new Vector(left.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < left.Length; i++)
            {
                if (right[i] == 0)
                    throw new DivideByZeroException($"Division by zero at index {i}");
                resultSpan[i] = left[i] / right[i];
            }

            return result;
        }

        /// <summary>
        /// Vector projection of 'a' onto 'b'
        /// </summary>
        public static Vector Project(Vector a, Vector b)
        {
            double dotProduct = a.DotProduct(b);
            double magnitudeSquared = b.MagnitudeSquared();

            if (magnitudeSquared == 0)
                throw new ArgumentException("Cannot project onto zero vector");

            return Vector.Multiply(b, dotProduct / magnitudeSquared);
        }

        /// <summary>
        /// Vector rejection of 'a' from 'b' (component perpendicular to b)
        /// </summary>
        public static Vector Reject(Vector a, Vector b)
        {
            return a - Project(a, b);
        }

        /// <summary>
        /// Linear interpolation between two vectors
        /// </summary>
        public static Vector Lerp(Vector from, Vector to, double t)
        {
            if (from.Length != to.Length)
                throw new ArgumentException("Vectors must have the same length");

            if (t < 0 || t > 1)
                throw new ArgumentOutOfRangeException(nameof(t), "Interpolation factor must be between 0 and 1");

            var result = new Vector(from.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < from.Length; i++)
            {
                resultSpan[i] = from[i] + t * (to[i] - from[i]);
            }

            return result;
        }

        /// <summary>
        /// Clamp vector elements to specified range
        /// </summary>
        public static Vector Clamp(Vector vector, double min, double max)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                double value = vector[i];
                resultSpan[i] = value < min ? min : (value > max ? max : value);
            }

            return result;
        }

        /// <summary>
        /// Normalize vector to unit length
        /// </summary>
        public static Vector Normalize(Vector vector)
        {
            double magnitude = vector.Magnitude();
            if (magnitude == 0)
                throw new ArgumentException("Cannot normalize zero vector");

            return Vector.Multiply(vector, 1.0 / magnitude);
        }

        /// <summary>
        /// Get absolute values of vector elements
        /// </summary>
        public static Vector Abs(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Abs(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get square root of vector elements
        /// </summary>
        public static Vector Sqrt(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Sqrt(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get exponential of vector elements
        /// </summary>
        public static Vector Exp(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Exp(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get natural logarithm of vector elements
        /// </summary>
        public static Vector Log(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Log(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get base-10 logarithm of vector elements
        /// </summary>
        public static Vector Log10(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Log10(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Raise vector elements to specified power
        /// </summary>
        public static Vector Pow(Vector vector, double exponent)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Pow(vector[i], exponent);
            }

            return result;
        }

        /// <summary>
        /// Get sine of vector elements
        /// </summary>
        public static Vector Sin(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Sin(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get cosine of vector elements
        /// </summary>
        public static Vector Cos(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Cos(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get tangent of vector elements
        /// </summary>
        public static Vector Tan(Vector vector)
        {
            var result = new Vector(vector.Length);
            var resultSpan = result.AsSpan();

            for (int i = 0; i < vector.Length; i++)
            {
                resultSpan[i] = Math.Tan(vector[i]);
            }

            return result;
        }

        /// <summary>
        /// Get minimum element value
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Min(Vector vector)
        {
            double min = vector[0];
            for (int i = 1; i < vector.Length; i++)
            {
                if (vector[i] < min)
                    min = vector[i];
            }
            return min;
        }

        /// <summary>
        /// Get maximum element value
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Max(Vector vector)
        {
            double max = vector[0];
            for (int i = 1; i < vector.Length; i++)
            {
                if (vector[i] > max)
                    max = vector[i];
            }
            return max;
        }

        /// <summary>
        /// Get sum of all elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sum(Vector vector)
        {
            double sum = 0.0;
            for (int i = 0; i < vector.Length; i++)
            {
                sum += vector[i];
            }
            return sum;
        }

        /// <summary>
        /// Get product of all elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Product(Vector vector)
        {
            double product = 1.0;
            for (int i = 0; i < vector.Length; i++)
            {
                product *= vector[i];
            }
            return product;
        }

        /// <summary>
        /// Get arithmetic mean of elements
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Mean(Vector vector)
        {
            return Sum(vector) / vector.Length;
        }

        /// <summary>
        /// Check if all elements satisfy a condition
        /// </summary>
        public static bool All(Vector vector, Predicate<double> predicate)
        {
            for (int i = 0; i < vector.Length; i++)
            {
                if (!predicate(vector[i]))
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Check if any element satisfies a condition
        /// </summary>
        public static bool Any(Vector vector, Predicate<double> predicate)
        {
            for (int i = 0; i < vector.Length; i++)
            {
                if (predicate(vector[i]))
                    return true;
            }
            return false;
        }

        /// <summary>
        /// Count elements that satisfy a condition
        /// </summary>
        public static int Count(Vector vector, Predicate<double> predicate)
        {
            int count = 0;
            for (int i = 0; i < vector.Length; i++)
            {
                if (predicate(vector[i]))
                    count++;
            }
            return count;
        }
    }
}