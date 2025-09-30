using System;

namespace ZiggyMath.Core
{
    /// <summary>
    /// Complex number structure
    /// </summary>
    public struct Complex : IEquatable<Complex>
    {
        public double Real { get; }
        public double Imaginary { get; }
        public double Magnitude => Math.Sqrt(Real * Real + Imaginary * Imaginary);
        public double Phase => Math.Atan2(Imaginary, Real);

        /// <summary>
        /// Create a complex number
        /// </summary>
        public Complex(double real, double imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        /// <summary>
        /// Create a real number (imaginary part = 0)
        /// </summary>
        public Complex(double real) : this(real, 0) { }

        /// <summary>
        /// Create complex number from polar coordinates
        /// </summary>
        public static Complex FromPolar(double magnitude, double phase)
        {
            return new Complex(
                magnitude * Math.Cos(phase),
                magnitude * Math.Sin(phase)
            );
        }

        /// <summary>
        /// Addition
        /// </summary>
        public static Complex operator +(Complex left, Complex right)
        {
            return new Complex(
                left.Real + right.Real,
                left.Imaginary + right.Imaginary
            );
        }

        /// <summary>
        /// Subtraction
        /// </summary>
        public static Complex operator -(Complex left, Complex right)
        {
            return new Complex(
                left.Real - right.Real,
                left.Imaginary - right.Imaginary
            );
        }

        /// <summary>
        /// Multiplication
        /// </summary>
        public static Complex operator *(Complex left, Complex right)
        {
            return new Complex(
                left.Real * right.Real - left.Imaginary * right.Imaginary,
                left.Real * right.Imaginary + left.Imaginary * right.Real
            );
        }

        /// <summary>
        /// Division
        /// </summary>
        public static Complex operator /(Complex left, Complex right)
        {
            double denominator = right.Real * right.Real + right.Imaginary * right.Imaginary;
            if (denominator == 0)
                throw new DivideByZeroException("Cannot divide by zero complex number");

            return new Complex(
                (left.Real * right.Real + left.Imaginary * right.Imaginary) / denominator,
                (left.Imaginary * right.Real - left.Real * right.Imaginary) / denominator
            );
        }

        /// <summary>
        /// Unary plus
        /// </summary>
        public static Complex operator +(Complex complex)
        {
            return complex;
        }

        /// <summary>
        /// Unary minus
        /// </summary>
        public static Complex operator -(Complex complex)
        {
            return new Complex(-complex.Real, -complex.Imaginary);
        }

        /// <summary>
        /// Equality comparison
        /// </summary>
        public static bool operator ==(Complex left, Complex right)
        {
            return left.Real == right.Real && left.Imaginary == right.Imaginary;
        }

        /// <summary>
        /// Inequality comparison
        /// </summary>
        public static bool operator !=(Complex left, Complex right)
        {
            return !(left == right);
        }

        /// <summary>
        /// Complex conjugate
        /// </summary>
        public Complex Conjugate()
        {
            return new Complex(Real, -Imaginary);
        }

        /// <summary>
        /// Complex reciprocal
        /// </summary>
        public Complex Reciprocal()
        {
            double denominator = Real * Real + Imaginary * Imaginary;
            if (denominator == 0)
                throw new DivideByZeroException("Cannot take reciprocal of zero complex number");

            return new Complex(Real / denominator, -Imaginary / denominator);
        }

        /// <summary>
        /// Check equality with another complex number
        /// </summary>
        public bool Equals(Complex other)
        {
            return this == other;
        }

        /// <summary>
        /// Check equality with object
        /// </summary>
        public override bool Equals(object obj)
        {
            return obj is Complex complex && Equals(complex);
        }

        /// <summary>
        /// Get hash code
        /// </summary>
        public override int GetHashCode()
        {
            return HashCode.Combine(Real, Imaginary);
        }

        /// <summary>
        /// Get string representation
        /// </summary>
        public override string ToString()
        {
            if (Imaginary == 0)
                return $"{Real}";
            else if (Real == 0)
                return $"{Imaginary}i";
            else if (Imaginary < 0)
                return $"{Real}{Imaginary}i";
            else
                return $"{Real}+{Imaginary}i";
        }
    }
}