using System;
using ZiggyMath.Core;

namespace ZiggyMath.Calculus
{
    /// <summary>
    /// Numerical integration methods
    /// </summary>
    public static class NumericalIntegration
    {
        public delegate double Function(double x);
        public delegate double Function2D(double x, double y);

        /// <summary>
        /// Trapezoidal rule for 1D integration
        /// </summary>
        public static double TrapezoidalRule(Function f, double a, double b, int intervals)
        {
            if (intervals <= 0)
                throw new ArgumentException("Number of intervals must be positive");

            double h = (b - a) / intervals;
            double sum = 0.5 * (f(a) + f(b));

            for (int i = 1; i < intervals; i++)
            {
                double x = a + i * h;
                sum += f(x);
            }

            return sum * h;
        }

        /// <summary>
        /// Simpson's rule for 1D integration
        /// </summary>
        public static double SimpsonRule(Function f, double a, double b, int intervals)
        {
            if (intervals <= 0)
                throw new ArgumentException("Number of intervals must be positive");
            if (intervals % 2 != 0)
                throw new ArgumentException("Number of intervals must be even for Simpson's rule");

            double h = (b - a) / intervals;
            double sum = f(a) + f(b);

            for (int i = 1; i < intervals; i++)
            {
                double x = a + i * h;
                if (i % 2 == 0)
                    sum += 2 * f(x);
                else
                    sum += 4 * f(x);
            }

            return sum * h / 3.0;
        }

        /// <summary>
        /// Adaptive quadrature for 1D integration
        /// </summary>
        public static double AdaptiveQuadrature(Function f, double a, double b, double tolerance)
        {
            return AdaptiveQuadratureRecursive(f, a, b, tolerance, 0);
        }

        /// <summary>
        /// Trapezoidal rule for 2D integration
        /// </summary>
        public static double TrapezoidalRule2D(Function2D f, double ax, double bx, double ay, double by, int nx, int ny)
        {
            if (nx <= 0 || ny <= 0)
                throw new ArgumentException("Number of intervals must be positive");

            double hx = (bx - ax) / nx;
            double hy = (by - ay) / ny;

            double sum = 0.0;

            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    double x = ax + i * hx;
                    double y = ay + j * hy;
                    double weight = 1.0;

                    if (i == 0 || i == nx)
                        weight *= 0.5;
                    if (j == 0 || j == ny)
                        weight *= 0.5;

                    sum += weight * f(x, y);
                }
            }

            return sum * hx * hy;
        }

        /// <summary>
        /// Monte Carlo integration for 1D functions
        /// </summary>
        public static double MonteCarlo(Function f, double a, double b, int samples)
        {
            if (samples <= 0)
                throw new ArgumentException("Number of samples must be positive");

            Random random = new Random();
            double sum = 0.0;

            for (int i = 0; i < samples; i++)
            {
                double x = a + (b - a) * random.NextDouble();
                sum += f(x);
            }

            return (b - a) * sum / samples;
        }

        /// <summary>
        /// Monte Carlo integration for 2D functions
        /// </summary>
        public static double MonteCarlo2D(Function2D f, double ax, double bx, double ay, double by, int samples)
        {
            if (samples <= 0)
                throw new ArgumentException("Number of samples must be positive");

            Random random = new Random();
            double sum = 0.0;

            for (int i = 0; i < samples; i++)
            {
                double x = ax + (bx - ax) * random.NextDouble();
                double y = ay + (by - ay) * random.NextDouble();
                sum += f(x, y);
            }

            return (bx - ax) * (by - ay) * sum / samples;
        }

        /// <summary>
        /// Gaussian quadrature (2-point) for 1D integration
        /// </summary>
        public static double GaussianQuadrature(Function f, double a, double b)
        {
            // 2-point Gaussian quadrature weights and points
            double[] weights = { 1.0, 1.0 };
            double[] points = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };

            double sum = 0.0;
            double mid = (a + b) / 2.0;
            double halfLength = (b - a) / 2.0;

            for (int i = 0; i < 2; i++)
            {
                double x = mid + halfLength * points[i];
                sum += weights[i] * f(x);
            }

            return halfLength * sum;
        }

        /// <summary>
        /// Recursive helper for adaptive quadrature
        /// </summary>
        private static double AdaptiveQuadratureRecursive(Function f, double a, double b, double tolerance, int depth)
        {
            const int maxDepth = 20;

            if (depth > maxDepth)
            {
                return TrapezoidalRule(f, a, b, 1000);
            }

            double mid = (a + b) / 2.0;
            double left = TrapezoidalRule(f, a, mid, 1);
            double right = TrapezoidalRule(f, mid, b, 1);
            double whole = TrapezoidalRule(f, a, b, 2);

            double error = Math.Abs(left + right - whole);

            if (error < tolerance)
            {
                return left + right;
            }
            else
            {
                return AdaptiveQuadratureRecursive(f, a, mid, tolerance / 2, depth + 1) +
                       AdaptiveQuadratureRecursive(f, mid, b, tolerance / 2, depth + 1);
            }
        }
    }
}