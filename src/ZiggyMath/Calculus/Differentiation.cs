using System;
using ZiggyMath.Core;

namespace ZiggyMath.Calculus
{
    /// <summary>
    /// Numerical differentiation
    /// </summary>
    public static class Differentiation
    {
        public delegate double Function(double x);

        /// <summary>
        /// Forward difference approximation
        /// </summary>
        public static double ForwardDifference(Function f, double x, double h)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");

            return (f(x + h) - f(x)) / h;
        }

        /// <summary>
        /// Backward difference approximation
        /// </summary>
        public static double BackwardDifference(Function f, double x, double h)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");

            return (f(x) - f(x - h)) / h;
        }

        /// <summary>
        /// Central difference approximation
        /// </summary>
        public static double CentralDifference(Function f, double x, double h)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");

            return (f(x + h) - f(x - h)) / (2 * h);
        }

        /// <summary>
        /// Richardson extrapolation for improved accuracy
        /// </summary>
        public static double RichardsonExtrapolation(Function f, double x, double h)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");

            double h2 = h / 2;
            double f1 = CentralDifference(f, x, h);
            double f2 = CentralDifference(f, x, h2);

            return f2 + (f2 - f1) / 3.0;
        }

        /// <summary>
        /// Second derivative using central difference
        /// </summary>
        public static double SecondDerivative(Function f, double x, double h)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");

            return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
        }

        /// <summary>
        /// Nth derivative using finite differences
        /// </summary>
        public static double NthDerivative(Function f, double x, double h, int n)
        {
            if (h <= 0)
                throw new ArgumentException("Step size h must be positive");
            if (n < 1)
                throw new ArgumentException("Derivative order must be positive");

            if (n == 1)
                return CentralDifference(f, x, h);
            if (n == 2)
                return SecondDerivative(f, x, h);

            // Use recursive approach for higher derivatives
            return NthDerivativeRecursive(f, x, h, n);
        }

        /// <summary>
        /// Compute gradient of multivariate function
        /// </summary>
        public static Vector Gradient(Func<double[], double> f, double[] point)
        {
            if (point == null || point.Length == 0)
                throw new ArgumentException("Point cannot be null or empty");

            int n = point.Length;
            double h = 1e-8;
            Vector gradient = new Vector(n);

            for (int i = 0; i < n; i++)
            {
                // Create function for partial derivative
                Function partial = x =>
                {
                    double[] point_i = new double[n];
                    Array.Copy(point, point_i, n);
                    point_i[i] = x;
                    return f(point_i);
                };

                gradient[i] = CentralDifference(partial, point[i], h);
            }

            return gradient;
        }

        /// <summary>
        /// Compute Jacobian matrix of vector-valued function
        /// </summary>
        public static Matrix Jacobian(Func<double[], double>[] functions, double[] point)
        {
            if (functions == null || functions.Length == 0)
                throw new ArgumentException("Functions array cannot be null or empty");
            if (point == null || point.Length == 0)
                throw new ArgumentException("Point cannot be null or empty");

            int m = functions.Length; // Number of functions
            int n = point.Length;     // Number of variables

            Matrix jacobian = new Matrix(m, n);
            double h = 1e-8;

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    // Create function for partial derivative
                    Func<double[], double> partial = x =>
                    {
                        double[] point_ij = new double[n];
                        Array.Copy(point, point_ij, n);
                        point_ij[j] = x[0];
                        return functions[i](point_ij);
                    };

                    // Create 1D function for differentiation
                    Function func1D = x => partial(new double[] { x });

                    jacobian[i, j] = CentralDifference(func1D, point[j], h);
                }
            }

            return jacobian;
        }

        /// <summary>
        /// Compute Hessian matrix of multivariate function
        /// </summary>
        public static Matrix Hessian(Func<double[], double> f, double[] point)
        {
            if (point == null || point.Length == 0)
                throw new ArgumentException("Point cannot be null or empty");

            int n = point.Length;
            Matrix hessian = new Matrix(n);

            // Compute gradient at point
            Vector gradient = Gradient(f, point);

            double h = 1e-6;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    if (i == j)
                    {
                        // Diagonal elements
                        hessian[i, j] = SecondDerivative(x => f(new double[] { x }), point[i], h);
                    }
                    else
                    {
                        // Off-diagonal elements using mixed partial derivatives
                        Func<double[], double> mixed = x =>
                        {
                            double[] point_ij = new double[n];
                            Array.Copy(point, point_ij, n);
                            point_ij[i] = x[0];
                            point_ij[j] = x[1];
                            return f(point_ij);
                        };

                        Function func2D = x =>
                        {
                            double[] args = { point[i], point[j] };
                            args[i] = x;
                            return mixed(args);
                        };

                        hessian[i, j] = CentralDifference(func2D, point[i], h);
                        hessian[j, i] = hessian[i, j]; // Symmetric
                    }
                }
            }

            return hessian;
        }

        /// <summary>
        /// Recursive helper for Nth derivative
        /// </summary>
        private static double NthDerivativeRecursive(Function f, double x, double h, int n)
        {
            if (n == 1)
                return CentralDifference(f, x, h);
            if (n == 2)
                return SecondDerivative(f, x, h);

            // Use central difference formula for higher derivatives
            double h2 = h * h;
            double f_plus = NthDerivativeRecursive(f, x + h, h, n - 1);
            double f_minus = NthDerivativeRecursive(f, x - h, h, n - 1);

            return (f_plus - f_minus) / (2 * h);
        }
    }
}