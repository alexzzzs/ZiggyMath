using System;
using ZiggyMath.Core;
using ZiggyMath.LinearAlgebra;

namespace ZiggyMath.Calculus
{
    /// <summary>
    /// Optimization algorithms
    /// </summary>
    public static class Optimization
    {
        public delegate double ObjectiveFunction(double[] parameters);
        public delegate double[] GradientFunction(double[] parameters);

        /// <summary>
        /// Golden section search for 1D optimization
        /// </summary>
        public static double GoldenSectionSearch(ObjectiveFunction f, double a, double b, double tolerance)
        {
            const double goldenRatio = 1.618033988749;
            const int maxIterations = 1000;

            double x1 = b - (b - a) / goldenRatio;
            double x2 = a + (b - a) / goldenRatio;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                if (Math.Abs(b - a) < tolerance)
                    return (a + b) / 2.0;

                if (f(new double[] { x1 }) < f(new double[] { x2 }))
                {
                    b = x2;
                    x2 = x1;
                    x1 = b - (b - a) / goldenRatio;
                }
                else
                {
                    a = x1;
                    x1 = x2;
                    x2 = a + (b - a) / goldenRatio;
                }
            }

            return (a + b) / 2.0;
        }

        /// <summary>
        /// Brent's method for 1D optimization
        /// </summary>
        public static double BrentMethod(ObjectiveFunction f, double a, double b, double tolerance)
        {
            const int maxIterations = 1000;
            const double goldenRatio = 1.618033988749;

            double x = (a + b) / 2.0;
            double w = x;
            double v = x;

            double fx = f(new double[] { x });
            double fw = fx;
            double fv = fx;

            double d = double.MaxValue;
            double e = 0.0;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double m = 0.5 * (a + b);
                double tol1 = tolerance * Math.Abs(x) + 1e-10;
                double tol2 = 2.0 * tol1;

                if (Math.Abs(x - m) <= tol2 - 0.5 * (b - a))
                    return x;

                if (Math.Abs(e) > tol1)
                {
                    // Parabolic fit
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);

                    if (q > 0)
                        p = -p;
                    else
                        q = -q;

                    double etemp = e;
                    e = d;

                    if (Math.Abs(p) >= Math.Abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                    {
                        e = (x >= m) ? a - x : b - x;
                        d = goldenRatio * e;
                    }
                    else
                    {
                        d = p / q;
                        double u_new = x + d;

                        if (u_new - a < tol2 || b - u_new < tol2)
                            d = (m >= x) ? tol1 : -tol1;
                    }
                }
                else
                {
                    e = (x >= m) ? a - x : b - x;
                    d = goldenRatio * e;
                }

                double u = x + (Math.Abs(d) >= tol1 ? d : (d >= 0 ? tol1 : -tol1));
                double fu = f(new double[] { u });

                if (fu <= fx)
                {
                    if (u >= x)
                        a = x;
                    else
                        b = x;

                    v = w; w = x; x = u;
                    fv = fw; fw = fx; fx = fu;
                }
                else
                {
                    if (u < x)
                        a = u;
                    else
                        b = u;

                    if (fu <= fw || w == x)
                    {
                        v = w; w = u;
                        fv = fw; fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u;
                        fv = fu;
                    }
                }
            }

            return x;
        }

        /// <summary>
        /// Gradient descent optimization
        /// </summary>
        public static double[] GradientDescent(ObjectiveFunction f, GradientFunction gradient,
                                             double[] initialGuess, double learningRate, int maxIterations)
        {
            if (initialGuess == null || initialGuess.Length == 0)
                throw new ArgumentException("Initial guess cannot be null or empty");

            double[] x = new double[initialGuess.Length];
            Array.Copy(initialGuess, x, initialGuess.Length);

            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double[] grad = gradient(x);
                double gradientNorm = 0.0;

                for (int i = 0; i < x.Length; i++)
                {
                    gradientNorm += grad[i] * grad[i];
                }
                gradientNorm = Math.Sqrt(gradientNorm);

                if (gradientNorm < tolerance)
                    break;

                // Update parameters
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] -= learningRate * grad[i];
                }
            }

            return x;
        }

        /// <summary>
        /// Newton's method for optimization
        /// </summary>
        public static double[] NewtonMethod(ObjectiveFunction f, GradientFunction gradient,
                                          Matrix hessian, double[] initialGuess, int maxIterations)
        {
            if (initialGuess == null || initialGuess.Length == 0)
                throw new ArgumentException("Initial guess cannot be null or empty");

            double[] x = new double[initialGuess.Length];
            Array.Copy(initialGuess, x, initialGuess.Length);

            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double[] grad = gradient(x);
                // For Newton's method, we would need the Hessian matrix
                // For now, skip this optimization step if Hessian is not available
                break;

                double gradientNorm = 0.0;
                for (int i = 0; i < x.Length; i++)
                {
                    gradientNorm += grad[i] * grad[i];
                }
                gradientNorm = Math.Sqrt(gradientNorm);

                if (gradientNorm < tolerance)
                    break;

                // Newton's method implementation would go here
                // For now, this is a placeholder
                break;
            }

            return x;
        }

        /// <summary>
        /// Conjugate gradient optimization
        /// </summary>
        public static double[] ConjugateGradient(ObjectiveFunction f, GradientFunction gradient,
                                              double[] initialGuess, int maxIterations)
        {
            if (initialGuess == null || initialGuess.Length == 0)
                throw new ArgumentException("Initial guess cannot be null or empty");

            int n = initialGuess.Length;
            double[] x = new double[n];
            Array.Copy(initialGuess, x, n);

            double[] g = gradient(x);
            double[] d = new double[n];

            // Initial direction is negative gradient
            for (int i = 0; i < n; i++)
            {
                d[i] = -g[i];
            }

            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double gradientNorm = 0.0;
                for (int i = 0; i < n; i++)
                {
                    gradientNorm += g[i] * g[i];
                }
                gradientNorm = Math.Sqrt(gradientNorm);

                if (gradientNorm < tolerance)
                    break;

                // Line search (simplified)
                double alpha = 1.0;
                double[] x_new = new double[n];
                for (int i = 0; i < n; i++)
                {
                    x_new[i] = x[i] + alpha * d[i];
                }

                double[] g_new = gradient(x_new);

                // Update direction using conjugate gradient formula
                double beta_numerator = 0.0;
                double beta_denominator = 0.0;

                for (int i = 0; i < n; i++)
                {
                    beta_numerator += g_new[i] * g_new[i];
                    beta_denominator += g[i] * g[i];
                }

                double beta = (beta_denominator != 0) ? beta_numerator / beta_denominator : 0.0;

                // Update direction
                for (int i = 0; i < n; i++)
                {
                    d[i] = -g_new[i] + beta * d[i];
                }

                // Update position and gradient
                Array.Copy(x_new, x, n);
                Array.Copy(g_new, g, n);
            }

            return x;
        }

        /// <summary>
        /// Helper method to compute Hessian matrix
        /// </summary>
        private static Matrix ComputeHessian(Func<double[], double> f, double[] point)
        {
            return Differentiation.Hessian(f, point);
        }
    }
}