using System;
using ZiggyMath.Core;
using ZiggyMath.LinearAlgebra;

namespace ZiggyMath.Statistics
{
    /// <summary>
    /// Regression analysis
    /// </summary>
    public static class Regression
    {
        /// <summary>
        /// Linear regression: y = mx + b
        /// </summary>
        public static (double slope, double intercept, double rSquared) LinearRegression(
            ReadOnlySpan<double> x, ReadOnlySpan<double> y)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (x.Length < 2)
                throw new ArgumentException("At least 2 data points required");

            int n = x.Length;
            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;

            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumXX += x[i] * x[i];
                sumYY += y[i] * y[i];
            }

            double meanX = sumX / n;
            double meanY = sumY / n;

            double slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
            double intercept = meanY - slope * meanX;

            // Calculate R-squared
            double ssRes = 0, ssTot = 0;
            for (int i = 0; i < n; i++)
            {
                double predicted = slope * x[i] + intercept;
                ssRes += (y[i] - predicted) * (y[i] - predicted);
                ssTot += (y[i] - meanY) * (y[i] - meanY);
            }

            double rSquared = ssTot != 0 ? 1 - ssRes / ssTot : 0;

            return (slope, intercept, rSquared);
        }

        /// <summary>
        /// Polynomial regression
        /// </summary>
        public static double[] PolynomialRegression(ReadOnlySpan<double> x, ReadOnlySpan<double> y, int degree)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (degree < 1)
                throw new ArgumentException("Degree must be at least 1");
            if (x.Length <= degree)
                throw new ArgumentException("Not enough data points for the given degree");

            int n = x.Length;
            int m = degree + 1;

            // Create Vandermonde matrix
            Matrix X = new Matrix(n, m);
            Vector Y = new Vector(n);

            for (int i = 0; i < n; i++)
            {
                Y[i] = y[i];
                double xi = 1.0;

                for (int j = 0; j < m; j++)
                {
                    X[i, j] = xi;
                    xi *= x[i];
                }
            }

            // Solve normal equations: (X^T * X) * coefficients = X^T * Y
            Matrix Xt = X.Transpose();
            Matrix XtX = Matrix.Multiply(Xt, X);
            Vector XtY = Matrix.Multiply(Xt, Y);

            Vector coefficients;
            GaussianElimination.Solve(XtX, XtY, out coefficients);

            double[] result = new double[m];
            for (int i = 0; i < m; i++)
            {
                result[i] = coefficients[i];
            }

            return result;
        }

        /// <summary>
        /// Multiple linear regression: Y = X * B
        /// </summary>
        public static (double[] coefficients, double rSquared) MultipleLinearRegression(
            Matrix X, Vector y)
        {
            if (X.Rows != y.Length)
                throw new ArgumentException("Matrix rows must equal vector length");

            int n = X.Rows;
            int p = X.Columns; // Number of predictors

            if (n <= p)
                throw new ArgumentException("Not enough observations for the number of predictors");

            // Solve normal equations: (X^T * X) * B = X^T * y
            Matrix Xt = X.Transpose();
            Matrix XtX = Matrix.Multiply(Xt, X);
            Vector XtY = Matrix.Multiply(Xt, y);

            Vector coefficients;
            GaussianElimination.Solve(XtX, XtY, out coefficients);

            // Calculate R-squared
            Vector yPred = Matrix.Multiply(X, coefficients);
            double ssRes = 0, ssTot = 0;
            double yMean = DescriptiveStats.Mean(y.AsReadOnlySpan());

            for (int i = 0; i < n; i++)
            {
                ssRes += (y[i] - yPred[i]) * (y[i] - yPred[i]);
                ssTot += (y[i] - yMean) * (y[i] - yMean);
            }

            double rSquared = ssTot != 0 ? 1 - ssRes / ssTot : 0;

            double[] result = new double[coefficients.Length];
            for (int i = 0; i < coefficients.Length; i++)
            {
                result[i] = coefficients[i];
            }

            return (result, rSquared);
        }

        /// <summary>
        /// Nonlinear regression using Gauss-Newton method
        /// </summary>
        public static double[] NonlinearRegression(Func<double[], double[], double> model,
                                                double[] x, double[] y, double[] initialGuess)
        {
            if (x == null || y == null || x.Length != y.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (initialGuess == null || initialGuess.Length == 0)
                throw new ArgumentException("Initial guess cannot be null or empty");

            const int maxIterations = 100;
            const double tolerance = 1e-10;

            double[] parameters = new double[initialGuess.Length];
            Array.Copy(initialGuess, parameters, initialGuess.Length);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // Compute residuals and Jacobian
                double[] residuals = new double[x.Length];
                double[,] jacobian = new double[x.Length, parameters.Length];

                for (int i = 0; i < x.Length; i++)
                {
                    double predicted = model(parameters, new double[] { x[i] });
                    residuals[i] = y[i] - predicted;

                    // Numerical differentiation for Jacobian
                    double h = 1e-8;
                    for (int j = 0; j < parameters.Length; j++)
                    {
                        double[] paramsPlus = new double[parameters.Length];
                        double[] paramsMinus = new double[parameters.Length];
                        Array.Copy(parameters, paramsPlus, parameters.Length);
                        Array.Copy(parameters, paramsMinus, parameters.Length);

                        paramsPlus[j] += h;
                        paramsMinus[j] -= h;

                        double fPlus = model(paramsPlus, new double[] { x[i] });
                        double fMinus = model(paramsMinus, new double[] { x[i] });

                        jacobian[i, j] = (fPlus - fMinus) / (2 * h);
                    }
                }

                // Solve normal equations: J^T * J * delta = J^T * residuals
                Matrix J = new Matrix(jacobian);
                Matrix Jt = J.Transpose();
                Matrix JtJ = Matrix.Multiply(Jt, J);
                Vector JtR = Matrix.Multiply(Jt, new Vector(residuals));

                Vector delta;
                GaussianElimination.Solve(JtJ, JtR, out delta);

                // Update parameters
                for (int j = 0; j < parameters.Length; j++)
                {
                    parameters[j] += delta[j];
                }

                // Check convergence
                double maxDelta = 0;
                for (int j = 0; j < parameters.Length; j++)
                {
                    if (Math.Abs(delta[j]) > maxDelta)
                        maxDelta = Math.Abs(delta[j]);
                }

                if (maxDelta < tolerance)
                    break;
            }

            return parameters;
        }

        /// <summary>
        /// Exponential regression: y = a * exp(b * x)
        /// </summary>
        public static (double a, double b, double rSquared) ExponentialRegression(
            ReadOnlySpan<double> x, ReadOnlySpan<double> y)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (x.Length < 2)
                throw new ArgumentException("At least 2 data points required");

            // Transform to linear: ln(y) = ln(a) + b * x
            double[] logY = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                if (y[i] <= 0)
                    throw new ArgumentException("All Y values must be positive for exponential regression");
                logY[i] = Math.Log(y[i]);
            }

            var (slope, intercept, rSquared) = LinearRegression(x, logY);

            double a = Math.Exp(intercept);
            double b = slope;

            return (a, b, rSquared);
        }

        /// <summary>
        /// Power regression: y = a * x^b
        /// </summary>
        public static (double a, double b, double rSquared) PowerRegression(
            ReadOnlySpan<double> x, ReadOnlySpan<double> y)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (x.Length < 2)
                throw new ArgumentException("At least 2 data points required");

            // Transform to linear: ln(y) = ln(a) + b * ln(x)
            double[] logX = new double[x.Length];
            double[] logY = new double[y.Length];

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] <= 0 || y[i] <= 0)
                    throw new ArgumentException("All X and Y values must be positive for power regression");
                logX[i] = Math.Log(x[i]);
                logY[i] = Math.Log(y[i]);
            }

            var (slope, intercept, rSquared) = LinearRegression(logX, logY);

            double a = Math.Exp(intercept);
            double b = slope;

            return (a, b, rSquared);
        }

        /// <summary>
        /// Logistic regression (simplified 2-class)
        /// </summary>
        public static (double[] coefficients, double logLikelihood) LogisticRegression(
            Matrix X, double[] y)
        {
            if (X.Rows != y.Length)
                throw new ArgumentException("Matrix rows must equal Y array length");
            if (X.Rows < 2)
                throw new ArgumentException("At least 2 observations required");

            int n = X.Rows;
            int p = X.Columns + 1; // +1 for intercept

            // Add intercept column
            Matrix XWithIntercept = new Matrix(n, p);
            for (int i = 0; i < n; i++)
            {
                XWithIntercept[i, 0] = 1.0; // Intercept
                for (int j = 1; j < p; j++)
                {
                    XWithIntercept[i, j] = X[i, j - 1];
                }
            }

            // Initial coefficients (zeros)
            double[] coefficients = new double[p];

            const int maxIterations = 100;
            const double tolerance = 1e-10;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double logLikelihood = 0.0;
                double[] gradient = new double[p];
                double[,] hessian = new double[p, p];

                for (int i = 0; i < n; i++)
                {
                    // Linear combination
                    double eta = 0.0;
                    for (int j = 0; j < p; j++)
                    {
                        eta += coefficients[j] * XWithIntercept[i, j];
                    }

                    double p_i = 1.0 / (1.0 + Math.Exp(-eta));
                    double residual = y[i] - p_i;

                    logLikelihood += y[i] * Math.Log(p_i + 1e-10) + (1 - y[i]) * Math.Log(1 - p_i + 1e-10);

                    // Gradient
                    for (int j = 0; j < p; j++)
                    {
                        gradient[j] += residual * XWithIntercept[i, j];
                    }

                    // Hessian (diagonal elements only for simplicity)
                    for (int j = 0; j < p; j++)
                    {
                        hessian[j, j] += p_i * (1 - p_i) * XWithIntercept[i, j] * XWithIntercept[i, j];
                    }
                }

                // Newton step
                double[] delta = new double[p];
                for (int j = 0; j < p; j++)
                {
                    if (hessian[j, j] != 0)
                    {
                        delta[j] = gradient[j] / hessian[j, j];
                        coefficients[j] += delta[j];
                    }
                }

                // Check convergence
                double maxDelta = 0;
                for (int j = 0; j < p; j++)
                {
                    if (Math.Abs(delta[j]) > maxDelta)
                        maxDelta = Math.Abs(delta[j]);
                }

                if (maxDelta < tolerance)
                    break;
            }

            return (coefficients, 0.0); // TODO: Return actual log likelihood
        }

        /// <summary>
        /// Evaluate polynomial at given x value
        /// </summary>
        public static double EvaluatePolynomial(double[] coefficients, double x)
        {
            double result = 0.0;
            double xPower = 1.0;

            for (int i = 0; i < coefficients.Length; i++)
            {
                result += coefficients[i] * xPower;
                xPower *= x;
            }

            return result;
        }

        /// <summary>
        /// Compute coefficient of determination (R-squared) for any regression
        /// </summary>
        public static double RSquared(ReadOnlySpan<double> observed, ReadOnlySpan<double> predicted)
        {
            if (observed.Length != predicted.Length)
                throw new ArgumentException("Observed and predicted arrays must have the same length");

            double meanObserved = DescriptiveStats.Mean(observed);
            double ssRes = 0, ssTot = 0;

            for (int i = 0; i < observed.Length; i++)
            {
                ssRes += (observed[i] - predicted[i]) * (observed[i] - predicted[i]);
                ssTot += (observed[i] - meanObserved) * (observed[i] - meanObserved);
            }

            return ssTot != 0 ? 1 - ssRes / ssTot : 0;
        }
    }
}