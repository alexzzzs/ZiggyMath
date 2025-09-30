using System;
using ZiggyMath.Core;

namespace ZiggyMath.Tests.UnitTests
{
    /// <summary>
    /// Unit tests for Vector class
    /// </summary>
    public static class VectorTests
    {
        public static int RunAllTests()
        {
            int passed = 0;
            int total = 0;

            Console.WriteLine("Running Vector tests...");

            total++; if (TestConstructor_WithSize_CreatesValidVector()) { passed++; } else { Console.WriteLine("FAILED: Constructor_WithSize_CreatesValidVector"); }
            total++; if (TestConstructor_WithData_CreatesValidVector()) { passed++; } else { Console.WriteLine("FAILED: Constructor_WithData_CreatesValidVector"); }
            total++; if (TestIndexer_GetAndSet_WorksCorrectly()) { passed++; } else { Console.WriteLine("FAILED: Indexer_GetAndSet_WorksCorrectly"); }
            total++; if (TestIndexer_OutOfBounds_ThrowsException()) { passed++; } else { Console.WriteLine("FAILED: Indexer_OutOfBounds_ThrowsException"); }
            total++; if (TestAdd_TwoVectors_ReturnsCorrectSum()) { passed++; } else { Console.WriteLine("FAILED: Add_TwoVectors_ReturnsCorrectSum"); }
            total++; if (TestSubtract_TwoVectors_ReturnsCorrectDifference()) { passed++; } else { Console.WriteLine("FAILED: Subtract_TwoVectors_ReturnsCorrectDifference"); }
            total++; if (TestMultiply_VectorAndScalar_ReturnsCorrectResult()) { passed++; } else { Console.WriteLine("FAILED: Multiply_VectorAndScalar_ReturnsCorrectResult"); }
            total++; if (TestDotProduct_TwoVectors_ReturnsCorrectResult()) { passed++; } else { Console.WriteLine("FAILED: DotProduct_TwoVectors_ReturnsCorrectResult"); }
            total++; if (TestMagnitude_UnitVector_ReturnsOne()) { passed++; } else { Console.WriteLine("FAILED: Magnitude_UnitVector_ReturnsOne"); }
            total++; if (TestNormalized_Vector_ReturnsUnitVector()) { passed++; } else { Console.WriteLine("FAILED: Normalized_Vector_ReturnsUnitVector"); }

            Console.WriteLine($"Vector tests: {passed}/{total} passed");
            return passed == total ? 0 : 1;
        }

        private static bool TestConstructor_WithSize_CreatesValidVector()
        {
            try
            {
                using var vector = new Vector(5);
                return vector.Length == 5 && vector.IsValid && vector.SizeInBytes == 40;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestConstructor_WithData_CreatesValidVector()
        {
            try
            {
                double[] data = { 1.0, 2.0, 3.0, 4.0, 5.0 };
                using var vector = new Vector(data);

                if (vector.Length != 5 || !vector.IsValid)
                    return false;

                for (int i = 0; i < 5; i++)
                {
                    if (vector[i] != data[i])
                        return false;
                }

                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestIndexer_GetAndSet_WorksCorrectly()
        {
            try
            {
                using var vector = new Vector(3);
                vector[0] = 10.0;
                vector[1] = 20.0;
                vector[2] = 30.0;

                return vector[0] == 10.0 && vector[1] == 20.0 && vector[2] == 30.0;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestIndexer_OutOfBounds_ThrowsException()
        {
            try
            {
                using var vector = new Vector(3);

                // These should throw exceptions
                var temp1 = vector[-1];
                var temp2 = vector[3];
                var temp3 = vector[10];

                return false; // Should not reach here
            }
            catch (IndexOutOfRangeException)
            {
                return true; // Expected exception
            }
            catch (Exception)
            {
                return false; // Wrong exception type
            }
        }

        private static bool TestAdd_TwoVectors_ReturnsCorrectSum()
        {
            try
            {
                using var v1 = new Vector(new double[] { 1, 2, 3 });
                using var v2 = new Vector(new double[] { 4, 5, 6 });
                using var result = Vector.Add(v1, v2);

                return result[0] == 5.0 && result[1] == 7.0 && result[2] == 9.0;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestSubtract_TwoVectors_ReturnsCorrectDifference()
        {
            try
            {
                using var v1 = new Vector(new double[] { 5, 7, 9 });
                using var v2 = new Vector(new double[] { 1, 2, 3 });
                using var result = Vector.Subtract(v1, v2);

                return result[0] == 4.0 && result[1] == 5.0 && result[2] == 6.0;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestMultiply_VectorAndScalar_ReturnsCorrectResult()
        {
            try
            {
                using var vector = new Vector(new double[] { 1, 2, 3 });
                using var result = Vector.Multiply(vector, 3.0);

                return result[0] == 3.0 && result[1] == 6.0 && result[2] == 9.0;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestDotProduct_TwoVectors_ReturnsCorrectResult()
        {
            try
            {
                using var v1 = new Vector(new double[] { 1, 2, 3 });
                using var v2 = new Vector(new double[] { 4, 5, 6 });
                double result = v1.DotProduct(v2);

                return Math.Abs(result - 32.0) < 1e-10; // 1*4 + 2*5 + 3*6 = 32
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestMagnitude_UnitVector_ReturnsOne()
        {
            try
            {
                using var vector = new Vector(new double[] { 1, 0, 0 });
                double magnitude = vector.Magnitude();

                return Math.Abs(magnitude - 1.0) < 1e-10;
            }
            catch (Exception)
            {
                return false;
            }
        }

        private static bool TestNormalized_Vector_ReturnsUnitVector()
        {
            try
            {
                using var vector = new Vector(new double[] { 3, 4 });
                using var normalized = vector.Normalized();

                double magnitude = normalized.Magnitude();
                bool correctMagnitude = Math.Abs(magnitude - 1.0) < 1e-10;
                bool correctValues = Math.Abs(normalized[0] - 0.6) < 1e-10 && Math.Abs(normalized[1] - 0.8) < 1e-10;

                return correctMagnitude && correctValues;
            }
            catch (Exception)
            {
                return false;
            }
        }
    }
}