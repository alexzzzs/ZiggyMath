using System;

namespace ZiggyMath.Tests
{
    class Program
    {
        static void Main(string[] args)
        {
            // Run comprehensive test suite
            int result = TestRunner.RunAllTests();
            Environment.Exit(result);
        }
    }
}