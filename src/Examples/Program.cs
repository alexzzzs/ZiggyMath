using System;
using Examples;

namespace Examples
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("🚀 ZiggyMath with ZiggyAlloc Integration Demo");
            Console.WriteLine("==============================================\n");

            try
            {
                // Run consolidated examples organized by functionality
                ZiggyMathExamples.RunAllExamples();

                Console.WriteLine("\n" + "=".PadRight(50, '='));
                Console.WriteLine("🎯 Performance Benchmarks");
                Console.WriteLine("=".PadRight(50, '=') + "\n");

                // Run performance benchmarks (including new allocator benchmarks)
                // PerformanceBenchmarks.RunBenchmarks();

                Console.WriteLine("\n" + "=".PadRight(50, '='));
                Console.WriteLine("🔍 Debug Integration Examples");
                Console.WriteLine("=".PadRight(50, '=') + "\n");

                // Run debug integration examples
                // DebugIntegrationExamples.RunDebugExamples();

                Console.WriteLine("\n" + "=".PadRight(50, '='));
                Console.WriteLine("✅ All examples completed successfully!");
                Console.WriteLine("=".PadRight(50, '='));
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\n❌ Error running examples: {ex.Message}");
                Console.WriteLine($"Stack trace: {ex.StackTrace}");
            }

            Console.WriteLine("\nPress any key to exit...");
            Console.ReadKey();
        }
    }
}