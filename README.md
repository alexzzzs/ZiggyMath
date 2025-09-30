# ZiggyMath

High-performance mathematics library with SIMD acceleration and zero GC pressure.

## Quick Start

```csharp
using ZiggyMath;
using ZiggyMath.Core;

// Create vectors
using var v1 = new Vector(1000);
using var v2 = new Vector(1000);

// Initialize and compute
for (int i = 0; i < 1000; i++)
{
    v1[i] = Math.Sin(i * 0.01);
    v2[i] = Math.Cos(i * 0.01);
}

double dotProduct = v1.DotProduct(v2);

// Matrix operations
using var matrix = new Matrix(100, 100);
using var result = Matrix.Multiply(matrix, matrix.Transpose());
```

## Performance

ZiggyMath delivers significant performance improvements:

- **5-29x faster memory operations** through SIMD
- **Zero GC pressure** for real-time applications
- **Hardware acceleration** with AVX2 support
- **Optimized memory access patterns** for cache efficiency

## Requirements

- .NET 8.0 or later
- x64 or ARM64 processor
- Unsafe code enabled (configured automatically)

## Installation

```bash
dotnet add package ZiggyMath
```

Or add to your .csproj:

```xml
<PackageReference Include="ZiggyMath" Version="1.0.0" />
```

## Dependency Injection Setup

```csharp
public void ConfigureServices(IServiceCollection services)
{
    // Register default allocator
    services.AddSingleton<IUnmanagedMemoryAllocator>(ZiggyMath.DefaultAllocator);

    // Or use custom allocator for specific workloads
    services.AddSingleton<IUnmanagedMemoryAllocator>(provider =>
        new HybridAllocator(new SystemMemoryAllocator()));
}
```

## Configuration Options

```csharp
public static class MathConfiguration
{
    public static IUnmanagedMemoryAllocator GetAllocatorForWorkload(MathWorkload workload)
    {
        return workload switch
        {
            MathWorkload.LinearAlgebra => new SlabAllocator(new SystemMemoryAllocator()),
            MathWorkload.SignalProcessing => new LargeBlockAllocator(new SystemMemoryAllocator()),
            MathWorkload.SmallComputations => new UnmanagedMemoryPool(new SystemMemoryAllocator()),
            _ => ZiggyMath.DefaultAllocator
        };
    }
}

public enum MathWorkload
{
    LinearAlgebra,
    SignalProcessing,
    SmallComputations,
    GeneralPurpose
}
```

## Acknowledgments

ZiggyMath is built on top of [ZiggyAlloc](https://github.com/alexzzzs/ziggyalloc), a high-performance memory management library that provides the foundation for ZiggyMath's zero-GC pressure design and advanced allocation strategies.

## Documentation

For detailed API documentation, see the [docs](docs/) directory.

## License

MIT License - see LICENSE file for details.