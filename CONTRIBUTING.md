# Contributing to ZiggyMath

Thank you for your interest in contributing to ZiggyMath! This document provides guidelines for contributing to this high-performance mathematics library.

## 🚀 Quick Start

1. **Fork** the repository on GitHub
2. **Clone** your fork locally
3. **Create** a feature branch from `main`
4. **Make** your changes
5. **Test** thoroughly
6. **Submit** a pull request

## 📋 Development Process

### 1. Setting Up Development Environment

```bash
# Clone your fork
git clone https://github.com/yourusername/ZiggyMath.git
cd ZiggyMath

# Install dependencies
dotnet restore

# Run tests to ensure everything works
dotnet test

# Run examples
cd src/Examples
dotnet run
```

### 2. Code Style Guidelines

#### C# Coding Standards
- Use **4 spaces** for indentation (no tabs)
- Follow **PascalCase** for class and method names
- Use **camelCase** for local variables and parameters
- Add **XML documentation** for all public APIs
- Use **async/await** for asynchronous operations
- Prefer **expression-bodied members** where appropriate

#### Example:
```csharp
/// <summary>
/// Compute the dot product of two vectors with SIMD optimization.
/// </summary>
/// <param name="left">First vector</param>
/// <param name="right">Second vector</param>
/// <returns>Dot product result</returns>
public static double DotProduct(this Vector left, Vector right)
    => VectorizedOperations.Sum(left.AsSpan(), right.AsSpan());
```

### 3. Testing Requirements

- **Unit tests** for all new functionality
- **Integration tests** for algorithm changes
- **Performance tests** for optimization changes
- **Edge case testing** for error conditions

```csharp
[TestMethod]
public void VectorAddition_ProducesCorrectResults()
{
    // Arrange
    using var context = new MathComputationContext(MathWorkload.SmallComputations);
    using var a = new Vector(new[] { 1.0, 2.0, 3.0 }, context.Allocator);
    using var b = new Vector(new[] { 4.0, 5.0, 6.0 }, context.Allocator);

    // Act
    using var result = a + b;

    // Assert
    Assert.AreEqual(5.0, result[0]);
    Assert.AreEqual(7.0, result[1]);
    Assert.AreEqual(9.0, result[2]);
}
```

### 4. Performance Requirements

- **Benchmark** new algorithms against existing implementations
- **Memory efficiency** must be maintained or improved
- **SIMD utilization** for compute-intensive operations
- **Allocation optimization** using ZiggyAlloc patterns

## 🔧 Types of Contributions

### 🐛 Bug Fixes
- **Reproduce** the issue first
- **Write a test** that demonstrates the fix
- **Implement** the minimal fix
- **Verify** no regressions

### ✨ New Features
- **Design discussion** via GitHub Issues
- **Implementation plan** with milestones
- **Comprehensive testing** including edge cases
- **Documentation updates** in API.md
- **Example code** in the Examples project

### 📚 Documentation
- **API documentation** for new public methods
- **Usage examples** for complex features
- **Performance guidance** for optimization
- **Migration guides** for breaking changes

### 🔍 Code Reviews
- **Clarity**: Code should be self-documenting
- **Performance**: Consider algorithmic complexity
- **Safety**: Memory safety and error handling
- **Testing**: Adequate test coverage

## 📝 Commit Guidelines

### Commit Message Format
```
<type>(<scope>): <description>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Test additions/modifications
- `perf`: Performance improvements
- `refactor`: Code restructuring
- `ci`: CI/CD changes

**Examples:**
```
feat(fft): Add 2D FFT implementation for image processing

- Implement 2D Cooley-Tukey algorithm
- Add comprehensive test coverage
- Update API documentation with examples

Closes #123
```

## 🚨 Pull Request Process

### Before Submitting

1. **Run all tests** locally
2. **Update documentation** if needed
3. **Add examples** for new features
4. **Check code style** consistency
5. **Performance validation** for algorithmic changes

### PR Template

```markdown
## Description
Brief description of changes

## Motivation
Why this change is needed

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests pass
- [ ] Performance benchmarks included
- [ ] Edge cases covered

## Documentation
- [ ] API.md updated
- [ ] Examples added
- [ ] README changes if needed

## Performance Impact
- [ ] Performance improved
- [ ] No performance regression
- [ ] Benchmarks included

## Related Issues
Closes #123
```

## 🏗️ Architecture Guidelines

### Memory Management
- **Always use `using` statements** for ZiggyMath objects
- **Prefer `MathComputationContext`** for workload optimization
- **Choose appropriate allocators** based on usage patterns
- **Avoid unnecessary allocations** in hot paths

### Algorithm Design
- **Consider cache locality** for large data structures
- **Use SIMD operations** for element-wise computations
- **Implement proper error handling** for edge cases
- **Document algorithmic complexity** in comments

### API Design
- **Consistent naming** across the library
- **Intuitive parameter order** (input, output, options)
- **Comprehensive XML documentation**
- **Example usage** in doc comments

## 📞 Getting Help

### Communication Channels
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: Design discussions and Q&A
- **Email**: For sensitive security issues

### Development Support
- **Documentation**: Check docs/ directory first
- **Examples**: Review src/Examples/ for usage patterns
- **Tests**: Study existing tests for implementation patterns

## 🎯 Quality Standards

### Code Quality
- **Zero warnings** in Release builds
- **Comprehensive error handling**
- **Input validation** for public APIs
- **Memory safety** with ZiggyAlloc patterns

### Performance Standards
- **Benchmark comparisons** for optimizations
- **Memory efficiency** tracking
- **SIMD utilization** where applicable
- **Scalability testing** for large inputs

### Testing Standards
- **Unit test coverage** > 90%
- **Integration tests** for complex workflows
- **Performance regression tests**
- **Cross-platform compatibility**

## 🚀 Recognition

Contributors will be recognized in:
- **README.md** contributors section
- **Release notes** for significant contributions
- **GitHub repository** contributor statistics

---

Thank you for contributing to ZiggyMath! Your efforts help make high-performance mathematical computing accessible to developers worldwide. 🎉