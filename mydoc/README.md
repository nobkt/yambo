# Yambo Documentation - mydoc

This directory contains specialized documentation for advanced computational techniques in Yambo.

## Contents

### Memory-Efficient BSE Calculation with Full Band Coverage
**File**: `memory_efficient_bse_full_bands.md`

**メモリ効率的な全バンドBSE計算の詳細理論と実装アルゴリズム**

This comprehensive technical document provides a rigorous mathematical framework and implementation algorithms for performing Bethe-Salpeter Equation (BSE) calculations with full band consideration, while maintaining memory requirements within practical limits (< 2TB).

**Key Features:**
- Complete mathematical formulation of BSE theory
- Memory-efficient on-the-fly kernel construction algorithms
- Haydock iterative diagonalization methods
- Advanced techniques including:
  - Band chunking
  - k-point partitioning
  - Tensor decomposition
  - Adaptive band selection
- Parallel implementation strategies
- Numerical benchmarks and validation
- No heuristic approximations or fallback methods

**Sections:**
1. Problem formulation and memory bottleneck analysis
2. Theoretical background (BSE, exchange/correlation kernels)
3. On-the-fly BSE kernel construction
4. Memory reduction techniques
5. Advanced algorithms (Chebyshev expansion, tensor decomposition)
6. Parallelization strategies
7. Numerical stability and convergence
8. Computational complexity analysis
9. Implementation details
10. Benchmarks and performance evaluation
11. Validation and verification
12. Extensions and applications
13. Implementation checklist and recommended settings
14. Summary and future outlook
15. Appendices (symbols, flowcharts, sample inputs)

**Memory Reduction:** ~1000x improvement (from PB-scale to GB-scale for typical systems)

**Target Audience:** Advanced users and developers working on excited-state calculations requiring full band consideration.

---

For questions or contributions, please refer to the main Yambo documentation and support channels at https://www.yambo-code.eu/
