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

### BSE Calculation with Coupling Terms and Full Band Coverage
**Files**: 
- `bse_full_bands_coupling_specification.md` - Detailed specification
- `bse_full_bands_coupling_design.md` - Detailed design and implementation

**BSE結合項を含む全バンド計算の詳細仕様書および詳細設計書**

These documents extend the memory-efficient BSE algorithms to handle the coupling mode (BSEmod="coupling"), where both resonant and anti-resonant (coupling) terms are included in the BSE Hamiltonian, while considering all bands.

**Key Features:**
- Complete mathematical formulation for coupling Hamiltonian
- Pseudo-Hermitian biorthogonal Lanczos method
- On-the-fly kernel construction for both resonant and coupling blocks
- Memory-efficient implementation maintaining rigorous mathematical treatment
- No heuristic approximations or fallback methods
- Detailed implementation algorithms with pseudocode and Fortran examples

**Sections (Specification):**
1. Mathematical formulation of coupling BSE Hamiltonian
2. Memory bottleneck analysis (coupling vs. resonant-only)
3. On-the-fly kernel construction strategy
4. Modified Haydock method for pseudo-Hermitian systems
5. Biorthogonal Lanczos iteration
6. Symmetry exploitation (time-reversal, spatial, pseudo-Hermitian)
7. Numerical stability and convergence
8. Performance metrics and benchmarks
9. Validation and quality assurance
10. Parameter recommendations

**Sections (Design):**
1. Module architecture and data flow
2. Data structure definitions (Fortran)
3. Detailed algorithms with pseudocode
4. Kernel action implementation (resonant and coupling blocks)
5. Memory management strategy
6. Parallelization schemes (MPI/OpenMP/hybrid)
7. Optimization techniques (FFT, caching, BLAS)
8. Testing and benchmarking procedures
9. Implementation checklist

**Memory Improvement:** Similar ~10^7× reduction as resonant-only case (PB-scale → GB-scale)

**Target Audience:** Developers implementing coupling BSE calculations with full band consideration.

---

---

### PR#3 Implementation Documentation (実装済み機能のドキュメント)

PR#3では、上記の仕様と設計に基づいて、メモリ効率的な全バンドBSEソルバーが実装されました。以下のドキュメントで実装の詳細と使用方法を説明しています。

**Files**:
- `PR3_feature_details.md` - 実装された機能の詳細
- `PR3_usage_guide.md` - 機能の使用方法ガイド
- `PR3_tutorial_inputs.md` - チュートリアル用入力ファイル集
- `PR3_installation_guide.md` - インストール手順書
- `PR3_test_examples.md` - 動作確認テスト例題

**PR#3で実装されたコンポーネント:**

| ファイル | 説明 |
|---------|------|
| `src/bse/K_Haydock_full_bands.F` | 双直交Lanczosソルバー本体 |
| `src/bse/K_dot_product_metric.F` | 計量内積・再直交化ルーチン |
| `src/modules/mod_BS_OnTheFly_kernel.F` | オンザフライカーネルモジュール |
| `src/modules/mod_BS_solvers.F` | ソルバーフラグの追加 |
| `src/bse/K_solvers.F` | ソルバーディスパッチの追加 |

**Key Implementation Features:**
- Biorthogonal Lanczos iteration for pseudo-Hermitian coupling Hamiltonian
- Metric operator Θ = diag(I, -I) for BSE coupling case
- Numerical stability checks (breakdown detection, imaginary part tolerance)
- Support for optics, kerr, magnons, and dichroism calculations
- No heuristic approximations or fallback methods

**Memory Reduction Achieved:** ~10^7× (PB-scale → GB-scale)

---

For questions or contributions, please refer to the main Yambo documentation and support channels at https://www.yambo-code.eu/
