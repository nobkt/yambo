# BSE全バンド計算における結合項を含む詳細設計書
## Detailed Design for Full-Band BSE Calculation with Coupling Terms

---

## 1. アーキテクチャ概要

### 1.1 モジュール構成

本設計では、以下の階層的モジュール構成を採用する：

```
BSE_Coupling_FullBands_Driver
├── Initialization_Module
│   ├── Parameter_Setup
│   ├── Memory_Allocation
│   └── Symmetry_Analysis
├── WaveFunction_Manager
│   ├── WF_Cache
│   ├── WF_Loader
│   └── WF_Prefetch
├── Kernel_Action_Module
│   ├── Resonant_Block
│   │   ├── Exchange_Resonant
│   │   └── Correlation_Resonant
│   └── Coupling_Block
│       ├── Exchange_Coupling
│       └── Correlation_Coupling
├── Lanczos_Solver
│   ├── Biorthogonal_Iteration
│   ├── Reorthogonalization
│   └── Convergence_Check
├── Response_Function_Calculator
│   ├── Tridiagonal_Solver
│   ├── Green_Function
│   └── Dielectric_Function
└── Output_Module
    ├── Spectrum_Writer
    └── Diagnostics
```

### 1.2 データフロー

```
[入力] → [初期化] → [Lanczos反復]
                         ↓
                    [カーネル作用]
                    ├─ [波動関数ロード]
                    ├─ [行列要素計算]
                    └─ [カーネル適用]
                         ↓
                    [ベクトル更新]
                         ↓
                    [収束判定] ─No→ [次の反復]
                         ↓ Yes
                    [応答関数計算]
                         ↓
                    [出力]
```

---

## 2. データ構造定義

(詳細は仕様書参照)

本設計書では、実装レベルの詳細を提供します。

---

## 3. 実装アルゴリズム

### 3.1 メインドライバのフローチャート

```
START
  ↓
[パラメータ読み込み]
  ↓
[メモリ割り当て]
  ↓
[波動関数キャッシュ初期化]
  ↓
[初期ベクトル設定]
  ↓
[反復カウンタ n = 0]
  ↓
[ループ開始: n < n_max AND NOT converged]
  ↓
  [カーネル作用計算]
  [H |ψ_n⟩ および H† |φ_n⟩]
    ↓
    ├─ [共鳴ブロック作用]
    │   ├─ 対角項
    │   ├─ 交換項（共鳴）
    │   └─ 相関項（共鳴）
    └─ [結合ブロック作用]
        ├─ 交換項（結合）
        └─ 相関項（結合）
  ↓
  [Lanczos係数計算]
  [a_n, b_{n+1}]
  ↓
  [ベクトル更新]
  [|ψ_{n+1}⟩, |φ_{n+1}⟩]
  ↓
  [双直交性チェック]
  ├─ OK → 次へ
  └─ NG → [再直交化]
  ↓
  [収束判定]
  ├─ 収束 → [ループ終了]
  └─ 非収束 → [n = n+1, ループ継続]
  ↓
[三重対角行列対角化]
  ↓
[誘電関数計算]
  ↓
[結果出力]
  ↓
[メモリ解放]
  ↓
END
```

---

## 4. 疑似コード

### 4.1 メインルーチン

```
PROCEDURE BSE_Coupling_FullBands_Main()
  
  // 1. 初期化
  CALL Read_Input_Parameters(params)
  CALL Allocate_Memory(params)
  CALL Initialize_WF_Cache(wf_cache, params)
  CALL Initialize_Kernel_Workspace(kws, params)
  CALL Setup_MPI_Distribution(params)
  
  // 2. 初期ベクトル
  CALL Setup_Initial_Vectors(psi_0, phi_0, params)
  psi_n = psi_0
  phi_n = phi_0
  n = 0
  converged = FALSE
  
  // 3. Lanczos反復
  WHILE (n < n_max AND NOT converged) DO
    
    // 3.1 カーネル作用
    CALL Apply_Hamiltonian_Coupling(psi_n, H_psi, kws, wf_cache)
    CALL Apply_Hamiltonian_Adjoint_Coupling(phi_n, H_phi, kws, wf_cache)
    
    // 3.2 Lanczos係数
    a[n] = Inner_Product_Metric(phi_n, H_psi)
    
    // 3.3 直交化
    psi_tilde = H_psi - a[n] * psi_n
    phi_tilde = H_phi - CONJ(a[n]) * phi_n
    
    IF (n > 0) THEN
      psi_tilde = psi_tilde - b[n] * psi_{n-1}
      phi_tilde = phi_tilde - b[n] * phi_{n-1}
    END IF
    
    // 3.4 規格化
    b_squared = Inner_Product_Metric(phi_tilde, psi_tilde)
    
    IF (REAL(b_squared) <= 0) THEN
      PRINT "Lanczos breakdown"
      converged = TRUE
      CONTINUE
    END IF
    
    b[n+1] = SQRT(ABS(REAL(b_squared)))
    psi_{n+1} = psi_tilde / b[n+1]
    phi_{n+1} = phi_tilde / b[n+1]
    
    // 3.5 再直交化（必要に応じて）
    ortho_error = Check_Biorthogonality(psi_{n+1}, phi_{n+1})
    IF (ortho_error > threshold) THEN
      CALL Reorthogonalize(psi_{n+1}, phi_{n+1})
    END IF
    
    // 3.6 収束チェック
    IF (n >= n_min) THEN
      CALL Calculate_Spectrum_Partial(a, b, n, spectrum)
      delta = Spectrum_Change(spectrum_current, spectrum_previous)
      IF (delta < tolerance) THEN
        converged = TRUE
      END IF
    END IF
    
    // 更新
    psi_{n-1} = psi_n
    phi_{n-1} = phi_n
    psi_n = psi_{n+1}
    phi_n = phi_{n+1}
    n = n + 1
    
  END WHILE
  
  // 4. 最終計算
  CALL Diagonalize_Tridiagonal(a, b, n, eigenvalues, eigenvectors)
  CALL Calculate_Dielectric_Function(eigenvalues, eigenvectors, spectrum)
  
  // 5. 出力
  CALL Write_Spectrum(spectrum, params)
  CALL Write_Diagnostics(n, convergence_history)
  
  // 6. クリーンアップ
  CALL Deallocate_Memory()
  CALL Finalize_MPI()
  
END PROCEDURE
```

### 4.2 カーネル作用（結合ハミルトニアン）

```
PROCEDURE Apply_Hamiltonian_Coupling(psi_in, H_psi_out, kws, wf_cache)
  
  INPUT:
    psi_in = {psi_plus, psi_minus}  // 結合ベクトル
    kws: Kernel workspace
    wf_cache: Wavefunction cache
  
  OUTPUT:
    H_psi_out = {H_psi_plus, H_psi_minus}
  
  // 初期化
  H_psi_plus = 0
  H_psi_minus = 0
  
  // 結合ハミルトニアンの作用:
  // [H^res   C    ] [ψ+]   [H^res ψ+ + C ψ-    ]
  // [-C*   -H^res*] [ψ-] = [-C* ψ+ - H^res* ψ-]
  
  // (1) H^res ψ+ の計算
  CALL Apply_Resonant_Block(psi_in.psi_plus, temp1, kws, wf_cache, conjugate=FALSE)
  H_psi_plus = H_psi_plus + temp1
  
  // (2) C ψ- の計算
  CALL Apply_Coupling_Block(psi_in.psi_minus, temp2, kws, wf_cache, conjugate=FALSE)
  H_psi_plus = H_psi_plus + temp2
  
  // (3) -C* ψ+ の計算
  CALL Apply_Coupling_Block(psi_in.psi_plus, temp3, kws, wf_cache, conjugate=TRUE)
  H_psi_minus = H_psi_minus - temp3
  
  // (4) -H^res* ψ- の計算
  CALL Apply_Resonant_Block(psi_in.psi_minus, temp4, kws, wf_cache, conjugate=TRUE)
  H_psi_minus = H_psi_minus - temp4
  
END PROCEDURE
```

### 4.3 共鳴ブロックの作用

```
PROCEDURE Apply_Resonant_Block(psi, H_psi, kws, wf_cache, conjugate)
  
  H_psi = 0
  
  // (A) 対角項: (E_c - E_v) ψ
  FOR ik = 1 TO n_k DO
    FOR ic = 1 TO n_c DO
      FOR iv = 1 TO n_v DO
        E_diff = E_conduction[ic,ik] - E_valence[iv,ik]
        IF (conjugate) THEN
          H_psi[iv,ic,ik] += CONJ(E_diff * psi[iv,ic,ik])
        ELSE
          H_psi[iv,ic,ik] += E_diff * psi[iv,ic,ik]
        END IF
      END FOR
    END FOR
  END FOR
  
  // (B) 交換項
  CALL Apply_Exchange_Resonant(psi, H_psi, kws, wf_cache, conjugate)
  
  // (C) 相関項
  CALL Apply_Correlation_Resonant(psi, H_psi, kws, wf_cache, conjugate)
  
END PROCEDURE
```

### 4.4 交換カーネル（共鳴）の詳細

```
PROCEDURE Apply_Exchange_Resonant(psi, H_psi, kws, wf_cache, conjugate)
  
  // k'ループ（MPI並列）
  FOR ikp = my_k_start TO my_k_end DO
    
    // k'の波動関数をロード
    CALL Load_Wavefunctions(ikp, all_bands, wf_cache)
    
    // kループ（OpenMP並列）
    PARALLEL FOR ik = 1 TO n_k DO
      
      // kの波動関数をロード
      CALL Load_Wavefunctions(ik, all_bands, wf_cache)
      
      // q = k - k'
      q = k_points[ik] - k_points[ikp]
      
      // バンドループ
      FOR iv = 1 TO n_v DO
        FOR ic = 1 TO n_c DO
          FOR ivp = 1 TO n_v DO
            FOR icp = 1 TO n_c DO
              
              K_ex = 0
              
              // G-vectorループ
              FOR iG = 1 TO n_G_exch DO
                
                q_G = q + G_vectors[iG]
                v_coul = 4 * PI / |q_G|^2
                
                // 遷移行列要素
                M_vvp = Matrix_Element(iv, ivp, ik, ikp, iG, wf_cache)
                M_ccp = Matrix_Element(ic, icp, ik, ikp, iG, wf_cache)
                
                IF (conjugate) THEN
                  K_ex += -2 * v_coul * M_vvp * CONJ(M_ccp)
                ELSE
                  K_ex += -2 * v_coul * M_vvp * CONJ(M_ccp)
                END IF
                
              END FOR
              
              K_ex = K_ex / Volume
              
              // カーネル作用
              IF (conjugate) THEN
                H_psi[iv,ic,ik] += CONJ(K_ex * psi[ivp,icp,ikp])
              ELSE
                H_psi[iv,ic,ik] += K_ex * psi[ivp,icp,ikp]
              END IF
              
            END FOR
          END FOR
        END FOR
      END FOR
      
    END PARALLEL FOR
    
  END FOR
  
  // MPI Allreduce
  CALL MPI_Allreduce_Sum(H_psi)
  
END PROCEDURE
```

### 4.5 結合ブロックの作用

```
PROCEDURE Apply_Coupling_Block(psi, H_psi, kws, wf_cache, conjugate)
  
  // 結合ブロックには対角項なし
  
  // (A) 交換項（結合）
  CALL Apply_Exchange_Coupling(psi, H_psi, kws, wf_cache, conjugate)
  
  // (B) 相関項（結合）
  CALL Apply_Correlation_Coupling(psi, H_psi, kws, wf_cache, conjugate)
  
END PROCEDURE
```

### 4.6 交換カーネル（結合）の詳細

```
PROCEDURE Apply_Exchange_Coupling(psi, H_psi, kws, wf_cache, conjugate)
  
  // 共鳴との違い: M_{vc'} と M_{v'c} の組み合わせ
  
  FOR ikp = my_k_start TO my_k_end DO
    CALL Load_Wavefunctions(ikp, all_bands, wf_cache)
    
    PARALLEL FOR ik = 1 TO n_k DO
      CALL Load_Wavefunctions(ik, all_bands, wf_cache)
      q = k_points[ik] - k_points[ikp]
      
      FOR iv = 1 TO n_v DO
        FOR ic = 1 TO n_c DO
          FOR ivp = 1 TO n_v DO
            FOR icp = 1 TO n_c DO
              
              K_ex_cpl = 0
              
              FOR iG = 1 TO n_G_exch DO
                q_G = q + G_vectors[iG]
                v_coul = 4 * PI / |q_G|^2
                
                // 注意: v-c' と v'-c の組み合わせ
                M_vcp = Matrix_Element(iv, icp, ik, ikp, iG, wf_cache)
                M_vpc = Matrix_Element(ivp, ic, ikp, ik, iG, wf_cache)
                
                IF (conjugate) THEN
                  K_ex_cpl += -2 * v_coul * CONJ(M_vcp * M_vpc)
                ELSE
                  K_ex_cpl += -2 * v_coul * M_vcp * M_vpc
                END IF
              END FOR
              
              K_ex_cpl = K_ex_cpl / Volume
              
              IF (conjugate) THEN
                H_psi[iv,ic,ik] += CONJ(K_ex_cpl * psi[ivp,icp,ikp])
              ELSE
                H_psi[iv,ic,ik] += K_ex_cpl * psi[ivp,icp,ikp]
              END IF
              
            END FOR
          END FOR
        END FOR
      END FOR
      
    END PARALLEL FOR
  END FOR
  
  CALL MPI_Allreduce_Sum(H_psi)
  
END PROCEDURE
```

---

## 5. データ構造実装（Fortran）

### 5.1 型定義

```fortran
! Lanczos結合ベクトル
type :: lanczos_vector_coupling
  complex(DP), allocatable :: psi_plus(:,:,:,:)   ! (n_v, n_c, n_k, n_spin)
  complex(DP), allocatable :: psi_minus(:,:,:,:)  ! (n_v, n_c, n_k, n_spin)
  real(DP) :: norm
  integer :: iteration
  logical :: is_right_vector
end type

! Lanczos状態
type :: lanczos_state
  type(lanczos_vector_coupling) :: psi_n, psi_nm1
  type(lanczos_vector_coupling) :: phi_n, phi_nm1
  complex(DP), allocatable :: a(:)
  real(DP), allocatable :: b(:)
  integer :: n_iter, n_iter_max
  logical :: converged
  real(DP) :: conv_epsilon
end type

! 波動関数キャッシュ
type :: wavefunction_cache
  complex(DP), allocatable :: wf(:,:,:)  ! (n_G, n_bands, n_k_cached)
  integer, allocatable :: k_list(:)
  integer, allocatable :: band_list(:)
  integer, allocatable :: access_time(:)
  integer :: n_hits, n_misses
  integer :: current_size, max_size
end type

! カーネルワークスペース
type :: kernel_workspace
  complex(DP), allocatable :: H_psi_plus(:,:,:,:)
  complex(DP), allocatable :: H_psi_minus(:,:,:,:)
  complex(DP), allocatable :: W(:,:,:)  ! (n_G, n_G, n_q)
  integer :: my_k_start, my_k_end
end type
```

---

## 6. メモリ管理戦略

### 6.1 メモリ割り当て計画

```
総メモリ = Lanczosベクトル + 波動関数キャッシュ + 遮蔽相互作用 + ワークスペース

Lanczosベクトル:
  4セット × 2成分 × n_k × n_v × n_c × 16 bytes
  = 8 × n_k × n_bands^2 × 16 bytes

波動関数キャッシュ:
  n_cache × n_k × n_bands × n_G × 16 bytes

遮蔽相互作用:
  n_q × n_G^2 × 16 bytes

ワークスペース:
  2 × 2 × n_k × n_v × n_c × 16 bytes
  = 4 × n_k × n_bands^2 × 16 bytes
```

### 6.2 メモリ削減テクニック

1. **遅延割り当て**: 必要になるまで配列を割り当てない
2. **早期解放**: 不要になった配列は即座に解放
3. **再利用**: 一時配列を複数の目的で再利用
4. **ストライド**: メモリアクセスパターンを最適化

---

## 7. 並列化実装

### 7.1 MPIレベル並列化

```
// k点分散
n_k_per_rank = n_k / n_ranks
my_k_start = rank * n_k_per_rank + 1
my_k_end = (rank + 1) * n_k_per_rank

// カーネル作用後の集約
CALL MPI_Allreduce(H_psi_local, H_psi_global, size, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD)
```

### 7.2 OpenMPレベル並列化

```fortran
!$OMP PARALLEL DO PRIVATE(ik, iv, ic, K_ex) SCHEDULE(DYNAMIC)
do ik = 1, n_k
  do iv = 1, n_v
    do ic = 1, n_c
      ! カーネル計算
    end do
  end do
end do
!$OMP END PARALLEL DO
```

### 7.3 ハイブリッド並列化

```
MPI: k'ループを分散
  OpenMP: kループを並列化
    SIMD: G-vectorループをベクトル化
```

---

## 8. 最適化技術

### 8.1 キャッシュ効率化

- **ループ順序**: 最速変化次元を内側に
- **ブロッキング**: キャッシュラインサイズに合わせてブロック化
- **プリフェッチ**: 次のデータを事前にロード

### 8.2 計算カーネルの最適化

- **BLAS/LAPACK利用**: 行列-ベクトル積はBLASで
- **FFT最適化**: FFTW等の高速ライブラリを使用
- **ベクトル化**: コンパイラの自動ベクトル化を促進

---

## 9. テスト戦略

### 9.1 単体テスト

- 各サブルーチンの独立テスト
- 既知の解析解との比較
- 対称性の検証

### 9.2 統合テスト

- 小規模系での全体動作確認
- 従来手法との比較
- 収束性の検証

### 9.3 性能テスト

- スケーラビリティ測定
- メモリ使用量の監視
- ボトルネックの特定

---

## 10. 実装チェックリスト

- [ ] データ構造の定義
- [ ] メモリ割り当て/解放ルーチン
- [ ] 波動関数キャッシュ管理
- [ ] 遷移行列要素計算（FFT版）
- [ ] 対角項の計算
- [ ] 交換カーネル（共鳴）
- [ ] 交換カーネル（結合）
- [ ] 相関カーネル（共鳴）
- [ ] 相関カーネル（結合）
- [ ] 双直交Lanczos反復
- [ ] 再直交化ルーチン
- [ ] 計量演算子付き内積
- [ ] 収束判定
- [ ] 三重対角行列対角化
- [ ] 誘電関数計算
- [ ] MPI並列化
- [ ] OpenMP並列化
- [ ] 入出力ルーチン
- [ ] エラーハンドリング
- [ ] テストスイート

---

**文書バージョン**: 1.0  
**作成日**: 2025-10-18  
**言語**: 日本語  
**文書形式**: Markdown  
**対応仕様書**: `bse_full_bands_coupling_specification.md` v1.0
