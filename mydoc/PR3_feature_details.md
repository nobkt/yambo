# PR#3 実装機能詳細ドキュメント
# Feature Details: Memory-Efficient Full-Band BSE Solver with Coupling Mode Support

---

## 1. 概要 (Overview)

PR#3では、Yamboコードに**メモリ効率的な全バンドBSE計算機能**を実装しました。この機能により、BSE計算において結合項（coupling terms）を含む場合（`BSEmod="coupling"`）でも、全てのバンドを考慮した計算が可能になります。

### 1.1 実装の目的

従来のBSE計算では、カーネル行列を完全に保存する必要があり、全バンド計算は以下の問題がありました：

- **メモリ要求量**: O(N_k² × N_bands⁴) → 数PB規模
- **結合項を含む場合**: メモリ要求量がさらに4倍に増加

本実装では、双直交（biorthogonal）Lanczos法とオンザフライカーネル構築を用いることで：

- **メモリ要求量**: O(N_k × N_bands² + N_G²) → 数十GB規模
- **削減率**: 約10⁷倍（PBスケール → GBスケール）

---

## 2. 実装されたコンポーネント

### 2.1 新規ソルバー: `K_Haydock_full_bands.F`

**ファイル位置**: `src/bse/K_Haydock_full_bands.F`

**機能**:
- 擬エルミート（pseudo-Hermitian）結合ハミルトニアンに対する双直交Lanczos反復
- 計量演算子 Θ = diag(I, -I) を用いた内積計算（既存の`K_dot_product`の`i_kind=1`を利用）
- `l_BS_ares_from_res=.TRUE.`（対称性によりA(n)=0）と`.FALSE.`の両方のケースに対応
- 数値安定性チェック（breakdown検出、虚数部許容度チェック）

**アルゴリズムの数学的基礎**:

結合BSEハミルトニアンの構造：
```
H_BSE = [ H^res      C     ]
        [ -C*    -H^res*   ]
```

このハミルトニアンは擬エルミート性を満たします：
```
H_BSE = Θ H_BSE† Θ⁻¹
```
ここで Θ = diag(I, -I) は計量演算子です。

### 2.2 新規モジュール: `mod_BS_OnTheFly_kernel.F`

**ファイル位置**: `src/modules/mod_BS_OnTheFly_kernel.F`

**機能**:
- LRU（Least Recently Used）方式の波動関数キャッシュ
- メモリ効率的なワークスペース管理
- メモリ使用量見積もりユーティリティ
- 既存の`K_multiply_by_V`インフラストラクチャとの統合

**提供される型と機能**:

```fortran
! 波動関数キャッシュエントリ
type :: WF_cache_entry
  complex(SP), allocatable :: wf(:,:)      ! (n_G, n_bands)
  integer                  :: ik           ! k点インデックス
  integer                  :: access_time  ! LRU用アクセス時間
  logical                  :: valid
end type WF_cache_entry

! オンザフライカーネルワークスペース
type :: OnTheFly_workspace_t
  complex(SP), allocatable :: M_vvp(:,:,:)    ! 行列要素 M_{v,v'}^G
  complex(SP), allocatable :: M_ccp(:,:,:)    ! 行列要素 M_{c,c'}^G
  complex(SP), allocatable :: M_vcp(:,:,:)    ! 行列要素 M_{v,c'}^G (結合用)
  complex(SP), allocatable :: M_vpc(:,:,:)    ! 行列要素 M_{v',c}^G (結合用)
  complex(SP), allocatable :: W_q(:,:)        ! 遮蔽相互作用 W_{GG'}(q)
  integer                  :: n_G_exch        ! 交換用G-vector数
  integer                  :: n_G_corr        ! 相関用G-vector数
  logical                  :: allocated
end type OnTheFly_workspace_t
```

### 2.3 計量内積ユーティリティ: `K_dot_product_metric.F`

**ファイル位置**: `src/bse/K_dot_product_metric.F`

**機能**:
- 擬エルミート系の計量内積計算
- 双直交性チェックと再直交化サポート

**提供されるサブルーチン**:

| サブルーチン | 説明 |
|------------|------|
| `K_dot_product_metric` | 計量演算子Θを用いた内積 <V\|Θ\|W> の計算 |
| `K_check_biorthogonality` | 双直交性条件 <φ_m\|ψ_n> = δ_mn の検証 |
| `K_reorthogonalize_biorthogonal` | 修正Gram-Schmidt法による再直交化 |

### 2.4 統合: `K_solvers.F`と`mod_BS_solvers.F`の変更

**変更点**:
- `mod_BS_solvers.F`: `run_Haydock_full_bands`フラグを追加
- `K_solvers.F`: 新しいソルバーへのディスパッチを追加（optics/kerr/magnons/dichroism）

---

## 3. アルゴリズムの詳細

### 3.1 双直交Lanczos反復

**初期化**:
1. 初期ベクトル |ψ₀⟩ を双極子遷移で設定
2. 左ベクトル |φ₀⟩ = Θ|ψ₀⟩ で初期化
3. 初期係数 C₁ = ⟨φ₀|Θ|ψ₀⟩ を計算

**反復** (n = 1, 2, ...):
```
Step 1: A_n = ⟨φ_n|H|ψ_n⟩
        (l_BS_ares_from_res=.TRUE.の場合、対称性によりA_n = 0)

Step 2: 直交化
        |ψ̃_{n+1}⟩ = H|ψ_n⟩ - A_n|ψ_n⟩ - B_n|ψ_{n-1}⟩

Step 3: B_{n+1}² = ⟨ψ_n|Θ|H|ψ_n⟩
        (数値安定性チェック: B² > 0, Im(B²)/Re(B²) < ε)

Step 4: 規格化
        |ψ_{n+1}⟩ = |ψ̃_{n+1}⟩ / B_{n+1}
        |φ_{n+1}⟩ = |ψ_{n+1}⟩  (擬エルミート性を利用)

Step 5: C_{n+1} = ⟨φ₀|Θ|ψ_n⟩

Step 6: 収束判定
        誘電関数の変化が閾値以下なら終了
```

### 3.2 応答関数の計算

三重対角行列からGreen関数を連分数展開で計算：

```
G(ω) = 1 / (ω - a₀ - b₁² / (ω - a₁ - b₂² / (ω - a₂ - ...)))
```

誘電関数：
```
ε_M(ω) = 1 - (8π/V) × Σ_λ |A_λ^(0)|² / (E_λ - ω - iη)
```

### 3.3 数値安定性

**監視項目**:

| パラメータ | 閾値 | 説明 |
|-----------|------|------|
| `ORTHO_TOL` | 10⁻¹² | 直交性許容誤差 |
| `ORTHO_WARN` | 10⁻¹⁰ | 直交性警告閾値 |
| `B_MIN` | 10⁻¹⁴ | 最小規格化係数（収束判定） |
| `IMAG_TOL` | 10⁻⁶ | B²の虚数部許容相対誤差 |

---

## 4. メモリ使用量

### 4.1 メモリ見積もり式

```
M_total = M_Lanczos + M_cache + M_workspace

M_Lanczos = 8 × N_k × N_bands² × 16 bytes
           (4セット × 2成分のLanczosベクトル)

M_cache = N_cache × N_bands × N_G × 16 bytes
         (波動関数キャッシュ)

M_workspace = 4 × (N_bands/2)² × N_G × 16 bytes
             (行列要素計算用ワークスペース)
```

### 4.2 数値例

典型的なパラメータ (N_k=1000, N_bands=500, N_G=5000, N_cache=100):
```
M_Lanczos   = 32 GB
M_cache     = 8 GB  
M_workspace = 40 GB
-----------------
M_total     ≈ 80 GB
```

従来手法との比較:
```
従来手法: 4000 PB (結合項含む全バンド)
提案手法: 80 GB
削減率: ≈ 5 × 10⁷ 倍
```

---

## 5. 対応する物理量

本実装は以下の物理量の計算に対応しています：

| 物理量 | OBS引数 | 説明 |
|--------|---------|------|
| 光学応答 | `"optics"` | 誘電関数、吸収スペクトル |
| カー効果 | `"kerr"` | 磁気光学効果 |
| マグノン | `"magnons"` | スピン励起 |
| 二色性 | `"dichroism"` | 円偏光二色性 |

---

## 6. 重要な注意事項

### 6.1 制限事項

- `BS_res_ares_n_mat==2`（2つのBSE行列を持つ場合）には対応していません
- 現時点では、完全なオンザフライ計算は`K_multiply_by_V`インフラストラクチャに委譲

### 6.2 ヒューリスティック処理・Fallbackの排除

**仕様に従い、本実装では以下を完全に排除しています**:

- ヒューリスティックなバンド選択やカットオフ
- 計算が困難な場合のFallback処理
- 近似的な対称性の適用

全ての処理は数学的に厳密な定式化に基づいています。

---

## 7. 関連ドキュメント

- [bse_full_bands_coupling_specification.md](bse_full_bands_coupling_specification.md) - 詳細仕様書
- [bse_full_bands_coupling_design.md](bse_full_bands_coupling_design.md) - 詳細設計書
- [memory_efficient_bse_full_bands.md](memory_efficient_bse_full_bands.md) - 理論的背景

---

**文書バージョン**: 1.0  
**作成日**: 2024-12-03  
**対応PR**: PR#3
