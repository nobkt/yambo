# PR#3 機能使用ガイド
# Usage Guide: Memory-Efficient Full-Band BSE Solver

---

## 1. はじめに (Introduction)

本ガイドでは、PR#3で実装されたメモリ効率的な全バンドBSEソルバーの使用方法を説明します。このソルバーは、BSE計算において結合項（coupling）を含む場合でも、全てのバンドを考慮した計算を効率的に実行できます。

---

## 2. 前提条件 (Prerequisites)

### 2.1 必要な事前計算

新機能を使用する前に、以下の計算が完了している必要があります：

1. **DFT計算** (Quantum ESPRESSO等)
   - 自己無撞着場計算（SCF）
   - バンド構造計算
   
2. **Yambo初期化**
   ```bash
   yambo -F
   ```

3. **GW計算**（オプション、準粒子補正を含める場合）
   ```bash
   yambo -p p -g n
   ```

4. **遮蔽相互作用の計算**
   ```bash
   yambo -b -o b -k sex
   ```

### 2.2 必要なデータベースファイル

- `SAVE/`ディレクトリ内：
  - `ns.db1` - システム情報
  - `ns.wf` - 波動関数
  - `ns.kb_pp*` - 擬ポテンシャル情報
  
- BSE用：
  - `ndb.BS_PAR_Q*` - BSEカーネル（事前計算が必要）

---

## 3. 基本的な使用方法 (Basic Usage)

### 3.1 入力変数の設定

新しいソルバーを使用するには、`run_Haydock_full_bands`フラグを有効にする必要があります。
これはYamboの入力ファイルで以下のように設定します：

```fortran
# BSE計算設定
optics                           # 光学応答計算を有効化
BSEmod= "coupling"               # 結合項を含むBSE計算
BSKmod= "SEX"                    # 遮蔽交換カーネル

# Haydockソルバー設定
BSSmod= "h"                      # Haydock反復法
BSHayTer                         # Haydock terminator（オプション）
BSEHay_full_bands= .TRUE.        # 全バンドソルバーを有効化

# 収束パラメータ
Haydock_treshold= 0.01           # 収束閾値（0.01 = 1%）
Haydock_iterMAX= 500             # 最大反復回数

# バンド範囲
BSEBands= 1 | 500                # 価電子バンド1～伝導バンド500

# エネルギー範囲
BEnRange= 0.0 | 10.0 eV          # 計算する周波数範囲
BEnSteps= 500                    # 周波数点数
BDmRange= 0.05 | 0.05 eV         # ブロードニング
```

### 3.2 実行コマンド

```bash
# 標準実行
yambo -F input_file -J output_dir -o b -y h

# MPI並列実行
mpirun -np 64 yambo -F input_file -J output_dir -o b -y h

# ハイブリッド並列実行
export OMP_NUM_THREADS=4
mpirun -np 16 yambo -F input_file -J output_dir -o b -y h
```

---

## 4. 詳細設定 (Advanced Configuration)

### 4.1 メモリ最適化

大規模系の計算では、以下のパラメータを調整してメモリ使用量を最適化できます：

```fortran
# 波動関数キャッシュサイズ（デフォルト: 100エントリ）
# キャッシュに保持する波動関数の最大数を指定
# メモリ使用量: N_cache × N_bands × N_G × 16 bytes
WFCacheSize= 50                  # キャッシュエントリ数を削減してメモリ節約

# バンドチャンキング（将来の拡張用）
BandChunkSize= 100               # 一度に処理するバンド数
```

### 4.2 収束制御

収束パラメータの詳細設定：

```fortran
# 収束閾値設定
Haydock_treshold= 0.001          # より厳密な収束（0.1%）
Haydock_iterMAX= 1000            # 反復回数上限を増加

# 出力頻度
Haydock_iterIO= 100              # 100反復ごとに中間結果出力
```

### 4.3 対称性の利用

対称性を活用してメモリと計算時間を削減：

```fortran
# 対称性設定
UseSymmetry                      # 空間対称性を利用
# TimeRevSymm                    # 時間反転対称性を利用
```

---

## 5. 出力ファイル (Output Files)

### 5.1 主要出力ファイル

計算完了後、以下のファイルが生成されます：

| ファイル名 | 内容 |
|-----------|------|
| `o-*.eps_q1_haydock_bse` | 誘電関数（実部・虚部） |
| `o-*.alpha_q1_haydock_bse` | 分極率 |
| `o-*.eel_q1_haydock_bse` | 電子エネルギー損失関数 |
| `r-*.report` | 計算レポート（収束履歴等） |
| `l-*.log` | 詳細ログ |

### 5.2 出力フォーマット

誘電関数ファイルの形式：
```
# E[eV]           Re[eps]       Im[eps]       Re[eps_o]     Im[eps_o]
  0.00000E+00     1.00000E+01   0.00000E+00   1.00000E+01   0.00000E+00
  0.02004E+00     1.00123E+01   0.00012E+00   1.00123E+01   0.00012E+00
  ...
```

### 5.3 収束履歴の確認

reportファイルで収束状況を確認：
```
[Haydock Full Bands] Iteration 50 Accuracy 5.23|0.01
[Haydock Full Bands] Iteration 100 Accuracy 1.12|0.01
[Haydock Full Bands] Iteration 150 Accuracy 0.45|0.01
[Haydock Full Bands] Iteration 175 Accuracy 0.008|0.01
Accuracy (reached)         0.008 [o/o]
Iterations                 175
```

---

## 6. トラブルシューティング (Troubleshooting)

### 6.1 よくある問題と解決策

#### 問題1: メモリ不足

**症状**: `YAMBO: Out of memory` エラー

**解決策**:
```fortran
# 波動関数キャッシュを削減
WFCacheSize= 20

# より多くのノードで並列実行
mpirun -np 128 yambo ...
```

#### 問題2: 収束しない

**症状**: 最大反復回数に達しても収束閾値を満たさない

**解決策**:
```fortran
# 収束閾値を緩和
Haydock_treshold= 0.05

# または反復回数を増加
Haydock_iterMAX= 2000

# terminatorを有効化
BSHayTer
```

#### 問題3: Lanczos breakdown

**症状**: `[Haydock Full Bands] Lanczos breakdown detected`

**原因**: BSE固有値に負の値が含まれている可能性

**解決策**:
- バンド範囲を見直す
- 遮蔽相互作用の計算パラメータを確認
- k点サンプリングを増やす

#### 問題4: 大きな虚数部の警告

**症状**: `[Haydock Full Bands] Large imaginary part in B^2`

**原因**: 数値精度の問題または不適切なパラメータ

**解決策**:
- G-vector数を増やす
- k点サンプリングを見直す

### 6.2 デバッグ情報の取得

詳細なデバッグ情報を得るには：

```bash
# 詳細ログを有効化
yambo -F input_file -J output_dir -o b -y h -V all
```

---

## 7. 性能チューニング (Performance Tuning)

### 7.1 並列化戦略

**推奨設定**:

| システムサイズ | MPIプロセス数 | OpenMPスレッド数 |
|--------------|--------------|-----------------|
| 小規模 (N_k < 100) | 8-16 | 4 |
| 中規模 (N_k ~ 1000) | 64-128 | 4 |
| 大規模 (N_k > 1000) | 256-1024 | 2-4 |

### 7.2 I/O最適化

大規模計算でのI/Oボトルネックを軽減：

```bash
# 高速ストレージを使用
export YAMBO_TMP=/scratch/$USER/yambo_tmp

# netCDF並列I/Oを有効化（対応している場合）
export YAMBO_NETCDF_PAR=1
```

### 7.3 メモリ使用量の見積もり

計算前にメモリ要求量を見積もる：

```
M_total ≈ 8 × N_k × N_bands² × 16 bytes (Lanczosベクトル)
        + N_cache × N_bands × N_G × 16 bytes (キャッシュ)
        + N_G² × 16 bytes (ワークスペース)
```

例: N_k=1000, N_bands=500, N_G=5000, N_cache=100
```
M_total ≈ 32 GB + 8 GB + 40 GB ≈ 80 GB
```

---

## 8. ベストプラクティス (Best Practices)

### 8.1 計算の流れ

推奨される計算手順：

1. **小規模テスト** (N_bands=50, 少数k点)
   - パラメータの妥当性確認
   - 収束挙動の確認

2. **中規模計算** (N_bands=200)
   - スペクトルの主要特徴の確認
   - メモリ使用量の確認

3. **本番計算** (全バンド)
   - 収束したパラメータで実行

### 8.2 収束の確認

収束を確認するためのチェックリスト：

- [ ] 誘電関数の形状が安定している
- [ ] ピーク位置が反復で変化しない
- [ ] 収束閾値を満たしている
- [ ] f-sum ruleがほぼ満たされている

### 8.3 結果の検証

結果の妥当性を確認する方法：

1. **和則の検証**
   ```python
   # f-sum rule: ∫ω×Im[ε(ω)]dω = π×ω_p²/2
   import numpy as np
   data = np.loadtxt('o-*.eps_q1_haydock_bse')
   omega = data[:,0]
   eps_im = data[:,2]
   integral = np.trapz(omega * eps_im, omega)
   print(f"f-sum rule integral: {integral}")
   ```

2. **既知の結果との比較**
   - 実験値との比較
   - 文献値との比較
   - 従来手法との比較

---

## 9. 参考文献 (References)

本実装に関連する理論的背景と詳細：

- [PR3_feature_details.md](PR3_feature_details.md) - 機能詳細
- [bse_full_bands_coupling_specification.md](bse_full_bands_coupling_specification.md) - 詳細仕様
- [bse_full_bands_coupling_design.md](bse_full_bands_coupling_design.md) - 設計書
- [memory_efficient_bse_full_bands.md](memory_efficient_bse_full_bands.md) - 理論背景

---

**文書バージョン**: 1.0  
**作成日**: 2024-12-03  
**対応PR**: PR#3
