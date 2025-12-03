# PR#3 チュートリアル入力ファイル集
# Tutorial Input Files for Memory-Efficient Full-Band BSE Calculation

---

## 概要

本ドキュメントでは、PR#3で実装されたメモリ効率的な全バンドBSEソルバーを使用するためのサンプル入力ファイルを提供します。

---

## 1. シリコン（バルク半導体）の例

### 1.1 基本入力ファイル: `Si_bse_full_bands_basic.in`

```fortran
#
# Silicon - Basic Full-Band BSE Calculation with Coupling
# メモリ効率的な全バンドBSE計算（結合項含む）
#
# 計算タイプ: 光学応答
optics                           # 光学応答計算
bse                              # BSE計算を有効化

# BSEモード設定
BSEmod= "coupling"               # 結合項を含む（反共鳴項も考慮）
BSKmod= "SEX"                    # 遮蔽交換カーネル

# ソルバー設定
BSSmod= "h"                      # Haydock反復法
BSHayTer                         # Haydock terminator有効化

# 全バンドソルバーを使用
% BSEHay_full_bands
.TRUE.
%

# バンド範囲（全バンド）
BSEBands=  1 | 100               # バンド1から100まで

# エネルギー範囲
BEnRange=  0.0  | 15.0  eV       # 0〜15 eV
BEnSteps=  500                   # 周波数点数

# ブロードニング
BDmRange=  0.05 | 0.05 eV        # 50 meV ローレンツ幅

# 収束パラメータ
Haydock_treshold= 0.02           # 2% 収束閾値
Haydock_iterMAX= 500             # 最大500反復
Haydock_iterIO= 50               # 50反復ごとに出力

# q点（光学極限）
% BLongDir
1.0 | 0.0 | 0.0                  # x方向偏光
%

# ゲージ
Gauge= "length"                  # 長さゲージ
```

### 1.2 高精度入力ファイル: `Si_bse_full_bands_accurate.in`

```fortran
#
# Silicon - High-Accuracy Full-Band BSE Calculation
# 高精度全バンドBSE計算
#
optics
bse

# BSEモード設定
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲（拡張）
BSEBands=  1 | 200               # より多くのバンドを含む

# エネルギー範囲
BEnRange=  0.0  | 20.0  eV
BEnSteps=  1000                  # 高密度周波数グリッド

# ブロードニング（小さめ）
BDmRange=  0.01 | 0.01 eV        # 10 meV

# 厳密な収束パラメータ
Haydock_treshold= 0.005          # 0.5% 収束閾値
Haydock_iterMAX= 1000            # 最大1000反復

# 複数偏光方向
% BLongDir
1.0 | 0.0 | 0.0                  # x方向
%

# 対称性
UseSymmetry
```

---

## 2. GaAs（III-V族半導体）の例

### 2.1 標準入力ファイル: `GaAs_bse_full_bands.in`

```fortran
#
# GaAs - Full-Band BSE Calculation with Coupling
# GaAs全バンドBSE計算
#
optics
bse

# BSEモード設定
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 80                # 価電子帯と伝導帯

# エネルギー範囲（GaAsのバンドギャップ付近）
BEnRange=  0.0  | 8.0  eV
BEnSteps=  400

# ブロードニング
BDmRange=  0.03 | 0.03 eV

# 収束パラメータ
Haydock_treshold= 0.02
Haydock_iterMAX= 500

# 偏光方向
% BLongDir
1.0 | 1.0 | 1.0                  # [111]方向
%

Gauge= "length"
```

---

## 3. hBN（2次元絶縁体）の例

### 3.1 標準入力ファイル: `hBN_bse_full_bands.in`

```fortran
#
# Hexagonal Boron Nitride - Full-Band BSE Calculation
# hBN全バンドBSE計算（2次元系）
#
optics
bse

# BSEモード設定
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 50

# エネルギー範囲（hBNの大きなギャップに対応）
BEnRange=  0.0  | 12.0  eV
BEnSteps=  600

# ブロードニング
BDmRange=  0.05 | 0.05 eV

# 収束パラメータ
Haydock_treshold= 0.02
Haydock_iterMAX= 500

# 面内偏光
% BLongDir
1.0 | 0.0 | 0.0                  # 面内x方向
%

# 2次元系用クーロン切断（事前設定が必要）
# CUTGeo= "box z"
# CUTBox= 0.0 | 0.0 | 30.0 bohr

Gauge= "length"
```

---

## 4. 磁気光学効果の例（Kerr効果）

### 4.1 カー効果入力ファイル: `kerr_effect.in`

```fortran
#
# Magneto-Optical Kerr Effect Calculation
# 磁気光学カー効果計算
#
kerr                             # カー効果計算
bse

# BSEモード設定
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 100

# エネルギー範囲
BEnRange=  0.0  | 10.0  eV
BEnSteps=  500

# ブロードニング
BDmRange=  0.05 | 0.05 eV

# 収束パラメータ
Haydock_treshold= 0.02
Haydock_iterMAX= 500

# 偏光方向（カー効果用）
% BLongDir
1.0 | 0.0 | 0.0
%
```

---

## 5. 二色性の例

### 5.1 円偏光二色性入力ファイル: `dichroism.in`

```fortran
#
# Circular Dichroism Calculation
# 円偏光二色性計算
#
dichroism                        # 二色性計算
bse

# BSEモード設定
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 100

# エネルギー範囲
BEnRange=  0.0  | 10.0  eV
BEnSteps=  500

# ブロードニング
BDmRange=  0.05 | 0.05 eV

# 収束パラメータ
Haydock_treshold= 0.02
Haydock_iterMAX= 500

# 偏光方向
% BLongDir
1.0 | 0.0 | 0.0
%
```

---

## 6. メモリ最適化設定の例

### 6.1 低メモリ設定: `low_memory.in`

```fortran
#
# Low Memory Configuration for Large Systems
# 大規模系向け低メモリ設定
#
optics
bse

BSEmod= "coupling"
BSKmod= "SEX"
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 500

# エネルギー範囲
BEnRange=  0.0  | 10.0  eV
BEnSteps=  200                   # 周波数点数を減らす

# ブロードニング
BDmRange=  0.1 | 0.1 eV          # 大きめのブロードニング

# 収束パラメータ（緩め）
Haydock_treshold= 0.05           # 5% 収束閾値
Haydock_iterMAX= 300             # 反復回数制限

# メモリ最適化
% WFCacheSize
20                               # キャッシュサイズを小さく
%

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

---

## 7. 並列計算設定の例

### 7.1 大規模並列設定: `parallel_config.in`

```fortran
#
# Parallel Configuration for HPC Systems
# HPCシステム向け並列設定
#
optics
bse

BSEmod= "coupling"
BSKmod= "SEX"
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# バンド範囲
BSEBands=  1 | 500

# エネルギー範囲
BEnRange=  0.0  | 15.0  eV
BEnSteps=  500

# ブロードニング
BDmRange=  0.05 | 0.05 eV

# 収束パラメータ
Haydock_treshold= 0.02
Haydock_iterMAX= 500
Haydock_iterIO= 20               # 頻繁な出力

# 偏光方向
% BLongDir
1.0 | 0.0 | 0.0
%

# 並列設定（ジョブスクリプトで設定）
# mpirun -np 128 yambo -F parallel_config.in -J output
# export OMP_NUM_THREADS=4

Gauge= "length"
UseSymmetry
```

---

## 8. 使用方法

### 8.1 入力ファイルの準備

1. 上記テンプレートをコピーして編集
2. システムに応じてバンド範囲とエネルギー範囲を調整
3. 利用可能なメモリに応じてキャッシュサイズを調整

### 8.2 実行コマンド

```bash
# 基本実行
yambo -F Si_bse_full_bands_basic.in -J Si_BSE_output -o b -y h

# 並列実行
mpirun -np 64 yambo -F Si_bse_full_bands_basic.in -J Si_BSE_output -o b -y h

# ハイブリッド並列実行
export OMP_NUM_THREADS=4
mpirun -np 16 yambo -F Si_bse_full_bands_basic.in -J Si_BSE_output -o b -y h
```

### 8.3 結果の確認

```bash
# 誘電関数の確認
gnuplot -e "plot 'o-Si_BSE_output.eps_q1_haydock_bse' u 1:3 w l"

# 収束履歴の確認
grep "Haydock Full Bands" r-Si_BSE_output.report
```

---

## 9. 注意事項

### 9.1 パラメータ選択のガイドライン

| パラメータ | 小規模系 | 中規模系 | 大規模系 |
|-----------|---------|---------|---------|
| BSEBands | 1\|50 | 1\|200 | 1\|500 |
| BEnSteps | 200 | 500 | 500 |
| Haydock_treshold | 0.01 | 0.02 | 0.05 |
| Haydock_iterMAX | 300 | 500 | 500 |
| WFCacheSize | 100 | 50 | 20 |

### 9.2 メモリ見積もり

```
メモリ ≈ 8 × N_k × N_bands² × 16 bytes + キャッシュ + ワークスペース
```

例：N_k=1000, N_bands=500の場合
```
≈ 8 × 1000 × 500² × 16 ≈ 32 GB (Lanczosベクトルのみ)
```

---

**文書バージョン**: 1.0  
**作成日**: 2024-12-03  
**対応PR**: PR#3
