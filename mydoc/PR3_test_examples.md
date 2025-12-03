# PR#3 テスト例題
# Test Examples for Memory-Efficient Full-Band BSE Solver Verification

---

## 概要

本ドキュメントでは、PR#3で実装されたメモリ効率的な全バンドBSEソルバーが正常に動作していることを確認するためのテスト例題を提供します。各テストは独立して実行可能であり、特定の機能の検証に焦点を当てています。

---

## 1. 基本動作確認テスト

### 1.1 テスト目的

新しいソルバー `K_Haydock_full_bands` が正常に呼び出され、基本的な反復が実行されることを確認します。

### 1.2 テスト入力ファイル: `test_basic_functionality.in`

```fortran
#
# Test 1: Basic Functionality Test
# 基本動作確認テスト
#
# 目的: ソルバーの起動と基本反復の確認
#
optics
bse

# BSEモード
BSEmod= "coupling"
BSKmod= "SEX"

# ソルバー設定
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

# 小さなテスト系
BSEBands=  1 | 10

# エネルギー範囲
BEnRange=  0.0 | 5.0 eV
BEnSteps= 100

# ブロードニング
BDmRange= 0.1 | 0.1 eV

# 収束パラメータ（緩め）
Haydock_treshold= 0.1
Haydock_iterMAX= 50

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 1.3 期待される結果

1. 計算が正常に開始される
2. `[Haydock Full Bands] Iteration N` メッセージが出力される
3. 最大反復回数以内に計算が終了する（収束または最大反復到達）
4. 誘電関数ファイルが生成される

### 1.4 確認コマンド

```bash
# テスト実行
yambo -F test_basic_functionality.in -J test1_output -o b -y h

# 結果確認
grep "Haydock Full Bands" r-test1_output.report
ls -la o-test1_output.eps_*

# 成功条件：
# - "Haydock Full Bands" 含むメッセージが存在
# - eps ファイルが生成されている
```

---

## 2. 収束性テスト

### 2.1 テスト目的

Lanczos反復が適切に収束し、指定した閾値を満たすことを確認します。

### 2.2 テスト入力ファイル: `test_convergence.in`

```fortran
#
# Test 2: Convergence Test
# 収束性テスト
#
# 目的: 収束閾値を満たすことの確認
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

BSEBands=  1 | 20

BEnRange=  0.0 | 8.0 eV
BEnSteps= 200

BDmRange= 0.05 | 0.05 eV

# 厳密な収束閾値
Haydock_treshold= 0.02
Haydock_iterMAX= 200
Haydock_iterIO= 10

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 2.3 期待される結果

1. 収束閾値 (0.02 = 2%) を満たして計算終了
2. 最終メッセージに `Accuracy (reached)` が閾値以下

### 2.4 確認コマンド

```bash
# テスト実行
yambo -F test_convergence.in -J test2_output -o b -y h

# 収束確認
tail -20 r-test2_output.report | grep -E "Accuracy|Iteration"

# 成功条件：
# - "Accuracy (reached)" の値が 0.02 以下
# - "Iterations" の値が Haydock_iterMAX 未満
```

---

## 3. 数値安定性テスト

### 3.1 テスト目的

長時間反復でも数値的に安定であり、breakdown や警告が発生しないことを確認します。

### 3.2 テスト入力ファイル: `test_numerical_stability.in`

```fortran
#
# Test 3: Numerical Stability Test
# 数値安定性テスト
#
# 目的: 長時間反復での安定性確認
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

BSEBands=  1 | 30

BEnRange=  0.0 | 10.0 eV
BEnSteps= 300

BDmRange= 0.05 | 0.05 eV

# 多くの反復を強制
Haydock_treshold= 0.001
Haydock_iterMAX= 300
Haydock_iterIO= 20

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 3.3 期待される結果

1. "Lanczos breakdown" 警告が発生しない
2. "Large imaginary part in B^2" 警告が発生しない
3. 計算が正常終了する

### 3.4 確認コマンド

```bash
# テスト実行
yambo -F test_numerical_stability.in -J test3_output -o b -y h

# 警告確認
grep -i "breakdown\|warning\|error" r-test3_output.report l-test3_output.log

# 成功条件：
# - "Lanczos breakdown" が存在しない
# - "Large imaginary part" が存在しないか、許容回数以下
# - 計算が正常終了
```

---

## 4. 結合項（Coupling）テスト

### 4.1 テスト目的

結合項を含むBSE計算が正しく実行されることを確認します。

### 4.2 テスト入力ファイル: `test_coupling_mode.in`

```fortran
#
# Test 4: Coupling Mode Test
# 結合項モードテスト
#
# 目的: BSEmod="coupling" での正しい動作確認
#
optics
bse

# 明示的に結合モードを指定
BSEmod= "coupling"
BSKmod= "SEX"
BSSmod= "h"
BSHayTer

% BSEHay_full_bands
.TRUE.
%

BSEBands=  1 | 25

BEnRange=  0.0 | 8.0 eV
BEnSteps= 200

BDmRange= 0.05 | 0.05 eV

Haydock_treshold= 0.02
Haydock_iterMAX= 200

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 4.3 期待される結果

1. 計算が正常に完了
2. 誘電関数が物理的に妥当な形状（正の虚部など）

### 4.4 確認コマンド

```bash
# テスト実行
yambo -F test_coupling_mode.in -J test4_output -o b -y h

# 結果確認
python3 << 'EOF'
import numpy as np
data = np.loadtxt('o-test4_output.eps_q1_haydock_bse', comments='#')
# 虚部が非負であることを確認
assert np.all(data[:,2] >= -1e-6), "Im[eps] should be non-negative"
print("Test 4 PASSED: Im[eps] is non-negative")
EOF
```

---

## 5. メモリ使用量テスト

### 5.1 テスト目的

大規模バンド数でもメモリ使用量が許容範囲内であることを確認します。

### 5.2 テスト入力ファイル: `test_memory_efficiency.in`

```fortran
#
# Test 5: Memory Efficiency Test
# メモリ効率テスト
#
# 目的: メモリ使用量が従来手法より大幅に削減されていることの確認
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

# 大きなバンド数
BSEBands=  1 | 100

BEnRange=  0.0 | 10.0 eV
BEnSteps= 100

BDmRange= 0.1 | 0.1 eV

# 収束パラメータ（緩め）
Haydock_treshold= 0.1
Haydock_iterMAX= 30

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 5.3 期待される結果

1. 計算が正常に完了（メモリ不足エラーなし）
2. メモリ使用量が理論値の範囲内

### 5.4 確認コマンド

```bash
# メモリ監視付きで実行
/usr/bin/time -v yambo -F test_memory_efficiency.in -J test5_output -o b -y h 2>&1 | tee test5_memory.log

# メモリ使用量確認
grep "Maximum resident set size" test5_memory.log

# 成功条件：
# - 計算が正常終了
# - メモリ使用量が100 GB以下（典型的な系の場合）
```

---

## 6. 物理量検証テスト

### 6.1 テスト目的

計算された誘電関数が物理的に妥当であることを確認します（和則など）。

### 6.2 テスト入力ファイル: `test_physical_validity.in`

```fortran
#
# Test 6: Physical Validity Test
# 物理的妥当性テスト
#
# 目的: f-sum rule等の和則の検証
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

BSEBands=  1 | 50

# 広いエネルギー範囲（和則検証用）
BEnRange=  0.0 | 30.0 eV
BEnSteps= 1000

BDmRange= 0.05 | 0.05 eV

Haydock_treshold= 0.01
Haydock_iterMAX= 300

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 6.3 期待される結果

1. f-sum rule が近似的に満たされる
2. Kramers-Kronig関係が満たされる

### 6.4 確認コマンド

```bash
# テスト実行
yambo -F test_physical_validity.in -J test6_output -o b -y h

# 和則検証
python3 << 'EOF'
import numpy as np

# データ読み込み
data = np.loadtxt('o-test6_output.eps_q1_haydock_bse', comments='#')
omega = data[:,0]
eps_im = data[:,2]

# f-sum rule: ∫ω×Im[ε(ω)]dω = π×ω_p²/2
# 数値積分
integrand = omega * eps_im
f_sum = np.trapz(integrand, omega)

print(f"f-sum rule integral: {f_sum:.4f}")
print("(Compare with π×ω_p²/2 for your system)")

# 虚部が正であることを確認
negative_count = np.sum(eps_im < -1e-6)
print(f"Negative Im[eps] count: {negative_count}")

if negative_count == 0:
    print("Test 6 PASSED: Physical consistency verified")
else:
    print(f"Test 6 WARNING: {negative_count} negative Im[eps] values")
EOF
```

---

## 7. 並列実行テスト

### 7.1 テスト目的

MPI並列およびOpenMP並列で正常に動作することを確認します。

### 7.2 テスト入力ファイル: `test_parallel_execution.in`

```fortran
#
# Test 7: Parallel Execution Test
# 並列実行テスト
#
# 目的: 並列計算での正常動作確認
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

BSEBands=  1 | 30

BEnRange=  0.0 | 10.0 eV
BEnSteps= 200

BDmRange= 0.05 | 0.05 eV

Haydock_treshold= 0.02
Haydock_iterMAX= 100

% BLongDir
1.0 | 0.0 | 0.0
%

Gauge= "length"
```

### 7.3 期待される結果

1. 並列実行が正常に完了
2. シリアル実行と同一の結果（許容誤差内）

### 7.4 確認コマンド

```bash
# シリアル実行
yambo -F test_parallel_execution.in -J test7_serial -o b -y h

# 並列実行
export OMP_NUM_THREADS=4
mpirun -np 4 yambo -F test_parallel_execution.in -J test7_parallel -o b -y h

# 結果比較
python3 << 'EOF'
import numpy as np

serial = np.loadtxt('o-test7_serial.eps_q1_haydock_bse', comments='#')
parallel = np.loadtxt('o-test7_parallel.eps_q1_haydock_bse', comments='#')

# 相対誤差計算
diff = np.abs(serial[:,2] - parallel[:,2])
max_diff = np.max(diff)
rel_diff = max_diff / (np.max(np.abs(serial[:,2])) + 1e-10)

print(f"Maximum absolute difference: {max_diff:.2e}")
print(f"Relative difference: {rel_diff:.2e}")

if rel_diff < 1e-6:
    print("Test 7 PASSED: Parallel results match serial")
else:
    print("Test 7 WARNING: Results differ significantly")
EOF
```

---

## 8. テスト実行スクリプト

### 8.1 一括テスト実行スクリプト: `run_all_tests.sh`

```bash
#!/bin/bash
#
# PR#3 Implementation Verification Tests
# 全テスト一括実行スクリプト
#

set -e  # エラー時に停止

echo "========================================"
echo "PR#3 Implementation Verification Tests"
echo "========================================"
echo ""

# テストディレクトリ作成
TEST_DIR="PR3_tests_$(date +%Y%m%d_%H%M%S)"
mkdir -p $TEST_DIR
cd $TEST_DIR

# テスト結果カウンター
PASSED=0
FAILED=0
WARNINGS=0

# テスト関数
run_test() {
    local test_name=$1
    local input_file=$2
    local check_command=$3
    
    echo "----------------------------------------"
    echo "Running: $test_name"
    echo "----------------------------------------"
    
    # 入力ファイル作成（heredocで記載）
    cat > $input_file
    
    # 実行
    if yambo -F $input_file -J ${test_name}_output -o b -y h > /dev/null 2>&1; then
        echo "  Execution: COMPLETED"
        
        # 結果チェック
        if eval "$check_command"; then
            echo "  Result: PASSED"
            ((PASSED++))
        else
            echo "  Result: FAILED"
            ((FAILED++))
        fi
    else
        echo "  Execution: FAILED"
        ((FAILED++))
    fi
    echo ""
}

# テストサマリー出力
echo ""
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Passed:   $PASSED"
echo "Failed:   $FAILED"
echo "Warnings: $WARNINGS"
echo "========================================"

if [ $FAILED -eq 0 ]; then
    echo "ALL TESTS PASSED!"
    exit 0
else
    echo "SOME TESTS FAILED!"
    exit 1
fi
```

---

## 9. 期待される出力例

### 9.1 正常終了時のレポート出力例

```
 <---> [08] Full-Band Haydock Solver (Biorthogonal) in the optics basis @q1
 <---> Accuracy (requested)       0.02000 [o/o]
 <---> [Haydock Full Bands] Iteration 1
 <---> [Haydock Full Bands] Iteration 10 Accuracy 45.32|0.02
 <---> [Haydock Full Bands] Iteration 50 Accuracy 5.12|0.02
 <---> [Haydock Full Bands] Iteration 100 Accuracy 0.89|0.02
 <---> [Haydock Full Bands] Iteration 125 Accuracy 0.015|0.02
 <---> Accuracy (reached)         0.01500 [o/o]
 <---> Iterations                 125
```

### 9.2 収束グラフの生成

```bash
# 収束履歴の抽出とプロット
grep "Haydock Full Bands.*Accuracy" r-test_output.report | \
  awk '{print $5, $6}' | tr '|' ' ' > convergence.dat

gnuplot << 'EOF'
set terminal png
set output 'convergence.png'
set xlabel 'Iteration'
set ylabel 'Accuracy [%]'
set logscale y
plot 'convergence.dat' using 1:2 with lines title 'Accuracy', \
     'convergence.dat' using 1:3 with lines title 'Threshold'
EOF
```

---

## 10. 注意事項

### 10.1 テスト前の準備

- SAVEディレクトリに必要なデータベースファイルが存在すること
- 遮蔽相互作用（W）が事前計算されていること
- 十分なメモリとディスク容量があること

### 10.2 テスト結果の解釈

- 収束閾値は系やパラメータに依存
- 小さな警告は必ずしもエラーではない
- 物理的妥当性の確認が最も重要

### 10.3 トラブルシューティング

問題が発生した場合は以下を確認：
1. 入力ファイルの構文
2. 必要なデータベースの存在
3. メモリ使用量
4. 詳細ログ（l-*.log）

---

**文書バージョン**: 1.0  
**作成日**: 2024-12-03  
**対応PR**: PR#3
