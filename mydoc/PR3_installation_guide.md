# PR#3 インストール手順書
# Installation Guide for Memory-Efficient Full-Band BSE Solver

---

## 1. 概要 (Overview)

本ドキュメントでは、PR#3で実装されたメモリ効率的な全バンドBSEソルバーを含むYamboのインストール手順を説明します。

---

## 2. システム要件 (System Requirements)

### 2.1 ハードウェア要件

| 項目 | 最小要件 | 推奨要件 |
|------|---------|---------|
| CPU | 4コア | 16コア以上 |
| メモリ | 16 GB | 64 GB以上 |
| ディスク | 50 GB | 200 GB以上 |
| ネットワーク | - | 高速Infiniband (並列計算時) |

### 2.2 ソフトウェア要件

**必須**:
- Fortran 2003以降対応コンパイラ (gfortran 4.9+, ifort 15+, nvfortran)
- C/C++コンパイラ (gcc, icc)
- MPI実装 (OpenMPI, MPICH, Intel MPI)
- BLAS/LAPACK ライブラリ
- FFTW3 または同等のFFTライブラリ
- netCDF (netcdf-c, netcdf-fortran)
- LibXC (交換相関汎関数ライブラリ)

**推奨**:
- ScaLAPACK (大規模並列計算用)
- SLEPc/PETSc (代替ソルバー用)
- HDF5 (データI/O最適化)
- OpenMP対応

### 2.3 依存ライブラリのバージョン

| ライブラリ | 最小バージョン | 推奨バージョン |
|-----------|--------------|--------------|
| GCC/gfortran | 4.9 | 10+ |
| OpenMPI | 2.0 | 4.0+ |
| FFTW3 | 3.3 | 3.3.8+ |
| netCDF | 4.4 | 4.7+ |
| LibXC | 3.0 | 5.0+ |
| ScaLAPACK | 2.0 | 2.1+ |

---

## 3. インストール手順 (Installation Steps)

### 3.1 ソースコードの取得

```bash
# Gitからクローン
git clone https://github.com/nobkt/yambo.git
cd yambo

# PR#3のブランチに切り替え（または既にマージ済みの場合はmasterを使用）
git checkout copilot/optimize-bse-calculation
# または
git checkout create-bse-calculation-specification  # PR#3のベースブランチ
```

### 3.2 依存ライブラリのインストール

#### Ubuntu/Debian系の場合

```bash
# 基本コンパイラとビルドツール
sudo apt update
sudo apt install -y build-essential gfortran g++ make autoconf automake libtool

# MPI
sudo apt install -y libopenmpi-dev openmpi-bin

# 数学ライブラリ
sudo apt install -y libblas-dev liblapack-dev libfftw3-dev

# netCDF
sudo apt install -y libnetcdf-dev libnetcdff-dev

# LibXC
sudo apt install -y libxc-dev

# ScaLAPACK（オプション）
sudo apt install -y libscalapack-openmpi-dev
```

#### CentOS/RHEL系の場合

```bash
# EPELリポジトリを有効化
sudo yum install -y epel-release

# 基本パッケージ
sudo yum install -y gcc gcc-gfortran gcc-c++ make autoconf automake libtool

# OpenMPI
sudo yum install -y openmpi openmpi-devel
module load mpi/openmpi-x86_64

# 数学ライブラリ
sudo yum install -y blas-devel lapack-devel fftw-devel

# netCDF
sudo yum install -y netcdf-devel netcdf-fortran-devel

# ScaLAPACK（オプション）
sudo yum install -y scalapack-openmpi-devel
```

#### HPCクラスタ環境の場合

```bash
# Intel環境の例
module load intel/2021
module load intelmpi/2021
module load mkl/2021
module load netcdf/4.7
module load libxc/5.0
module load fftw/3.3.8
```

### 3.3 Yamboの設定

```bash
cd yambo

# 設定スクリプトを実行
./configure --help  # オプション一覧を確認

# 基本設定（OpenMPI + gfortran）
./configure \
  --enable-open-mp \
  --enable-mpi \
  --with-mpi-libs="-lmpi" \
  --with-blas-libs="-lblas" \
  --with-lapack-libs="-llapack" \
  --with-fft-libs="-lfftw3" \
  --with-netcdf-path="/usr" \
  --with-libxc-path="/usr"

# Intel環境での設定例
./configure \
  FC=ifort \
  CC=icc \
  --enable-open-mp \
  --enable-mpi \
  --with-mpi-libs="-lmpi" \
  --with-blas-libs="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential" \
  --with-lapack-libs="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential" \
  --with-fft-includedir="$MKLROOT/include/fftw" \
  --with-fft-libs="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential" \
  --with-netcdf-path="$NETCDF_PATH" \
  --with-libxc-path="$LIBXC_PATH"
```

### 3.4 コンパイル

```bash
# ビルド（並列コンパイル）
make -j8

# または特定のターゲットのみビルド
make yambo -j8      # yambo実行ファイルのみ
make ypp -j8        # 後処理ツール

# インストール（オプション）
make install PREFIX=/path/to/install
```

### 3.5 ビルドの確認

```bash
# 実行ファイルの確認
ls -la bin/
# yambo, ypp, p2y などが存在することを確認

# バージョン確認
./bin/yambo -h

# 機能確認
./bin/yambo -features
# BSE関連機能が有効になっていることを確認
```

---

## 4. 確認テスト (Verification Tests)

### 4.1 基本テスト

```bash
# テストディレクトリに移動（テストデータが必要）
cd tests

# 基本的なBSE計算テスト
./run_tests.sh bse

# 全バンドBSEテスト
./run_tests.sh bse_full_bands
```

### 4.2 インタラクティブテスト

```bash
# インタラクティブモードでの確認
./bin/yambo

# 以下のメッセージが表示されれば正常
# <---> Initialization
# <---> Checking input file
```

---

## 5. 環境変数の設定 (Environment Variables)

### 5.1 推奨設定

```bash
# ~/.bashrc または ~/.bash_profile に追加

# Yamboのパス
export YAMBO_HOME=/path/to/yambo
export PATH=$YAMBO_HOME/bin:$PATH

# OpenMPスレッド数
export OMP_NUM_THREADS=4
export OMP_STACKSIZE=512M

# MPI設定
export OMPI_MCA_btl=^openib  # InfiniBandの問題を回避

# 一時ディレクトリ（高速ストレージを推奨）
export YAMBO_TMP=/scratch/$USER/yambo_tmp
mkdir -p $YAMBO_TMP

# メモリ制限（オプション）
ulimit -s unlimited  # スタックサイズ
```

### 5.2 HPCジョブスクリプト例

#### SLURM

```bash
#!/bin/bash
#SBATCH --job-name=yambo_bse
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=128G

module load openmpi
module load netcdf
module load libxc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=512M

cd $SLURM_SUBMIT_DIR

srun yambo -F input.in -J output -o b -y h
```

#### PBS/Torque

```bash
#!/bin/bash
#PBS -N yambo_bse
#PBS -l nodes=4:ppn=16
#PBS -l walltime=24:00:00
#PBS -l mem=128gb

module load openmpi
module load netcdf
module load libxc

export OMP_NUM_THREADS=4
export OMP_STACKSIZE=512M

cd $PBS_O_WORKDIR

mpirun -np 64 yambo -F input.in -J output -o b -y h
```

---

## 6. トラブルシューティング (Troubleshooting)

### 6.1 コンパイルエラー

#### 問題: netCDFが見つからない

```
Error: Cannot find netCDF library
```

**解決策**:
```bash
# netCDFのパスを明示的に指定
./configure --with-netcdf-path=/path/to/netcdf \
            --with-netcdf-libs="-lnetcdff -lnetcdf"
```

#### 問題: Fortranコンパイラエラー

```
Error: Unrecognized command-line option
```

**解決策**:
```bash
# コンパイラフラグを確認
./configure FC=gfortran FCFLAGS="-O2 -fallow-argument-mismatch"
```

### 6.2 実行時エラー

#### 問題: メモリ不足

```
YAMBO: Out of memory
```

**解決策**:
```bash
# スタックサイズを増加
ulimit -s unlimited

# OpenMPスタックを増加
export OMP_STACKSIZE=1G

# MPI並列数を増やしてメモリを分散
mpirun -np 128 yambo ...
```

#### 問題: MPI通信エラー

```
MPI_ABORT was invoked
```

**解決策**:
```bash
# OpenMPI設定を調整
export OMPI_MCA_btl="self,tcp"
export OMPI_MCA_mpi_cuda_support=0

# Intel MPIの場合
export I_MPI_FABRICS=shm:tcp
```

### 6.3 機能関連の問題

#### 問題: BSE計算が開始しない

**確認事項**:
1. SAVEディレクトリが存在するか確認
2. 必要なデータベースファイル（ns.db1, ns.wf等）が存在するか確認
3. 入力ファイルの構文エラーがないか確認

```bash
# 入力ファイルの確認
yambo -F input.in -I  # 初期化モード
```

---

## 7. アップデート手順 (Update Procedure)

### 7.1 既存インストールのアップデート

```bash
cd yambo

# 変更を保存
git stash

# 最新版を取得
git fetch origin
git checkout main
git pull

# 変更を復元
git stash pop

# 再ビルド
make clean
make -j8
```

### 7.2 設定の再適用

```bash
# 設定をリセットする場合
./configure --help  # オプション確認
./configure [以前と同じオプション]
make -j8
```

---

## 8. 付録：コンパイルオプション一覧

### 8.1 主要オプション

| オプション | 説明 |
|-----------|------|
| `--enable-mpi` | MPI並列を有効化 |
| `--enable-open-mp` | OpenMP並列を有効化 |
| `--with-blas-libs` | BLASライブラリパス |
| `--with-lapack-libs` | LAPACKライブラリパス |
| `--with-fft-libs` | FFTライブラリパス |
| `--with-netcdf-path` | netCDFインストールパス |
| `--with-libxc-path` | LibXCインストールパス |
| `--with-scalapack-libs` | ScaLAPACKライブラリパス |

### 8.2 最適化オプション

| オプション | 説明 |
|-----------|------|
| `--enable-time-profile` | タイミング情報を有効化 |
| `--enable-memory-profile` | メモリ使用量追跡を有効化 |
| `--enable-dp` | 倍精度演算（デフォルト） |
| `--enable-etsf-io` | ETSF I/Oサポート |

---

**文書バージョン**: 1.0  
**作成日**: 2024-12-03  
**対応PR**: PR#3
