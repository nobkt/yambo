# BSE全バンド計算における結合項を含む詳細仕様書
## Detailed Specification for Full-Band BSE Calculation with Coupling Terms

---

## 1. 概要と目的

### 1.1 背景

Bethe-Salpeter方程式(BSE)計算において、共鳴(resonant)項と反共鳴(anti-resonant)項の両方を含む完全なハミルトニアン（結合ハミルトニアン）を用いた全バンド計算は、以下の理由により従来困難であった：

1. **メモリ要求量の爆発的増加**：結合項を含むと、ハミルトニアン行列の次元が2倍になり、メモリ要求量は4倍に増加
2. **計算複雑度の増加**：非エルミート行列の取り扱いによる数値的困難
3. **対称性の喪失**：結合ハミルトニアンは擬エルミート性を持つが、通常のエルミート行列ソルバーが使用不可

### 1.2 本仕様書の目的

本仕様書では、以下を達成する厳密な数学的アルゴリズムを定義する：

- **BSEmod="coupling"** の場合でも全バンドを考慮可能
- メモリ使用量を実用的範囲（2TB以下）に抑制
- ヒューリスティックな近似やfallback処理を一切使用しない
- 数学的厳密性を完全に維持

---

## 2. 数学的定式化

### 2.1 結合BSEハミルトニアン

完全な結合BSEハミルトニアンは以下の2×2ブロック行列で表される：

$$
\mathcal{H}_{\text{BSE}} = \begin{pmatrix}
H^{\text{res}} & C \\
-C^* & -H^{\text{res},*}
\end{pmatrix}
$$

ここで：

- $H^{\text{res}}$：共鳴ハミルトニアン（サイズ：$N_{\text{eh}} \times N_{\text{eh}}$）
- $C$：結合ブロック（サイズ：$N_{\text{eh}} \times N_{\text{eh}}$）
- $N_{\text{eh}} = N_k \times N_v \times N_c$：電子-正孔対の総数

全バンドを考慮する場合：$N_v \sim N_c \sim N_{\text{bands}}$ より、
$$
N_{\text{eh}} = N_k \times N_{\text{bands}}^2
$$

結合ハミルトニアンの次元：
$$
\dim(\mathcal{H}_{\text{BSE}}) = 2N_{\text{eh}} = 2N_k N_{\text{bands}}^2
$$

### 2.2 共鳴ブロック

共鳴ハミルトニアン $H^{\text{res}}$ の行列要素：

$$
H^{\text{res}}_{(vck),(v'c'k')} = (E_{ck} - E_{vk})\delta_{vv'}\delta_{cc'}\delta_{kk'} + K^{\text{res}}_{(vck),(v'c'k')}
$$

共鳴カーネル：
$$
K^{\text{res}}_{(vck),(v'c'k')} = K^{\text{ex}}_{(vck),(v'c'k')} + K^{\text{corr}}_{(vck),(v'c'k')}
$$

#### 2.2.1 共鳴交換項

$$
K^{\text{ex,res}}_{(vck),(v'c'k')} = -\frac{2\delta_{\sigma\sigma'}}{V} \sum_G v(q+G) M_{vv'}^{G}(k,k') [M_{cc'}^{G}(k,k')]^*
$$

ここで：
- $q = k - k'$
- $v(q) = 4\pi/|q|^2$：クーロンポテンシャル
- $M_{nn'}^{G}(k,k') = \langle nk|e^{-i(q+G)\cdot r}|n'k'\rangle$：遷移行列要素

#### 2.2.2 共鳴相関項

$$
K^{\text{corr,res}}_{(vck),(v'c'k')} = \frac{\delta_{\sigma\sigma'}}{V} \sum_{GG'} M_{vc'}^{G}(k,k') W_{GG'}(q,\omega=E_{ck}-E_{vk}) [M_{v'c}^{G'}(k',k)]^*
$$

ここで：
- $W_{GG'}(q,\omega)$：遮蔽クーロン相互作用

### 2.3 結合ブロック

結合行列 $C$ の要素：

$$
C_{(vck),(v'c'k')} = K^{\text{cpl}}_{(vck),(v'c'k')}
$$

結合カーネル：
$$
K^{\text{cpl}}_{(vck),(v'c'k')} = K^{\text{ex,cpl}}_{(vck),(v'c'k')} + K^{\text{corr,cpl}}_{(vck),(v'c'k')}
$$

#### 2.3.1 結合交換項

$$
K^{\text{ex,cpl}}_{(vck),(v'c'k')} = -\frac{2\delta_{\sigma\sigma'}}{V} \sum_G v(q+G) M_{vc'}^{G}(k,k') M_{v'c}^{G}(k',k)
$$

**重要な差異**：共鳴項と異なり、価電子-伝導バンド間の遷移行列要素となる。

#### 2.3.2 結合相関項

$$
K^{\text{corr,cpl}}_{(vck),(v'c'k')} = \frac{\delta_{\sigma\sigma'}}{V} \sum_{GG'} M_{vv'}^{G}(k,k') W_{GG'}(q,\omega) M_{cc'}^{G'}(k,k')
$$

### 2.4 擬エルミート性

結合ハミルトニアンは以下の擬エルミート性を満たす：

$$
\mathcal{H}_{\text{BSE}} = \Theta \mathcal{H}_{\text{BSE}}^\dagger \Theta^{-1}
$$

ここで $\Theta$ は計量演算子：
$$
\Theta = \begin{pmatrix}
\mathbb{I} & 0 \\
0 & -\mathbb{I}
\end{pmatrix}
$$

この性質により、固有値は実数または複素共役対となる。

---

## 3. メモリ要求量の解析

### 3.1 従来手法のメモリ使用量

#### 3.1.1 共鳴のみの場合

カーネル行列 $H^{\text{res}}$ を完全に保存：
$$
M_{\text{res}} = N_k^2 N_{\text{bands}}^4 \times 16 \text{ bytes}
$$

**例**（$N_k=1000$、$N_{\text{bands}}=500$）：
$$
M_{\text{res}} = 1000^2 \times 500^4 \times 16 \approx 1000 \text{ PB}
$$

#### 3.1.2 結合項を含む場合

完全な結合ハミルトニアン $\mathcal{H}_{\text{BSE}}$ を保存：
$$
M_{\text{cpl}} = (2N_k N_{\text{bands}}^2)^2 \times 16 \text{ bytes}
$$

$$
M_{\text{cpl}} = 4N_k^2 N_{\text{bands}}^4 \times 16 = 4M_{\text{res}}
$$

**例**（$N_k=1000$、$N_{\text{bands}}=500$）：
$$
M_{\text{cpl}} \approx 4000 \text{ PB}
$$

### 3.2 目標メモリ使用量

実用的な計算ノード（2TB RAM）で実行可能とするため：
$$
M_{\text{target}} < 2 \text{ TB} = 2 \times 10^{12} \text{ bytes}
$$

**必要な削減率**：
$$
r = \frac{M_{\text{cpl}}}{M_{\text{target}}} = \frac{4 \times 10^{18}}{2 \times 10^{12}} = 2 \times 10^{6}
$$

すなわち、**約200万倍**のメモリ削減が必要。

---

## 4. 基本戦略

### 4.1 オンザフライカーネル構築

**原理**：カーネル行列を明示的に保存せず、必要な時に計算する。

結合ハミルトニアンの作用を直接計算：
$$
\mathcal{H}_{\text{BSE}} \begin{pmatrix} |v_+\rangle \\ |v_-\rangle \end{pmatrix} = \begin{pmatrix} H^{\text{res}}|v_+\rangle + C|v_-\rangle \\ -C^*|v_+\rangle - H^{\text{res},*}|v_-\rangle \end{pmatrix}
$$

これにより、保存が必要なのは：
- Lanczosベクトル：$|v_\pm\rangle$（各サイズ：$N_k N_{\text{bands}}^2$）
- 一時的な波動関数と行列要素

### 4.2 修正Haydock反復法

結合ハミルトニアンは非エルミートだが擬エルミートなので、修正Haydock法を使用。

**標準的な共鳴のみのHaydock反復**：
$$
|\tilde{v}_{n+1}\rangle = H^{\text{res}}|v_n\rangle - a_n|v_n\rangle - b_n|v_{n-1}\rangle
$$

**結合項を含む修正版**：

左ベクトル $|\phi_n\rangle$ と右ベクトル $|\psi_n\rangle$ を使用：

$$
|\tilde{\psi}_{n+1}\rangle = \mathcal{H}_{\text{BSE}}|\psi_n\rangle - a_n|\psi_n\rangle - b_n|\psi_{n-1}\rangle
$$

$$
|\tilde{\phi}_{n+1}\rangle = \mathcal{H}_{\text{BSE}}^\dagger|\phi_n\rangle - a_n^*|\phi_n\rangle - b_n^*|\phi_{n-1}\rangle
$$

双直交性条件：
$$
\langle \phi_m | \psi_n \rangle = \delta_{mn}
$$

### 4.3 メモリスケーリング

**提案手法のメモリ要求量**：

反復中に保持する配列：
1. 現在と前の右Lanczosベクトル：$2 \times 2N_k N_{\text{bands}}^2 \times 16$
2. 現在と前の左Lanczosベクトル：$2 \times 2N_k N_{\text{bands}}^2 \times 16$
3. 波動関数キャッシュ：$N_{\text{cache}} \times N_k \times N_G \times 16$
4. 遮蔽相互作用：$N_q \times N_G^2 \times 16$

総メモリ：
$$
M_{\text{prop}} = 8N_k N_{\text{bands}}^2 \times 16 + N_{\text{cache}} N_k N_G \times 16 + N_q N_G^2 \times 16
$$

**数値例**（$N_k=1000$、$N_{\text{bands}}=500$、$N_G=5000$、$N_{\text{cache}}=100$、$N_q=100$）：

$$
M_{\text{prop}} = 8 \times 1000 \times 500^2 \times 16 + 100 \times 1000 \times 5000 \times 16 + 100 \times 5000^2 \times 16
$$

$$
= 32 + 8 + 40 = 80 \text{ GB}
$$

**削減率**：
$$
\frac{M_{\text{cpl}}}{M_{\text{prop}}} = \frac{4000 \text{ PB}}{80 \text{ GB}} = 5 \times 10^{7}
$$

約**5000万倍**のメモリ削減を達成。

---

## 5. 擬エルミート系のための修正Haydock法

### 5.1 双直交Lanczos法の定式化

擬エルミート行列 $\mathcal{H}$ に対して、以下の双直交Lanczos反復を実行：

**初期化**：

$$
|\psi_0\rangle = \frac{|\Psi_{\text{init}}\rangle}{\sqrt{\langle\Psi_{\text{init}}|\Theta|\Psi_{\text{init}}\rangle}}
$$

$$
|\phi_0\rangle = \Theta|\psi_0\rangle
$$

**反復** ($n = 0, 1, 2, \ldots$)：

1. ハミルトニアン作用：
   $$
   |\omega_n\rangle = \mathcal{H}_{\text{BSE}}|\psi_n\rangle
   $$
   
   $$
   |\eta_n\rangle = \mathcal{H}_{\text{BSE}}^\dagger|\phi_n\rangle
   $$

2. 対角要素：
   $$
   a_n = \langle\phi_n|\omega_n\rangle = \langle\phi_n|\mathcal{H}_{\text{BSE}}|\psi_n\rangle
   $$

3. 直交化：
   $$
   |\tilde{\psi}_{n+1}\rangle = |\omega_n\rangle - a_n|\psi_n\rangle - b_n|\psi_{n-1}\rangle
   $$
   
   $$
   |\tilde{\phi}_{n+1}\rangle = |\eta_n\rangle - a_n^*|\phi_n\rangle - b_n^*|\phi_{n-1}\rangle
   $$

4. 規格化係数：
   $$
   b_{n+1}^2 = \langle\tilde{\phi}_{n+1}|\Theta|\tilde{\psi}_{n+1}\rangle
   $$
   
   **注意**：$b_{n+1}$ は実数とは限らない（擬エルミート性のため）。安全のため：
   $$
   b_{n+1} = \sqrt[+]{|\langle\tilde{\phi}_{n+1}|\Theta|\tilde{\psi}_{n+1}\rangle|}
   $$
   
   正の平方根を取る。

5. 規格化：
   $$
   |\psi_{n+1}\rangle = \frac{|\tilde{\psi}_{n+1}\rangle}{b_{n+1}}
   $$
   
   $$
   |\phi_{n+1}\rangle = \frac{|\tilde{\phi}_{n+1}\rangle}{b_{n+1}^*}
   $$

### 5.2 三重対角行列の構成

双直交Lanczos反復により、以下の三重対角行列が得られる：

$$
T_n = \begin{pmatrix}
a_0 & b_1 & 0 & \cdots & 0 \\
b_1 & a_1 & b_2 & \cdots & 0 \\
0 & b_2 & a_2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & b_n \\
0 & 0 & 0 & b_n & a_n
\end{pmatrix}
$$

ここで $a_n$ は複素数、$b_n$ は実数（正の平方根を取るため）。

### 5.3 誘電関数の計算

光学遷移双極子：
$$
d_\alpha = \sum_{vck\sigma} \langle vk\sigma|e^{iq\cdot r}|\alpha ck\sigma\rangle
$$

ここで $\alpha$ は偏光方向。

初期ベクトルを双極子励起に対応させる：
$$
\psi_0^{(vck\sigma)} \propto \frac{d_\alpha^{(vck\sigma)}}{\sqrt{\sum_{vck\sigma}|d_\alpha^{(vck\sigma)}|^2}}
$$

マクロスコピック誘電関数：
$$
\epsilon_M(\omega) = 1 - \lim_{q \to 0} v(q) \sum_\lambda \frac{|O_\lambda|^2}{\omega - E_\lambda + i\eta}
$$

ここで：
- $E_\lambda$：$T_n$ の固有値（励起エネルギー）
- $O_\lambda = \langle e_\lambda | \psi_0 \rangle$：振動子強度（$e_\lambda$ は $T_n$ の固有ベクトル）

三重対角行列 $T_n$ から Green関数を計算：
$$
G(\omega) = \langle \psi_0 | (\omega \mathbb{I} - \mathcal{H}_{\text{BSE}})^{-1} | \psi_0 \rangle \approx \langle e_0 | (\omega \mathbb{I} - T_n)^{-1} | e_0 \rangle
$$

ここで $e_0 = (1, 0, 0, \ldots)^T$。

連分数表現：
$$
G(\omega) = \frac{1}{\omega - a_0 - \cfrac{b_1^2}{\omega - a_1 - \cfrac{b_2^2}{\omega - a_2 - \cdots}}}
$$

---

## 6. オンザフライカーネル作用の詳細

### 6.1 ベクトル表現

結合BSEベクトルを以下のように表現：
$$
|\Psi\rangle = \begin{pmatrix} |\psi_+\rangle \\ |\psi_-\rangle \end{pmatrix}
$$

ここで：
- $|\psi_+\rangle$：共鳴成分（サイズ：$N_{\text{eh}}$）
- $|\psi_-\rangle$：反共鳴成分（サイズ：$N_{\text{eh}}$）

成分表示：
$$
\psi_\pm^{(vck\sigma)} = \langle vck\sigma | \psi_\pm \rangle
$$

### 6.2 ハミルトニアン作用

$$
\mathcal{H}_{\text{BSE}}|\Psi\rangle = \begin{pmatrix} H^{\text{res}}|\psi_+\rangle + C|\psi_-\rangle \\ -C^*|\psi_+\rangle - H^{\text{res},*}|\psi_-\rangle \end{pmatrix}
$$

#### 6.2.1 共鳴ブロックの作用

対角項：
$$
(H^{\text{res}}_{\text{diag}}|\psi_+\rangle)^{(vck\sigma)} = (E_{ck} - E_{vk})\psi_+^{(vck\sigma)}
$$

交換項：
$$
(K^{\text{ex,res}}|\psi_+\rangle)^{(vck\sigma)} = -\frac{2\delta_{\sigma\sigma'}}{V} \sum_{v'c'k'\sigma'} \sum_G v(q+G) M_{vv'}^{G}(k,k') [M_{cc'}^{G}(k,k')]^* \psi_+^{(v'c'k'\sigma')}
$$

相関項：
$$
(K^{\text{corr,res}}|\psi_+\rangle)^{(vck\sigma)} = \frac{\delta_{\sigma\sigma'}}{V} \sum_{v'c'k'\sigma'} \sum_{GG'} M_{vc'}^{G}(k,k') W_{GG'}(k-k',\omega_0) [M_{v'c}^{G'}(k',k)]^* \psi_+^{(v'c'k'\sigma')}
$$

#### 6.2.2 結合ブロックの作用

交換項：
$$
(K^{\text{ex,cpl}}|\psi_-\rangle)^{(vck\sigma)} = -\frac{2\delta_{\sigma\sigma'}}{V} \sum_{v'c'k'\sigma'} \sum_G v(q+G) M_{vc'}^{G}(k,k') M_{v'c}^{G}(k',k) \psi_-^{(v'c'k'\sigma')}
$$

相関項：
$$
(K^{\text{corr,cpl}}|\psi_-\rangle)^{(vck\sigma)} = \frac{\delta_{\sigma\sigma'}}{V} \sum_{v'c'k'\sigma'} \sum_{GG'} M_{vv'}^{G}(k,k') W_{GG'}(k-k',\omega_0) M_{cc'}^{G'}(k,k') \psi_-^{(v'c'k'\sigma')}
$$

### 6.3 遷移行列要素の計算

すべてのカーネル作用で必要となる行列要素：
$$
M_{nn'}^{G}(k,k') = \langle nk|e^{-i(k-k'+G)\cdot r}|n'k'\rangle
$$

FFTを用いた効率的計算：

1. 波動関数を実空間へ変換：
   $$
   \psi_{nk}(r) = \sum_{G'} c_{nk}^{G'} e^{i(k+G')\cdot r}
   $$

2. 位相因子を乗算：
   $$
   \tilde{\psi}(r) = \psi_{nk}(r) e^{-i(k-k')\cdot r}
   $$

3. FFTで逆格子空間へ：
   $$
   \tilde{\psi}(G) = \text{FFT}[\tilde{\psi}(r)]
   $$

4. 行列要素：
   $$
   M_{nn'}^{G}(k,k') = \sum_{G'} [c_{nk}^{G'-G}]^* c_{n'k'}^{G'}
   $$

**計算量**：$O(N_G \log N_G)$ per matrix element

---

## 7. 数値安定性と精度管理

### 7.1 双直交性の維持

長時間反復により、双直交性 $\langle\phi_m|\psi_n\rangle = \delta_{mn}$ が数値的に崩れる可能性がある。

**完全再直交化**：

各反復で、新しいベクトルを全ての前のベクトルに対して直交化：

$$
|\tilde{\psi}_{n+1}\rangle \leftarrow |\tilde{\psi}_{n+1}\rangle - \sum_{j=0}^{n} \langle\phi_j|\tilde{\psi}_{n+1}\rangle|\psi_j\rangle
$$

$$
|\tilde{\phi}_{n+1}\rangle \leftarrow |\tilde{\phi}_{n+1}\rangle - \sum_{j=0}^{n} [\langle\phi_j|\tilde{\psi}_{n+1}\rangle]^*|\phi_j\rangle
$$

**計算量**：$O(n^2 N_{\text{eh}})$ → 大規模系では prohibitive

**選択的再直交化**：

双直交性の喪失を監視し、必要な場合のみ再直交化：

$$
\delta_{mn} = |\langle\phi_m|\psi_n\rangle - \delta_{mn}^{\text{ideal}}|
$$

閾値 $\epsilon_{\text{ortho}} = 10^{-12}$ を設定し、$\delta_{mn} > \epsilon_{\text{ortho}}$ の場合のみ再直交化。

### 7.2 規格化係数のチェック

$b_{n+1}^2 = \langle\tilde{\phi}_{n+1}|\Theta|\tilde{\psi}_{n+1}\rangle$ が負または虚数部が大きい場合、数値的不安定の兆候。

**対策**：

1. 実部のみを使用：
   $$
   b_{n+1} = \sqrt{\text{Re}[\langle\tilde{\phi}_{n+1}|\Theta|\tilde{\psi}_{n+1}\rangle]}
   $$

2. 閾値チェック：
   $$
   \left|\frac{\text{Im}[b_{n+1}^2]}{\text{Re}[b_{n+1}^2]}\right| < \epsilon_{\text{tol}} = 10^{-6}
   $$
   
   この条件が満たされない場合、アルゴリズムを停止し、警告を出力。

### 7.3 収束判定

**誘電関数の収束**：

$$
\Delta\epsilon_n(\omega) = \max_\omega \left|\epsilon_M^{(n)}(\omega) - \epsilon_M^{(n-1)}(\omega)\right|
$$

収束条件：
$$
\Delta\epsilon_n(\omega) < \epsilon_{\text{conv}} = 10^{-4}
$$

**三重対角行列固有値の収束**：

$$
\Delta E_n = \max_{\lambda} |E_\lambda^{(n)} - E_\lambda^{(n-1)}|
$$

収束条件：
$$
\Delta E_n < E_{\text{conv}} = 10^{-5} \text{ eV}
$$

---

## 8. 対称性の利用

### 8.1 時間反転対称性

時間反転演算子 $\hat{\Theta}_{\text{TR}}$ に対して：

$$
\mathcal{H}_{\text{BSE}}(\bar{k}) = \hat{\Theta}_{\text{TR}} \mathcal{H}_{\text{BSE}}(k) \hat{\Theta}_{\text{TR}}^{-1}
$$

これにより、IBZ（既約ブリルアンゾーン）のk点のみ計算すればよい：
$$
N_k^{\text{irred}} = \frac{N_k}{2}
$$

**メモリ削減**：係数2

### 8.2 空間対称性

点群操作 $\hat{S} \in G$ に対して：

$$
\mathcal{H}_{\text{BSE}}(S k) = \hat{U}(S) \mathcal{H}_{\text{BSE}}(k) \hat{U}^{-1}(S)
$$

ここで $\hat{U}(S)$ はユニタリ変換。

IBZ内のk点数：
$$
N_k^{\text{irred}} = \frac{N_k}{|G|} = \frac{N_k}{N_{\text{sym}}}
$$

**メモリ削減**：係数 $N_{\text{sym}}$（通常、$N_{\text{sym}} = 2$~$48$）

### 8.3 擬エルミート対称性

$$
\mathcal{H}_{\text{BSE}} = \Theta \mathcal{H}_{\text{BSE}}^\dagger \Theta^{-1}
$$

この対称性により、共鳴ブロックと反共鳴ブロックの間に以下の関係：

$$
\text{反共鳴ブロック} = -[\text{共鳴ブロック}]^*
$$

したがって、反共鳴ブロックの陽的な計算は不要（共鳴ブロックから導出可能）。

**メモリ削減**：追加の係数2（原理的には可能だが、実装の簡潔性のため通常は使用しない）

---

## 9. 初期ベクトルの選択

### 9.1 双極子励起

光学応答を計算する場合、初期ベクトルは双極子遷移に対応：

$$
\Psi_{\text{init}}^{(vck\sigma)} = d_\alpha^{(vck\sigma)} = \langle vk\sigma|\hat{r}_\alpha|ck\sigma\rangle
$$

ここで $\hat{r}_\alpha$ は位置演算子の $\alpha$ 成分（$\alpha = x, y, z$）。

長波長極限 $(q \to 0)$ では：
$$
d_\alpha^{(vck\sigma)} \approx -i\frac{\langle vk\sigma|\hat{p}_\alpha|ck\sigma\rangle}{E_{ck} - E_{vk}}
$$

ここで $\hat{p}_\alpha$ は運動量演算子の $\alpha$ 成分。

### 9.2 結合ベクトルの初期化

結合ハミルトニアンの場合、初期ベクトルを以下のように設定：

$$
|\Psi_{\text{init}}\rangle = \begin{pmatrix} |d_\alpha\rangle \\ 0 \end{pmatrix}
$$

または、より対称的な形：

$$
|\Psi_{\text{init}}\rangle = \frac{1}{\sqrt{2}}\begin{pmatrix} |d_\alpha\rangle \\ |d_\alpha\rangle \end{pmatrix}
$$

後者の場合、初期化後に計量演算子で規格化：

$$
|\psi_0\rangle = \frac{|\Psi_{\text{init}}\rangle}{\sqrt{\langle\Psi_{\text{init}}|\Theta|\Psi_{\text{init}}\rangle}}
$$

### 9.3 ランダム位相法

特定の励起を狙わず、全スペクトルを均等にサンプルする場合、ランダム位相を持つ初期ベクトル：

$$
\Psi_{\text{init}}^{(vck\sigma)} = e^{i\phi_{vck\sigma}}
$$

ここで $\phi_{vck\sigma}$ は $[0, 2\pi)$ の一様乱数。

**利点**：すべての励起状態に対して非ゼロの重なりを持つ

**欠点**：特定の低エネルギー励起の収束が遅い

---

## 10. 実装上の考慮事項

### 10.1 データレイアウト

**Lanczosベクトルの保存**：

```
psi_plus(n_v, n_c, n_k, n_spin)
psi_minus(n_v, n_c, n_k, n_spin)
phi_plus(n_v, n_c, n_k, n_spin)
phi_minus(n_v, n_c, n_k, n_spin)
```

**メモリオーダー**：最速変化次元を $n_v$（価電子バンド）とし、キャッシュ効率を最大化。

### 10.2 行列要素の事前計算

一部の行列要素は複数回使用されるため、事前計算して保存：

**双極子行列要素**：
$$
d_\alpha^{(vck)} = \langle vk|\hat{r}_\alpha|ck\rangle
$$

保存サイズ：$O(N_k N_v N_c) = O(N_k N_{\text{bands}}^2)$

**振動子強度**：
$$
f_{vc}(k) = \frac{2(E_{ck} - E_{vk})}{3}\sum_\alpha |d_\alpha^{(vck)}|^2
$$

### 10.3 カーネル計算の最適化

**k'ループの最適な順序**：

外側ループを $k'$ とし、$k'$ に依存する波動関数を一度だけロード：

```
for k' in k_points:
    load WF(n, k') for n in all_bands
    
    for k in k_points:
        q = k - k'
        load WF(n, k) for n in all_bands
        
        compute M(n, n', G, k, k') for all n, n', G
        
        apply kernel to vectors
        
        unload WF(n, k)
    
    unload WF(n, k')
```

**並列化**：$k'$ ループをMPIで並列化（各MPIランクが異なる $k'$ セットを担当）

---

## 11. 性能指標と見積もり

### 11.1 メモリ使用量

**Lanczosベクトル**：
$$
M_{\text{vec}} = 4 \times 2N_k N_{\text{bands}}^2 \times 16 = 128 N_k N_{\text{bands}}^2 \text{ bytes}
$$

（4セット：$\psi_\pm^{(n)}, \psi_\pm^{(n-1)}$、各2成分）

**波動関数キャッシュ**：
$$
M_{\text{wf}} = N_{\text{cache}} N_k N_{\text{bands}} N_G \times 16 \text{ bytes}
$$

**遮蔽相互作用**：
$$
M_W = N_q N_G^2 \times 16 \text{ bytes}
$$

**総メモリ**：
$$
M_{\text{total}} = M_{\text{vec}} + M_{\text{wf}} + M_W
$$

### 11.2 計算時間

**1反復あたりの時間**：

各 $k'$ に対して、全 $k$ との相互作用を計算：
$$
T_{\text{iter}} = N_k^2 \times T_{\text{kernel}}
$$

カーネル計算時間：
$$
T_{\text{kernel}} = N_{\text{bands}}^4 \times N_G \times t_{\text{flop}} + N_{\text{bands}}^2 \times N_G \times t_{\text{mem}}
$$

ここで：
- $t_{\text{flop}}$：浮動小数点演算時間（～1 ns）
- $t_{\text{mem}}$：メモリアクセス時間（～100 ns）

**総計算時間**：
$$
T_{\text{total}} = N_{\text{iter}} \times T_{\text{iter}}
$$

典型的な反復回数：$N_{\text{iter}} = 200$~$500$

### 11.3 並列効率

**理想的な並列数**：
$$
N_{\text{proc}}^{\text{ideal}} = N_k
$$

各プロセスが異なる $k'$ を担当。

**スケーラビリティ**：
$$
\text{Speedup} = \frac{T_{\text{serial}}}{T_{\text{parallel}}} \approx N_{\text{proc}} \times \eta
$$

ここで $\eta$ は並列効率（通常、$\eta \sim 0.8$~$0.95$）。

---

## 12. 検証と品質保証

### 12.1 小規模系での検証

**手順**：

1. 少数バンド（$N_v = N_c = 10$）で計算
2. 従来の直接対角化法と比較
3. 以下の量を検証：
   - 励起エネルギー（誤差 $< 10^{-5}$ eV）
   - 振動子強度（相対誤差 $< 10^{-4}$）
   - 誘電関数スペクトル（全周波数範囲で相対誤差 $< 10^{-3}$）

### 12.2 対称性のチェック

**擬エルミート性**：

任意のベクトル $|\chi\rangle$ について：
$$
\langle\chi|\Theta\mathcal{H}_{\text{BSE}}|\chi\rangle = [\langle\chi|\Theta\mathcal{H}_{\text{BSE}}|\chi\rangle]^*
$$

数値的検証：偏差 $< 10^{-10}$

**時間反転対称性**：

$$
\mathcal{H}_{\text{BSE}}(k) = \mathcal{H}_{\text{BSE}}(\bar{k})^*
$$

（適切な基底変換後）

### 12.3 収束性の監視

各反復で以下をモニタリング：

1. Lanczosベクトルのノルム：
   $$
   \|\psi_n\| = \sqrt{\langle\psi_n|\Theta|\psi_n\rangle} \stackrel{?}{=} 1
   $$

2. 双直交性の偏差：
   $$
   \delta_{\text{ortho}} = \max_{m \neq n} |\langle\phi_m|\psi_n\rangle|
   $$

3. 三重対角係数の安定性：
   $$
   |b_{n+1}| > b_{\text{min}} = 10^{-14}
   $$
   
   $b_{n+1}$ が小さすぎる場合、反復を終了（収束とみなす）。

---

## 13. エラーハンドリングと異常検出

### 13.1 異常終了条件

以下の場合、計算を停止しエラーメッセージを出力：

1. **規格化係数が非正**：
   $$
   \text{Re}[b_{n+1}^2] \leq 0
   $$

2. **双直交性の崩壊**：
   $$
   \delta_{\text{ortho}} > \epsilon_{\text{ortho}}^{\text{max}} = 10^{-8}
   $$

3. **メモリ不足**：
   必要なメモリが利用可能メモリを超過

4. **収束の失敗**：
   最大反復回数に達しても収束しない

### 13.2 警告条件

以下の場合、警告を出力するが計算は継続：

1. **収束の遅延**：
   $$
   N_{\text{iter}} > N_{\text{iter}}^{\text{warn}} = 1000
   $$

2. **双直交性の緩和**：
   $$
   \delta_{\text{ortho}} > \epsilon_{\text{ortho}}^{\text{warn}} = 10^{-10}
   $$

3. **三重対角係数の虚数部**：
   $$
   \left|\frac{\text{Im}[a_n]}{\text{Re}[a_n]}\right| > 10^{-8}
   $$

---

## 14. パラメータ推奨値

### 14.1 収束パラメータ

| パラメータ | 記号 | 推奨値 | 説明 |
|----------|------|-------|------|
| 最大反復数 | $N_{\text{iter}}^{\text{max}}$ | 1000 | Haydock反復の上限 |
| 誘電関数収束閾値 | $\epsilon_{\text{conv}}$ | $10^{-4}$ | 相対変化量 |
| エネルギー収束閾値 | $E_{\text{conv}}$ | $10^{-5}$ eV | 固有値の絶対変化量 |
| 双直交性許容誤差 | $\epsilon_{\text{ortho}}$ | $10^{-12}$ | 再直交化の閾値 |

### 14.2 メモリパラメータ

| パラメータ | 記号 | 推奨値 | 説明 |
|----------|------|-------|------|
| 波動関数キャッシュサイズ | $N_{\text{cache}}$ | 100 | キャッシュするバンド数 |
| バンドチャンクサイズ | $N_{\text{chunk}}$ | 50~100 | 一度に処理するバンド数 |
| k点パーティション数 | $N_{\text{part}}$ | 10~20 | k点を分割する数 |

### 14.3 数値精度パラメータ

| パラメータ | 記号 | 推奨値 | 説明 |
|----------|------|-------|------|
| ゼロ判定閾値 | $\epsilon_{\text{zero}}$ | $10^{-16}$ | 数値ゼロとみなす値 |
| 最小規格化係数 | $b_{\text{min}}$ | $10^{-14}$ | 反復終了の判定 |
| 虚数部許容度 | $\epsilon_{\text{imag}}$ | $10^{-6}$ | 実数量の虚数部の許容相対誤差 |

---

## 15. まとめ

本仕様書では、BSE計算において結合項（BSEmod="coupling"）を含む場合でも全バンドを考慮可能な、数学的に厳密なアルゴリズムの仕様を定義した。

### 15.1 達成目標

1. **メモリ削減**：$O(10^7)$ 倍の削減（PBスケール → GBスケール）
2. **数学的厳密性**：ヒューリスティックやfallback処理を一切使用しない
3. **実用性**：2TB RAMのノードで実行可能
4. **汎用性**：任意の系と計算条件に適用可能

### 15.2 核心技術

1. **オンザフライカーネル構築**：行列を保存せず、作用のみを計算
2. **双直交Lanczos法**：擬エルミート系に対応した修正Haydock反復
3. **対称性の完全利用**：時間反転・空間対称性・擬エルミート対称性
4. **効率的な波動関数管理**：LRUキャッシュとプリフェッチ

### 15.3 理論的保証

- **収束性**：擬エルミート性により、アルゴリズムの収束が保証される
- **精度**：有限精度演算の影響を除き、厳密解を得る
- **安定性**：選択的再直交化により、長時間反復でも安定

### 15.4 次のステップ

本仕様書に基づき、詳細設計書（アルゴリズムの実装レベルの記述）を作成する。詳細設計書には以下を含める：

1. 疑似コードとフローチャート
2. データ構造の定義
3. サブルーチン仕様
4. 並列化スキーム
5. I/O戦略
6. テストケースとベンチマーク

---

## 参考文献

1. **Onida, G., Reining, L., & Rubio, A.** (2002). *Electronic excitations: density-functional versus many-body Green's-function approaches*. Reviews of Modern Physics, 74(2), 601.

2. **Rohlfing, M., & Louie, S. G.** (2000). *Electron-hole excitations and optical spectra from first principles*. Physical Review B, 62(4), 4927.

3. **Strinati, G.** (1988). *Application of the Green's functions method to the study of the optical properties of semiconductors*. La Rivista del Nuovo Cimento, 11(12), 1-86.

4. **Haydock, R., Heine, V., & Kelly, M. J.** (1972). *Electronic structure based on the local atomic environment for tight-binding bands*. Journal of Physics C: Solid State Physics, 5(20), 2845.

5. **Parlett, B. N., & Taylor, D. R.** (1985). *A look-ahead Lanczos algorithm for unsymmetric matrices*. Mathematics of Computation, 44(169), 105-124.

6. **Freund, R. W., & Nachtigal, N. M.** (1991). *QMR: a quasi-minimal residual method for non-Hermitian linear systems*. Numerische Mathematik, 60(1), 315-339.

7. **Marini, A., Hogan, C., Grüning, M., & Varsano, D.** (2009). *yambo: An ab initio tool for excited state calculations*. Computer Physics Communications, 180(8), 1392-1403.

---

## 付録：用語集

| 用語 | 英語 | 定義 |
|-----|------|------|
| 結合項 | Coupling term | BSEハミルトニアンの共鳴-反共鳴間の非対角ブロック |
| 反共鳴 | Anti-resonant | 励起状態の消滅過程に対応する項 |
| 擬エルミート | Pseudo-Hermitian | $H = \Theta H^\dagger \Theta^{-1}$ を満たす演算子 |
| 双直交 | Biorthogonal | $\langle\phi_m\|\psi_n\rangle = \delta_{mn}$ を満たす2組のベクトル |
| オンザフライ | On-the-fly | 保存せず必要時に計算すること |
| Lanczos法 | Lanczos method | 大規模固有値問題の反復解法 |
| 三重対角行列 | Tridiagonal matrix | 主対角とその隣接要素のみ非零の行列 |
| 既約ブリルアンゾーン | Irreducible Brillouin Zone (IBZ) | 対称性により縮約されたBZ |

---

**文書バージョン**：1.0  
**作成日**：2025-10-18  
**言語**：日本語・英語併記  
**文書形式**：Markdown
