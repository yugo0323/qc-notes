# 第 1 回　量子化学で解きたい方程式：Born–Oppenheimer 近似と電子のハミルトニアン

!!! abstract "この回の道筋"
    量子化学の計算プログラムが解いているのは、分子全体のシュレーディンガー方程式ではありません。核を止めたときの「電子のシュレーディンガー方程式」です。この回では、その方程式がどこから来るのか（Born–Oppenheimer 近似）と、その答えから何が分かるのか（ポテンシャルエネルギー曲面）を説明します。

    | 節 | 内容 |
    |---|---|
    | 1 | 分子のハミルトニアン：5 つの項 |
    | 2 | 原子単位：式から定数を消す |
    | 3 | Born–Oppenheimer 近似：核と電子を分ける |
    | 4 | ポテンシャルエネルギー曲面：H₂ の計算例 |
    | 5 | Born–Oppenheimer 近似が悪くなるとき |

    **前に読んでおく回**：なし（この連載の最初の回です）

## 1. 分子のハミルトニアン：5 つの項

分子は、いくつかの原子核と電子でできています。
量子力学では、分子の状態は波動関数 $\Psi$ で表され、時間によらない状態はシュレーディンガー方程式

$$
\hat H \Psi = E \Psi
$$

を満たします。
$\hat H$ はハミルトニアン（Hamiltonian）と呼ばれ、系の全エネルギーを表す演算子です。
$E$ はその状態のエネルギーです。

分子のハミルトニアンを書くために、記号を決めます。

- 電子に $i, j$ の番号を付け、位置を $\mathbf r_i$ とします。電子の数を $N$ とします。
- 原子核に $A, B$ の番号を付け、位置を $\mathbf R_A$、質量を $M_A$、原子番号を $Z_A$ とします。
- $r_{iA} = |\mathbf r_i - \mathbf R_A|$、$r_{ij} = |\mathbf r_i - \mathbf r_j|$、$R_{AB} = |\mathbf R_A - \mathbf R_B|$ と書きます。
- 電子の質量を $m_\mathrm e$、電気素量を $e$、真空の誘電率を $\varepsilon_0$ とします。

相対論の効果と、スピンが関わる磁気的な相互作用を無視すると、分子のハミルトニアンは次の 5 つの項の和です。

$$
\hat H =
\underbrace{-\sum_{i}\frac{\hbar^2}{2m_\mathrm e}\nabla_i^2}_{\hat T_\mathrm e}
\underbrace{-\sum_{A}\frac{\hbar^2}{2M_A}\nabla_A^2}_{\hat T_\mathrm N}
\underbrace{-\sum_{i}\sum_{A}\frac{Z_A e^2}{4\pi\varepsilon_0\, r_{iA}}}_{\hat V_\mathrm{eN}}
\underbrace{+\sum_{i<j}\frac{e^2}{4\pi\varepsilon_0\, r_{ij}}}_{\hat V_\mathrm{ee}}
\underbrace{+\sum_{A<B}\frac{Z_A Z_B e^2}{4\pi\varepsilon_0\, R_{AB}}}_{\hat V_\mathrm{NN}}
$$

各項の意味は次のとおりです。

| 項 | 意味 | 符号 |
|---|---|---|
| $\hat T_\mathrm e$ | 電子の運動エネルギー | 正 |
| $\hat T_\mathrm N$ | 核の運動エネルギー | 正 |
| $\hat V_\mathrm{eN}$ | 電子と核の引き合い（クーロン引力） | 負 |
| $\hat V_\mathrm{ee}$ | 電子どうしの反発 | 正 |
| $\hat V_\mathrm{NN}$ | 核どうしの反発 | 正 |

$\nabla_i^2$ は電子 $i$ の座標についてのラプラシアン（2 階微分の和）、$\nabla_A^2$ は核 $A$ の座標についてのラプラシアンです。
和の $i<j$ は、同じ組を 2 回数えないための書き方です。

この方程式は、水素原子のように電子が 1 個で核が 1 個なら厳密に解けます。
しかし電子が 2 個以上になると、$\hat V_\mathrm{ee}$ がすべての電子の座標を結び付けるため、厳密には解けません。
さらに核も動くので、変数の数は $3N + 3M$（$M$ は核の数）に増えます。
そこで、まず変数を減らす工夫が必要になります。それが第 3 節の Born–Oppenheimer 近似です。
その前に、式を見やすくするために単位を決めます。

**ここまでで分かったこと**：分子のハミルトニアンは、電子と核の運動エネルギー、電子と核の引力、電子どうしと核どうしの反発の 5 項でできています。

## 2. 原子単位：式から定数を消す

上の式には $\hbar$、$m_\mathrm e$、$e$、$4\pi\varepsilon_0$ が何度も出てきます。
これらをすべて 1 とする単位系を、原子単位（atomic units）と呼びます。
量子化学の論文やプログラムの出力は、ほとんどがこの単位で書かれています。

原子単位では、長さとエネルギーの単位が次のように決まります。

| 量 | 原子単位での名前 | SI などでの値 |
|---|---|---|
| 長さ | ボーア（bohr）、$a_0$ | 0.529177 Å |
| エネルギー | ハートリー（hartree）、$E_\mathrm h$ | 27.211 eV ＝ 627.5 kcal/mol ＝ 2625.5 kJ/mol |
| 質量 | 電子の質量 $m_\mathrm e$ | — |
| 電荷 | 電気素量 $e$ | — |

数値は CODATA の推奨値 [Tiesinga 2021] によります。
たとえば水素原子の基底状態のエネルギーは、原子単位でちょうど $-\tfrac12\,E_\mathrm h$ です（核の質量を無限大とした場合）。

原子単位では、分子のハミルトニアンは次のように短くなります。

$$
\hat H =
-\sum_{i}\frac12\nabla_i^2
-\sum_{A}\frac{1}{2M_A}\nabla_A^2
-\sum_{i}\sum_{A}\frac{Z_A}{r_{iA}}
+\sum_{i<j}\frac{1}{r_{ij}}
+\sum_{A<B}\frac{Z_A Z_B}{R_{AB}}
$$

ここで $M_A$ は、電子の質量を 1 とした核の質量です。
陽子の質量は電子の約 1836 倍です [Tiesinga 2021]。
したがって、いちばん軽い水素の核でも $1/(2M_A)$ は $1/(2\cdot 1836)$ ほどしかありません。
この「核はとても重い」という事実が、次の節の出発点です。

**ここまでで分かったこと**：原子単位では $\hbar = m_\mathrm e = e = 4\pi\varepsilon_0 = 1$ とし、エネルギーをハートリー、長さをボーアで測ります。この連載でも、断りがなければ原子単位を使います。

## 3. Born–Oppenheimer 近似：核と電子を分ける

### 3.1 考え方

核は電子より 1800 倍以上重いので、同じ程度の力を受けても、核はずっとゆっくり動きます。
電子から見ると、核はほとんど止まっているように見えます。
そこで、次の 2 段階に分けて問題を解きます。

1. 核の位置 $\mathbf R = (\mathbf R_1, \mathbf R_2, \dots)$ を固定し、その中で電子の状態を求める。
2. 求めた電子のエネルギーを、核が感じるポテンシャルとして使い、核の運動を考える。

この分け方を Born–Oppenheimer 近似と呼びます。
Born と Oppenheimer は 1927 年に、電子と核の質量比の 4 乗根 $(m_\mathrm e/M)^{1/4}$ を小さなパラメータとして分子のエネルギーを展開し、この分け方の根拠を示しました [Born 1927]。
教科書では、Szabo–Ostlund の 2.1.2 節に分かりやすい説明があります [Szabo 1996]。

### 3.2 電子のシュレーディンガー方程式

核を止めると、核の運動エネルギー $\hat T_\mathrm N$ を考えなくてよくなります。
残りの項のうち、電子の座標を含むものを集めたのが、電子のハミルトニアンです。

$$
\hat H_\mathrm{el}
= -\sum_{i}\frac12\nabla_i^2
-\sum_{i}\sum_{A}\frac{Z_A}{r_{iA}}
+\sum_{i<j}\frac{1}{r_{ij}}
$$

核どうしの反発 $V_\mathrm{NN} = \sum_{A<B} Z_A Z_B / R_{AB}$ は、核を止めればただの数です。
そこで、電子の方程式を次のように書きます。

$$
\hat H_\mathrm{el}\,\psi_k(\mathbf r;\mathbf R) = E_k^\mathrm{el}(\mathbf R)\,\psi_k(\mathbf r;\mathbf R)
$$

$\mathbf r$ はすべての電子の座標をまとめたもの、$k = 0, 1, 2, \dots$ は電子の状態の番号（0 が基底状態）です。
「$;\mathbf R$」は、核の位置が変数ではなく、外から与える値（パラメータ）であることを表します。
核の位置を変えるたびに、この方程式を解き直します。

量子化学の計算プログラムが解いているのは、この電子のシュレーディンガー方程式です。
プログラムが「全エネルギー」として出力するのは、電子のエネルギーに核の反発を足した値

$$
E_k(\mathbf R) = E_k^\mathrm{el}(\mathbf R) + V_\mathrm{NN}(\mathbf R)
$$

です。
この連載の第 3 回以降で扱う手法（Hartree–Fock 法、電子相関の方法、密度汎関数理論）は、どれもこの方程式を近似的に解く方法です。

### 3.3 核の運動

次に、核の運動を考えます。
分子全体の波動関数を、電子の状態 $\psi_k$ と核の波動関数 $\chi(\mathbf R)$ の積

$$
\Psi(\mathbf r, \mathbf R) \approx \chi(\mathbf R)\,\psi_k(\mathbf r;\mathbf R)
$$

で近似します。
導出は付録 A に回しますが、核の波動関数は次の方程式を満たします。

$$
\left[\,-\sum_A \frac{1}{2M_A}\nabla_A^2 + E_k(\mathbf R)\right]\chi(\mathbf R) = E\,\chi(\mathbf R)
$$

この式は、ポテンシャル $E_k(\mathbf R)$ の中を核が動く、というシュレーディンガー方程式の形をしています。
つまり、電子のエネルギー $E_k(\mathbf R)$ が、核にとってのポテンシャルの役目をします。
分子の振動や化学反応は、このポテンシャルの上での核の運動として理解できます。

付録 A で示すように、この近似では、核の運動によって電子の状態が $\psi_k$ から別の $\psi_l$ に移る効果（非断熱結合）を捨てています。
捨てた項には $1/M_A$ がかかっているため、核が重いほど近似がよくなります。

**ここまでで分かったこと**：核を止めて電子の方程式を解き、得られたエネルギー $E_k(\mathbf R)$ を核のポテンシャルとして使います。これが Born–Oppenheimer 近似で、量子化学の計算の土台です。

## 4. ポテンシャルエネルギー曲面：H₂ の計算例

核の位置 $\mathbf R$ の関数としての $E_k(\mathbf R)$ を、ポテンシャルエネルギー曲面（potential energy surface、PES）と呼びます。
核が $M$ 個ある分子では、並進と回転の 6 個（直線分子では 5 個）を除いた $3M-6$ 個の変数の関数です。
変数が 1 個の場合は、ポテンシャルエネルギー曲線と呼ぶこともあります。

ポテンシャルエネルギー曲面から、化学で使う多くの量が分かります。

- **平衡構造**：曲面の極小点。
- **振動数**：極小点のまわりの曲がり具合（2 階微分）。
- **反応のエネルギー**：反応物と生成物の極小点のエネルギー差。
- **活性化エネルギー**：反応物の極小点と、反応経路の上でいちばん高い点（遷移状態、1 次の鞍点）のエネルギー差。

いちばん簡単な例として、H₂ 分子のポテンシャルエネルギー曲線を計算してみます。
変数は H–H の距離 $R$ だけです。
距離ごとに電子の方程式を解きます。
解き方には、使った基底関数（第 6 回で説明します）の範囲で厳密な解を与える完全 CI（full CI、[第 8 回](08-configuration-interaction.md)で説明します）を使います。
計算には PySCF [Sun 2020] を使い、基底関数は cc-pVDZ [Dunning 1989] とします。

```python
import numpy as np
from pyscf import gto, scf, fci

def e_h2(r):
    # r は H–H 距離（Å）。全エネルギー E(R)（ハートリー）を返す
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {r}", basis="cc-pvdz", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    return fci.FCI(mf).kernel()[0]

R = np.round(np.arange(0.40, 4.001, 0.02), 2)
E = np.array([e_h2(r) for r in R])
i = np.argmin(E)

# 遠く離れた 2 個の H 原子のエネルギー
h = gto.M(atom="H 0 0 0", basis="cc-pvdz", spin=1, verbose=0)
e_h = scf.UHF(h).run().e_tot
De = 2 * e_h - E[i]
print(f"最小の点      R = {R[i]:.2f} Å   E = {E[i]:.5f} Eh")
print(f"H 原子 2 個    2E(H) = {2*e_h:.5f} Eh")
print(f"解離エネルギー De = {De:.5f} Eh = {De*27.211386:.2f} eV")
```

出力は次のとおりです。

```text
最小の点      R = 0.76 Å   E = -1.16367 Eh
H 原子 2 個    2E(H) = -0.99856 Eh
解離エネルギー De = 0.16512 Eh = 4.49 eV
```

曲線を図にすると、次のようになります。

![H₂ のポテンシャルエネルギー曲線（FCI/cc-pVDZ）。横軸は H–H 距離、縦軸は全エネルギー。破線は 2 個の H 原子に分かれたときのエネルギー。](img/01-h2-pes.png)

この図から、次のことが読み取れます。

- 距離が短いと、核どうしの反発 $V_\mathrm{NN} = 1/R$ が大きくなり、エネルギーが急に上がります。
- 距離が長いと、エネルギーは 2 個の H 原子のエネルギーの和に近づきます。
- その間に極小があり、これが H₂ の平衡距離（この計算では 0.76 Å）です。極小の深さ $D_\mathrm e$ が、結合を切るのに必要なエネルギー（零点振動を含まない値）です。

H 原子のエネルギーは厳密には $-0.5\,E_\mathrm h$ なので、2 個で $-1\,E_\mathrm h$ のはずです。
上の計算では $-0.99856\,E_\mathrm h$ と、少し高く出ました。
この差は、電子の方程式の解き方ではなく、基底関数が有限であることから来ます。
基底関数の選び方は第 6 回で扱います。

**ここまでで分かったこと**：核の位置を変えながら電子の方程式を解くと、ポテンシャルエネルギー曲面が得られます。平衡構造、振動、反応のエネルギーは、この曲面から読み取ります。

## 5. Born–Oppenheimer 近似が悪くなるとき

第 3 節で捨てた非断熱結合は、付録 A の式 (A.4) のように、電子の状態 $\psi_k$ を核の座標で微分した量を含みます。
付録 A で示すように、この結合は 2 つの電子状態のエネルギー差 $E_l - E_k$ に反比例します。差が小さいほど、結合は大きくなります。
したがって、Born–Oppenheimer 近似が悪くなるのは、主に次の場合です。

- 2 つの電子状態のエネルギーが近づく、または交わる場所（円錐交差など）。光化学反応や、励起状態からの緩和がこれに当たります。
- 核が軽く、速く動く場合（水素の移動など）。ただし多くの場合、近似はまだよく成り立ちます。

基底状態の平衡構造の近くで、励起状態とのエネルギー差が十分大きければ、Born–Oppenheimer 近似はたいへんよい近似です。
この連載の第 I〜IV 部は、すべてこの近似の中で話を進めます。
近似が破れる場合は、[第 22 回](22-conical-intersections.md)で扱います。

!!! note "「全エネルギー」の意味に注意"
    計算プログラムが出す全エネルギーは、静止した核の位置での電子のエネルギーと核の反発の和です。零点振動エネルギーや熱の寄与は含みません。反応熱などを実験と比べるときは、ポテンシャルエネルギー曲面の上で振動を別に計算して足します。

**ここまでで分かったこと**：Born–Oppenheimer 近似は、電子状態どうしのエネルギーが近いところで悪くなります。それ以外では、量子化学の計算の土台として信頼できます。

## まとめ

- 分子のハミルトニアンは、電子と核の運動エネルギー、電子と核の引力、電子どうしと核どうしの反発の 5 項でできています。
- 原子単位では $\hbar = m_\mathrm e = e = 4\pi\varepsilon_0 = 1$ とし、1 $E_\mathrm h$ = 27.211 eV = 627.5 kcal/mol です。
- Born–Oppenheimer 近似では、核を止めて電子のシュレーディンガー方程式を解きます。量子化学の計算プログラムが解いているのはこの方程式です。
- 電子のエネルギーを核の位置の関数として見たものがポテンシャルエネルギー曲面で、構造・振動・反応はここから読み取ります。
- 電子状態のエネルギーが近づくところでは、この近似が悪くなります。

## 付録 A　核の方程式と非断熱結合の導出

分子全体のシュレーディンガー方程式 $\hat H\Psi = E\Psi$ から出発します。
ハミルトニアンを $\hat H = \hat T_\mathrm N + \hat H_\mathrm{el} + V_\mathrm{NN}$ と分けます。
各 $\mathbf R$ で、電子の状態 $\{\psi_k(\mathbf r;\mathbf R)\}$ は規格直交な完全系になるとします。
すると、全体の波動関数はいつでも

$$
\Psi(\mathbf r,\mathbf R) = \sum_l \chi_l(\mathbf R)\,\psi_l(\mathbf r;\mathbf R) \tag{A.1}
$$

と展開できます（Born–Huang の展開 [Born 1954]）。
これを $\hat H\Psi = E\Psi$ に入れます。
$\hat H_\mathrm{el} + V_\mathrm{NN}$ は $\psi_l$ にかかると $E_l(\mathbf R)$ を出します。
$\hat T_\mathrm N$ は $\chi_l$ と $\psi_l$ の両方に微分がかかるので、積の微分

$$
\nabla_A^2(\chi_l\psi_l) = (\nabla_A^2\chi_l)\psi_l + 2(\nabla_A\chi_l)\cdot(\nabla_A\psi_l) + \chi_l\nabla_A^2\psi_l \tag{A.2}
$$

を使います。
左から $\psi_k^*$ をかけて電子の座標で積分し、$\langle\psi_k|\psi_l\rangle = \delta_{kl}$ を使うと、次の式が得られます。

$$
\left[\hat T_\mathrm N + E_k(\mathbf R)\right]\chi_k + \sum_l \Lambda_{kl}\,\chi_l = E\,\chi_k \tag{A.3}
$$

$$
\Lambda_{kl} = -\sum_A \frac{1}{2M_A}\Big[\,2\,\langle\psi_k|\nabla_A\psi_l\rangle\cdot\nabla_A + \langle\psi_k|\nabla_A^2\psi_l\rangle\Big] \tag{A.4}
$$

ここで $\langle\cdots\rangle$ は電子の座標だけについての積分です。
ここまでは近似を使っていません。

$\Lambda_{kl}$ には $1/M_A$ がかかっていて、核が重いほど小さくなります。
$\Lambda_{kl}$ をすべて捨てると、式 (A.3) は異なる $k$ の間で結び付かなくなり、本文の核の方程式になります。
これが Born–Oppenheimer 近似です。
対角項 $\Lambda_{kk}$ だけを残す近似は、断熱近似と呼ばれ、$\Lambda_{kk}$ の寄与は断熱補正（diagonal Born–Oppenheimer correction）と呼ばれます。

最後に、非対角の結合の大きさを見積もります。
$\hat H_\mathrm{el}\psi_l = E_l^\mathrm{el}\psi_l$ の両辺を $\mathbf R_A$ で微分し、左から $\psi_k^*$（$k \neq l$）をかけて積分すると、

$$
\langle\psi_k|\nabla_A\psi_l\rangle = \frac{\langle\psi_k|(\nabla_A\hat H_\mathrm{el})|\psi_l\rangle}{E_l^\mathrm{el} - E_k^\mathrm{el}} \qquad (k\neq l) \tag{A.5}
$$

が得られます。
分母がエネルギー差なので、2 つの状態が近づくと結合が大きくなります。
本文の第 5 節は、この式にもとづいています。

## 参考文献

- M. Born, R. Oppenheimer, "Zur Quantentheorie der Molekeln", Ann. Phys. 389, 457 (1927). DOI: [10.1002/andp.19273892002](https://doi.org/10.1002/andp.19273892002)
- A. Szabo, N. S. Ostlund, "Modern Quantum Chemistry" (Dover, 1996). ISBN: 9780486691862
- M. Born, K. Huang, "Dynamical Theory of Crystal Lattices" (Oxford University Press, 1954). ISBN: 9780198503699
- E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor, "CODATA recommended values of the fundamental physical constants: 2018", Rev. Mod. Phys. 93, 025010 (2021). DOI: [10.1103/RevModPhys.93.025010](https://doi.org/10.1103/RevModPhys.93.025010)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
- T. H. Dunning Jr., "Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through neon and hydrogen", J. Chem. Phys. 90, 1007 (1989). DOI: [10.1063/1.456153](https://doi.org/10.1063/1.456153)
