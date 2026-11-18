# 第 3 回　Hartree–Fock 法のおさらい

!!! abstract "この回の道筋"
    Hartree–Fock 法は、波動関数を 1 個の Slater 決定子に限ったときに、いちばんエネルギーが低くなる軌道を求める方法です。この回では、変分原理から Hartree–Fock 方程式を導き、基底関数で展開した Roothaan–Hall 方程式を、自己無撞着場（SCF）の繰り返しで解く手順を確かめます。最後に、軌道エネルギーの意味（Koopmans の定理）を、H₂O の計算で確かめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 変分原理：エネルギーの上限 |
    | 2 | 決定子 1 個のエネルギー |
    | 3 | Fock 演算子：クーロンと交換 |
    | 4 | Roothaan–Hall 方程式 |
    | 5 | SCF の手順：H₂O を自分で解く |
    | 6 | 軌道エネルギーと Koopmans の定理 |

    **前に読んでおく回**：量子化学の理論 [第 1 回](01-born-oppenheimer.md)、[第 2 回](02-determinants.md)

!!! tip "物理的なイメージ"
    Hartree–Fock 法が求めているのは、各電子が「ほかのすべての電子がつくる平均の場」（平均のクーロン反発と、交換の効果）の中を独立に動くとしたときの、最もよい軌道です。ところが、その平均の場は、求めたい軌道そのものから決まります。そこで、軌道を仮に決める、平均の場をつくる、その中で軌道を解き直す、という手順を、軌道が変わらなくなるまで繰り返します。これが自己無撞着場（SCF）の考え方です。電子どうしが瞬間ごとに避け合う動きは、平均にならされて失われます。これが第 7 回で扱う電子相関です。

## 1. 変分原理：エネルギーの上限

第 1 回で見たように、量子化学では電子のシュレーディンガー方程式 $\hat H\Psi = E\Psi$ を解きます。
電子が 2 個以上あると、これは厳密には解けません。
そこで、近似の波動関数を作り、その「よさ」をエネルギーで測ります。
その根拠になるのが変分原理です。

**変分原理**：規格化された任意の波動関数 $\tilde\Psi$ について、

$$
\tilde E = \langle\tilde\Psi|\hat H|\tilde\Psi\rangle \ \ge\ E_0
$$

が成り立ちます。$E_0$ は基底状態の厳密なエネルギーです。

理由は次のとおりです。
$\hat H$ の固有関数 $\Psi_n$（エネルギー $E_n \ge E_0$）は完全系をなすので、$\tilde\Psi = \sum_n c_n\Psi_n$ と展開できます。
すると $\tilde E = \sum_n |c_n|^2 E_n \ge E_0\sum_n|c_n|^2 = E_0$ となります。

したがって、波動関数にいくつかの調節できる部分（パラメータ）を持たせ、$\tilde E$ がいちばん低くなるように選べば、その範囲で最もよい近似が得られます。
Hartree–Fock 法では、波動関数を 1 個の Slater 決定子に限り、その中の軌道を調節します。

**ここまでで分かったこと**：近似の波動関数のエネルギーは、厳密なエネルギーより必ず高くなります。エネルギーが低いほど、よい近似です。

## 2. 決定子 1 個のエネルギー

以下では、閉殻の分子を考えます。
閉殻とは、$N$ 個の電子が $N/2$ 個の空間軌道 $\phi_1, \dots, \phi_{N/2}$ に、$\alpha$ と $\beta$ の対で入っている状態です。
各空間軌道に同じ形の $\alpha$ と $\beta$ の軌道を使う方法を、制限 Hartree–Fock（restricted Hartree–Fock、RHF）と呼びます。
空間軌道は実数で、規格直交であるとします。

第 2 回の付録で 2 電子の決定子のエネルギーを求めました。
同じ計算を $N$ 電子に広げると、閉殻の決定子 $\Phi_0$ のエネルギーは次のようになります（原子単位）。

$$
E_\mathrm{HF} = \langle\Phi_0|\hat H_\mathrm{el}|\Phi_0\rangle
= 2\sum_{i}^{N/2} h_{ii} + \sum_{i}^{N/2}\sum_{j}^{N/2}\big(2J_{ij} - K_{ij}\big)
$$

記号は第 2 回と同じです。

- $h_{ii} = \langle\phi_i|\hat h|\phi_i\rangle$。$\hat h = -\frac12\nabla^2 - \sum_A Z_A/r_A$ は、電子 1 個の運動エネルギーと核からの引力です。
- $J_{ij} = (ii|jj)$ はクーロン積分、$K_{ij} = (ij|ji)$ は交換積分です。$(ij|kl)$ は化学者の表記です。

各項の数え方は次のとおりです。
1 電子の項は、各空間軌道に 2 個の電子が入るので 2 倍です。
2 電子の項では、軌道 $i$ と $j$ の電子の組が 4 組（$\alpha\alpha$、$\alpha\beta$、$\beta\alpha$、$\beta\beta$）あります。
4 組すべてにクーロンの反発があり、交換は同じスピンの 2 組だけに出ます。
$i, j$ の両方について和をとると組を 2 回数えるので、全体を半分にして $2J_{ij} - K_{ij}$ になります。
$i = j$ の項は $2J_{ii} - K_{ii} = J_{ii}$ で、同じ軌道に入った $\alpha$ と $\beta$ の電子の反発です。

**ここまでで分かったこと**：閉殻の決定子のエネルギーは、1 電子の積分 $h_{ii}$ と、クーロン積分 $J_{ij}$・交換積分 $K_{ij}$ で書けます。

## 3. Fock 演算子：クーロンと交換

### 3.1 Hartree–Fock 方程式

前の節の $E_\mathrm{HF}$ を、軌道が規格直交であるという条件のもとで最小にします。
導出は付録 A に回しますが、最小になる軌道は次の方程式を満たします。

$$
\hat f\,\phi_i = \varepsilon_i\,\phi_i
$$

これを Hartree–Fock 方程式と呼びます。
$\hat f$ は Fock 演算子（Fock operator）、$\varepsilon_i$ は軌道エネルギー（orbital energy）です。
閉殻の場合、Fock 演算子は次の形です。

$$
\hat f = \hat h + \sum_{j}^{N/2}\big(2\hat J_j - \hat K_j\big)
$$

$\hat J_j$ はクーロン演算子、$\hat K_j$ は交換演算子で、任意の関数 $\phi$ に次のように作用します。

$$
\hat J_j\,\phi(\mathbf r_1) = \left[\int \frac{|\phi_j(\mathbf r_2)|^2}{r_{12}}\,d\mathbf r_2\right]\phi(\mathbf r_1)
$$

$$
\hat K_j\,\phi(\mathbf r_1) = \left[\int \frac{\phi_j^*(\mathbf r_2)\,\phi(\mathbf r_2)}{r_{12}}\,d\mathbf r_2\right]\phi_j(\mathbf r_1)
$$

### 3.2 方程式の意味

クーロン演算子 $\hat J_j$ は、軌道 $\phi_j$ の電子の雲がつくる静電ポテンシャルです。
したがって $\sum_j 2\hat J_j$ は、すべての電子の雲がつくる平均の静電ポテンシャルを表します。
Hartree–Fock 法では、各電子が、ほかの電子の瞬間の位置ではなく、この平均のポテンシャルの中を動くと考えます。
このため、Hartree–Fock 法は平均場（mean-field）の方法と呼ばれます。

交換演算子 $\hat K_j$ は、作用する関数 $\phi$ を積分の中に入れてしまうので、ある点でのポテンシャルの値として書けません（非局所的な演算子です）。
これは第 2 回で見た交換の効果、つまり同じスピンの電子が避け合う効果を表します。

和の中には $j = i$ の項も含まれます。
これは、電子が自分自身の雲と反発しているように見えます。
しかし $\langle\phi_i|\hat J_i|\phi_i\rangle = J_{ii} = K_{ii} = \langle\phi_i|\hat K_i|\phi_i\rangle$ なので、同じスピンの自分自身との反発は、交換の項でちょうど打ち消されます。
Hartree–Fock 法には、このような自己相互作用の誤差がありません。

### 3.3 自己無撞着場

Fock 演算子は、占有軌道 $\phi_j$ を使って作られています。
つまり、方程式の答え（軌道）が、方程式そのものに入っています。
このため、Hartree–Fock 方程式は一度では解けず、次の繰り返しで解きます。

1. 軌道を推測する。
2. その軌道で Fock 演算子を作る。
3. Hartree–Fock 方程式を解いて、新しい軌道を得る。
4. 軌道が変わらなくなるまで 2〜3 を繰り返す。

軌道がつくる場と、その場の中の軌道が一致した状態を、自己無撞着場（self-consistent field、SCF）と呼びます。
この考え方は Hartree が 1928 年に導入し [Hartree 1928]、Fock と Slater が 1930 年に反対称化した波動関数（決定子）に広げました [Fock 1930; Slater 1930]。

**ここまでで分かったこと**：Hartree–Fock の軌道は、Fock 演算子の固有関数です。Fock 演算子は占有軌道に依存するので、繰り返し計算（SCF）で解きます。

## 4. Roothaan–Hall 方程式

Hartree–Fock 方程式は微分と積分を含む方程式で、分子ではそのままでは解けません。
そこで、軌道を、あらかじめ決めた $K$ 個の関数 $g_\mu(\mathbf r)$ の一次結合で表します。

$$
\phi_i(\mathbf r) = \sum_{\mu=1}^{K} C_{\mu i}\,g_\mu(\mathbf r)
$$

$g_\mu$ を基底関数（basis function）と呼び、ふつうは各原子を中心とする原子軌道のような関数を使います（第 6 回で詳しく扱います）。
これは、学部で学んだ LCAO（原子軌道の一次結合）の考え方そのものです。

この展開を Hartree–Fock 方程式に入れ、左から $g_\mu^*$ をかけて積分すると、係数 $C_{\mu i}$ についての行列の方程式

$$
\mathbf F\mathbf C = \mathbf S\mathbf C\boldsymbol\varepsilon
$$

になります。
これを Roothaan–Hall 方程式と呼びます。
1951 年に Roothaan と Hall がそれぞれ示しました [Roothaan 1951; Hall 1951]。
各記号は $K\times K$ の行列です。

- $\mathbf S$：重なり行列。$S_{\mu\nu} = \int g_\mu^*\,g_\nu\,d\mathbf r$。基底関数は互いに直交しないので、単位行列になりません。
- $\mathbf C$：各列が 1 つの軌道の係数。
- $\boldsymbol\varepsilon$：軌道エネルギーを並べた対角行列。
- $\mathbf F$：Fock 行列。$F_{\mu\nu} = \int g_\mu^*\,\hat f\,g_\nu\,d\mathbf r$。

基底関数を実数とすると、閉殻の Fock 行列は次のように書けます。

$$
F_{\mu\nu} = h_{\mu\nu} + \sum_{\lambda\sigma} D_{\lambda\sigma}\Big[(\mu\nu|\lambda\sigma) - \tfrac12(\mu\lambda|\nu\sigma)\Big],
\qquad
D_{\lambda\sigma} = 2\sum_{i}^{N/2} C_{\lambda i}C_{\sigma i}
$$

$D_{\lambda\sigma}$ は密度行列で、電子密度を $\rho(\mathbf r) = \sum_{\lambda\sigma}D_{\lambda\sigma}g_\lambda(\mathbf r)g_\sigma(\mathbf r)$ と表します。
$(\mu\nu|\lambda\sigma)$ は、基底関数についての電子反発積分（化学者の表記）です。
括弧の第 1 項がクーロン、第 2 項が交換の寄与です。
エネルギーは、同じ記号で次のように書けます。

$$
E_\mathrm{HF} = \frac12\sum_{\mu\nu}D_{\mu\nu}\big(h_{\mu\nu} + F_{\mu\nu}\big) + V_\mathrm{NN}
$$

$V_\mathrm{NN}$ は核どうしの反発です。
基底関数を有限にしたので、問題は「行列を作り、一般化固有値問題を解く」ことに変わりました。
得られる $K$ 個の軌道のうち、エネルギーの低い $N/2$ 個を占有軌道、残りを仮想軌道（空軌道）と呼びます。

**ここまでで分かったこと**：軌道を基底関数で展開すると、Hartree–Fock 方程式は行列の方程式 $\mathbf F\mathbf C = \mathbf S\mathbf C\boldsymbol\varepsilon$ になります。$\mathbf F$ は密度行列 $\mathbf D$ を通して $\mathbf C$ に依存します。

## 5. SCF の手順：H₂O を自分で解く

第 3.3 節の繰り返しを、行列の言葉で書き直すと次のようになります。

1. 積分 $\mathbf S$、$\mathbf h$、$(\mu\nu|\lambda\sigma)$ を計算する。
2. 最初の密度行列 $\mathbf D$ を推測する。
3. $\mathbf D$ から Fock 行列 $\mathbf F$ を作り、エネルギーを計算する。
4. $\mathbf F\mathbf C = \mathbf S\mathbf C\boldsymbol\varepsilon$ を解き、エネルギーの低い $N/2$ 個の軌道から新しい $\mathbf D$ を作る。
5. エネルギー（と $\mathbf D$）が変わらなくなるまで 3〜4 を繰り返す。

この手順を、H₂O について numpy で書いてみます。
構造は O–H 距離 0.958 Å、H–O–H 角 104.5° とし、基底関数は cc-pVDZ [Dunning 1989] を使います。
積分の計算だけを PySCF [Sun 2020] に任せます。
最初の推測には、電子間の反発を無視した $\mathbf F = \mathbf h$ の解を使います。

```python
import numpy as np
import scipy.linalg
from pyscf import gto, scf

r, th = 0.958, np.deg2rad(104.5)   # O–H 距離（Å）と H–O–H 角
geom = [["O", (0, 0, 0)],
        ["H", (0,  r*np.sin(th/2), r*np.cos(th/2))],
        ["H", (0, -r*np.sin(th/2), r*np.cos(th/2))]]
mol = gto.M(atom=geom, basis="cc-pvdz", unit="Angstrom", verbose=0)

S = mol.intor("int1e_ovlp")                              # 重なり行列 S
Hcore = mol.intor("int1e_kin") + mol.intor("int1e_nuc")  # 1 電子の部分 h
eri = mol.intor("int2e")                                 # (μν|λσ)（化学者の表記）
nocc = mol.nelectron // 2
Enuc = mol.energy_nuc()
print("基底関数の数 K =", mol.nao, "  占有軌道の数 =", nocc)

def density(F):
    eps, C = scipy.linalg.eigh(F, S)        # FC = SCε を解く
    return 2 * C[:, :nocc] @ C[:, :nocc].T

D = density(Hcore)                          # 最初の推測：電子間反発を無視
E_old = 0.0
for it in range(1, 100):
    J = np.einsum("pqrs,rs->pq", eri, D)
    K = np.einsum("prqs,rs->pq", eri, D)
    F = Hcore + J - 0.5 * K                 # Fock 行列
    E = 0.5 * np.sum(D * (Hcore + F)) + Enuc
    if it <= 5 or it % 5 == 0 or abs(E - E_old) < 1e-10:
        print(f"{it:3d}  E = {E:.10f}  ΔE = {E-E_old:.1e}")
    if abs(E - E_old) < 1e-10:
        break
    E_old = E
    D = density(F)

mf = scf.RHF(mol).run()
print(f"自作の SCF  : {E:.8f} Eh")
print(f"PySCF の RHF: {mf.e_tot:.8f} Eh（{mf.cycles} 回で収束）")
```

出力は次のとおりです（途中の行は一部だけ表示しています）。

```text
基底関数の数 K = 24   占有軌道の数 = 5
  1  E = -68.8742008387  ΔE = -6.9e+01
  2  E = -69.9500047151  ΔE = -1.1e+00
  3  E = -73.3054234635  ΔE = -3.4e+00
  4  E = -73.4268430572  ΔE = -1.2e-01
  5  E = -74.6944568108  ΔE = -1.3e+00
 10  E = -76.0199831017  ΔE = -1.3e-02
 15  E = -76.0267330389  ΔE = -5.6e-05
 20  E = -76.0267606207  ΔE = -2.3e-07
 25  E = -76.0267607333  ΔE = -9.3e-10
 28  E = -76.0267607338  ΔE = -3.4e-11
自作の SCF  : -76.02676073 Eh
PySCF の RHF: -76.02676073 Eh（8 回で収束）
```

短いコードの繰り返しで、PySCF と同じエネルギー $-76.02676\,E_\mathrm h$ が得られました。
エネルギーは繰り返しのたびに下がり、28 回で $10^{-10}\,E_\mathrm h$ の範囲で変わらなくなりました。

PySCF は同じ答えに 8 回で到達しています。
これは、これまでの何回かの Fock 行列と誤差を組み合わせて次の Fock 行列を作る、DIIS という加速法を使っているためです [Pulay 1980; Pulay 1982]。
上の単純な繰り返しは、分子によっては振動して収束しないこともあります。
実際のプログラムは、よい初期推測と DIIS などの加速法を組み合わせています。

**ここまでで分かったこと**：SCF は「密度行列から Fock 行列を作り、対角化して新しい密度行列を作る」ことの繰り返しです。H₂O でも数十行のコードで解けます。

## 6. 軌道エネルギーと Koopmans の定理

### 6.1 H₂O の軌道エネルギー

前の節の計算で得た軌道エネルギーを表示します。

```python
for i, e in enumerate(mf.mo_energy[:7]):
    tag = "占有" if i < nocc else "空  "
    print(f"  {i+1}  {tag}  {e:10.5f}  {e*27.211386:9.3f}")
```

```text
  1  占有   -20.55057   -559.209
  2  占有    -1.33632    -36.363
  3  占有    -0.69887    -19.017
  4  占有    -0.56646    -15.414
  5  占有    -0.49310    -13.418
  6  空       0.18544      5.046
  7  空       0.25614      6.970
```

左の列がハートリー、右の列が eV です。
5 個の占有軌道は、対称性の記号でエネルギーの低い順に 1a₁、2a₁、1b₂、3a₁、1b₁ です。
いちばん低い 1a₁ は、ほぼ酸素の 1s 軌道（内殻）です。
最高被占軌道（HOMO）の 1b₁ は、分子面に垂直な酸素の 2p 軌道で、非共有電子対に当たります。

### 6.2 軌道エネルギーの意味

軌道エネルギーは、Fock 演算子の期待値です。

$$
\varepsilon_i = \langle\phi_i|\hat f|\phi_i\rangle = h_{ii} + \sum_{j}^{N/2}\big(2J_{ij} - K_{ij}\big)
$$

これは、軌道 $\phi_i$ の電子 1 個の運動エネルギーと核からの引力に、ほかのすべての電子との反発を足したものです。
第 2 節の $E_\mathrm{HF}$ と比べると、

$$
E_\mathrm{HF} - V_\mathrm{NN} = \sum_{i}^{N/2}\big(h_{ii} + \varepsilon_i\big) \neq 2\sum_{i}^{N/2}\varepsilon_i
$$

です。
軌道エネルギーを単純に足すと、電子どうしの反発を 2 回数えてしまうためです。
全エネルギーは、軌道エネルギーの和ではありません。

### 6.3 Koopmans の定理

それでは、軌道エネルギーそのものは何を表すのでしょうか。
答えを与えるのが、1934 年に Koopmans が示した Koopmans の定理です [Koopmans 1934]。教科書では Szabo–Ostlund の 3.3.1 節で詳しく扱われています [Szabo 1996]。
導出は付録 B に回しますが、結果は次のとおりです。

**Koopmans の定理**：$N$ 電子の Hartree–Fock の決定子から、占有軌道 $\phi_i$ の電子を 1 個取り除いた決定子を作ります。このとき、ほかの軌道は変えずにそのまま使います。そのエネルギーの差は

$$
E_{N-1} - E_N = -\varepsilon_i
$$

です。

つまり、占有軌道のエネルギーの符号を変えたものは、イオン化エネルギーの近似になります。
H₂O では、HOMO から $-\varepsilon_\mathrm{HOMO} = 13.42$ eV が得られます。
実験のイオン化エネルギーは、NIST Chemistry WebBook の評価値で 12.621 eV です（関連リンク）。
近似としては、それなりによい値です。

この近似には、2 つの誤差があります。

1. **軌道の緩和を無視している**：電子が 1 個抜けると、残りの電子の軌道は形を変えて、エネルギーを下げます。Koopmans の定理はこれを無視するので、イオン化エネルギーを大きめに見積もります。
2. **電子相関を無視している**：Hartree–Fock 法そのものの誤差です。ふつうは、中性分子の方が電子が多いぶん相関エネルギーが大きく、これを入れるとイオン化エネルギーは大きくなります。

2 つの誤差は逆向きなので、一部が打ち消し合います。
これを確かめるために、陽イオン H₂O⁺ を別に Hartree–Fock 法で解き、エネルギー差（ΔSCF 法）を求めます。
陽イオンは電子の数が奇数で閉殻ではないので、$\alpha$ と $\beta$ の電子に別々の空間軌道を許す非制限 Hartree–Fock（unrestricted Hartree–Fock、UHF）法 [Pople 1954] を使います。

```python
neu = gto.M(atom=geom, basis="cc-pvdz", verbose=0)                     # geom は前の節と同じ構造
cat = gto.M(atom=geom, basis="cc-pvdz", charge=1, spin=1, verbose=0)  # H2O+（2 重項）
e0 = scf.RHF(neu).run()
e1 = scf.UHF(cat).run()
print(f"Koopmans  : -ε(HOMO)        = {-e0.mo_energy[4]*27.211386:.2f} eV")
print(f"ΔSCF      : E(H2O+) - E(H2O) = {(e1.e_tot-e0.e_tot)*27.211386:.2f} eV")
```

```text
Koopmans  : -ε(HOMO)        = 13.42 eV
ΔSCF      : E(H2O+) - E(H2O) = 10.74 eV
```

陽イオンの軌道を緩和させると、イオン化エネルギーは 10.74 eV まで下がり、実験値の 12.62 eV を大きく下回りました。
緩和だけを入れると、相関の誤差だけが残るためです。
Koopmans の定理がそれなりによい値を与えたのは、2 つの誤差が打ち消し合ったためだと分かります。

!!! warning "仮想軌道のエネルギーには注意"
    Koopmans の定理を仮想軌道に使うと、電子親和力の近似が得られるように見えます。しかし、仮想軌道の場合は緩和と相関の誤差が同じ向きに働くので、誤差が大きくなります [Szabo 1996]。また、仮想軌道のエネルギーは基底関数によって大きく変わります。上の H₂O の LUMO のエネルギーが正なのは、そのまま「電子が付かない」ことを意味するわけではありません。

**ここまでで分かったこと**：占有軌道のエネルギーの符号を変えたものは、イオン化エネルギーの近似です（Koopmans の定理）。軌道の緩和と電子相関の誤差が打ち消し合うので、それなりによい値になります。

## まとめ

- 変分原理により、近似の波動関数のエネルギーは厳密なエネルギーより高くなります。Hartree–Fock 法は、決定子 1 個の中でエネルギーが最も低くなる軌道を選びます。
- Hartree–Fock の軌道は Fock 演算子 $\hat f = \hat h + \sum_j(2\hat J_j - \hat K_j)$ の固有関数で、各電子はほかの電子の平均の場の中を動きます。
- 軌道を基底関数で展開すると、Roothaan–Hall 方程式 $\mathbf F\mathbf C = \mathbf S\mathbf C\boldsymbol\varepsilon$ になり、密度行列を更新する SCF の繰り返しで解けます。
- 全エネルギーは軌道エネルギーの和ではありません。占有軌道のエネルギーの符号を変えたものが、イオン化エネルギーの近似になります（Koopmans の定理）。

## 付録 A　Hartree–Fock 方程式の導出

規格直交の条件 $\langle\phi_i|\phi_j\rangle = \delta_{ij}$ のもとで $E_\mathrm{HF}$ を最小にするため、ラグランジュの未定乗数 $\lambda_{ji}$ を使って

$$
\mathcal L = E_\mathrm{HF} - 2\sum_{ij}^{N/2}\lambda_{ji}\big(\langle\phi_i|\phi_j\rangle - \delta_{ij}\big)
$$

を考えます（係数 2 は後の式を見やすくするためです）。
軌道を $\phi_i \to \phi_i + \delta\phi_i$ と少し変えたときの $\mathcal L$ の変化を、$\delta\phi_i$ の 1 次まで求めます。

1 電子の項からは $2\sum_i\langle\delta\phi_i|\hat h|\phi_i\rangle$ と、その複素共役が出ます。
2 電子の項 $\sum_{ij}(2J_{ij} - K_{ij})$ では、$\delta\phi_i$ が $i$ の側にも $j$ の側にも現れます。
$J_{ij} = J_{ji}$、$K_{ij} = K_{ji}$ を使って 2 つをまとめると、$2\sum_{ij}\langle\delta\phi_i|(2\hat J_j - \hat K_j)|\phi_i\rangle$ と、その複素共役になります。
したがって、

$$
\delta\mathcal L = 2\sum_i\Big\langle\delta\phi_i\Big|\,\hat f\phi_i - \sum_j\lambda_{ji}\phi_j\Big\rangle + \text{(複素共役)}
$$

です。
任意の $\delta\phi_i$ で $\delta\mathcal L = 0$ となるためには、

$$
\hat f\,\phi_i = \sum_j^{N/2}\lambda_{ji}\,\phi_j
$$

が必要です。
$\boldsymbol\lambda$ はエルミート行列なので、ユニタリ変換 $\mathbf U$ で対角化できます。
占有軌道どうしをユニタリ変換で混ぜても、決定子は位相因子しか変わらず（第 2 回の性質 4）、$\hat f$ も変わりません。
そこで $\boldsymbol\lambda$ を対角にする軌道を選ぶと、本文の $\hat f\phi_i = \varepsilon_i\phi_i$ になります。
この軌道を正準軌道（canonical orbital）と呼びます。

## 付録 B　Koopmans の定理の導出

スピン軌道 $\chi_a$ で書くと、$N$ 電子の決定子のエネルギーは

$$
E_N = \sum_a^{N} h_{aa} + \frac12\sum_{a}^{N}\sum_{b}^{N}\big([aa|bb] - [ab|ba]\big)
$$

です。$[ab|cd]$ はスピン軌道についての化学者の表記の積分で、$a = b$ の項は 0 になります。
ここから、軌道を変えずに $\chi_c$ を取り除いた決定子のエネルギー $E_{N-1}$ を引くと、$\chi_c$ を含む項だけが残ります。

$$
E_N - E_{N-1} = h_{cc} + \sum_{b}^{N}\big([cc|bb] - [cb|bc]\big)
$$

右辺は、スピン軌道 $\chi_c$ の軌道エネルギー $\varepsilon_c = \langle\chi_c|\hat f|\chi_c\rangle$ そのものです。
したがって $E_{N-1} - E_N = -\varepsilon_c$ となります。
閉殻の場合、スピンを積分すると本文の $\varepsilon_i = h_{ii} + \sum_j(2J_{ij} - K_{ij})$ に一致します。

## 参考文献

- A. Szabo, N. S. Ostlund, "Modern Quantum Chemistry" (Dover, 1996). ISBN: 9780486691862
- D. R. Hartree, "The Wave Mechanics of an Atom with a Non-Coulomb Central Field. Part I. Theory and Methods", Math. Proc. Camb. Philos. Soc. 24, 89 (1928). DOI: [10.1017/S0305004100011919](https://doi.org/10.1017/S0305004100011919)
- V. Fock, "Näherungsmethode zur Lösung des quantenmechanischen Mehrkörperproblems", Z. Phys. 61, 126 (1930). DOI: [10.1007/BF01340294](https://doi.org/10.1007/BF01340294)
- J. C. Slater, "Note on Hartree's Method", Phys. Rev. 35, 210 (1930). DOI: [10.1103/PhysRev.35.210.2](https://doi.org/10.1103/PhysRev.35.210.2)
- C. C. J. Roothaan, "New Developments in Molecular Orbital Theory", Rev. Mod. Phys. 23, 69 (1951). DOI: [10.1103/RevModPhys.23.69](https://doi.org/10.1103/RevModPhys.23.69)
- G. G. Hall, "The molecular orbital theory of chemical valency VIII. A method of calculating ionization potentials", Proc. R. Soc. Lond. A 205, 541 (1951). DOI: [10.1098/rspa.1951.0048](https://doi.org/10.1098/rspa.1951.0048)
- J. A. Pople, R. K. Nesbet, "Self-Consistent Orbitals for Radicals", J. Chem. Phys. 22, 571 (1954). DOI: [10.1063/1.1740120](https://doi.org/10.1063/1.1740120)
- T. Koopmans, "Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms", Physica 1, 104 (1934). DOI: [10.1016/S0031-8914(34)90011-2](https://doi.org/10.1016/S0031-8914(34)90011-2)
- P. Pulay, "Convergence acceleration of iterative sequences. The case of SCF iteration", Chem. Phys. Lett. 73, 393 (1980). DOI: [10.1016/0009-2614(80)80396-4](https://doi.org/10.1016/0009-2614(80)80396-4)
- P. Pulay, "Improved SCF convergence acceleration", J. Comput. Chem. 3, 556 (1982). DOI: [10.1002/jcc.540030413](https://doi.org/10.1002/jcc.540030413)
- T. H. Dunning Jr., "Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through neon and hydrogen", J. Chem. Phys. 90, 1007 (1989). DOI: [10.1063/1.456153](https://doi.org/10.1063/1.456153)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- NIST Chemistry WebBook, Water（イオン化エネルギーの評価値）：<https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=20>
