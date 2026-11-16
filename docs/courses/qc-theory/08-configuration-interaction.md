# 第 8 回　配置間相互作用（CI）

!!! abstract "この回の道筋"
    HF の決定子に、電子を励起させた決定子（励起配置）を足し合わせ、その係数を変分法で決める方法が配置間相互作用（CI）です。この回では、励起配置の記号を決め、Brillouin の定理から「2 電子励起が最も効く」ことを導きます。そのうえで CIS、CISD、FCI を順に見て、CISD の大きな弱点であるサイズ無矛盾性の欠如と、その補正（Davidson 補正）を計算で確かめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 励起配置と CI の波動関数 |
    | 2 | Brillouin の定理と、相関エネルギーを決める配置 |
    | 3 | CIS、CISD、FCI |
    | 4 | 計算例：H₂O の CISD と FCI |
    | 5 | サイズ無矛盾性の欠如 |
    | 6 | Davidson 補正 |

    **前に読んでおく回**：[量子化学の理論 第 7 回](07-electron-correlation.md)、[第 3 回](03-hartree-fock.md)

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書きます。

## 1. 励起配置と CI の波動関数

HF 法で得たスピン軌道を、電子が入っているもの（占有軌道）と、空いているもの（仮想軌道）に分けます。記号は次のように決めます。

- $i, j, k, \dots$：占有スピン軌道
- $a, b, c, \dots$：仮想スピン軌道
- $\Phi_0$：HF の決定子

$\Phi_0$ から、占有軌道 $i$ の電子を仮想軌道 $a$ に移した決定子を $\Phi_i^a$ と書き、**1 電子励起配置**と呼びます。同じように、$i, j$ の 2 個の電子を $a, b$ に移した決定子 $\Phi_{ij}^{ab}$ を **2 電子励起配置**と呼びます。3 電子、4 電子の励起配置も同じように作れます。

これらの決定子は、どれも同じ HF 軌道からできていて、互いに直交します。しかも、0 電子から $N$ 電子までの励起配置をすべて集めると、その基底関数で作れる $N$ 電子の波動関数を、どれでも表せます。そこで、厳密な波動関数をこれらの和で表します。

$$
\Psi = c_0 \Phi_0 + \sum_{i,a} c_i^a \Phi_i^a + \sum_{i<j,\,a<b} c_{ij}^{ab} \Phi_{ij}^{ab} + \sum_{i<j<k,\,a<b<c} c_{ijk}^{abc} \Phi_{ijk}^{abc} + \cdots
$$

係数 $c$ は、エネルギーの期待値が最小になるように決めます（変分法）。これは、決定子を基底としてハミルトニアンの行列 $\mathbf H$（要素 $H_{IJ} = \langle \Phi_I | \hat H | \Phi_J \rangle$）をつくり、固有値問題

$$
\mathbf H \mathbf c = E \mathbf c
$$

を解くことと同じです。最も小さい固有値が基底状態、その次からが励起状態のエネルギーです。この方法を**配置間相互作用**（configuration interaction, CI）と呼びます。

すべての励起配置を入れた CI を**完全 CI**（full CI, FCI）と呼びます。FCI は、その基底関数の中で Schrödinger 方程式を厳密に解いたことになります。ただし後で見るように、決定子の数が爆発的に増えるので、FCI は小さな系でしか解けません。そこで、励起の段数で打ち切ります。

**ここまでで分かったこと**：CI は、HF の決定子と励起配置の線形結合で波動関数を表し、係数をハミルトニアンの行列の対角化で決める方法です。

## 2. Brillouin の定理と、相関エネルギーを決める配置

どの励起配置が大事かを考えるために、行列要素 $\langle \Phi_0 | \hat H | \Phi_I \rangle$ を調べます。必要になるのは、決定子どうしのハミルトニアンの行列要素を与える **Slater–Condon 則**です。その結果だけを書きます（導出は[第 4 回](04-second-quantization.md)や Szabo–Ostlund の 2 章にあります）。

- 二つの決定子が **3 個以上**のスピン軌道で違うと、行列要素は 0 です。ハミルトニアンが 2 電子演算子までしか含まないからです。
- 1 個だけ違う場合：$\langle \Phi_0 | \hat H | \Phi_i^a \rangle = \langle i | \hat h | a \rangle + \sum_j \langle ij \| aj \rangle$
- 2 個違う場合：$\langle \Phi_0 | \hat H | \Phi_{ij}^{ab} \rangle = \langle ij \| ab \rangle$

ここで $\hat h$ は 1 電子ハミルトニアン（運動エネルギーと核の引力）です。2 電子積分は**物理学者の表記**で書き、反対称化したものを $\langle pq \| rs \rangle = \langle pq | rs \rangle - \langle pq | sr \rangle$ とします。

$$
\langle pq | rs \rangle = \int\!\!\int \phi_p^*(\mathbf x_1)\, \phi_q^*(\mathbf x_2)\, \frac{1}{r_{12}}\, \phi_r(\mathbf x_1)\, \phi_s(\mathbf x_2)\, d\mathbf x_1 d\mathbf x_2
$$

1 個違う場合の式の右辺は、Fock 演算子 $\hat f$ の行列要素 $f_{ia} = \langle i | \hat f | a \rangle$ そのものです。HF 方程式が解けていると、占有軌道と仮想軌道の間の Fock 行列の要素は 0 になります。したがって

$$
\langle \Phi_0 | \hat H | \Phi_i^a \rangle = f_{ia} = 0
$$

です。これを **Brillouin の定理**と呼びます。HF 法の解が、軌道を少し回転させてもエネルギーが変わらない点（停留点）であることを、言い換えたものです。

PySCF で確かめます。HF の解の Fock 行列を分子軌道の基底に直し、占有と仮想の間の要素の最大値を見ます。

```python
import numpy as np
from pyscf import gto, scf

mol = gto.M(atom="O 0 0 0; H 0 0.759 0.588; H 0 -0.759 0.588", basis="6-31g", verbose=0)
mf = scf.RHF(mol).run()
C = mf.mo_coeff
F = C.T @ mf.get_fock() @ C          # MO 基底での Fock 行列
nocc = mol.nelectron // 2
print("max |F_ia| =", f"{np.abs(F[:nocc, nocc:]).max():.1e}")
```

```
max |F_ia| = 1.9e-07
```

0 でないのは、SCF の収束の判定が有限だからです。収束を厳しくすれば、もっと小さくなります。

この定理から、相関エネルギーについて大事なことが言えます。$\Psi$ を厳密な固有関数とし、$\langle \Phi_0 | \Psi \rangle = 1$ となるように規格化します（中間規格化、$c_0 = 1$）。$\hat H \Psi = E \Psi$ の左から $\Phi_0$ をかけると、$E = \langle \Phi_0 | \hat H | \Psi \rangle$ です。右辺に $\Psi$ の展開を入れると、1 電子励起の項は Brillouin の定理で消えます。3 電子以上の励起の項は Slater–Condon 則で消えます。残るのは次の式です。

$$
E = E_\mathrm{HF} + \sum_{i<j,\,a<b} c_{ij}^{ab} \langle ij \| ab \rangle
$$

つまり、**相関エネルギーは 2 電子励起配置の係数だけで決まります**。もちろん、その係数 $c_{ij}^{ab}$ の値は、ほかの励起配置との混ざり合いを通して、1 電子・3 電子・4 電子励起の影響も受けます。

**ここまでで分かったこと**：Brillouin の定理により、HF の決定子は 1 電子励起配置と直接には混ざりません。相関エネルギーは 2 電子励起配置の係数で表されるので、2 電子励起が最も重要です。

## 3. CIS、CISD、FCI

### CIS

$\Phi_0$ と 1 電子励起配置だけを使う CI を **CIS**（CI singles）と呼びます。Brillouin の定理により、$\Phi_0$ と 1 電子励起配置の間の行列要素は 0 です。そのため、CIS の基底状態のエネルギーは HF のエネルギーと同じで、相関は入りません。

一方、1 電子励起配置どうしの行列を対角化すると、励起状態のエネルギーと波動関数が得られます。CIS は、励起状態を調べる最も簡単な方法として使われてきました [Foresman 1992]。

### CISD

基底状態の相関を入れる最も簡単な方法は、2 電子励起配置を加えることです。$\Phi_0$、1 電子励起、2 電子励起の配置を使う CI を **CISD**（CI singles and doubles）と呼びます。1 電子励起配置は、$\Phi_0$ とは直接混ざりませんが、2 電子励起配置と混ざるので、間接的にエネルギーに効きます。

CISD の計算量を見積もります。占有軌道の数を $o$、仮想軌道の数を $v$ とすると、2 電子励起配置の数はおよそ $o^2 v^2$ です。1 回の反復で最も重い計算は、$o^2 v^4$ に比例します。$o$ と $v$ はどちらも系の大きさ $N$ に比例するので、CISD の計算量は $N^6$ に比例します。

### FCI

FCI の決定子の数は、空間軌道の数を $n$、α 電子と β 電子の数を $N_\alpha$、$N_\beta$ とすると

$$
\binom{n}{N_\alpha} \binom{n}{N_\beta}
$$

です。たとえば、次の節で使う H₂O の 6-31G 基底では $n = 13$、$N_\alpha = N_\beta = 5$ なので、$\binom{13}{5}^2 = 1287^2 \approx 1.7 \times 10^6$ 個です。基底を cc-pVDZ（$n = 24$）にするだけで、$\binom{24}{5}^2 = 42504^2 \approx 1.8 \times 10^9$ 個に増えます。軌道と電子を増やすと、決定子の数は指数関数的に増えます。

**ここまでで分かったこと**：CIS は基底状態の相関を入れません。CISD は 2 電子励起までを入れて相関の大部分を取り込み、計算量は $N^6$ です。FCI は厳密ですが、決定子の数が指数関数的に増えます。

## 4. 計算例：H₂O の CISD と FCI

H₂O を 6-31G 基底で計算し、CISD と FCI の相関エネルギーを比べます。平衡に近い構造（O–H 0.96 Å、H–O–H 104.5°）から、二つの O–H 結合を同時に 1.5 倍、2 倍に伸ばした構造も計算します。最後の二つの列（CISD+Q と $c_0^2$）は 6 節で使います。

```python
import numpy as np
from pyscf import gto, scf, ci, fci

def h2o(scale):
    """平衡構造（O–H 0.96 Å、HOH 104.5°）の O–H 結合を scale 倍に伸ばす。"""
    r, th = 0.96 * scale, np.deg2rad(104.5) / 2
    return f"O 0 0 0; H 0 {r*np.sin(th):.6f} {r*np.cos(th):.6f}; H 0 {-r*np.sin(th):.6f} {r*np.cos(th):.6f}"

print("scale   E_c(CISD)  E_c(CISD+Q)  E_c(FCI)    c0^2")
for scale in [1.0, 1.5, 2.0]:
    mol = gto.M(atom=h2o(scale), basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    myci = ci.CISD(mf).run()
    c0 = myci.ci[0]                               # HF 決定子の係数（規格化済み）
    e_q = (1 - c0**2) * myci.e_corr               # Davidson 補正
    e_fci = fci.FCI(mf).kernel()[0]
    print(f"{scale:4.1f} {myci.e_corr:11.5f} {myci.e_corr + e_q:11.5f}"
          f" {e_fci - mf.e_tot:11.5f} {c0**2:7.3f}")
```

```
scale   E_c(CISD)  E_c(CISD+Q)  E_c(FCI)    c0^2
 1.0    -0.13030    -0.13545    -0.13712   0.960
 1.5    -0.17712    -0.19474    -0.19702   0.900
 2.0    -0.23926    -0.29421    -0.29305   0.770
```

この計算は、FCI の部分に 2 分ほどかかります（CISD は数秒です）。

CISD が取り込んだ相関エネルギーの割合は、平衡構造で $0.13030 / 0.13712 = 95\%$ です。結合を 1.5 倍に伸ばすと 90%、2 倍では 82% に下がります。結合を伸ばすと、第 7 回で見た静的相関が強くなります。HF の決定子の重み $c_0^2$ も 0.96 から 0.77 に下がり、2 電子励起までで打ち切る近似が悪くなります。

**ここまでで分かったこと**：平衡構造の H₂O では、CISD は相関エネルギーの 95% を取り込みます。結合を伸ばすと、その割合は下がります。

## 5. サイズ無矛盾性の欠如

CISD には、もっと根本的な弱点があります。それを見るために、十分に離れていて相互作用しない二つの分子 A と B を考えます。正しいエネルギーは、それぞれのエネルギーの和 $E(\mathrm A) + E(\mathrm B)$ のはずです。手法がこの関係を満たすことを、**サイズ無矛盾性**（size consistency）と呼びます。

CISD はこれを満たしません。理由は次のとおりです。A と B が離れていると、厳密な波動関数は、それぞれの波動関数の積 $\Psi_\mathrm A \Psi_\mathrm B$ です。$\Psi_\mathrm A$ と $\Psi_\mathrm B$ はどちらも 2 電子励起配置を含みます。その積には、「A で 2 電子励起、同時に B でも 2 電子励起」という配置が現れます。これは全体から見ると 4 電子励起配置です。CISD は 4 電子励起配置を含まないので、A+B 全体を CISD で計算すると、この寄与を取りこぼします。

実際に計算します。平衡距離の H₂ を、互いに 50 Å 離して $n$ 個並べ、全体を CISD で計算します。電子が 2 個の H₂ 1 個では、CISD は FCI と同じです。したがって、正しい答えは「1 分子あたりの相関エネルギーが $n$ によらず一定」です。

```python
from pyscf import gto, scf, ci

print("n   E_c(CISD)/n   E_c(CISD+Q)/n")
for n in [1, 2, 4, 8]:
    # H2（R = 0.74 Å）を n 個、互いに 50 Å 離して並べる
    atoms = "; ".join(f"H {50*k} 0 0; H {50*k} 0 0.74" for k in range(n))
    mol = gto.M(atom=atoms, basis="cc-pvdz", verbose=0)
    mf = scf.RHF(mol).run()
    myci = ci.CISD(mf).run()
    e_q = (1 - myci.ci[0]**2) * myci.e_corr
    print(f"{n}  {myci.e_corr / n:12.5f}  {(myci.e_corr + e_q) / n:12.5f}")
```

```
n   E_c(CISD)/n   E_c(CISD+Q)/n
1      -0.03467      -0.03526
2      -0.03411      -0.03520
4      -0.03309      -0.03500
8      -0.03138      -0.03445
```

2 列目の CISD の値を見ます。$n = 1$ の $-0.03467\ E_\mathrm{h}$ は正しい値（第 7 回の FCI の値と同じ）です。ところが $n$ が増えると、1 分子あたりの相関エネルギーの絶対値は小さくなり、$n = 8$ では 1 割ほど減っています。分子の数を増やすほど、取りこぼす割合が増えます。最小基底の H₂ を $n$ 個並べたモデルでは、$n$ が大きいとき、CISD の全相関エネルギーが $n$ ではなく $\sqrt n$ に比例します。導出は付録 A に回します（Szabo–Ostlund の 4 章にも同じ例があります）。

この欠点は、大きな分子ほど深刻です。たとえば反応エネルギーを計算するとき、大きさの違う分子どうしを比べると、誤差が打ち消し合いません。このため現在では、CISD を単独で使うことはあまりなく、サイズ無矛盾性を満たす結合クラスター法（第 10 回）が使われています。

**ここまでで分かったこと**：CISD は、離れた分子の同時励起（全体としての 4 電子励起）を含まないため、サイズ無矛盾性を満たしません。系が大きくなるほど、取り込める相関の割合が減ります。

## 6. Davidson 補正

CISD で取りこぼした 4 電子励起などの寄与を、簡単に見積もる方法があります。CISD の波動関数を規格化したときの HF の決定子の係数を $c_0$、CISD の相関エネルギーを $\Delta E_\mathrm{CISD}$ として、

$$
\Delta E_\mathrm{Q} = (1 - c_0^2)\, \Delta E_\mathrm{CISD}
$$

を足します。これを **Davidson 補正**と呼び、補正したものを CISD+Q と書きます [Langhoff 1974]。Q は 4 電子励起（quadruples）の頭文字です。HF の決定子の重みが小さいほど、つまり励起配置の重みが大きいほど、補正が大きくなります。

4 節と 5 節の表の CISD+Q の列を見ます。

- H₂O（4 節）：平衡構造では、FCI との差が $0.0068$ から $0.0017\ E_\mathrm{h}$ に縮みます。1.5 倍の構造でも同じように改善します。2 倍の構造では、CISD+Q（$-0.29421$）が FCI（$-0.29305$）を下回り、補正が大きすぎます。
- H₂ を $n$ 個並べた系（5 節）：1 分子あたりの値の $n$ による変化が、CISD より小さくなります。ただし $n = 8$ でもまだ一定ではありません。また $n = 1$ では、CISD がすでに厳密なのに補正が加わり、正しい値より低くなります。

Davidson 補正は、こうした簡単な見積もりの式です。変分法ではないので、エネルギーが FCI より低くなることもあります。サイズ無矛盾性も、厳密には回復しません。CISD を使うときの手軽な改善策ですが、限界があることを知って使ってください。

**ここまでで分かったこと**：Davidson 補正は、CISD のエネルギーを手軽に改善し、サイズ無矛盾性の欠如を小さくします。ただし近似的な補正で、過大になることもあります。

## まとめ

- CI は、HF の決定子と励起配置の線形結合で波動関数を表し、係数を対角化で決める変分法です。
- Brillouin の定理により、HF の決定子は 1 電子励起配置と直接には混ざらず、相関エネルギーは 2 電子励起の係数で決まります。
- CIS は基底状態の相関を入れず、CISD は相関の大部分を入れます（計算量 $N^6$）。FCI は厳密ですが、決定子の数が指数関数的に増えます。
- CISD はサイズ無矛盾性を満たさず、系が大きいほど相関を取りこぼします。
- Davidson 補正（CISD+Q）はこれを手軽に改善しますが、近似的で、過大になることもあります。

## 付録 A　最小基底の H₂ を $n$ 個並べた系の CISD

各 H₂ に最小基底を使うと、軌道は $\sigma_g$ と $\sigma_u$ の 2 個だけです。H₂ 1 個の FCI の波動関数は $\sigma_g^2$（HF）と $\sigma_u^2$（2 電子励起）の組み合わせになります（1 電子励起は対称性で混ざりません）。

$n$ 個の H₂ が十分に離れているとき、CISD の空間は次の $n + 1$ 個の決定子でできます。

- $\Phi_0$：すべての H₂ が $\sigma_g^2$ にある HF の決定子
- $\Phi_k$（$k = 1, \dots, n$）：$k$ 番目の H₂ だけが $\sigma_u^2$ に励起した決定子

全体の 1 電子励起は、1 個の H₂ の中の励起（対称性が合わない）か、離れた分子の間の励起（軌道の重なりがない）です。どちらも、ここに挙げた決定子との行列要素が 0 なので、混ざりません。ハミルトニアンの行列要素は、$E_\mathrm{HF}$ を基準にして次のとおりです。

- $\langle \Phi_0 | \hat H - E_\mathrm{HF} | \Phi_k \rangle = K$（H₂ 1 個の中の 2 電子積分で、$k$ によりません）
- $\langle \Phi_k | \hat H - E_\mathrm{HF} | \Phi_k \rangle = 2\Delta$（1 個の H₂ を励起するのに要するエネルギー、$\Delta > 0$）
- $\langle \Phi_k | \hat H | \Phi_l \rangle = 0$（$k \ne l$。二つの決定子は 4 個のスピン軌道で違うので、Slater–Condon 則により 0）

基底状態ではすべての $\Phi_k$ の係数が等しいので、$\Psi = \Phi_0 + c \sum_k \Phi_k$ とおけます。相関エネルギーを $E_\mathrm c$ とすると、固有値問題は

$$
E_\mathrm c = n K c, \qquad K + 2\Delta\, c = E_\mathrm c\, c
$$

です。第 2 式から $c = K / (E_\mathrm c - 2\Delta)$ を第 1 式に入れると、$E_\mathrm c^2 - 2\Delta E_\mathrm c - n K^2 = 0$ となります。負の解をとって

$$
E_\mathrm c = \Delta - \sqrt{\Delta^2 + n K^2}
$$

を得ます。$n = 1$ ではこれが H₂ 1 個の厳密な相関エネルギー $\varepsilon$ です。正しい答えは $n\varepsilon$ ですが、$n$ が大きいと $E_\mathrm c \approx -\sqrt n\, |K|$ となり、$\sqrt n$ にしか比例しません。

## 参考文献

- A. Szabo, N. S. Ostlund, "Modern Quantum Chemistry" (Dover, 1996). ISBN: 9780486691862
- T. Helgaker, P. Jørgensen, J. Olsen, "Molecular Electronic-Structure Theory" (Wiley, 2000). ISBN: 9780471967552
- S. R. Langhoff, E. R. Davidson, "Configuration interaction calculations on the nitrogen molecule", Int. J. Quantum Chem. 8, 61 (1974). DOI: [10.1002/qua.560080106](https://doi.org/10.1002/qua.560080106)
- J. B. Foresman, M. Head-Gordon, J. A. Pople, M. J. Frisch, "Toward a systematic molecular orbital theory for excited states", J. Phys. Chem. 96, 135 (1992). DOI: [10.1021/j100180a030](https://doi.org/10.1021/j100180a030)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
