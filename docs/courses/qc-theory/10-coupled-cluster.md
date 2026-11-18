# 第 10 回　結合クラスター法：CCSD(T)

!!! abstract "この回の道筋"
    結合クラスター法（CC 法）は、励起演算子の指数関数で波動関数を表す方法です。指数関数の形のおかげで、CISD の弱点だったサイズ無矛盾性が自然に満たされます。この回では、CC 法の考え方から CCSD と CCSD(T) までを順に見て、CCSD(T) が「ゴールドスタンダード」と呼ばれる理由と、それが破綻する場面を計算で確かめます。最後に、破綻を事前に見分ける T1・D1 診断を紹介します。

    | 節 | 内容 |
    |---|---|
    | 1 | 指数型の波動演算子 |
    | 2 | サイズ無矛盾性 |
    | 3 | CC の方程式と CCSD |
    | 4 | (T)：3 電子励起の摂動補正 |
    | 5 | 計算コスト |
    | 6 | 「ゴールドスタンダード」と、その限界 |
    | 7 | T1・D1 診断 |

    **前に読んでおく回**：[量子化学の理論 第 8 回](08-configuration-interaction.md)、[第 9 回](09-mp2.md)、[第 4 回](04-second-quantization.md)

!!! tip "物理的なイメージ"
    結合クラスター法が求めているのは、電子対の避け合い（2 電子励起）を、分子のあちこちで同時に、独立に起こせる形で表した波動関数です。励起演算子の指数関数 $e^{\hat T}$ を展開すると、2 電子励起を 2 か所で同時に起こした 4 電子励起なども、積の形で自動的に入ります。離れた 2 つの分子では、それぞれの中の励起の積がそのまま現れるので、エネルギーが足し算になります（サイズ無矛盾）。(T) は、3 個の電子が同時に避け合う効果を、摂動で安く見積もった補正です。

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書き、1 m$E_\mathrm{h}$ = 0.001 $E_\mathrm{h}$ ≈ 0.63 kcal/mol です。記号は第 8 回と同じで、$i, j, k$ は占有スピン軌道、$a, b, c$ は仮想スピン軌道、$\Phi_0$ は HF の決定子です。2 電子積分は、物理学者の表記で反対称化した $\langle pq \| rs \rangle = \langle pq | rs \rangle - \langle pq | sr \rangle$ を使います。

## 1. 指数型の波動演算子

CC 法では、励起配置を**励起演算子**でつくります。第二量子化（[第 4 回](04-second-quantization.md)）の生成演算子 $\hat a_a^\dagger$ と消滅演算子 $\hat a_i$ を使うと、1 電子励起と 2 電子励起の演算子は次のように書けます。

$$
\hat T_1 = \sum_{ia} t_i^a\, \hat a_a^\dagger \hat a_i, \qquad
\hat T_2 = \frac14 \sum_{ijab} t_{ij}^{ab}\, \hat a_a^\dagger \hat a_b^\dagger \hat a_j \hat a_i
$$

$\hat T_1 \Phi_0$ は 1 電子励起配置の和、$\hat T_2 \Phi_0$ は 2 電子励起配置の和です。係数 $t$ を**振幅**（amplitude）と呼びます。3 電子励起の $\hat T_3$ なども同じように定義します。

CI 法では、波動関数を $(1 + \hat T_1 + \hat T_2 + \cdots)\Phi_0$ の形の**線形結合**で表しました。CC 法では、これを**指数関数**に置き換えます [Čížek 1966]。

$$
\Psi_\mathrm{CC} = e^{\hat T} \Phi_0, \qquad \hat T = \hat T_1 + \hat T_2 + \cdots
$$

$e^{\hat T}$ を**波動演算子**と呼びます。CC 法の全体像は、総説 [Bartlett 2007] にまとまっています。$\hat T$ のすべての励起を入れれば、CC 法は FCI と同じ答えになります。違いが出るのは、$\hat T$ を途中で打ち切ったときです。

たとえば $\hat T = \hat T_2$ で打ち切ったとします。指数関数を展開すると

$$
e^{\hat T_2} \Phi_0 = \Big(1 + \hat T_2 + \tfrac12 \hat T_2^2 + \tfrac16 \hat T_2^3 + \cdots\Big) \Phi_0
$$

となります。$\tfrac12 \hat T_2^2 \Phi_0$ は 4 電子励起配置です。ただし、その係数は独立な変数ではなく、2 電子励起の振幅の積 $t_{ij}^{ab} t_{kl}^{cd}$ で決まります。つまり CC 法は、「2 組の電子対が独立に、同時に励起する」寄与を、新しい変数を増やさずに取り込みます。このように振幅の積で表される高次の励起を、**非連結**（disconnected）な励起と呼びます。

**ここまでで分かったこと**：CC 法は、励起演算子 $\hat T$ の指数関数で波動関数を表します。$\hat T$ を打ち切っても、振幅の積の形で高次の励起が入ります。

## 2. サイズ無矛盾性

第 8 回で見たように、CISD は、離れた二つの分子 A と B の「同時の 2 電子励起」（全体で 4 電子励起）を含まないので、サイズ無矛盾ではありませんでした。CC 法ではこの問題が起きません。

A と B が十分に離れているとき、$\hat T$ は A の中の励起 $\hat T_\mathrm A$ と B の中の励起 $\hat T_\mathrm B$ の和になります。二つは別々の軌道に作用するので交換し、指数関数は積に分かれます。

$$
e^{\hat T_\mathrm A + \hat T_\mathrm B} \Phi_0 = e^{\hat T_\mathrm A} e^{\hat T_\mathrm B} \Phi_0
$$

これは、A と B の CC の波動関数の積です。$e^{\hat T_\mathrm A}$ と $e^{\hat T_\mathrm B}$ の展開の積に、A と B の同時の励起が自動的に入っています。3 節で述べるエネルギーの式も A と B の和に分かれるので、$E(\mathrm{A+B}) = E(\mathrm A) + E(\mathrm B)$ が成り立ちます。

第 8 回と同じ計算（H₂ を 50 Å 離して $n$ 個並べる）を、CCSD（3 節）で行います。

```python
from pyscf import gto, scf, cc

print("n   E_c(CCSD)/n")
for n in [1, 2, 4, 8]:
    atoms = "; ".join(f"H {50*k} 0 0; H {50*k} 0 0.74" for k in range(n))
    mol = gto.M(atom=atoms, basis="cc-pvdz", verbose=0)
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    print(f"{n}  {mycc.e_corr / n:12.5f}")
```

```
n   E_c(CCSD)/n
1      -0.03467
2      -0.03467
4      -0.03467
8      -0.03467
```

CISD では、1 分子あたりの相関エネルギーが $n = 8$ で $-0.03138$ まで減りました。CCSD では、$n$ によらず正しい値 $-0.03467$（H₂ 1 個の FCI の値）のままです。

**ここまでで分かったこと**：指数関数は和を積に変えるので、CC 法は打ち切っても（CCSD でも）サイズ無矛盾です。

## 3. CC の方程式と CCSD

### 振幅の決め方

CI 法では、変分法で係数を決めました。CC 法で同じことをすると、エネルギーの式が無限に続く項を含み、実用になりません。そこで CC 法では、**射影**で振幅を決めます。

Schrödinger 方程式 $\hat H e^{\hat T} \Phi_0 = E e^{\hat T} \Phi_0$ の両辺に、左から $e^{-\hat T}$ をかけます。

$$
\bar H \Phi_0 = E \Phi_0, \qquad \bar H = e^{-\hat T} \hat H e^{\hat T}
$$

$\bar H$ を**相似変換したハミルトニアン**と呼びます。この式に、左から $\Phi_0$、1 電子励起配置、2 電子励起配置をかけると、次の式が得られます。

$$
E = \langle \Phi_0 | \bar H | \Phi_0 \rangle, \qquad
0 = \langle \Phi_i^a | \bar H | \Phi_0 \rangle, \qquad
0 = \langle \Phi_{ij}^{ab} | \bar H | \Phi_0 \rangle
$$

後の二つが、振幅 $t_i^a$ と $t_{ij}^{ab}$ を決める連立方程式（**振幅方程式**）です。振幅について非線形なので、反復して解きます。

$\bar H$ は交換子の級数（Baker–Campbell–Hausdorff 展開）

$$
\bar H = \hat H + [\hat H, \hat T] + \frac{1}{2!}[[\hat H, \hat T], \hat T] + \frac{1}{3!}[[[\hat H, \hat T], \hat T], \hat T] + \frac{1}{4!}[[[[\hat H, \hat T], \hat T], \hat T], \hat T]
$$

で書けます。この級数は、4 重の交換子で**ちょうど終わります**。$\hat H$ が 2 電子演算子までしか含まないためです（理由は[第 4 回](04-second-quantization.md)の第二量子化の規則から出ます。詳しくは [Crawford 2000] を参照してください）。このため、振幅方程式は有限の項で書けます。

エネルギーの式は特に簡単です。HF の軌道を使うと（$f_{ia} = 0$）、

$$
E_\mathrm{CC} = E_\mathrm{HF} + \sum_{i<j} \sum_{a<b} \langle ij \| ab \rangle \big( t_{ij}^{ab} + t_i^a t_j^b - t_i^b t_j^a \big)
$$

となります。第 8 回で、厳密な相関エネルギーは 2 電子励起の係数で決まることを見ました。CC 法では、その係数が $t_{ij}^{ab} + t_i^a t_j^b - t_i^b t_j^a$ に置き換わっています。

CC 法のエネルギーは、変分法で求めたものではありません。そのため、CI 法と違って、FCI より低いエネルギーが出ることもあります（6 節で実例を見ます）。

### CCSD

$\hat T = \hat T_1 + \hat T_2$ で打ち切った CC 法を **CCSD**（CC singles and doubles）と呼びます [Purvis 1982]。CCSD は、変数の数は CISD と同じですが、指数関数の展開を通して $\hat T_2^2$（非連結な 4 電子励起）や $\hat T_1 \hat T_2$（非連結な 3 電子励起）などを含みます。

**ここまでで分かったこと**：CC 法の振幅は、相似変換したハミルトニアンへの射影で決めます。方程式は有限の項で書けますが、エネルギーは変分的ではありません。CCSD は $\hat T_1$ と $\hat T_2$ までを使う CC 法です。

## 4. (T)：3 電子励起の摂動補正

CCSD が含まない最も大きな寄与は、**連結した** 3 電子励起、つまり $\hat T_3$ です。$\hat T_3$ まで含めた CCSDT は正確ですが、計算量が $N^8$ に比例し、小さな分子にしか使えません（5 節）。

そこで、$\hat T_3$ の寄与を摂動論で見積もり、CCSD のエネルギーに足します。これが **CCSD(T)** です [Raghavachari 1989]。考え方は次のとおりです。

1. まず CCSD を解き、$t_i^a$ と $t_{ij}^{ab}$ を得ます。
2. 3 電子励起の振幅を、MP 摂動論（第 9 回）の形の式で、CCSD の振幅から**一度だけ**見積もります。分母は第 9 回と同じく軌道エネルギーの差 $\varepsilon_i + \varepsilon_j + \varepsilon_k - \varepsilon_a - \varepsilon_b - \varepsilon_c$ です。
3. その振幅から、エネルギーの補正を計算します。補正は、摂動論の 4 次で現れる 3 電子励起の項と、5 次で現れる 1 電子励起と 3 電子励起の結合の項からなります。

括弧の (T) は、3 電子励起（triples）を反復せずに摂動で入れたことを表します。

**ここまでで分かったこと**：CCSD(T) は、CCSD の振幅を使って、連結した 3 電子励起の寄与を摂動論で一度だけ見積もり、足す方法です。

## 5. 計算コスト

占有軌道の数を $o$、仮想軌道の数を $v$ とします。系の大きさ $N$ に対して $o, v \propto N$ です。

| 手法 | 最も重い計算 | 次数 | 反復 |
|---|---|---|---|
| MP2 | 積分の変換 | $N^5$ | なし |
| CCSD | $o^2 v^4$ | $N^6$ | あり |
| (T) | $o^3 v^4$ | $N^7$ | なし（1 回） |
| CCSDT | $o^3 v^5$ | $N^8$ | あり |

CCSD(T) の計算量は $N^7$ です。系の大きさが 2 倍になると、計算時間はおよそ $2^7 = 128$ 倍になります。また、必要なメモリやディスクも $o^2 v^2$ 以上で増えます。このため、通常の CCSD(T) を使えるのは、比較的小さな分子に限られます。大きな分子に CCSD(T) を使うための工夫が、第 11 回の局所相関法です。

**ここまでで分かったこと**：CCSD は $N^6$ の反復計算、(T) は $N^7$ の一度きりの計算です。

## 6. 「ゴールドスタンダード」と、その限界

### 計算で比べる

第 8 回と同じ H₂O（6-31G 基底）で、CCSD と CCSD(T) を FCI と比べます。結合の長さを平衡の 1.0 倍から 2.5 倍まで変えます。最後の二つの列（T1、D1）は 7 節で使います。

```python
import numpy as np
from pyscf import gto, scf, cc, fci

def h2o(scale):
    """平衡構造（O–H 0.96 Å、HOH 104.5°）の O–H 結合を scale 倍に伸ばす。"""
    r, th = 0.96 * scale, np.deg2rad(104.5) / 2
    return f"O 0 0 0; H 0 {r*np.sin(th):.6f} {r*np.cos(th):.6f}; H 0 {-r*np.sin(th):.6f} {r*np.cos(th):.6f}"

print("scale  CCSD-FCI  CCSD(T)-FCI    T1     D1   (単位 mEh)")
for scale in [1.0, 1.5, 2.0, 2.5]:
    mol = gto.M(atom=h2o(scale), basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    e_t = mycc.ccsd_t()
    e_fci = fci.FCI(mf).kernel()[0]
    d_ccsd = (mycc.e_tot - e_fci) * 1000
    d_t = (mycc.e_tot + e_t - e_fci) * 1000
    print(f"{scale:4.1f} {d_ccsd:9.2f} {d_t:11.2f}"
          f"  {mycc.get_t1_diagnostic():6.3f} {mycc.get_d1_diagnostic():6.3f}")
```

```
scale  CCSD-FCI  CCSD(T)-FCI    T1     D1   (単位 mEh)
 1.0      1.53        0.53   0.006  0.013
 1.5      5.60        1.10   0.019  0.051
 2.0      9.92       -5.05   0.035  0.093
 2.5     -4.67      -44.63   0.052  0.120
```

FCI の計算に時間がかかるので、全体で数分かかります（T1、D1 の列は、T1 と D1 の診断の値で、単位はありません）。

### 平衡構造では非常に正確

平衡構造では、CCSD(T) と FCI の差は 0.53 m$E_\mathrm{h}$（約 0.3 kcal/mol）です。第 8 回の CISD（6.8 m$E_\mathrm{h}$）や、第 9 回の MP2（8.1 m$E_\mathrm{h}$）よりはるかに小さく、CCSD（1.53 m$E_\mathrm{h}$）と比べても (T) で誤差が 3 分の 1 になっています。

このように、1 個の決定子が波動関数の大部分を占める分子（単参照性の強い分子）では、CCSD(T) は同じ基底関数の FCI に非常に近い値を与えます。しかもサイズ無矛盾で、$N^7$ の計算量で済みます。このため CCSD(T) は、ほかの手法の精度を測る基準として広く使われ、「ゴールドスタンダード」と呼ばれることがあります。

ただし、ここで比べたのは「同じ基底関数の中での FCI との差」です。実験値に近い値を得るには、基底関数も十分に大きくする（完全基底の極限に外挿する）必要があります。小さな基底関数での CCSD(T) は、手法としては正確でも、答えとしては不正確です。

### 結合を伸ばすと破綻する

結合を伸ばすと、様子が変わります。

- 1.5 倍では、CCSD(T) の誤差は 1.1 m$E_\mathrm{h}$ で、まだ良い精度です。
- 2.0 倍では、CCSD(T) のエネルギーが FCI より 5 m$E_\mathrm{h}$ **低く**なります。変分的でないので、こうしたことが起こります。
- 2.5 倍では、CCSD(T) は FCI より 45 m$E_\mathrm{h}$（約 28 kcal/mol）も低くなり、まったく信頼できません。

原因は、第 7 回の静的相関です。結合を伸ばすと HF の決定子の重みが下がり（第 8 回の $c_0^2$ を参照）、振幅が大きくなります。CCSD はまだある程度対応できますが、(T) の補正は摂動論なので、第 9 回の MP 級数と同じく、軌道エネルギーの差が小さい状況で破綻します。結合の解離、ビラジカル、多くの遷移金属錯体のように**多参照性**（複数の配置が大きな重みで混ざる性質）が強い系では、CCSD(T) を「ゴールドスタンダード」として信じてはいけません。

**ここまでで分かったこと**：単参照性の強い分子では、CCSD(T) は FCI に非常に近く、基準として使われます。多参照性が強くなると、(T) は破綻し、FCI より大きく下にずれることもあります。

## 7. T1・D1 診断

CCSD(T) が信頼できるかどうかを、FCI と比べずに見分ける目安があります。代表的なものが、CCSD の 1 電子励起の振幅 $t_i^a$ を使う二つの診断です。ここでは閉殻の分子で、$t_i^a$ を空間軌道の振幅とします。

**T1 診断**は、$t_i^a$ を並べたベクトルの長さを、相関させた電子の数 $N$ の平方根で割ったものです [Lee 1989]（書誌は末尾の「関連リンク」にあります）。

$$
\mathcal T_1 = \sqrt{\frac{\sum_{ia} (t_i^a)^2}{N}}
$$

多参照性が強いと、HF の軌道が良い出発点でなくなり、それを補うために 1 電子励起の振幅が大きくなります。Lee と Taylor は、$\mathcal T_1$ が 0.02 を超えると、単参照の手法の信頼性が下がり、多参照の手法が必要になる可能性が高いとしました。

**D1 診断**は、$t_i^a$ を行列と見たときの最大の特異値（$\mathbf t \mathbf t^\mathsf T$ の最大固有値の平方根）です [Janssen 1998]。$\mathcal T_1$ は全体の平均なので、分子の一部だけで振幅が大きいと、大きな分子では薄まってしまいます。D1 は最大の値を見るので、この薄まりが起きにくい量です。

6 節の表で、二つの診断の変化を見ます。$\mathcal T_1$ は平衡構造の 0.006 から、1.5 倍で 0.019、2.0 倍で 0.035 に増えます。0.02 を超えたあたりから、CCSD(T) が FCI からずれ始めていることが分かります。D1 も同じように増えています。

ただし、診断の目安の値は、系の種類によって違います。たとえば 3d 遷移金属を含む分子については、$\mathcal T_1 > 0.05$、D1 $> 0.15$ などを、多参照性が強いことの目安とする提案があります [Jiang 2012]。同じ論文では、これらの値を下回っていても多参照性が強い分子（Cr₂ など）があることも指摘されています。診断は「危ないかもしれない」という警告として使い、診断の値が小さいことを、正しさの保証とは考えないでください。

**ここまでで分かったこと**：T1・D1 診断は、CCSD の 1 電子励起の振幅から多参照性を見積もる目安です。目安の値は系の種類で違い、値が小さくても安全とは限りません。

## まとめ

- CC 法は $e^{\hat T}\Phi_0$ の形で波動関数を表し、$\hat T$ を打ち切ってもサイズ無矛盾です。
- 振幅は相似変換したハミルトニアンへの射影で決め、エネルギーは変分的ではありません。
- CCSD(T) は、CCSD（$N^6$）に、連結した 3 電子励起の摂動補正（$N^7$）を足した方法です。
- 単参照性の強い分子では FCI に非常に近く「ゴールドスタンダード」と呼ばれますが、多参照性が強いと破綻します。
- T1・D1 診断は多参照性の目安で、警告として使います。

## 参考文献

- T. Helgaker, P. Jørgensen, J. Olsen, "Molecular Electronic-Structure Theory" (Wiley, 2000). ISBN: 9780471967552
- J. Čížek, "On the Correlation Problem in Atomic and Molecular Systems. Calculation of Wavefunction Components in Ursell-Type Expansion Using Quantum-Field Theoretical Methods", J. Chem. Phys. 45, 4256 (1966). DOI: [10.1063/1.1727484](https://doi.org/10.1063/1.1727484)
- G. D. Purvis III, R. J. Bartlett, "A full coupled-cluster singles and doubles model: The inclusion of disconnected triples", J. Chem. Phys. 76, 1910 (1982). DOI: [10.1063/1.443164](https://doi.org/10.1063/1.443164)
- K. Raghavachari, G. W. Trucks, J. A. Pople, M. Head-Gordon, "A fifth-order perturbation comparison of electron correlation theories", Chem. Phys. Lett. 157, 479 (1989). DOI: [10.1016/S0009-2614(89)87395-6](https://doi.org/10.1016/S0009-2614(89)87395-6)
- T. D. Crawford, H. F. Schaefer III, "An Introduction to Coupled Cluster Theory for Computational Chemists", Rev. Comput. Chem. 14, 33 (2000). DOI: [10.1002/9780470125915.ch2](https://doi.org/10.1002/9780470125915.ch2)
- R. J. Bartlett, M. Musiał, "Coupled-cluster theory in quantum chemistry", Rev. Mod. Phys. 79, 291 (2007). DOI: [10.1103/RevModPhys.79.291](https://doi.org/10.1103/RevModPhys.79.291)
- C. L. Janssen, I. M. B. Nielsen, "New diagnostics for coupled-cluster and Møller–Plesset perturbation theory", Chem. Phys. Lett. 290, 423 (1998). DOI: [10.1016/S0009-2614(98)00504-1](https://doi.org/10.1016/S0009-2614(98)00504-1)
- W. Jiang, N. J. DeYonker, A. K. Wilson, "Multireference Character for 3d Transition-Metal-Containing Molecules", J. Chem. Theory Comput. 8, 460 (2012). DOI: [10.1021/ct2006852](https://doi.org/10.1021/ct2006852)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- T. J. Lee, P. R. Taylor, "A diagnostic for determining the quality of single-reference electron correlation methods", Int. J. Quantum Chem. Symp. 23, 199 (1989)（NASA 技術報告サーバーの記録）：https://ntrs.nasa.gov/citations/19900002238
