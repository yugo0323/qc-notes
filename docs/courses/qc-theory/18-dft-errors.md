# 第 18 回　DFT の誤差：自己相互作用誤差と静的相関

!!! abstract "この回の道筋"
    近似の汎関数の誤差の多くは、2 つの「厳密な条件」を破ることから説明できます。電子の数が分数のときにエネルギーが直線になる条件と、スピンが分数のときにエネルギーが一定になる条件です。この回では、1 電子の自己相互作用誤差から始めて、この 2 つの条件と、それを破ったときに起こる失敗（H₂⁺ と H₂ の解離）を、計算で確かめます。最後に、誤差を「汎関数のせい」と「密度のせい」に分ける考え方（DC-DFT）を紹介します。

    | 節 | 内容 |
    |---|---|
    | 1 | 1 電子の自己相互作用誤差 |
    | 2 | 分数電荷とエネルギーの直線性 |
    | 3 | 多電子の自己相互作用誤差と非局在化誤差 |
    | 4 | 計算例：H₂⁺ の解離 |
    | 5 | 分数スピンと静的相関誤差 |
    | 6 | 計算例：H₂ の解離 |
    | 7 | 密度駆動誤差と汎関数駆動誤差 |

    **前に読んでおく回**：量子化学の理論 [第 7 回](07-electron-correlation.md)、[第 15 回](15-kohn-sham.md)、[第 16 回](16-functionals.md)

!!! tip "物理的なイメージ"
    この回で調べるのは、近似の汎関数がもつ「電子が自分自身と反発してしまう」誤り（自己相互作用誤差）とその仲間です。電子 1 個の系でも、近似の汎関数では、クーロン反発と交換のエネルギーが正しく打ち消し合わず、電子が自分自身の分布と反発します。その結果、エネルギーは電子の数について下に凸に曲がり、電子をいくつかの場所に広げた状態を好みすぎます（H₂⁺ の解離の失敗）。また、スピンが分数の状態を正しく扱えないと、H₂ の解離のように静的相関の強い状態で、エネルギーが高くなりすぎます。

この回では原子単位を使います。エネルギーの換算は 1 $E_\mathrm h$ = 27.211 eV = 627.5 kcal/mol です。
計算はすべて PySCF [Sun 2020] で、基底関数は aug-cc-pVQZ です。

## 1. 1 電子の自己相互作用誤差

[第 15 回](15-kohn-sham.md)で見たように、Kohn–Sham 法のエネルギーには、古典的な静電エネルギー

$$
J[\rho] = \frac12\iint\frac{\rho(\mathbf r_1)\rho(\mathbf r_2)}{r_{12}}\,d\mathbf r_1d\mathbf r_2
$$

が入っています。
$J$ は、電子が自分自身の雲と反発する分も含みます。
電子が 1 個の系（H 原子や H₂⁺）では、電子どうしの反発は本当は 0 です。
したがって、厳密な汎関数は、どんな 1 電子の密度についても

$$
E_\mathrm{xc}[\rho] = -J[\rho]
$$

を満たします。
自分自身との反発を、交換相関エネルギーがちょうど打ち消すということです。
Hartree–Fock 法では、第 3 回で見たとおり、交換の項がこれを自動的に打ち消します。

ところが、LDA や GGA の $E_\mathrm{xc}$ は、一般にこの条件を満たしません。
残った $J + E_\mathrm{xc}$ を、1 電子の自己相互作用誤差（self-interaction error、SIE）と呼びます。
Perdew と Zunger は、この誤差を軌道ごとに差し引く補正（自己相互作用補正、SIC）を提案しました [Perdew 1981]。

電子がちょうど 1 個の H 原子の全エネルギーだけを見ると、この誤差は目立たないことがあります。
誤差がはっきり現れるのは、電子が分数個あるように見える状況です。
次の節で、その意味を説明します。

**ここまでで分かったこと**：厳密な汎関数では、1 電子の系で $E_\mathrm{xc} = -J$ となり、自分自身との反発が打ち消されます。近似の汎関数ではこれが破れ、自己相互作用誤差が残ります。

## 2. 分数電荷とエネルギーの直線性

### 2.1 直線性の条件

電子の数 $N$ を分数に広げます。
分数の電子数は、たとえば、遠く離れた 2 つの原子の間で電子 1 個を分け合っているとき、それぞれの原子について現れます。
Perdew、Parr、Levy、Balduz は、電子の数が $N_0$ と $N_0 + 1$ の状態を $1-\omega$ と $\omega$ の割合で混ぜた統計的な混合（アンサンブル）を考えました [Perdew 1982]。
その結果、厳密なエネルギーは、整数の電子数の間で直線になります。

$$
E(N_0 + \omega) = (1 - \omega)\,E(N_0) + \omega\,E(N_0 + 1),\qquad 0 \le \omega \le 1
$$

直線の傾きは、$N_0 - 1 < N < N_0$ では $-I$（$I$ はイオン化エネルギー）、$N_0 < N < N_0 + 1$ では $-A$（$A$ は電子親和力）です。
傾きは整数の $N_0$ で不連続に変わります。

[第 15 回](15-kohn-sham.md)の Janak の定理によれば、$\partial E/\partial n_\mathrm{HOMO} = \varepsilon_\mathrm{HOMO}$ でした。
エネルギーが直線なら、傾きは一定なので、$\varepsilon_\mathrm{HOMO}$ は占有数によらず $-I$ に等しくなります。
[第 15 回](15-kohn-sham.md)の $-\varepsilon_\mathrm{HOMO} = I$ は、この直線性から理解できます。

### 2.2 H 原子の E(N) を計算する

H 原子で、1s の α 電子の占有数 $n$ を 0 から 1 まで変えて、エネルギーを計算します。
[第 15 回](15-kohn-sham.md)と同じように、占有数を固定して SCF を解きます。
PySCF は電子 1 個の系を特別な近道で解くので、ここでは通常の SCF の手順（`scf.hf.kernel`）を呼びます。

```python
import numpy as np
from pyscf import gto, scf, dft

basis = "aug-cc-pvqz"
methods = ["HF", "PBE", "B3LYP", "WB97X-V"]
h = gto.M(atom="H 0 0 0", basis=basis, spin=1, verbose=0)

def e_occ(method, na, nb=0.0):
    """1s の α、β の占有数を na、nb に固定してエネルギーを求める。"""
    mf = scf.UHF(h) if method == "HF" else dft.UKS(h)
    if method != "HF":
        mf.xc = method
    def get_occ(mo_energy=None, mo_coeff=None):
        occ = np.zeros((2, h.nao))
        occ[0, 0], occ[1, 0] = na, nb
        return occ
    mf.get_occ = get_occ
    return scf.hf.kernel(mf)[1]       # 1 電子系の近道を使わず、SCF で解く

print("n     " + "".join(f"{m:>10s}" for m in methods))
for n in [0.25, 0.5, 0.75, 1.0]:
    print(f"{n:4.2f} " + "".join(f"{e_occ(m, n):10.5f}" for m in methods))
```

```text
n             HF       PBE     B3LYP   WB97X-V
0.25   -0.12499  -0.16142  -0.15535  -0.13473
0.50   -0.24997  -0.30319  -0.29432  -0.26618
0.75   -0.37496  -0.41611  -0.41018  -0.38847
1.00   -0.49995  -0.49993  -0.50239  -0.50074
```

$n = 0$ ではエネルギーは 0 です。
電子が 1 個以下の H 原子の厳密なエネルギーは、$E(n) = -0.5\,n$ の直線です。

- **HF**：1 電子の系では自己相互作用がないので、直線になります。
- **PBE、B3LYP**：直線より下に曲がっています（下に凸、convex）。$n = 1/2$ では、両端を結んだ直線より、PBE で $0.053\,E_\mathrm h$（1.45 eV）、B3LYP で $0.043\,E_\mathrm h$ 低くなります。
- **ωB97X-V**：長距離の交換を厳密な交換にしているので、曲がりは $0.016\,E_\mathrm h$ まで小さくなりますが、消えてはいません。

PBE は、$n = 1$ の H 原子ではほぼ正しいのに、分数の電子では大きく誤ります。
分数の電子の雲は、自分自身との反発 $J$ が $n^2$ に比例して小さくなります。
近似の $E_\mathrm{xc}$ は、それを正しい割合で打ち消さないので、エネルギーが曲がります。
[第 15 回](15-kohn-sham.md)の H₂O の計算で、PBE の $E(n)$ が曲がっていたのも同じ理由です。

**ここまでで分かったこと**：厳密なエネルギーは、整数の電子数の間で直線です。PBE や B3LYP のエネルギーは下に凸に曲がり、分数の電子をもつ状態を低く見積もりすぎます。

## 3. 多電子の自己相互作用誤差と非局在化誤差

前の節の H 原子は 1 電子の系でした。
電子が多い系でも、同じ曲がりが問題になります。
Mori-Sánchez、Cohen、Yang は、1 電子の自己相互作用誤差を除いた汎関数でも、多電子の系では分数電荷のエネルギーが直線にならないことを示し、これを多電子の自己相互作用誤差と呼びました [Mori-Sánchez 2006]。
直線性からのずれが、自己相互作用誤差を多電子に広げた定義になります。

エネルギーが下に凸だと、電子を 1 か所に置くより、いくつかの場所に分数ずつ広げた方がエネルギーが低くなります。
そのため、近似の汎関数は、電子を広げすぎる傾向を持ちます。
これを非局在化誤差（delocalization error）と呼びます [Mori-Sánchez 2008]。
反対に、HF は上に凸（concave）で、電子を 1 か所に集めすぎる局在化誤差を持ちます [Mori-Sánchez 2008]。

非局在化誤差は、次のような失敗の原因になります。

- **分数の電荷への解離**：半局所汎関数では、中性の分子 XY を引き離すと、X$^{+q}$ と Y$^{-q}$ のように分数の電荷を帯びた断片に分かれることがあります。たとえば NaCl は Na$^{+0.4}$ と Cl$^{-0.4}$ に解離し得ます [Ruzsinszky 2006]。
- **バンドギャップの過小評価**：固体のバンドギャップの問題も、分数電荷のエネルギーの曲がりから説明されます [Mori-Sánchez 2008]。

Cohen、Mori-Sánchez、Yang は、この非局在化誤差と、第 5 節で扱う静的相関誤差が、DFT の多くの目立った失敗の原因であり、分数電荷と分数スピンの見方で理解できると述べています [Cohen 2008a]。

**ここまでで分かったこと**：分数電荷のエネルギーが下に凸だと、電子を広げすぎる非局在化誤差が生じます。半局所汎関数と多くの混成汎関数はこの誤差を持ち、HF は逆の局在化誤差を持ちます。

## 4. 計算例：H₂⁺ の解離

非局在化誤差がいちばん単純に現れるのが、H₂⁺ の解離です。
H₂⁺ は電子が 1 個なので、HF はこの基底関数の中での厳密解です。
原子核の間の距離 $R$ を大きくすると、正しい答えは H 原子 + H⁺ で、エネルギーは $-0.5\,E_\mathrm h$ に近づきます。

ここでは、分子の対称性を保った解（電子が両方の核に半分ずつある解）を求めます。
離れた距離では、DFT の SCF が電子を片方に寄せた解に落ちることがあるので、`symmetry=True` で対称性を課します。

```python
from pyscf import gto, scf, dft

def make_mf(mol, method):
    if method == "HF":
        return scf.UHF(mol)
    mf = dft.UKS(mol); mf.xc = method
    return mf

print("R/Å  " + "".join(f"{m:>10s}" for m in methods))
for R in [0.74, 1.06, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0]:
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {R}", basis=basis, charge=1, spin=1,
                symmetry=True, verbose=0)
    print(f"{R:4.2f} " + "".join(f"{make_mf(mol, m).run().e_tot:10.5f}" for m in methods))
```

```text
R/Å          HF       PBE     B3LYP   WB97X-V
0.74   -0.56961  -0.56960  -0.57184  -0.56940
1.06   -0.60253  -0.60927  -0.61008  -0.60741
2.00   -0.55220  -0.58297  -0.57937  -0.57084
3.00   -0.51516  -0.57127  -0.56296  -0.54331
4.00   -0.50355  -0.57497  -0.56398  -0.53505
6.00   -0.50019  -0.58439  -0.57109  -0.53255
8.00   -0.49999  -0.58985  -0.57543  -0.53240
10.00   -0.49997  -0.59315  -0.57806  -0.53237
```

HF（厳密解）は、$R$ とともに $-0.5\,E_\mathrm h$ に近づきます。
一方、PBE と B3LYP は、離れるほどかえってエネルギーが下がり、10 Å では厳密解より $0.09\,E_\mathrm h$（PBE、2.5 eV）や $0.08\,E_\mathrm h$（B3LYP）も低くなります。
平衡距離 1.06 Å でも、PBE は $0.007\,E_\mathrm h$（4.2 kcal/mol）低すぎます。

離れた H₂⁺ は、電子を半分ずつ持つ 2 個の H 原子です。
したがって、エネルギーは前の節の $E(1/2)$ の 2 倍に、半分の電荷どうしの静電的な相互作用を足したものになるはずです。
核の間の距離を $R$ とすると、核どうしの反発 $1/R$、各電子雲と相手の核の引力 $2\times(-\frac12\cdot\frac1R)$、2 つの電子雲の反発 $\frac14\cdot\frac1R$ を足して、合計 $+\frac1{4R}$ です。
PBE で $R = 10$ Å（18.90 bohr）を計算すると、

$$
2E_\mathrm{PBE}(1/2) + \frac{1}{4R} = 2\times(-0.30319) + 0.01323 = -0.59315\,E_\mathrm h
$$

となり、表の値と一致します。
H₂⁺ の解離の誤差は、H 原子の $E(n)$ の曲がりそのものです。
曲がりの小さい ωB97X-V では、誤差も小さくなっています。

なお、PBE では、電子を片方の核に寄せた解（正しい H + H⁺）より、両方に半分ずつ広げた解の方がエネルギーが低くなります。
対称性を課さなかった最初の試みでは、6 Å 以上で SCF が寄せた解に落ち、$-0.4993\,E_\mathrm h$ 付近の値を与えました。
これは、PBE のエネルギーの最小ではなく、たまたま SCF が止まった別の解です。
汎関数の最小が誤った（広がりすぎた）密度を選ぶことが、非局在化誤差の本質です。

**ここまでで分かったこと**：H₂⁺ の解離では、PBE や B3LYP は電子を両方の核に半分ずつ広げた状態を低く見積もりすぎ、エネルギーが下がりすぎます。その大きさは、分数電荷のエネルギーの曲がりで定量的に説明できます。

## 5. 分数スピンと静的相関誤差

次に、電子の数は整数のまま、スピンが分数になる状況を考えます。
H 原子で、α 電子 1/2 個と β 電子 1/2 個が 1s に入った状態です。
これは、α スピンの H 原子と β スピンの H 原子を半分ずつ混ぜたアンサンブルで、どちらも同じエネルギーの基底状態です。
Cohen、Mori-Sánchez、Yang は、このような分数スピンの状態について、厳密な汎関数のエネルギーは、混ぜたもとの状態のエネルギーと同じで一定であることを証明しました [Cohen 2008b]。
H 原子なら、$-0.5\,E_\mathrm h$ のままです。
これを一定性の条件（constancy condition）と呼びます。

分数スピンの状態は、静的相関（第 7 回）の強い系で自然に現れます。
H₂ を解離させると、正しい状態は 2 個の H 原子で、それぞれのスピンは α とも β とも決まりません。
スピンを制限した（RKS の）計算では、密度は各原子で α と β が半分ずつ、つまり分数スピンの H 原子 2 個になります。
したがって、解離の極限の RKS のエネルギーは、分数スピンの H 原子のエネルギーの 2 倍です。

近似の汎関数は、一定性の条件を大きく破ります [Cohen 2008b]。
その誤差を静的相関誤差（static correlation error）と呼びます。
結合の解離や、Mott 絶縁体のバンド構造などの強相関の系で、大きな誤差の原因になります。

**ここまでで分かったこと**：厳密な汎関数では、分数スピンの状態のエネルギーは一定です。近似の汎関数はこれを破り、結合の解離などで静的相関誤差を生みます。

## 6. 計算例：H₂ の解離

分数スピンの H 原子のエネルギーと、スピンを制限した H₂ の解離曲線を計算します。
`e_occ` は第 2.2 節の関数です。
離れた距離では、σg 軌道と σu 軌道のエネルギーがほぼ等しくなり、SCF が収束しにくくなります。
そこで対称性を使い、2 個の電子を σg 軌道に入れることを指定します（`irrep_nelec`）。

```python
from pyscf import fci

print("H 原子     " + "".join(f"{m:>10s}" for m in methods))
print("(1, 0)     " + "".join(f"{e_occ(m, 1.0, 0.0):10.5f}" for m in methods))
print("(1/2, 1/2) " + "".join(f"{e_occ(m, 0.5, 0.5):10.5f}" for m in methods))

def rks(mol, method):
    mf = scf.RHF(mol) if method == "HF" else dft.RKS(mol)
    if method != "HF":
        mf.xc = method
    mf.irrep_nelec = {"A1g": 2}       # 2 個の電子を σg 軌道に入れる
    return mf.run()

print("H2 (R/Å)  " + "".join(f"{m:>10s}" for m in methods) + "       FCI")
for R in [0.74, 2.0, 4.0, 10.0]:
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {R}", basis=basis, symmetry=True, verbose=0)
    es = [rks(mol, m).e_tot for m in methods]
    mol_c1 = gto.M(atom=f"H 0 0 0; H 0 0 {R}", basis=basis, verbose=0)   # FCI は対称性なしで
    e_fci = fci.FCI(scf.RHF(mol_c1).run()).kernel()[0]
    print(f"{R:5.2f}     " + "".join(f"{e:10.5f}" for e in es) + f"{e_fci:10.5f}")
```

```text
H 原子             HF       PBE     B3LYP   WB97X-V
(1, 0)       -0.49995  -0.49993  -0.50239  -0.50074
(1/2, 1/2)   -0.35764  -0.45888  -0.44805  -0.39315
H2 (R/Å)          HF       PBE     B3LYP   WB97X-V       FCI
 0.74       -1.13348  -1.16653  -1.18054  -1.17308  -1.17386
 2.00       -0.92633  -0.99858  -1.00189  -0.98326  -1.02193
 4.00       -0.79245  -0.92390  -0.91582  -0.86133  -0.99999
10.00       -0.74177  -0.91775  -0.90138  -0.81277  -0.99990
```

まず H 原子です。
α と β が 1/2 ずつの状態のエネルギーは、厳密には $-0.5\,E_\mathrm h$ のはずです。
ところが、どの方法でもエネルギーが高くなり、PBE で $0.041\,E_\mathrm h$（1.1 eV）、B3LYP で $0.052\,E_\mathrm h$、ωB97X-V で $0.107\,E_\mathrm h$、HF で $0.142\,E_\mathrm h$ も高すぎます。
この表では、厳密な交換を多く含む方法ほど、誤差が大きくなっています。

次に H₂ です。
FCI は、この基底関数の中での厳密解で、10 Å では H 原子 2 個のエネルギー $-1.0\,E_\mathrm h$ に近づきます。
スピンを制限した計算は、どれもそこまで下がりません。
PBE の 10 Å の値 $-0.91775\,E_\mathrm h$ は、分数スピンの H 原子のエネルギーの 2 倍 $2\times(-0.45888) = -0.91776\,E_\mathrm h$ と一致します。
解離の極限での誤差は、分数スピンの H 原子の誤差そのものです。
B3LYP と ωB97X-V も、ほぼ同じ関係になっています。
ただし、厳密な交換は 2 つの原子にまたがる軌道を通して離れた点を結ぶので、一致は PBE ほど正確ではありません。

HF の誤差は、第 7 回で見たとおり、イオン型の配置が半分混ざることによるものです。
DFT の誤差はそれより小さいものの、PBE でも $0.082\,E_\mathrm h$（2.2 eV）残ります。
平衡距離の付近（0.74 Å）では、DFT は FCI に近い値を与えています。
静的相関誤差は、結合が切れて静的相関が強くなるほど大きくなります。

第 7 回で見たように、スピンの対称性を破った解（UKS）を使えば、解離の極限のエネルギーは正しくなります。
しかし、それは各原子に α か β の電子を整数個ずつ置くことで、分数スピンの状態を避けているだけです。
そのかわり、スピン汚染が生じます。
スピンを保ったまま正しいエネルギーを与えるには、一定性の条件を満たす汎関数が必要です。Cohen らは、この条件に近づけることを汎関数の開発の大きな目標にすべきだと述べています [Cohen 2008b]。

**ここまでで分かったこと**：スピンを制限した H₂ の解離の極限のエネルギーは、分数スピンの H 原子のエネルギーの 2 倍です。近似の汎関数は分数スピンのエネルギーを高く見積もるので、解離の極限で大きな誤差（静的相関誤差）が出ます。

## 7. 密度駆動誤差と汎関数駆動誤差

### 7.1 誤差を 2 つに分ける

近似の汎関数 $\tilde E$ で SCF を解くと、エネルギーだけでなく密度 $\tilde\rho$ も近似になります。
Kim、Sim、Burke は、エネルギーの誤差を次の 2 つに分けました [Kim 2013]。

$$
\Delta E = \tilde E[\tilde\rho] - E[\rho]
= \underbrace{\tilde E[\rho] - E[\rho]}_{\Delta E_\mathrm F}
+ \underbrace{\tilde E[\tilde\rho] - \tilde E[\rho]}_{\Delta E_\mathrm D}
$$

$\rho$ は厳密な密度、$E[\rho]$ は厳密なエネルギーです。

- **汎関数駆動誤差** $\Delta E_\mathrm F$：正しい密度を使っても残る、汎関数そのものの誤差です。
- **密度駆動誤差** $\Delta E_\mathrm D$：密度が正しくないことから来る誤差です。$\tilde\rho$ は $\tilde E$ を最小にする密度なので、$\Delta E_\mathrm D \le 0$ です。

ふつうの計算では、汎関数駆動誤差が大部分を占めます。
しかし、密度駆動誤差が大きくなる「異常な」場合があります。
Kim らは、その例として、H⁻ などの陰イオンの電子親和力、引き離した NaCl（PBE では Na$^{0.4+}$ と Cl$^{0.4-}$ に分かれます）、溶液中のラジカルの構造を挙げています。
また、HOMO と LUMO の差（軌道のギャップ）が小さいことが、密度駆動誤差が大きいことの目安になるとしています。

### 7.2 DC-DFT

密度駆動誤差が大きい場合は、より正確な密度で汎関数を評価すれば誤差が減ります。
これを密度補正 DFT（density-corrected DFT、DC-DFT）と呼びます。
Kim らは、HF の密度で DFT のエネルギーを評価する方法（HF-DFT）を提案しました [Kim 2013]。
自己相互作用誤差で密度が広がりすぎる場合、HF の密度は、近似の汎関数の密度より正確なことが多いためです。
Sim、Song、Burke は、厳密な密度が分からなくても、密度駆動誤差がどれくらい重要かを見積もる方法を示しました [Sim 2018]。
また、密度の誤差が大きい場合に混成汎関数の厳密な交換の割合を調節することがありますが、DC-DFT を使えばその必要はないと述べています。

### 7.3 H₂⁺ の誤差を分ける

H₂⁺ は電子が 1 個なので、HF の密度がこの基底関数の中での厳密な密度です。
そこで、PBE の誤差を 2 つに分けてみます。
`ks.energy_tot(dm=...)` で、PBE のエネルギーを HF の密度行列で評価できます。

```python
print(" R/Å   E(PBE)    ΔE      ΔE_F     ΔE_D   （単位 Eh）")
for R in [1.06, 4.0, 10.0]:
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {R}", basis=basis, charge=1, spin=1,
                symmetry=True, verbose=0)
    hf = scf.UHF(mol).run()              # 1 電子系なので、HF の密度がこの基底での厳密な密度
    ks = dft.UKS(mol); ks.xc = "PBE"; ks.run()
    e_exact = hf.e_tot
    e_pbe_on_hf = ks.energy_tot(dm=hf.make_rdm1())   # PBE の汎関数を HF の密度で評価
    dE   = ks.e_tot - e_exact
    dE_F = e_pbe_on_hf - e_exact
    dE_D = ks.e_tot - e_pbe_on_hf
    print(f"{R:5.2f} {ks.e_tot:9.5f} {dE:8.5f} {dE_F:8.5f} {dE_D:8.5f}")
```

```text
 R/Å   E(PBE)    ΔE      ΔE_F     ΔE_D   （単位 Eh）
 1.06  -0.60927 -0.00673 -0.00559 -0.00114
 4.00  -0.57497 -0.07141 -0.06566 -0.00576
10.00  -0.59315 -0.09318 -0.08829 -0.00489
```

どの距離でも、誤差の 9 割前後は汎関数駆動誤差 $\Delta E_\mathrm F$ です。
対称性を保った H₂⁺ では、正しい密度も両方の核に半分ずつあるので、PBE の密度の形はそれほど間違っていません。
誤差の原因は、半分ずつの電子の雲に対する $E_\mathrm{xc}$ の値そのもの、つまり第 2 節の $E(n)$ の曲がりです。
このような場合は、HF の密度を使っても誤差はほとんど減りません。
密度駆動誤差が大きいのは、第 7.1 節の NaCl のように、近似の汎関数が密度そのものを誤った形にしてしまう場合です。

**ここまでで分かったこと**：エネルギーの誤差は、汎関数駆動誤差と密度駆動誤差に分けられます。ふつうは前者が大きく、後者が大きい異常な場合には、HF の密度などを使う DC-DFT で誤差を減らせます。

## まとめ

- 厳密な汎関数は、1 電子の系で $E_\mathrm{xc} = -J$ を満たします。近似の汎関数ではこれが破れ、自己相互作用誤差が残ります。
- 厳密なエネルギーは、分数電荷について直線になります。半局所汎関数や B3LYP のエネルギーは下に凸に曲がり、電子を広げすぎる非局在化誤差を生みます。H₂⁺ の解離での大きな誤差は、この曲がりで定量的に説明できます。
- 厳密なエネルギーは、分数スピンについて一定です。近似の汎関数はこれを破り、H₂ の解離のような静的相関の強い系で誤差（静的相関誤差）を生みます。
- エネルギーの誤差は、汎関数駆動誤差と密度駆動誤差に分けられます。密度駆動誤差が大きいときは、HF の密度などで評価する DC-DFT が有効です。

## 参考文献

- J. P. Perdew, A. Zunger, "Self-interaction correction to density-functional approximations for many-electron systems", Phys. Rev. B 23, 5048 (1981). DOI: [10.1103/PhysRevB.23.5048](https://doi.org/10.1103/PhysRevB.23.5048)
- J. P. Perdew, R. G. Parr, M. Levy, J. L. Balduz Jr., "Density-Functional Theory for Fractional Particle Number: Derivative Discontinuities of the Energy", Phys. Rev. Lett. 49, 1691 (1982). DOI: [10.1103/PhysRevLett.49.1691](https://doi.org/10.1103/PhysRevLett.49.1691)
- P. Mori-Sánchez, A. J. Cohen, W. Yang, "Many-electron self-interaction error in approximate density functionals", J. Chem. Phys. 125, 201102 (2006). DOI: [10.1063/1.2403848](https://doi.org/10.1063/1.2403848)
- P. Mori-Sánchez, A. J. Cohen, W. Yang, "Localization and Delocalization Errors in Density Functional Theory and Implications for Band-Gap Prediction", Phys. Rev. Lett. 100, 146401 (2008). DOI: [10.1103/PhysRevLett.100.146401](https://doi.org/10.1103/PhysRevLett.100.146401)
- A. Ruzsinszky, J. P. Perdew, G. I. Csonka, O. A. Vydrov, G. E. Scuseria, "Spurious fractional charge on dissociated atoms: Pervasive and resilient self-interaction error of common density functionals", J. Chem. Phys. 125, 194112 (2006). DOI: [10.1063/1.2387954](https://doi.org/10.1063/1.2387954)
- A. J. Cohen, P. Mori-Sánchez, W. Yang, "Insights into Current Limitations of Density Functional Theory", Science 321, 792 (2008). DOI: [10.1126/science.1158722](https://doi.org/10.1126/science.1158722)
- A. J. Cohen, P. Mori-Sánchez, W. Yang, "Fractional spins and static correlation error in density functional theory", J. Chem. Phys. 129, 121104 (2008). DOI: [10.1063/1.2987202](https://doi.org/10.1063/1.2987202)
- M.-C. Kim, E. Sim, K. Burke, "Understanding and Reducing Errors in Density Functional Calculations", Phys. Rev. Lett. 111, 073003 (2013). DOI: [10.1103/PhysRevLett.111.073003](https://doi.org/10.1103/PhysRevLett.111.073003)
- E. Sim, S. Song, K. Burke, "Quantifying Density Errors in DFT", J. Phys. Chem. Lett. 9, 6385 (2018). DOI: [10.1021/acs.jpclett.8b02855](https://doi.org/10.1021/acs.jpclett.8b02855)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
