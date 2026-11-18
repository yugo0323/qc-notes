# 第 17 回　分散力の補正

!!! abstract "この回の道筋"
    分散力（London 分散力）は、すべての原子と分子の間に働く弱い引力です。Hartree–Fock 法と半局所汎関数は、この力を原理的に表せません。この回では、その理由を確かめたうえで、実際に使われている補正（DFT-D2/D3/D4、非局所汎関数 VV10、多体の分散）の考え方と式を紹介します。計算例では、Ar₂ の相互作用エネルギーを比べます。

    | 節 | 内容 |
    |---|---|
    | 1 | 分散力とは |
    | 2 | 半局所汎関数が分散力を表せない理由 |
    | 3 | 計算例：Ar₂ の相互作用エネルギー |
    | 4 | 原子対の補正：DFT-D2、D3、D4 |
    | 5 | 非局所汎関数：VV10 |
    | 6 | 多体の分散：ATM 項と MBD |

    **前に読んでおく回**：量子化学の理論 [第 15 回](15-kohn-sham.md)、[第 16 回](16-functionals.md)

!!! tip "物理的なイメージ"
    分散力が表しているのは、離れた原子どうしの電子の「瞬間的なゆらぎ」が、互いにそろうことで生まれる引力です。一方の原子の電子が一瞬かたよると瞬間的な双極子ができ、その電場で相手の原子の電子もかたよって、引き合います。時間で平均すると双極子は 0 ですが、二つのゆらぎの相関は 0 になりません。各点の近くの密度しか見ない半局所汎関数や、平均の場で考える HF 法は、この離れた電子の相関を原理的に表せないので、外から補正を足します。

この回では原子単位を使います。エネルギーの換算は 1 $E_\mathrm h$ = 27.211 eV = 627.5 kcal/mol = 219474.6 cm$^{-1}$ です。

## 1. 分散力とは

Ar 原子は電荷を持たず、永久双極子もありません。
それでも、Ar 原子どうしは弱く引き合い、低温では液体や固体になります。
この引力が分散力（dispersion force）です。
London が量子力学で説明したので、London 分散力とも呼ばれます [London 1937]。

分散力は、次のように理解できます。
原子の電子の雲は、平均すれば球対称ですが、各瞬間には少しゆがんでおり、瞬間的な双極子があります。
この瞬間の双極子が、もう一方の原子に双極子を誘起します。
2 つの双極子のゆらぎが互いにそろうように相関すると、エネルギーが下がります。
これは、離れた 2 つの原子の電子の間の相関です（第 7 回の言葉では、電子相関の一部です）。

原子が十分に離れているとき、分散のエネルギーは原子間距離 $R$ の 6 乗に反比例します。

$$
E_\mathrm{disp} \approx -\frac{C_6}{R^6}
$$

$C_6$ を分散係数と呼びます。
電子を 1 個ずつばねにつないだ簡単な模型（Drude 模型）で、この形を導けます。
導出は付録 A に回しますが、同じ 2 個の原子では次の結果になります。

$$
C_6 = \frac34\,\omega_0\,\alpha^2
$$

$\alpha$ は原子の分極率、$\omega_0$ はばねの振動数（原子単位では $\hbar = 1$ なので、エネルギーと同じ単位）です。
分極率の大きい、つまり電子の雲がゆがみやすい原子ほど、分散力は強くなります。
次の項として $C_8/R^8$、$C_{10}/R^{10}$ などもあります。

分散のエネルギーは、1 対ごとには小さくても、大きな分子では多数の原子の対について足し合わさります。
そのため、分子結晶の安定性、大きな分子の立体配座、タンパク質と薬の結合などで、大きな寄与になります。

**ここまでで分かったこと**：分散力は、離れた原子の電子のゆらぎの相関から生じる引力で、遠くでは $-C_6/R^6$ の形です。分極率の大きい原子ほど強く働きます。

## 2. 半局所汎関数が分散力を表せない理由

Hartree–Fock 法には、電子相関がありません。
したがって、分散力はまったく入りません。

半局所汎関数（LDA、GGA、meta-GGA）は相関を含みますが、やはり分散力を表せません。
理由は、交換相関エネルギーの形にあります（[第 16 回](16-functionals.md)）。

$$
E_\mathrm{xc}^\mathrm{SL}[\rho] = \int e_\mathrm{xc}\big(\rho(\mathbf r), \nabla\rho(\mathbf r), \tau(\mathbf r)\big)\,d\mathbf r
$$

各点 $\mathbf r$ のエネルギー密度は、その点の近くの密度だけで決まります。
いま、2 つの原子 A と B が十分に離れていて、電子の雲が重ならないとします。
すると、空間の各点では、A か B の一方の密度しかありません。
積分は A の領域と B の領域に分かれ、

$$
E_\mathrm{xc}^\mathrm{SL}[\rho_\mathrm A + \rho_\mathrm B] = E_\mathrm{xc}^\mathrm{SL}[\rho_\mathrm A] + E_\mathrm{xc}^\mathrm{SL}[\rho_\mathrm B]
$$

となります。
交換相関エネルギーから、A と B の相互作用はまったく出てきません。
離れた電子の相関を表すには、点 $\mathbf r$ と、遠く離れた点 $\mathbf r'$ の情報を同時に使う必要があるためです。
混成汎関数で混ぜる厳密な交換も、相関ではないので、分散力を与えません。

平衡距離の付近では、2 つの原子の電子の雲が少し重なります。
この重なりの領域で、半局所汎関数は何らかの引力や斥力を与えます。
しかしそれは分散力を正しく表したものではなく、汎関数によって大きく違います。
平均場の方法と半局所汎関数が遠距離の電子相関（分散力）を含まないことは、総説でも強調されています [Grimme 2016]。

**ここまでで分かったこと**：半局所汎関数のエネルギーは各点の近くの密度だけで決まるので、重ならない 2 つの原子の間の相関（分散力）を表せません。HF にも分散力はありません。

## 3. 計算例：Ar₂ の相互作用エネルギー

Ar₂ の相互作用エネルギー

$$
\Delta E(R) = E(\mathrm{Ar_2}) - 2E(\mathrm{Ar})
$$

を、PySCF [Sun 2020] で計算します。
比べる方法は次のとおりです。

- HF
- PBE（GGA）、B3LYP（混成）：分散補正なし
- VV10：非局所相関を含む汎関数（第 5 節で扱います）[Vydrov 2010]
- ωB97X-V：VV10 を含む範囲分離混成汎関数（[第 16 回](16-functionals.md)）[Mardirossian 2014]
- CCSD(T)：分散力を含む高精度の方法（[第 10 回](10-coupled-cluster.md)）

基底関数は aug-cc-pVTZ です。
有限の基底関数では、二量体の中の各原子が相手の基底関数を借りてエネルギーを下げ、相互作用が見かけ上強くなります（基底関数重ね合わせ誤差、BSSE）。
これを除くため、原子のエネルギーも、相手の位置に基底関数だけを置いて計算します（counterpoise 補正）。

```python
from pyscf import gto, scf, dft, cc

basis = "aug-cc-pvtz"
cm = 219474.63          # 1 Eh を cm^-1 に

def energy(atom, method):
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    if method == "CCSD(T)":
        mf = scf.RHF(mol).run()
        nfz = 5 * sum(1 for a in mol._atom if a[0] == "Ar")   # Ar の 1s2s2p を凍結
        c = cc.CCSD(mf, frozen=nfz).run()
        return c.e_tot + c.ccsd_t()
    if method == "HF":
        return scf.RHF(mol).run().e_tot
    mf = dft.RKS(mol)
    mf.xc = method
    mf.grids.level = 4
    return mf.run().e_tot

def interaction(R, method):
    """counterpoise 補正した相互作用エネルギー（cm^-1）"""
    dimer = f"Ar 0 0 0; Ar 0 0 {R}"
    mono  = f"Ar 0 0 0; ghost-Ar 0 0 {R}"     # 相手の位置に基底関数だけを置く
    return (energy(dimer, method) - 2 * energy(mono, method)) * cm

Rs = [3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 5.0, 6.0]
methods = ["HF", "PBE", "B3LYP", "VV10", "WB97X-V", "CCSD(T)"]
print("R/Å  " + "".join(f"{m:>9s}" for m in methods))
for R in Rs:
    print(f"{R:4.1f} " + "".join(f"{interaction(R, m):9.1f}" for m in methods))
```

出力は次のとおりです（単位は cm$^{-1}$。負の値が引力）。

```text
R/Å         HF      PBE    B3LYP     VV10  WB97X-V  CCSD(T)
 3.2     705.6    370.2    526.6    157.6    225.5    288.7
 3.4     343.7    114.3    254.1    -38.0     -1.1     50.3
 3.6     165.9      5.1    130.6    -99.6    -82.5    -42.1
 3.8      79.5    -35.0     75.3   -105.5    -99.4    -69.3
 4.0      37.8    -44.3     49.7    -91.6    -90.6    -69.5
 4.4       8.5    -33.7     28.1    -56.9    -57.2    -49.1
 5.0       0.9    -13.2     12.2    -25.7    -24.3    -24.0
 6.0       0.0     -1.6      1.9     -8.1     -7.3     -7.7
```

![Ar₂ の相互作用エネルギー（counterpoise 補正、aug-cc-pVTZ）。赤の十字は、精密な計算による極小（3.762 Å、−99.35 cm⁻¹）。](img/17-ar2.png)

比べるための基準として、CCSD(T) を大きな基底関数の極限まで外挿し、さらに高い励起の効果も加えた精密な計算があります [Patkowski 2010]。
それによると、Ar₂ のポテンシャルの極小は 3.762 Å にあり、深さは 99.35 cm$^{-1}$（0.284 kcal/mol）です（図の赤の十字）。
上の CCSD(T) の極小（約 70 cm$^{-1}$）がこれより浅いのは、aug-cc-pVTZ では相関エネルギーを取りきれないためです。

表と図から、次のことが分かります。

1. **HF と B3LYP**：どの距離でも正（斥力）で、Ar₂ は結合しません。
2. **PBE**：4.0 Å 付近に 44 cm$^{-1}$ の極小があります。しかし、遠くでは急に 0 に近づきます。5.0 Å から 6.0 Å で、CCSD(T) は 24.0 から 7.7 cm$^{-1}$ へ約 3.1 分の 1 になり、$1/R^6$ から期待される $(6/5)^6 \approx 3.0$ 分の 1 とほぼ一致します。PBE は 13.2 から 1.6 cm$^{-1}$ へ 8 分の 1 になり、$1/R^6$ よりずっと速く消えます。PBE の引力は、第 2 節で述べた密度の重なりから来るもので、分散力ではありません。
3. **VV10 と ωB97X-V**：極小は 3.8 Å 付近で、深さは 100 cm$^{-1}$ 前後です。精密な値に近く、遠くでも CCSD(T) と同じように減ります。VV10 の論文でも、Ar₂ の結合をわずかに強く見積もることが報告されています [Vydrov 2010]。

分散補正のない汎関数どうしでも、PBE は弱く結合させ、B3LYP は斥力を与えます。
分散力を含まない汎関数の平衡付近の結果は、汎関数の選び方で大きく変わり、物理的な意味がありません。
分散力を正しく入れるには、汎関数とは別の補正か、非局所の項が必要です。

**ここまでで分かったこと**：Ar₂ では、HF と B3LYP は結合せず、PBE の弱い引力も $1/R^6$ の形ではありません。非局所相関 VV10 を含む汎関数は、精密な計算に近い結合を与えます。

## 4. 原子対の補正：DFT-D2、D3、D4

### 4.1 考え方

いちばん簡単な補正は、汎関数のエネルギーに、原子の対ごとの $-C_6/R^6$ を足すことです。

$$
E_\mathrm{DFT\text{-}D} = E_\mathrm{KS\text{-}DFT} + E_\mathrm{disp}
$$

このとき、2 つの問題があります。

1. 原子が近づくと $1/R^6$ は発散します。
2. 原子が近い距離では、汎関数自身が電子相関の一部をすでに含んでいます。そのまま足すと、相関を 2 回数えることになります。

そこで、近い距離で補正を弱める減衰関数（damping function）をかけます。
減衰関数がどの距離から効くかを決めるパラメータを、汎関数ごとに合わせます。
この方法を DFT-D と呼び、Grimme らが系統的に発展させてきました [Grimme 2004]。

### 4.2 DFT-D2

DFT-D2 は、次の形の補正です [Grimme 2006]。

$$
E_\mathrm{disp}^\mathrm{D2} = -s_6\sum_{A<B}\frac{C_6^{AB}}{R_{AB}^6}\,f_\mathrm{dmp}(R_{AB})
$$

$C_6^{AB}$ は、元素ごとに決めた $C_6$ を組み合わせて作ります。
$f_\mathrm{dmp}$ は、元素ごとの半径の和より近い距離で 0 に近づく関数です。
$s_6$ は汎関数ごとに 1 個だけ決める倍率で、BLYP、PBE、TPSS、B3LYP の値が報告されています。
元素の $C_6$ は、分子の環境によらず一定とされています。

### 4.3 DFT-D3

DFT-D3 は、D2 の考え方を大きく改良したものです [Grimme 2010]。
主な違いは次の 3 点です。

1. **$C_6$ を第一原理計算から求める**：さまざまな化学的環境にある原子について、時間依存 DFT で振動数に依存する分極率を計算し、そこから $C_6^{AB}$ を求めます [Caldeweyher 2017]。
2. **環境への依存**：同じ炭素でも、sp³ の炭素と sp² の炭素では $C_6$ が違います。D3 は、分子の構造から各原子の「配位数」（coordination number、何個の原子とどの程度結合しているかを表す、小数も取る数）を計算し、それに応じて $C_6$ を内挿します。
3. **$C_8$ 項**：$C_8^{AB}$ を、$C_6^{AB}$ から既知の関係式で求めて加えます。

2 体の項は次の形です（和は異なる原子の対について。$\frac12$ は同じ対を 2 回数えないため）。

$$
E^{(2)} = -\frac12\sum_{A\neq B}\sum_{n=6,8}s_n\frac{C_n^{AB}}{R_{AB}^n}\,f_{\mathrm d,n}(R_{AB})
$$

元の D3 では、減衰関数として、近い距離で 0 になる形を使います（ゼロ減衰）。

$$
f_{\mathrm d,n}(R) = \frac{1}{1 + 6\,\big(R/(s_{r,n}R_0^{AB})\big)^{-\alpha_n}}
$$

$R_0^{AB}$ は原子の対ごとの半径、$\alpha_6 = 14$、$\alpha_8 = 16$ です。
汎関数ごとに合わせるパラメータは、$s_8$ と $s_{r,6}$ の 2 個だけです。$s_{r,8} = 1$ に固定し、$s_6$ も 1 に固定します（二重混成汎関数では $s_6 < 1$）。

### 4.4 Becke–Johnson 減衰：D3(BJ)

ゼロ減衰には、近い距離で補正のエネルギーが急に 0 に戻るため、原子の間に見かけの斥力が生じるという欠点があります。
そこで、Becke と Johnson の分散のモデル [Becke 2005] で使われた形の減衰が D3 に取り入れられました [Grimme 2011]。

$$
E^{(2)} = -\frac12\sum_{A\neq B}\sum_{n=6,8}s_n\frac{C_n^{AB}}{R_{AB}^n + \big(a_1R_0^{AB} + a_2\big)^n},
\qquad R_0^{AB} = \sqrt{C_8^{AB}/C_6^{AB}}
$$

$R_{AB} \to 0$ でも分母は 0 にならず、補正は一定の値に近づきます。
汎関数ごとのパラメータは $s_8$、$a_1$、$a_2$ の 3 個です。
Grimme らは、ゼロ減衰と BJ 減衰を 12 個の汎関数で比べ、結果の差は小さいが、BJ 減衰には近距離での見かけの斥力がないことを示しました [Grimme 2011]。
また、HF に補正を加える場合は、BJ 減衰だけが推奨されています。

### 4.5 DFT-D4

DFT-D4 は D3 の後継です [Caldeweyher 2017; Caldeweyher 2019]。
D3 との主な違いは、原子の部分電荷を使うことです。
各原子の分極率を、配位数に加えて部分電荷に応じて変え、そこから $C_6$ をその場で計算します。
分子どうしの $C_6$ の再現では、D3 の平均相対誤差 4.7% に対して 3.8% に改善しています。
特に金属を含む系で、電荷への依存を入れた効果が大きいと報告されています。
60 個以上の汎関数について BJ 型の減衰のパラメータが決められています [Caldeweyher 2019]。

主要族の分子のベンチマーク GMTKN55 でも、熱化学や非共有結合相互作用を扱うには分散補正が必要だと結論されています [Goerigk 2017]。
総説によれば、現在の主な補正は、遠距離の分散エネルギーを典型的に 5% 程度の相対誤差で与えます [Grimme 2016]。

**ここまでで分かったこと**：DFT-D は、原子対ごとの $-C_n/R^n$ に減衰関数をかけて足す補正です。D3 は環境に応じた $C_6$ と $C_8$ を使い、BJ 減衰は近距離で一定の値に近づきます。D4 は部分電荷も使います。

## 5. 非局所汎関数：VV10

もう一つの方法は、交換相関汎関数そのものに、離れた 2 点の密度を結ぶ項を加えることです。
このような汎関数を、非局所相関汎関数と呼びます。
代表的な例に、vdW-DF [Dion 2004] と VV10 [Vydrov 2010] があります。
ここでは、ωB97X-V や ωB97M-V（[第 16 回](16-functionals.md)）にも使われている VV10 を説明します。

VV10 は、非局所相関エネルギーを次の形に書きます。

$$
E_\mathrm c^\mathrm{nl} = \frac12\iint\rho(\mathbf r)\,\Phi(\mathbf r, \mathbf r')\,\rho(\mathbf r')\,d\mathbf r\,d\mathbf r'
$$

$\Phi$ は、2 点の密度と密度の勾配から決まる関数（核、kernel）です。
2 点が遠く離れると、$\Phi$ は次の形になります。

$$
\Phi \to -\frac{3}{2\,\omega_0(\mathbf r)\,\omega_0(\mathbf r')\,\big[\omega_0(\mathbf r) + \omega_0(\mathbf r')\big]\,|\mathbf r - \mathbf r'|^6}
$$

$\omega_0(\mathbf r)$ は、その点の密度と勾配から決まる、局所的な励起のエネルギーの目安です。
この形のおかげで、離れた分子の間に $-C_6/R^6$ の引力が自然に出てきます。
付録 A の Drude 模型の $C_6 \propto \omega_0\alpha^2$ を、空間の各点の小さな振動子について足し合わせたものと考えることができます。

VV10 には、調節するパラメータが 2 個あります。

- $C$：$\omega_0$ の式の中のパラメータで、$C_6$ が正確になるように決めます。元の論文では $C = 0.0093$ です。
- $b$：近い距離で $1/R^6$ を減衰させる強さを決めます。半局所の部分（rPW86 の交換と PBE の相関）と組み合わせ、22 個の分子間相互作用（S22）の結合エネルギーに合わせて $b = 5.9$ と決められました。

VV10 は原子の種類や配位数の表を使わず、密度だけから分散を計算します。
そのかわり、$\mathbf r$ と $\mathbf r'$ の 2 重の積分が必要なので、DFT-D より計算が重くなります。

**ここまでで分かったこと**：VV10 は、離れた 2 点の密度を結ぶ非局所相関の項で、遠くでは $-C_6/R^6$ を与えます。パラメータは $C_6$ と結合エネルギーに合わせた 2 個です。

## 6. 多体の分散：ATM 項と MBD

ここまでの補正は、原子の対ごとの和でした。
しかし、分散のエネルギーは、厳密には対の和では書けません。
3 個以上の原子が同時に関わる寄与（多体の寄与）があるためです。

### 6.1 ATM 項

3 個の原子の間の最初の多体の項は、Axilrod と Teller が導いた 3 体の項です [Axilrod 1943]。
Axilrod–Teller–Muto（ATM）項とも呼ばれます。

$$
E^{(3)} = \sum_{A<B<C} C_9^{ABC}\,\frac{3\cos\theta_A\cos\theta_B\cos\theta_C + 1}{\big(R_{AB}R_{BC}R_{CA}\big)^3}
$$

$\theta_A$ などは、3 個の原子が作る三角形の内角です。
$C_9^{ABC} > 0$ なので、正三角形のような配置（角度の余弦がすべて正）では、この項は正、つまり斥力になります。
一直線の配置では、引力になります。

DFT-D3 では、$C_9^{ABC} \approx \sqrt{C_6^{AB}C_6^{BC}C_6^{CA}}$ と近似して、この項を加えることができます [Grimme 2010]。
グラフェンの層の間の結合では、3 体の項を入れると結合が約 10% 弱くなりました [Grimme 2010]。
D4 も ATM 項を含みます [Caldeweyher 2019]。

### 6.2 MBD

多体の効果をすべての次数で入れる方法の一つが、多体分散（many-body dispersion、MBD）です [Tkatchenko 2012]。
各原子を、分極率を持つ振動子（付録 A の Drude 模型の振動子）で表し、それらを双極子の相互作用でつないだ系のシュレーディンガー方程式を解きます。
連成した振動子の零点エネルギーの変化が、分散のエネルギーになります。
付録 A は 2 個の振動子の場合で、MBD はそれを多数の振動子に広げたものです。

MBD では、ある原子の分極率が、周りの原子の分極によって打ち消されたり強められたりする効果（遮蔽）も入ります。
各原子の分極率と $C_6$ は、分子の電子密度を原子に分け、自由な原子の精密な値を基準にして求めます [Tkatchenko 2009]。
この遮蔽と多体の分散のエネルギーは、小さな分子でも無視できず、生体分子の立体配座のエネルギーや分子結晶の結合では欠かせないことが示されています [Tkatchenko 2012]。

**ここまでで分かったこと**：分散のエネルギーには、原子の対の和では書けない多体の寄与があります。3 体の ATM 項は配置によって斥力にも引力にもなります。MBD は、連成した振動子のモデルで多体の効果をすべて含めます。

## まとめ

- 分散力は、離れた電子のゆらぎの相関から生じる引力で、遠くでは $-C_6/R^6$ の形です。
- HF と半局所汎関数は、重ならない密度の間の相関を表せないので、分散力を原理的に含みません。Ar₂ では、HF と B3LYP は結合せず、PBE の弱い引力も遠くで消えます。
- DFT-D は原子対の $-C_n/R^n$ に減衰関数をかけて足します。D3 は環境に応じた $C_6$ と $C_8$、BJ 減衰は近距離で有限、D4 は部分電荷も使います。
- VV10 は密度だけから分散を計算する非局所相関汎関数で、ωB97X-V や ωB97M-V に含まれています。
- 多体の分散は、ATM 項（3 体）や MBD（すべての次数）で入れます。

## 付録 A　Drude 模型による $C_6$ の導出

原子を、電荷 $q$ と質量 $m$ の粒子が、振動数 $\omega_0$ の等方的なばねで核につながれたものと考えます（Drude 模型）。
静的な分極率は $\alpha = q^2/(m\omega_0^2)$ です。
同じ 2 個の原子を $z$ 軸の上に距離 $R$ だけ離して置き、それぞれの粒子の変位を $\mathbf u_1$、$\mathbf u_2$ とします。
$R$ が変位より十分大きいとき、2 つの双極子 $q\mathbf u_1$、$q\mathbf u_2$ の相互作用は

$$
V = \frac{q^2}{R^3}\big(x_1x_2 + y_1y_2 - 2z_1z_2\big)
$$

です。
ハミルトニアンは $x$、$y$、$z$ の方向ごとに分かれ、各方向は結合した 2 個の調和振動子になります。
結合の強さを $V = \kappa\,u_1u_2$ と書くと、$z$ 方向では $\kappa = -2q^2/R^3$、$x$ と $y$ 方向では $\kappa = q^2/R^3$ です。
変数を $u_\pm = (u_1 \pm u_2)/\sqrt2$ に変えると、振動数が

$$
\omega_\pm = \omega_0\sqrt{1 \pm \frac{\kappa}{m\omega_0^2}}
$$

の 2 個の独立な振動子になります。
基底状態のエネルギーは零点エネルギーの和 $\frac12(\omega_+ + \omega_-)$ です（$\hbar = 1$）。
$x = \kappa/(m\omega_0^2)$ が小さいとき、$\sqrt{1+x} + \sqrt{1-x} \approx 2 - x^2/4$ なので、相互作用がないときからのエネルギーの変化は

$$
\Delta E = -\frac{\omega_0}{8}\,x^2
$$

です。
$x^2 = \kappa^2/(m\omega_0^2)^2 = \alpha^2\kappa^2/q^4$ を使い、3 方向を足すと、

$$
\Delta E = -\frac{\omega_0}{8}\cdot\frac{\alpha^2}{R^6}\big(4 + 1 + 1\big) = -\frac34\,\frac{\omega_0\,\alpha^2}{R^6}
$$

となります。
これが本文の $C_6 = \frac34\omega_0\alpha^2$ です。
古典的には、ばねの粒子は静止していて双極子はありません。
引力は、零点振動（量子的なゆらぎ）が相関することから生じています。

## 参考文献

- F. London, "The general theory of molecular forces", Trans. Faraday Soc. 33, 8 (1937). DOI: [10.1039/TF937330008B](https://doi.org/10.1039/TF937330008B)
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, "A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu", J. Chem. Phys. 132, 154104 (2010). DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- S. Grimme, A. Hansen, J. G. Brandenburg, C. Bannwarth, "Dispersion-Corrected Mean-Field Electronic Structure Methods", Chem. Rev. 116, 5105 (2016). DOI: [10.1021/acs.chemrev.5b00533](https://doi.org/10.1021/acs.chemrev.5b00533)
- K. Patkowski, K. Szalewicz, "Argon pair potential at basis set and excitation limits", J. Chem. Phys. 133, 094304 (2010). DOI: [10.1063/1.3478513](https://doi.org/10.1063/1.3478513)
- S. Grimme, "Accurate description of van der Waals complexes by density functional theory including empirical corrections", J. Comput. Chem. 25, 1463 (2004). DOI: [10.1002/jcc.20078](https://doi.org/10.1002/jcc.20078)
- S. Grimme, "Semiempirical GGA-type density functional constructed with a long-range dispersion correction", J. Comput. Chem. 27, 1787 (2006). DOI: [10.1002/jcc.20495](https://doi.org/10.1002/jcc.20495)
- A. D. Becke, E. R. Johnson, "A density-functional model of the dispersion interaction", J. Chem. Phys. 123, 154101 (2005). DOI: [10.1063/1.2065267](https://doi.org/10.1063/1.2065267)
- S. Grimme, S. Ehrlich, L. Goerigk, "Effect of the damping function in dispersion corrected density functional theory", J. Comput. Chem. 32, 1456 (2011). DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- E. Caldeweyher, C. Bannwarth, S. Grimme, "Extension of the D3 dispersion coefficient model", J. Chem. Phys. 147, 034112 (2017). DOI: [10.1063/1.4993215](https://doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth, S. Grimme, "A generally applicable atomic-charge dependent London dispersion correction", J. Chem. Phys. 150, 154122 (2019). DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
- M. Dion, H. Rydberg, E. Schröder, D. C. Langreth, B. I. Lundqvist, "Van der Waals Density Functional for General Geometries", Phys. Rev. Lett. 92, 246401 (2004). DOI: [10.1103/PhysRevLett.92.246401](https://doi.org/10.1103/PhysRevLett.92.246401)
- O. A. Vydrov, T. Van Voorhis, "Nonlocal van der Waals density functional: The simpler the better", J. Chem. Phys. 133, 244103 (2010). DOI: [10.1063/1.3521275](https://doi.org/10.1063/1.3521275)
- N. Mardirossian, M. Head-Gordon, "ωB97X-V: A 10-parameter, range-separated hybrid, generalized gradient approximation density functional with nonlocal correlation, designed by a survival-of-the-fittest strategy", Phys. Chem. Chem. Phys. 16, 9904 (2014). DOI: [10.1039/C3CP54374A](https://doi.org/10.1039/C3CP54374A)
- B. M. Axilrod, E. Teller, "Interaction of the van der Waals Type Between Three Atoms", J. Chem. Phys. 11, 299 (1943). DOI: [10.1063/1.1723844](https://doi.org/10.1063/1.1723844)
- A. Tkatchenko, M. Scheffler, "Accurate Molecular Van Der Waals Interactions from Ground-State Electron Density and Free-Atom Reference Data", Phys. Rev. Lett. 102, 073005 (2009). DOI: [10.1103/PhysRevLett.102.073005](https://doi.org/10.1103/PhysRevLett.102.073005)
- A. Tkatchenko, R. A. DiStasio Jr., R. Car, M. Scheffler, "Accurate and Efficient Method for Many-Body van der Waals Interactions", Phys. Rev. Lett. 108, 236402 (2012). DOI: [10.1103/PhysRevLett.108.236402](https://doi.org/10.1103/PhysRevLett.108.236402)
- L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, "A look at the density functional theory zoo with the advanced GMTKN55 database for general main group thermochemistry, kinetics and noncovalent interactions", Phys. Chem. Chem. Phys. 19, 32184 (2017). DOI: [10.1039/C7CP04913G](https://doi.org/10.1039/C7CP04913G)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- DFT-D3 プログラムの説明書（Grimme のグループ。D3 の 2 体・3 体の項の式）：<https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3>
- simple-dftd3（D3 の再実装）の説明書：<https://dftd3.readthedocs.io/>
- dftd4 の説明書：<https://dftd4.readthedocs.io/>
