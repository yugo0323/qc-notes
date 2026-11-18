# 第 16 回　汎関数の階段：LDA、GGA、meta-GGA、混成、範囲分離

!!! abstract "この回の道筋"
    Kohn–Sham 法で近似するのは、交換相関エネルギー $E_\mathrm{xc}$ だけです（[第 15 回](15-kohn-sham.md)）。この回では、$E_\mathrm{xc}$ の近似（汎関数）を、使う材料の多さで 5 段に並べた「ヤコブの梯子」に沿って紹介します。各段の代表的な汎関数が、理論の条件だけから作られたのか、実験値などに合わせて作られたのかも区別します。最後に、H₂O の原子化エネルギーの計算と、大きなベンチマークで分かっている傾向を比べます。

    | 節 | 内容 |
    |---|---|
    | 1 | ヤコブの梯子 |
    | 2 | 第 1 段：LDA |
    | 3 | 第 2 段：GGA |
    | 4 | 第 3 段：meta-GGA |
    | 5 | 第 4 段：混成汎関数と B3LYP |
    | 6 | 範囲分離混成汎関数 |
    | 7 | 第 5 段：二重混成汎関数 |
    | 8 | 計算例：H₂O の原子化エネルギー |
    | 9 | ベンチマークで分かっている傾向 |

    **前に読んでおく回**：量子化学の理論 [第 15 回](15-kohn-sham.md)

!!! tip "物理的なイメージ"
    汎関数が近似しているのは、交換相関エネルギー、つまり電子が互いを避けることでエネルギーがどれだけ下がるかです。ヤコブの梯子は、その見積もりに使う材料の多さの順です。各点の密度だけを見る LDA、まわりの密度の傾きも見る GGA、運動エネルギー密度も見る meta-GGA は、その点の近くの情報だけで決めます。混成汎関数は、軌道を通して離れた点の情報（厳密な交換）を混ぜ、二重混成汎関数は空軌道も使います。上の段ほど材料が多く、計算は重くなります。

この回では原子単位を使います。エネルギーの換算は 1 $E_\mathrm h$ = 27.211 eV = 627.5 kcal/mol です。

## 1. ヤコブの梯子

厳密な $E_\mathrm{xc}[\rho]$ の形は分かりません。
そこで、$E_\mathrm{xc}$ を、各点での簡単な量の関数として近似します。

$$
E_\mathrm{xc} \approx \int e_\mathrm{xc}\big(\rho(\mathbf r), \nabla\rho(\mathbf r), \dots\big)\,d\mathbf r
$$

$e_\mathrm{xc}$ は交換相関エネルギーの密度です。
材料（$e_\mathrm{xc}$ の変数）を増やすほど、近似の自由度が増え、計算は重くなります。

Perdew と Schmidt は、材料の多さで近似を段に分け、聖書の物語になぞらえて「ヤコブの梯子」と呼びました [Perdew 2001]。
梯子の上には、化学で必要な精度（目安として 1 kcal/mol）の「天国」があります。
段と材料は次のとおりです。

| 段 | 名前 | 加わる材料 |
|---|---|---|
| 1 | LDA（局所密度近似） | 密度 $\rho$ |
| 2 | GGA（一般化勾配近似） | 密度の勾配 $\nabla\rho$ |
| 3 | meta-GGA | 運動エネルギー密度 $\tau$ |
| 4 | 混成（hybrid） | 占有軌道（厳密な交換） |
| 5 | 二重混成（double hybrid）など | 仮想軌道 |

第 1〜3 段は、各点の近くの情報だけを使うので、半局所（semilocal）汎関数と呼ばれます。
第 4 段からは、軌道を通して空間の離れた点の情報が入ります。

梯子を上る方法は、大きく 2 通りあります。

- **非経験的（nonempirical）**：厳密な $E_\mathrm{xc}$ が満たすと分かっている条件（拘束条件、exact constraint）を満たすように式を決めます。
- **経験的（empirical）**：式の中にパラメータを置き、分子の実験値や精密な計算値に合うように決めます。

以下では、各汎関数がどちらで作られたかも書きます。

**ここまでで分かったこと**：汎関数は、使う材料（$\rho$、$\nabla\rho$、$\tau$、占有軌道、仮想軌道）で 5 段に分けられます。作り方には、理論の条件で決める方法と、データに合わせる方法があります。

## 2. 第 1 段：LDA

局所密度近似（LDA）は、各点の $e_\mathrm{xc}$ を、同じ密度の一様な電子ガスの値で置き換えます（[第 15 回](15-kohn-sham.md)）。
一様な電子ガスは、正の電荷が一様に広がった中を電子が動く模型です。

交換の部分は、一様な電子ガスで厳密に計算でき、[第 15 回](15-kohn-sham.md)の式になります。

$$
E_\mathrm x^\mathrm{LDA}[\rho] = -\frac34\Big(\frac3\pi\Big)^{1/3}\int\rho^{4/3}\,d\mathbf r
$$

この式を交換の近似として使う考え方は Slater にさかのぼるので、Slater の交換とも呼ばれます [Slater 1951]。

相関の部分は、一様な電子ガスでも式では解けません。
Ceperley と Alder は、量子モンテカルロ法でその相関エネルギーを精密に計算しました [Ceperley 1980]。
Vosko、Wilk、Nusair（VWN）は、その結果に合う式を作りました [Vosko 1980]。
Slater の交換と VWN の相関を組み合わせた LDA を、SVWN と呼びます。

LDA は、一様な電子ガスという極限では正しい近似です。
その意味で、非経験的な汎関数です。
一方、分子の密度は一様からほど遠く、LDA は結合エネルギーを大きく見積もりすぎます。
第 8 節の計算でも、H₂O の原子化エネルギーを 30 kcal/mol 以上大きく見積もります。

**ここまでで分かったこと**：LDA は一様な電子ガスの交換と相関を使います。相関は量子モンテカルロの結果に合わせた式（VWN など）です。分子の結合エネルギーは大きく見積もりすぎます。

## 3. 第 2 段：GGA

密度の変化の速さを入れるため、勾配 $\nabla\rho$ を材料に加えます。
勾配は、次の無次元の量にまとめて使います。

$$
s(\mathbf r) = \frac{|\nabla\rho|}{2(3\pi^2)^{1/3}\rho^{4/3}}
$$

$s$ は、その点の近くで密度がどれだけ急に変わるかを表します。一様な電子ガスでは $s = 0$ です。
交換の GGA は、LDA の交換に、$s$ の関数 $F_\mathrm x(s)$ をかけた形に書けます。

$$
E_\mathrm x^\mathrm{GGA} = \int e_\mathrm x^\mathrm{LDA}(\rho)\,F_\mathrm x(s)\,d\mathbf r,
\qquad F_\mathrm x(0) = 1
$$

$F_\mathrm x$ を増強因子（enhancement factor）と呼びます。
GGA の違いは、主に $F_\mathrm x$ と相関の式の選び方にあります。

代表的な GGA は次のとおりです。

- **B88（交換）**：Becke の交換汎関数です [Becke 1988]。遠くでの交換エネルギー密度の正しい振る舞いを満たすように式を作り、残る 1 個のパラメータを、原子の厳密な交換エネルギーに合わせて決めました。
- **LYP（相関）**：Lee、Yang、Parr の相関汎関数です [Lee 1988]。He 原子の相関エネルギーに合わせて作られた Colle–Salvetti の式を、密度の汎関数に書き直したものです。
- **BLYP**：B88 と LYP の組み合わせです。部分的に経験的な GGA です。
- **PBE**：Perdew、Burke、Ernzerhof の GGA です [Perdew 1996]。LDA の部分を除き、すべてのパラメータが基礎的な定数で、厳密な条件から決められています。非経験的な GGA の代表です。

**ここまでで分かったこと**：GGA は、密度の勾配から作った $s$ を使って LDA を補正します。B88 と LYP は原子のデータを一部使い、PBE は条件だけで作られています。

## 4. 第 3 段：meta-GGA

meta-GGA は、Kohn–Sham 軌道から作る運動エネルギー密度

$$
\tau(\mathbf r) = \frac12\sum_{i}^{\mathrm{occ}}|\nabla\phi_i(\mathbf r)|^2
$$

を材料に加えます（和は占有されたスピン軌道について）。
$\tau$ は軌道を使いますが、各点の値だけで決まるので、計算の手間は GGA と大きく変わりません。

$\tau$ を加える利点の一つは、「その点の電子が、1 個の軌道だけで表される領域か」を見分けられることです。
閉殻で 1 個の空間軌道 $\phi$ だけが占有されているとき、$\rho = 2\phi^2$ なので、

$$
\tau = |\nabla\phi|^2 = \frac{|\nabla\rho|^2}{8\rho} \equiv \tau_\mathrm W
$$

です（$\phi$ は実数としました）。
$\tau_\mathrm W$ は、$\rho$ と $\nabla\rho$ だけで書けます。
そこで $\tau$ と $\tau_\mathrm W$ を比べると、たとえば H 原子や H₂ の結合の付近のような 1 軌道の領域を見分けられます。
こうした領域では、自己相互作用（[第 15 回](15-kohn-sham.md)）を打ち消す条件を、汎関数に課すことができます。

代表的な meta-GGA は、どれも非経験的です。

- **TPSS**：Tao、Perdew、Staroverov、Scuseria の meta-GGA です [Tao 2003]。経験的なパラメータを使わずに、厳密な条件を満たすように作られています。
- **SCAN**：meta-GGA が満たしうる既知の 17 個の厳密な条件を、すべて満たすように作られました [Sun 2015]。さらに、希ガス原子など、「適切な基準系（appropriate norms）」で厳密かほぼ厳密になるようにしてあります。
- **r²SCAN**：SCAN は、数値的に扱いにくい（数値積分の精度を上げないと結果が安定しにくい）汎関数でした。r²SCAN は、数値的に扱いやすくした rSCAN を、厳密な条件を満たすように直したものです [Furness 2020]。

**ここまでで分かったこと**：meta-GGA は、運動エネルギー密度 $\tau$ を加えます。$\tau$ を使うと、1 個の軌道で表される領域を見分けられます。TPSS、SCAN、r²SCAN は非経験的に作られています。

## 5. 第 4 段：混成汎関数と B3LYP

### 5.1 厳密な交換を混ぜる理由

第 4 段では、Kohn–Sham 軌道で計算した Hartree–Fock の交換エネルギー $E_\mathrm x^\mathrm{HF}$ を混ぜます。
$E_\mathrm x^\mathrm{HF}$ を、ここでは厳密な交換と呼びます（[第 15 回](15-kohn-sham.md)で定義した交換エネルギー $E_\mathrm x$ そのものです）。

混ぜる理由は、電子どうしの相互作用の強さを $\lambda$ 倍した系を考えると分かります。
$\lambda$ を 0 から 1 まで変えるとき、密度は本物の $\rho$ のまま保つとします。
$\lambda = 0$ は Kohn–Sham の参照系で、$E_\mathrm{xc}$ は厳密な交換だけになります。
$\lambda = 1$ は本物の系です。
$E_\mathrm{xc}$ は、この間で $\lambda$ について積分した形に書けます（断熱接続、adiabatic connection）。
半局所汎関数は、$\lambda = 0$ の端で誤差が大きいので、その端を厳密な交換で置き換えると改善が見込めます。

Perdew、Ernzerhof、Burke は、この考え方から、混ぜる割合をおよそ $1/4$ とする根拠を示しました [Perdew 1996b]。

$$
E_\mathrm{xc}^\mathrm{hyb} = E_\mathrm{xc}^\mathrm{GGA} + \frac14\big(E_\mathrm x^\mathrm{HF} - E_\mathrm x^\mathrm{GGA}\big)
$$

GGA に PBE を使ったものが **PBE0** です [Adamo 1999]。
混ぜる割合も理論から決めたので、PBE0 は経験的なパラメータを持たない混成汎関数です。

混成汎関数の計算では、$E_\mathrm x^\mathrm{HF}$ の微分から、Hartree–Fock と同じ交換演算子 $\hat K$ が出ます。
Kohn–Sham 方程式の $v_\mathrm{xc}$ の一部が、非局所的な $\hat K$ に置き換わります。

### 5.2 B3PW91 と B3LYP

Becke は 1993 年に、3 個のパラメータを持つ混成汎関数を提案しました [Becke 1993]。

$$
E_\mathrm{xc} = E_\mathrm{xc}^\mathrm{LDA} + a_0\big(E_\mathrm x^\mathrm{HF} - E_\mathrm x^\mathrm{LDA}\big) + a_\mathrm x\,\Delta E_\mathrm x^\mathrm{B88} + a_\mathrm c\,\Delta E_\mathrm c^\mathrm{PW91}
$$

$\Delta E^\mathrm{B88}_\mathrm x$ と $\Delta E^\mathrm{PW91}_\mathrm c$ は、GGA の交換と相関の、LDA への補正の部分です。
3 個のパラメータは、56 個の原子化エネルギーなどの実験値に合わせて決められ、$a_0 = 0.20$、$a_\mathrm x = 0.72$、$a_\mathrm c = 0.81$ になりました。
この汎関数は B3PW91 と呼ばれます。

B3LYP は、この式の相関の部分を LYP に置き換えたものです [Stephens 1994]。
パラメータは決め直さずに、そのまま使われました [Casida 2026]。
LYP は、LDA の部分と勾配の補正の部分に分かれていません。
そこで B3LYP は、相関の LDA の部分の 19% を VWN で補い、次の形になっています。

$$
E_\mathrm{xc}^\mathrm{B3LYP} = 0.08\,E_\mathrm x^\mathrm{Slater} + 0.20\,E_\mathrm x^\mathrm{HF} + 0.72\,E_\mathrm x^\mathrm{B88} + 0.19\,E_\mathrm c^\mathrm{VWN} + 0.81\,E_\mathrm c^\mathrm{LYP}
$$

### 5.3 B3LYP の中の VWN の版

ここに、よく知られた混乱があります。
VWN の論文には、相関エネルギーの式が何通りか載っています [Vosko 1980]。
主に使われるのは、乱雑位相近似（RPA）の計算結果に合わせた式と、Ceperley–Alder の量子モンテカルロの結果に合わせた式（VWN5 と呼ばれます）です。
LDA として正しいのは後者です。
ところが Gaussian の B3LYP の VWN は前者で、Stephens らの論文はどちらを使ったかを明記していませんでした [Casida 2026]。
このため、プログラムによって B3LYP の中身が違うことになりました [Hertwig 1997]。

PySCF [Sun 2020] では、`B3LYP` は Gaussian と同じ RPA の版で、VWN5 を使う版は `B3LYP5` という名前です。
H₂O（構造は[第 15 回](15-kohn-sham.md)と同じ、aug-cc-pVTZ）で確かめます。
`geom` は[第 15 回](15-kohn-sham.md)のコードと同じです。

```python
from pyscf import gto, dft
mol = gto.M(atom=geom, basis="aug-cc-pvtz", verbose=0)
for xc in ["B3LYP",                                     # PySCF の既定
           ".2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN_RPA",
           "B3LYP5",
           ".2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN5"]:
    e = dft.RKS(mol, xc=xc).run().e_tot
    print(f"{e:.6f}  {xc}")
```

```text
-76.466157  B3LYP
-76.466157  .2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN_RPA
-76.429101  B3LYP5
-76.429101  .2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN5
```

2 つの版の全エネルギーは 0.037 $E_\mathrm h$（約 23 kcal/mol）も違います。
エネルギーの差（反応エネルギーなど）では、この違いの大部分が打ち消し合います。
第 8 節の原子化エネルギーでは、差は 0.6 kcal/mol でした。
それでも、別のプログラムの結果と比べるときは、どちらの版かを必ず確かめてください。

**ここまでで分かったこと**：混成汎関数は厳密な交換を混ぜます。PBE0 は割合 1/4 を理論から決め、B3LYP は 3 個のパラメータを実験値に合わせた B3PW91 から作られました。B3LYP の VWN には 2 つの版があり、プログラムによって違います。

## 6. 範囲分離混成汎関数

混成汎関数では、厳密な交換の割合は、電子間の距離によらず一定でした。
範囲分離（range separation）では、電子間のクーロン相互作用を、距離 $r_{12}$ で 2 つに分けます。

$$
\frac{1}{r_{12}} = \underbrace{\frac{\operatorname{erfc}(\omega r_{12})}{r_{12}}}_{\text{短距離}} + \underbrace{\frac{\operatorname{erf}(\omega r_{12})}{r_{12}}}_{\text{長距離}}
$$

$\operatorname{erf}$ は誤差関数、$\operatorname{erfc} = 1 - \operatorname{erf}$ です。
$\omega$（単位は bohr$^{-1}$）が、短距離と長距離の境目のおよその逆数を決めます。
そのうえで、交換の短距離の部分と長距離の部分で、厳密な交換の割合を変えます。

長距離の部分を 100% 厳密な交換にするものを、長距離補正（long-range corrected、LC）汎関数と呼びます。
Iikura、Tsuneda、Yanai、Hirao は、GGA の交換にこの補正を加える方法を示しました [Iikura 2001]。
半局所汎関数の交換には長距離の相互作用が欠けており、それが共役ポリエンの分極率の過大評価などの原因になることが、この論文で確かめられています。

代表的な範囲分離混成汎関数は、どれも経験的に作られています。

- **ωB97X**：Chai と Head-Gordon の汎関数です [Chai 2008]。B97 型の柔軟な式のパラメータを、熱化学などのデータに合わせました。短距離にも厳密な交換を少し混ぜます（短距離に厳密な交換を混ぜない版は ωB97 と呼ばれます）。対称なラジカル陽イオンの解離や、離れた分子間の電荷移動励起での、従来の混成汎関数の失敗を大きく減らします。
- **ωB97X-V**：パラメータは 10 個で、非局所相関（VV10、[第 17 回](17-dispersion.md)で扱います）を含みます [Mardirossian 2014]。式の項の組み合わせ 16,383 通りを訓練データに合わせ、試験データで最もよいものを選びました。
- **ωB97M-V**：meta-GGA を土台にした範囲分離混成汎関数で、VV10 を含みます [Mardirossian 2016]。約 $10\times10^9$ 通りの候補から、870 個の訓練データで 12 個のパラメータの式を選びました。

**ここまでで分かったこと**：範囲分離混成汎関数は、クーロン相互作用を距離で分け、短距離と長距離で厳密な交換の割合を変えます。ωB97 系は、多数のデータに合わせた経験的な汎関数です。

## 7. 第 5 段：二重混成汎関数

第 5 段では、仮想軌道も使います。
代表は、二重混成（double hybrid）汎関数です。
Grimme の B2PLYP は、混成汎関数の相関の一部を、MP2（[第 9 回](09-mp2.md)）と同じ形の 2 次摂動の相関に置き換えます [Grimme 2006]。

$$
E_\mathrm{xc}^\mathrm{B2PLYP} = (1 - a_\mathrm x)E_\mathrm x^\mathrm{B88} + a_\mathrm x E_\mathrm x^\mathrm{HF} + (1 - c)E_\mathrm c^\mathrm{LYP} + c\,E_\mathrm c^\mathrm{PT2}
$$

$E_\mathrm c^\mathrm{PT2}$ は、Kohn–Sham の軌道と軌道エネルギーで計算した MP2 の相関エネルギーです。
2 個のパラメータは生成熱のデータ（G2/97）に合わせて決められ、$a_\mathrm x = 0.53$、$c = 0.27$ になりました。
厳密な交換の割合が 53% と大きいことが特徴です。

計算は 2 段階です。
まず、$E_\mathrm c^\mathrm{PT2}$ を除いた部分で Kohn–Sham 方程式を解き、その軌道で $E_\mathrm c^\mathrm{PT2}$ を計算して足します。
計算の手間は MP2 と同じ程度に増えます。

**ここまでで分かったこと**：二重混成汎関数は、仮想軌道を使う 2 次摂動の相関を混ぜます。B2PLYP は、厳密な交換 53% と PT2 の相関 27% を、データに合わせて決めています。

## 8. 計算例：H₂O の原子化エネルギー

梯子の各段の汎関数で、H₂O の原子化エネルギー

$$
D_\mathrm e = E(\mathrm O) + 2E(\mathrm H) - E(\mathrm{H_2O})
$$

を計算します。
O 原子は三重項、H 原子は二重項として、非制限の方法（UKS）で解きます。
基準には、CCSD(T)（[第 10 回](10-coupled-cluster.md)で扱う高精度の方法）の値を使います。
CCSD(T) は基底関数への依存が大きいので、aug-cc-pVTZ と aug-cc-pVQZ の相関エネルギーから、基底関数が無限に大きい極限へ外挿しました（$X^{-3}$ 則 [Helgaker 1997]、O の 1s は相関させない）。
得られた基準値は $D_\mathrm e$ = 233.1 kcal/mol です。
DFT は基底関数への依存が小さいので、aug-cc-pVQZ だけで計算します。

```python
import numpy as np
from pyscf import gto, dft, mp

basis = "aug-cc-pvqz"
h2o = gto.M(atom=geom, basis=basis, verbose=0)     # geom は[第 15 回](15-kohn-sham.md)と同じ
o   = gto.M(atom="O 0 0 0", basis=basis, spin=2, verbose=0)
h   = gto.M(atom="H 0 0 0", basis=basis, spin=1, verbose=0)
kcal = 627.5095

def energy(mol, xc):
    """xc の汎関数で全エネルギーを返す。B2PLYP は PT2 の項を足す。"""
    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.grids.level = 4
    if xc == "B2PLYP":
        mf.xc = "0.53*HF + 0.47*B88, 0.73*LYP"
        mf.run()
        pt = mp.MP2(mf).run() if mol.spin == 0 else mp.UMP2(mf).run()
        return mf.e_tot + 0.27 * pt.e_corr
    mf.xc = xc
    if xc.upper().endswith("-V"):
        mf.nlcgrids.level = 1      # VV10 の積分の格子
    return mf.run().e_tot

def atomization(xc):
    return (energy(o, xc) + 2*energy(h, xc) - energy(h2o, xc)) * kcal

ref = 233.1   # CCSD(T)/CBS（本文を参照）
for xc in ["HF", "SVWN", "BLYP", "PBE", "TPSS", "SCAN", "R2SCAN",
           "B3LYP", "PBE0", "WB97X-V", "WB97M-V", "B2PLYP"]:
    d = atomization(xc)
    print(f"{xc:8s} {d:7.1f} {d-ref:+6.1f}")
```

```text
HF         155.9  -77.2
SVWN       266.5  +33.4
BLYP       232.9   -0.2
PBE        234.4   +1.3
TPSS       228.6   -4.5
SCAN       228.9   -4.2
R2SCAN     230.0   -3.1
B3LYP      231.2   -1.9
PBE0       227.2   -5.9
WB97X-V    230.3   -2.8
WB97M-V    232.1   -1.0
B2PLYP     231.5   -1.6
```

右の列が基準からの誤差（kcal/mol）です。
CCSD(T) の基準値は、この例のための参考値です。相関させていない O の 1s などの寄与が、まだ 1 kcal/mol ほどの不確かさとして残ります。

表から、次のことが読み取れます。

1. HF は、相関がないので結合を大きく弱く見積もります（第 7 回）。
2. LDA（SVWN）は、逆に 33 kcal/mol も強く見積もります。
3. GGA より上の段では、誤差はどれも 6 kcal/mol 以内に収まります。

ただし、段を上るほど誤差が小さくなる、という単純な順序にはなっていません。
この分子 1 個では、BLYP の誤差がいちばん小さくなりました。
分子 1 個、性質 1 種類の比較では、誤差の打ち消し合いの偶然が結果を左右します。
汎関数の良し悪しは、次の節のような大きなベンチマークで判断します。

**ここまでで分かったこと**：H₂O の原子化エネルギーでは、LDA は大きく強く結合を見積もり、GGA より上の段ではどれも数 kcal/mol の誤差です。1 個の分子だけでは、汎関数の順位は決められません。

## 9. ベンチマークで分かっている傾向

GMTKN55 は、主要族元素の熱化学、反応速度（障壁）、非共有結合相互作用を集めたベンチマークです [Goerigk 2017]。
55 個のデータセットに、1505 個の相対エネルギーが含まれます。
Goerigk らは、分散補正（[第 17 回](17-dispersion.md)）の有無を含む 217 通りの汎関数を評価し、83 個を詳しく調べました。
その論文の要旨は、次の傾向を述べています。

1. 熱化学や非共有結合相互作用を扱うときは、Minnesota 系の汎関数も含め、分散補正が必要です。
2. 最も信頼できるのは二重混成汎関数で、使えるなら使うべきです。特に DSD-BLYP-D3(BJ)、DSD-PBEP86-D3(BJ)、B2GPPLYP-D3(BJ) が挙げられています。
3. 混成汎関数では、ωB97X-V、M052X-D3(0)、ωB97X-D3 が最もよく、範囲分離しない混成汎関数では PW6B95-D3(BJ) が推奨されています。
4. meta-GGA では SCAN-D3(BJ) が推奨されます。それ以外の meta-GGA は、GGA の revPBE-D3(BJ)、B97-D3(BJ)、OLYP-D3(BJ) に負けています。
5. B3LYP のように広く使われている汎関数の多くは、推奨に入っていません。

つまり、平均としては梯子を上るほどよくなります。
ただし、同じ段の中でも汎関数によって差が大きく、段だけでは決まりません。
また、どの段でも分散補正を加えることが前提です。
「よく使われているから」という理由で汎関数を選ばず、扱う問題に近いベンチマークの結果を見て選ぶことが大切です。

**ここまでで分かったこと**：GMTKN55 では、平均として二重混成、混成、半局所の順によい結果が得られます。分散補正はほぼ必須で、B3LYP は推奨に入っていません。

## まとめ

- 汎関数は、使う材料で 5 段に分けられます（ヤコブの梯子）：LDA（$\rho$）、GGA（$\nabla\rho$）、meta-GGA（$\tau$）、混成（占有軌道）、二重混成（仮想軌道）。
- 非経験的な汎関数（PBE、TPSS、SCAN、r²SCAN、PBE0）は厳密な条件から、経験的な汎関数（B3LYP、ωB97X-V、ωB97M-V、B2PLYP）はデータに合わせて作られています。
- B3LYP の VWN には RPA の版と VWN5 の版があり、プログラムによって違うので、全エネルギーを比べるときは注意が必要です。
- 範囲分離混成汎関数は、長距離の交換を厳密な交換にして、半局所汎関数の長距離の欠点を補います。
- GMTKN55 では、二重混成と一部の範囲分離混成が最もよく、分散補正はほぼ必須です。

## 参考文献

- J. P. Perdew, K. Schmidt, "Jacob's ladder of density functional approximations for the exchange-correlation energy", AIP Conf. Proc. 577, 1 (2001). DOI: [10.1063/1.1390175](https://doi.org/10.1063/1.1390175)
- J. C. Slater, "A Simplification of the Hartree-Fock Method", Phys. Rev. 81, 385 (1951). DOI: [10.1103/PhysRev.81.385](https://doi.org/10.1103/PhysRev.81.385)
- D. M. Ceperley, B. J. Alder, "Ground State of the Electron Gas by a Stochastic Method", Phys. Rev. Lett. 45, 566 (1980). DOI: [10.1103/PhysRevLett.45.566](https://doi.org/10.1103/PhysRevLett.45.566)
- S. H. Vosko, L. Wilk, M. Nusair, "Accurate spin-dependent electron liquid correlation energies for local spin density calculations: a critical analysis", Can. J. Phys. 58, 1200 (1980). DOI: [10.1139/p80-159](https://doi.org/10.1139/p80-159)
- A. D. Becke, "Density-functional exchange-energy approximation with correct asymptotic behavior", Phys. Rev. A 38, 3098 (1988). DOI: [10.1103/PhysRevA.38.3098](https://doi.org/10.1103/PhysRevA.38.3098)
- C. Lee, W. Yang, R. G. Parr, "Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density", Phys. Rev. B 37, 785 (1988). DOI: [10.1103/PhysRevB.37.785](https://doi.org/10.1103/PhysRevB.37.785)
- J. P. Perdew, K. Burke, M. Ernzerhof, "Generalized Gradient Approximation Made Simple", Phys. Rev. Lett. 77, 3865 (1996). DOI: [10.1103/PhysRevLett.77.3865](https://doi.org/10.1103/PhysRevLett.77.3865)
- J. Tao, J. P. Perdew, V. N. Staroverov, G. E. Scuseria, "Climbing the Density Functional Ladder: Nonempirical Meta–Generalized Gradient Approximation Designed for Molecules and Solids", Phys. Rev. Lett. 91, 146401 (2003). DOI: [10.1103/PhysRevLett.91.146401](https://doi.org/10.1103/PhysRevLett.91.146401)
- J. Sun, A. Ruzsinszky, J. P. Perdew, "Strongly Constrained and Appropriately Normed Semilocal Density Functional", Phys. Rev. Lett. 115, 036402 (2015). DOI: [10.1103/PhysRevLett.115.036402](https://doi.org/10.1103/PhysRevLett.115.036402)
- J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, J. Sun, "Accurate and Numerically Efficient r2SCAN Meta-Generalized Gradient Approximation", J. Phys. Chem. Lett. 11, 8208 (2020). DOI: [10.1021/acs.jpclett.0c02405](https://doi.org/10.1021/acs.jpclett.0c02405)
- J. P. Perdew, M. Ernzerhof, K. Burke, "Rationale for mixing exact exchange with density functional approximations", J. Chem. Phys. 105, 9982 (1996). DOI: [10.1063/1.472933](https://doi.org/10.1063/1.472933)
- C. Adamo, V. Barone, "Toward reliable density functional methods without adjustable parameters: The PBE0 model", J. Chem. Phys. 110, 6158 (1999). DOI: [10.1063/1.478522](https://doi.org/10.1063/1.478522)
- A. D. Becke, "Density-functional thermochemistry. III. The role of exact exchange", J. Chem. Phys. 98, 5648 (1993). DOI: [10.1063/1.464913](https://doi.org/10.1063/1.464913)
- P. J. Stephens, F. J. Devlin, C. F. Chabalowski, M. J. Frisch, "Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra Using Density Functional Force Fields", J. Phys. Chem. 98, 11623 (1994). DOI: [10.1021/j100096a001](https://doi.org/10.1021/j100096a001)
- R. H. Hertwig, W. Koch, "On the parameterization of the local correlation functional. What is Becke-3-LYP?", Chem. Phys. Lett. 268, 345 (1997). DOI: [10.1016/S0009-2614(97)00207-8](https://doi.org/10.1016/S0009-2614(97)00207-8)
- M. E. Casida, "Axel Dieter Becke: Not just the “B” in “B3LYP”", arXiv:2607.21465 (2026).
- H. Iikura, T. Tsuneda, T. Yanai, K. Hirao, "A long-range correction scheme for generalized-gradient-approximation exchange functionals", J. Chem. Phys. 115, 3540 (2001). DOI: [10.1063/1.1383587](https://doi.org/10.1063/1.1383587)
- J.-D. Chai, M. Head-Gordon, "Systematic optimization of long-range corrected hybrid density functionals", J. Chem. Phys. 128, 084106 (2008). DOI: [10.1063/1.2834918](https://doi.org/10.1063/1.2834918)
- N. Mardirossian, M. Head-Gordon, "ωB97X-V: A 10-parameter, range-separated hybrid, generalized gradient approximation density functional with nonlocal correlation, designed by a survival-of-the-fittest strategy", Phys. Chem. Chem. Phys. 16, 9904 (2014). DOI: [10.1039/C3CP54374A](https://doi.org/10.1039/C3CP54374A)
- N. Mardirossian, M. Head-Gordon, "ωB97M-V: A combinatorially optimized, range-separated hybrid, meta-GGA density functional with VV10 nonlocal correlation", J. Chem. Phys. 144, 214110 (2016). DOI: [10.1063/1.4952647](https://doi.org/10.1063/1.4952647)
- S. Grimme, "Semiempirical hybrid density functional with perturbative second-order correlation", J. Chem. Phys. 124, 034108 (2006). DOI: [10.1063/1.2148954](https://doi.org/10.1063/1.2148954)
- T. Helgaker, W. Klopper, H. Koch, J. Noga, "Basis-set convergence of correlated calculations on water", J. Chem. Phys. 106, 9639 (1997). DOI: [10.1063/1.473863](https://doi.org/10.1063/1.473863)
- L. Goerigk, A. Hansen, C. Bauer, S. Ehrlich, A. Najibi, S. Grimme, "A look at the density functional theory zoo with the advanced GMTKN55 database for general main group thermochemistry, kinetics and noncovalent interactions", Phys. Chem. Chem. Phys. 19, 32184 (2017). DOI: [10.1039/C7CP04913G](https://doi.org/10.1039/C7CP04913G)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- Libxc（汎関数のライブラリ）の汎関数の一覧：<https://libxc.gitlab.io/functionals/>
- PySCF の DFT モジュールの説明：<https://pyscf.org/user/dft.html>
