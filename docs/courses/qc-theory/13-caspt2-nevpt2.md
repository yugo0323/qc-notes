# 第 13 回　多参照摂動論：CASPT2 と NEVPT2

!!! abstract "この回の道筋"
    CASSCF は静的相関を正しく扱いますが、動的相関が足りません。この回では、CASSCF の波動関数を出発点にして、動的相関を 2 次の摂動で足す方法を学びます。代表的な二つの方法、CASPT2 と NEVPT2 は、0 次のハミルトニアンの選び方が違います。その違いから、侵入者状態の問題、レベルシフト、IPEA シフトをめぐる議論が生まれます。最後に、結合を伸ばした H₂O で、MP2、CASSCF、CASPT2、NEVPT2 を FCI と比べます。

    | 節 | 内容 |
    |---|---|
    | 1 | CASSCF の上に摂動を足す考え方 |
    | 2 | 1 次の波動関数の空間：内部縮約 |
    | 3 | CASPT2：一般化 Fock 演算子を使う 0 次ハミルトニアン |
    | 4 | 侵入者状態とレベルシフト |
    | 5 | IPEA シフトとその議論 |
    | 6 | NEVPT2：Dyall のハミルトニアン |
    | 7 | 計算例：H₂O の二つの O–H 結合を伸ばす |
    | 8 | 多状態版：MS-CASPT2、XMS-CASPT2、QD-NEVPT2 |

    **前に読んでおく回**：[量子化学の理論 第 9 回](09-mp2.md)、[第 12 回](12-casscf.md)

!!! tip "物理的なイメージ"
    CASPT2 と NEVPT2 が求めているのは、CASSCF の波動関数を出発点にして、活性空間の外へ電子が一瞬跳び出す過程（動的相関）を、2 次の摂動で足した補正です。第 9 回の MP2 の考え方を、多配置の出発点に広げたものです。各項の分母は「跳び出しに必要なエネルギー」の見積もりで、その見積もり方（0 次のハミルトニアン）が二つの方法の違いです。見積もりが偶然 0 に近くなると補正が発散し（侵入者状態）、これを防ぐための工夫がレベルシフトです。

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書きます。1 $E_\mathrm{h}$ = 27.211 eV = 627.5 kcal/mol、1 m$E_\mathrm{h}$ = $10^{-3}\ E_\mathrm{h}$ です。2 電子積分は、空間軌道についての化学者の表記 $(pq|rs)$ を使います。軌道の添字は、不活性軌道を $i, j$、活性軌道を $t, u, v, w$、空軌道を $a, b$、どれでもよい軌道を $p, q, r, s$ で表します。$\hat E_{pq}$ は第 4 回の、スピンを変えずに電子を $q$ から $p$ へ移す演算子です。

## 1. CASSCF の上に摂動を足す考え方

第 9 回の MP2 は、HF の決定子を 0 次の波動関数とし、Fock 演算子の和を 0 次のハミルトニアン $\hat H_0$ にとる摂動論でした。2 次のエネルギーは、次の式で与えられました。

$$
E^{(2)} = \sum_{n \ne 0} \frac{|\langle \Phi_n | \hat V | \Phi_0 \rangle|^2}{E_0^{(0)} - E_n^{(0)}}
$$

第 9 回の 6 節で見たように、静的相関が強いと、分母 $E_0^{(0)} - E_n^{(0)}$ が小さくなり、MP2 は信頼できなくなります。HF の決定子が、そもそも良い出発点ではないからです。

そこで、出発点を CASSCF の波動関数 $\Psi_0$ に替えます。静的相関は $\Psi_0$ の中に入っているので、残りの動的相関だけを摂動で足せばよいはずです。この考え方の方法を、まとめて**多参照摂動論**（multireference perturbation theory）と呼びます。エネルギーは

$$
E = E_\mathrm{CASSCF} + E^{(2)}
$$

の形になります。

ただし、MP2 をそのまま広げることはできません。摂動論を使うには、$\Psi_0$ を固有関数にもつ $\hat H_0$ が必要です。HF の決定子は Fock 演算子の和の固有関数でしたが、CASSCF の波動関数は多くの配置の重ね合わせなので、簡単な 1 電子演算子の固有関数にはなりません。$\hat H_0$ をどう選ぶかが、多参照摂動論の中心の問題です。この選び方の違いが、3 節の CASPT2 と 6 節の NEVPT2 の違いです。

**ここまでで分かったこと**：多参照摂動論は、CASSCF の波動関数を出発点にして、動的相関を 2 次の摂動で足す方法です。その鍵は、CASSCF の波動関数を固有関数にもつ 0 次のハミルトニアンの選び方です。

## 2. 1 次の波動関数の空間：内部縮約

$\hat H_0$ の話の前に、摂動で足す配置（摂動の相手）をどう選ぶかを決めます。

$\hat H$ には 2 電子演算子までしか含まれないので、$\hat H |\Psi_0\rangle$ に現れるのは、$\Psi_0$ から電子を最大 2 個動かした配置だけです。動かし方は、どの種類の軌道からどの種類の軌道へ動かすかで分類できます。たとえば、

- 不活性軌道 2 個から空軌道 2 個へ（$\hat E_{ai} \hat E_{bj}$ 型）
- 活性軌道から空軌道へ、不活性軌道から活性軌道へ、などの組み合わせ

です。活性軌道の中だけでの動かし方は、すでに CAS の中に含まれているので除きます。

CAS の配置 1 個ずつに、こうした励起をかけた配置をすべて使うと、配置の数は、CAS の配置の数と、外側の励起の数の積になり、非常に大きくなります。そこで、CASPT2 でも NEVPT2 の多くの実装でも、**励起の演算子を $\Psi_0$ 全体にかけた関数**

$$
\hat E_{pq} \hat E_{rs} |\Psi_0\rangle
$$

を摂動の相手に使います。$\Psi_0$ の中の CAS の係数の比は、固定したままです。これを**内部縮約**（internal contraction）と呼びます [Andersson 1990; Andersson 1992]。摂動の相手の数は、CAS の配置の数によらず、軌道の添字の組の数で決まります。そのかわり、計算には $\Psi_0$ の 3 電子・4 電子の密度行列が必要になります（[第 14 回](14-dmrg-sci.md)で、これが大きな活性空間での壁になることを見ます）。

**ここまでで分かったこと**：多参照摂動論では、励起の演算子を CASSCF の波動関数全体にかけた関数を、摂動の相手に使います（内部縮約）。配置の数は抑えられますが、高次の密度行列が必要になります。

## 3. CASPT2：一般化 Fock 演算子を使う 0 次ハミルトニアン

Andersson、Malmqvist、Roos らは、CASSCF の波動関数を出発点とする 2 次の摂動論 **CASPT2** を 1990 年と 1992 年に発表しました [Andersson 1990; Andersson 1992]。CASPT2 の $\hat H_0$ は、MP2 の Fock 演算子を多配置に広げたものです。

まず、CASSCF の 1 電子密度行列 $\gamma_{rs}$（第 12 回）を使って、**一般化 Fock 演算子**を定義します。

$$
\hat F = \sum_{pq} f_{pq}\, \hat E_{pq}, \qquad
f_{pq} = h_{pq} + \sum_{rs} \gamma_{rs} \Big[ (pq|rs) - \frac12 (pr|qs) \Big]
$$

$\Psi_0$ が閉殻の決定子 1 個のときは、$\gamma$ は占有軌道で 2、それ以外で 0 なので、$f_{pq}$ は HF の Fock 行列そのものになります。

ところが、$\hat F$ そのものは $\Psi_0$ を固有関数にもちません。そこで、空間をいくつかの部分に分け、それぞれの中だけで $\hat F$ を働かせます。

$$
\hat H_0 = \hat P_0 \hat F \hat P_0 + \hat P_K \hat F \hat P_K + \hat P_{SD} \hat F \hat P_{SD} + \hat P_X \hat F \hat P_X
$$

ここで、$\hat P_0$ は $\Psi_0$ への射影演算子、$\hat P_K$ は CAS の中で $\Psi_0$ に直交する部分への射影、$\hat P_{SD}$ は 2 節の摂動の相手がつくる空間への射影、$\hat P_X$ は残りへの射影です。こうすると、$\Psi_0$ は $\hat H_0$ の固有関数になり、固有値は $E_0^{(0)} = \langle \Psi_0 | \hat F | \Psi_0 \rangle$ です。

この $\hat H_0$ の利点は、$\Psi_0$ が閉殻の決定子 1 個（活性軌道がない場合）のとき、CASPT2 が MP2 に一致することです。第 9 回の MP2 を、そのまま多配置に広げたものと考えられます。

一方、弱点もあります。$\hat F$ は 1 電子演算子なので、活性軌道の中の電子どうしの反発を、平均の形でしか含みません。この弱点が、4 節と 5 節で扱う問題（侵入者状態と、開殻の系での系統的な誤差）につながります。

**ここまでで分かったこと**：CASPT2 は、CASSCF の密度行列でつくった一般化 Fock 演算子を、空間ごとに区切って 0 次のハミルトニアンにします。単一の閉殻配置では MP2 に一致します。

## 4. 侵入者状態とレベルシフト

### 侵入者状態

2 次のエネルギーの分母 $E_0^{(0)} - E_n^{(0)}$ は、摂動の相手の 0 次のエネルギーが $E_0^{(0)}$ に近いと、0 に近づきます。分母が 0 に近い項は非常に大きくなり、エネルギーが発散したり、ポテンシャルエネルギー曲線にとげが現れたりします。このような摂動の相手を、**侵入者状態**（intruder state）と呼びます。

CASPT2 では、侵入者状態が起こりえます。$\hat H_0$ が 1 電子演算子なので、摂動の相手の 0 次のエネルギーは軌道のエネルギーの足し引きで決まり、活性軌道の電子の反発がきちんと入りません。そのため、本当はエネルギーの高い配置が、0 次では $\Psi_0$ と同じくらいのエネルギーに見えることがあるからです。

### 実数のレベルシフト

侵入者状態への対処として、Roos と Andersson は、$\hat H_0$ に定数 $\varepsilon$ を足して分母をずらす**レベルシフト**を提案しました [Roos 1995]。分母の 1 個の項で書くと、

$$
\frac{1}{E_0^{(0)} - E_n^{(0)}} \;\longrightarrow\; \frac{1}{E_0^{(0)} - E_n^{(0)} - \varepsilon}
$$

です（$\varepsilon > 0$）。基底状態では、分母は負なので、$\varepsilon$ を足すと分母の大きさが増え、0 に近づかなくなります。シフトのせいで変わった 2 次のエネルギーは、近似的に補正して元に戻します [Roos 1995]。彼らは、この方法で侵入者状態が取り除け、残りの相関エネルギーへの影響は小さいことを示しました。

### 虚数のレベルシフト

実数のシフトにも限界があります。励起状態では、分母が正にも負にもなりうるので、小さなシフトは特異点の位置を動かすだけです。大きなシフトは、別の場所に新しい発散を生むことがあります。Forsberg と Malmqvist はこの点を指摘し、**虚数のレベルシフト**を提案しました [Forsberg 1997]。分母に虚数 $i\varepsilon$ を足して、エネルギーの実部だけをとると、

$$
\mathrm{Re}\, \frac{1}{\Delta + i\varepsilon} = \frac{\Delta}{\Delta^2 + \varepsilon^2}, \qquad \Delta = E_0^{(0)} - E_n^{(0)}
$$

となります。右辺は、$\Delta$ が正でも負でも 0 でも有限です。$|\Delta|$ が $\varepsilon$ より十分に大きい項は、元の $1/\Delta$ とほとんど変わりません。Forsberg らは、N₂ と Cr₂ への応用で、この操作によるポテンシャルの歪みが小さいことを示しました [Forsberg 1997]。どちらのシフトを使う場合も、使った値を必ず報告し、結果がその値にどれだけ依存するかを確かめてください。

**ここまでで分かったこと**：分母が 0 に近い摂動の相手（侵入者状態）は、エネルギーを発散させます。CASPT2 では、実数または虚数のレベルシフトでこれを抑えます。

## 5. IPEA シフトとその議論

CASPT2 には、侵入者状態とは別に、系統的な誤差も知られています。Ghigo、Roos、Malmqvist は、元の CASPT2 が開殻の系の相関エネルギーを大きく見積もりすぎ、その結果、解離エネルギーや励起エネルギーが小さくなりすぎることを指摘しました [Ghigo 2004]。

原因は、活性軌道の対角の Fock 要素 $f_{tt}$ にあります。3 節の $f_{pq}$ の式から、$f_{tt}$ は、その軌道から電子を取り去るエネルギー（イオン化エネルギー，IP）と、電子を加えるエネルギー（電子親和力，EA）を、占有数で重みづけた平均に近い値になります [Ghigo 2004]。

$$
f_{tt} \approx -\frac12 \Big[ \gamma_{tt}\, \mathrm{IP}_t + (2 - \gamma_{tt})\, \mathrm{EA}_t \Big]
$$

占有数 $\gamma_{tt}$ が 2 なら $-\mathrm{IP}_t$、0 なら $-\mathrm{EA}_t$ になり、それぞれ正しい値です。ところが占有数が 1 に近い開殻の軌道では、IP と EA の平均になります。その軌道から電子を取り去る励起にも、加える励起にも、同じ平均の値が使われるので、分母が小さく見積もられ、相関エネルギーが大きくなりすぎます。

Ghigo らは、活性軌道から電子を取り去る励起と、活性軌道へ電子を加える励起とで、分母を別々にずらす修正を提案しました。ずらす大きさは、パラメータ $\varepsilon$ と活性軌道の占有数で決まります。これを **IPEA シフト**と呼び、標準の値は $\varepsilon = 0.25\ E_\mathrm{h}$ です。49 個の二原子分子の解離エネルギーで、誤差の根平均二乗が 0.22 eV（シフトなし）から 0.09 eV（$\varepsilon = 0.25\ E_\mathrm{h}$）に下がったと報告されています [Ghigo 2004; Zobel 2017 のまとめによる]。

一方、励起エネルギーについては議論があります。

- Zobel、Nogueira、González は、有機分子の励起エネルギーについて、IPEA シフトを使う根拠がないと結論しました [Zobel 2017]。文献に報告された、シフトなしの CASPT2 の励起エネルギー 356 個を参照値と比べると、平均の過小評価は 0.02 eV にすぎませんでした。13 個の二原子・三原子分子の 137 個の状態で FCI と比べても、過小評価は 0.05 eV でした。さらに、IPEA シフトによる変化は、系の大きさ（動的相関の大きさ）や基底関数に強く依存し、普遍的な補正にはなっていないことを示しました。
- Sarkar、Loos、Boggio-Pasqua、Jacquemin は、35 個の有機分子の 284 個の励起状態について、誤差 ±0.03 eV 程度の高精度な参照値と比べました [Sarkar 2022]。その結果、IPEA シフトを使った CASPT2 と、部分縮約版の NEVPT2（6 節）が、平均絶対誤差 0.11 eV と 0.13 eV で、どちらもかなり信頼できる値を与えました。

このように、IPEA シフトの是非は、性質（解離エネルギーか励起エネルギーか）、系の種類、参照値の選び方によって評価が分かれています。CASPT2 を使うときは、IPEA シフトの有無と値を明記し、比べたい文献と条件をそろえてください。

**ここまでで分かったこと**：IPEA シフトは、開殻の活性軌道で分母が小さく見積もられる CASPT2 の系統的な誤差を補正する経験的な修正です。解離エネルギーは改善しますが、励起エネルギーへの使用には賛否があります。

## 6. NEVPT2：Dyall のハミルトニアン

CASPT2 の問題の多くは、$\hat H_0$ が 1 電子演算子で、活性電子どうしの反発を平均でしか扱わないことから来ていました。そこで、**活性軌道の中では本物の 2 電子相互作用を残す**という選び方があります。

Dyall は、次の 0 次のハミルトニアンを提案しました [Dyall 1995]。

$$
\hat H_\mathrm{D} = \sum_i \varepsilon_i \hat E_{ii} + \sum_a \varepsilon_a \hat E_{aa} + C + \hat H_\mathrm{act}
$$

$$
\hat H_\mathrm{act} = \sum_{tu} h_{tu}^\mathrm{eff}\, \hat E_{tu} + \frac12 \sum_{tuvw} (tu|vw) \big( \hat E_{tu} \hat E_{vw} - \delta_{uv} \hat E_{tw} \big)
$$

各項の意味は次のとおりです。

- 不活性軌道と空軌道には、MP2 と同じく軌道エネルギー $\varepsilon_i$、$\varepsilon_a$ を使います。
- 活性軌道には、本物のハミルトニアンを活性軌道だけに制限したもの $\hat H_\mathrm{act}$ を使います。$h_{tu}^\mathrm{eff} = h_{tu} + \sum_i \big[ 2(tu|ii) - (ti|iu) \big]$ は、不活性軌道の電子がつくる平均の場を含めた 1 電子積分です。
- $C$ は、$\Psi_0$ の固有値を本物のエネルギーに合わせるための定数です。

CASSCF の波動関数は、活性空間の中での FCI なので、$\hat H_\mathrm{act}$ の固有関数です。したがって、$\Psi_0$ はそのまま $\hat H_\mathrm{D}$ の固有関数になり、CASPT2 のように射影演算子で区切る必要がありません。

Angeli、Cimiraglia、Malrieu らは、この Dyall のハミルトニアンを使った 2 次の摂動論 **NEVPT2**（$n$-electron valence state perturbation theory）を 2001 年に発表しました [Angeli 2001]。名前の「$n$ 電子の原子価状態」は、活性電子の状態を、2 電子相互作用を含めたまま扱うことを指しています。NEVPT2 には、次の性質があります [Angeli 2001]。

1. **侵入者状態がありません。** 摂動の相手の 0 次のエネルギーには、活性電子どうしの反発がきちんと入ります。そのため、分母が 0 に近づくことが起こらず、レベルシフトが要りません。
2. **サイズ無矛盾です。** 十分に離れた二つの分子の NEVPT2 のエネルギーは、それぞれのエネルギーの和になります（第 8 回のサイズ無矛盾性）。

NEVPT2 には、摂動の相手をどこまで縮約するかで、**強縮約**（strongly contracted, SC）と**部分縮約**（partially contracted, PC）の二つの版があります [Angeli 2001; Angeli 2002]。SC 版は、同じ種類の励起をまとめて 1 個の関数にするので、安く計算できます。PC 版は、2 節の内部縮約に近い形で、より柔軟です。

NEVPT2 の計算には、活性軌道の 4 電子の密度行列が必要です [Angeli 2002]。CASPT2 でも、4 電子の密度行列を一般化 Fock 行列と縮約した量が必要です。活性空間が大きくなると、この部分が計算の重さを決めます。

**ここまでで分かったこと**：NEVPT2 は、活性軌道の中で本物の 2 電子相互作用を残す Dyall のハミルトニアンを 0 次にとります。そのため侵入者状態がなく、レベルシフトが要りません。サイズ無矛盾性も満たします。

## 7. 計算例：H₂O の二つの O–H 結合を伸ばす

第 8 回と第 9 回で使った H₂O（6-31G 基底）を、もう一度計算します。平衡に近い構造と、二つの O–H 結合を 1.5 倍、2 倍に伸ばした構造です。第 9 回で見たとおり、2 倍に伸ばすと MP2 は大きく外れました。この系では FCI が計算できるので、各方法の誤差を正確に比べられます。

活性空間は、二つの O–H 結合の σ 軌道 2 個と σ\* 軌道 2 個に、電子 4 個を入れた CAS(4,4) にします。PySCF [Sun 2020] の `mrpt.NEVPT` は SC-NEVPT2、`mrpt.caspt2.CASPT2` は内部縮約の CASPT2 です。CASPT2 は、IPEA シフトなしと、$\varepsilon = 0.25\ E_\mathrm{h}$ の両方で計算します。どの方法でも、すべての電子を相関させています。

```python
import numpy as np
from pyscf import gto, scf, mp, mcscf, mrpt, fci

def h2o(scale):
    """平衡構造（O–H 0.96 Å、HOH 104.5°）の O–H 結合を scale 倍に伸ばす。"""
    r, th = 0.96 * scale, np.deg2rad(104.5) / 2
    return f"O 0 0 0; H 0 {r*np.sin(th):.6f} {r*np.cos(th):.6f}; H 0 {-r*np.sin(th):.6f} {r*np.cos(th):.6f}"

print("scale    HF      MP2   CASSCF   NEVPT2  CASPT2  CASPT2(IPEA)  （E - E_FCI、単位 mEh）")
for scale in [1.0, 1.5, 2.0]:
    mol = gto.M(atom=h2o(scale), basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    e_mp2 = mp.MP2(mf).run().e_tot
    mc = mcscf.CASSCF(mf, 4, 4).run()          # 活性空間：O–H の σ 2 個と σ* 2 個
    e_nev = mc.e_tot + mrpt.NEVPT(mc).kernel()  # SC-NEVPT2 の相関エネルギーを足す
    e_pt2 = mrpt.caspt2.CASPT2(mc).kernel()[0]
    e_ipea = mrpt.caspt2.CASPT2(mc, ipea_shift=0.25).kernel()[0]
    e_fci = fci.FCI(mf).kernel()[0]
    errs = [mf.e_tot, e_mp2, mc.e_tot, e_nev, e_pt2, e_ipea]
    print(f"{scale:4.1f} " + " ".join(f"{1000*(e - e_fci):7.2f}" for e in errs))
```

```
scale    HF      MP2   CASSCF   NEVPT2  CASPT2  CASPT2(IPEA)  （E - E_FCI、単位 mEh）
 1.0  137.12    8.07   83.85    6.54    7.44    7.80
 1.5  197.02   19.07   80.36    7.63    6.73    7.59
 2.0  293.05   45.74   69.89    9.75    4.91    6.43
```

全体で 2 分ほどかかり、そのほとんどは FCI の計算です。表の値は、各方法のエネルギーと FCI のエネルギーの差です。

表から、次のことが読み取れます。

- **MP2**：誤差が 8 m$E_\mathrm{h}$ から 46 m$E_\mathrm{h}$ へ、6 倍近くに増えます。第 9 回で見たとおりです。
- **CASSCF**：誤差は 70〜84 m$E_\mathrm{h}$ と、どの構造でも大きいままです。これは第 12 回で見た、動的相関の不足です。結合を伸ばすと誤差が少し減るのは、伸ばした結合の分の動的相関が減るからです。
- **NEVPT2 と CASPT2**：誤差はどの構造でも 5〜10 m$E_\mathrm{h}$ に収まり、構造による変化も 3 m$E_\mathrm{h}$（約 2 kcal/mol）程度です。CASSCF で静的相関を、摂動で動的相関を扱う二段構えが、うまく働いています。

化学で大事なのは、全エネルギーそのものより、構造の間のエネルギー差です。構造によって誤差がほぼ一定であれば、エネルギー差の誤差は小さくなります。MP2 では平衡構造と 2 倍の構造の間で誤差が 38 m$E_\mathrm{h}$（約 24 kcal/mol）変わるのに対し、NEVPT2 と CASPT2 では 3 m$E_\mathrm{h}$ 程度です。

ただし、これは小さな基底関数の 1 個の例にすぎません。IPEA シフトの有無や NEVPT2 と CASPT2 のどちらが良いかを、この表から一般的に結論することはできません（5 節のとおり、系と性質によって評価が分かれます）。

**ここまでで分かったこと**：結合を伸ばした H₂O で、CASPT2 と NEVPT2 の誤差は構造によらずほぼ一定でした。MP2 のように静的相関で破綻することも、CASSCF のように動的相関が大きく足りないこともありません。

## 8. 多状態版：MS-CASPT2、XMS-CASPT2、QD-NEVPT2

ここまでの方法は、CASSCF の状態 1 個ずつに摂動を足します。これを**単状態**（single-state）の摂動論と呼びます。二つの状態のエネルギーが近いと、単状態の方法では困ることがあります。

たとえば、イオン性の状態と共有結合性の状態が交わる回避交差の近くでは、CASSCF の段階での二つの状態の混ざり方が、動的相関を入れると変わります。CASSCF の状態を固定したまま別々に摂動を足すと、この混ざり方の変化を取り込めず、ポテンシャルエネルギー曲線がおかしくなります。

**MS-CASPT2**（multi-state CASPT2）は、この問題に対処します [Finley 1998]。考えたい状態 $\Psi_I$（第 12 回の状態平均 CASSCF で得たもの）について、それぞれの 1 次の波動関数 $\Psi_I^{(1)}$ を求めます。そのうえで、状態の間の行列

$$
H^\mathrm{eff}_{IJ} = \delta_{IJ} E_I^\mathrm{CAS} + \frac12 \Big( \langle \Psi_I | \hat H | \Psi_J^{(1)} \rangle + \langle \Psi_J | \hat H | \Psi_I^{(1)} \rangle \Big)
$$

を有効ハミルトニアンとして対角化します。$E_I^\mathrm{CAS}$ は状態 $I$ の CASSCF のエネルギーです。対角の要素は、単状態の CASPT2 のエネルギーと同じです。こうすると、摂動の後で状態どうしが混ざり直すことができます。

MS-CASPT2 にも問題がありました。各状態の $\hat H_0$ が別々に定義されているため、結果が、出発点の CASSCF の状態をどう組み合わせて選ぶかに依存しうるのです。Granovsky は、すべての状態に共通の 0 次のハミルトニアンを使い、この依存性をなくした XMCQDPT2 を提案しました [Granovsky 2011]。同じ考えを CASPT2 に取り入れたのが **XMS-CASPT2** で、回避交差や円錐交差の近くでより良いポテンシャルを与えることが示されています [Shiozaki 2011]。

NEVPT2 にも、同じ目的の多状態版 **QD-NEVPT2**（quasi-degenerate NEVPT2）があります [Angeli 2004]。LiF のイオン性と共有結合性の状態の回避交差や、エチレンの原子価状態とリュードベリ状態の混ざり合いで、その有用性が示されました。

**ここまでで分かったこと**：エネルギーの近い状態がある場合は、有効ハミルトニアンを対角化して状態を混ぜ直す多状態版を使います。CASPT2 では MS 版と XMS 版、NEVPT2 では QD 版があります。

## まとめ

- 多参照摂動論は、CASSCF の波動関数を出発点にし、動的相関を 2 次の摂動で足します。摂動の相手には、励起の演算子を波動関数全体にかけた内部縮約の関数を使います。
- CASPT2 は、一般化 Fock 演算子から 0 次のハミルトニアンをつくります。単一の閉殻配置では MP2 に一致しますが、侵入者状態が起こりうるので、実数や虚数のレベルシフトを使います。
- IPEA シフトは開殻の軌道での系統的な誤差を補正する経験的な修正で、解離エネルギーを改善しますが、励起エネルギーへの使用には賛否があります。
- NEVPT2 は、活性軌道に本物の 2 電子相互作用を残す Dyall のハミルトニアンを使い、侵入者状態がなく、サイズ無矛盾です。
- 状態どうしが近い場合は、MS/XMS-CASPT2 や QD-NEVPT2 などの多状態版を使います。

## 参考文献

- T. Helgaker, P. Jørgensen, J. Olsen, "Molecular Electronic-Structure Theory" (Wiley, 2000). ISBN: 9780471967552
- K. Andersson, P.-Å. Malmqvist, B. O. Roos, A. J. Sadlej, K. Wolinski, "Second-order perturbation theory with a CASSCF reference function", J. Phys. Chem. 94, 5483 (1990). DOI: [10.1021/j100377a012](https://doi.org/10.1021/j100377a012)
- K. Andersson, P.-Å. Malmqvist, B. O. Roos, "Second-order perturbation theory with a complete active space self-consistent field reference function", J. Chem. Phys. 96, 1218 (1992). DOI: [10.1063/1.462209](https://doi.org/10.1063/1.462209)
- B. O. Roos, K. Andersson, "Multiconfigurational perturbation theory with level shift — the Cr2 potential revisited", Chem. Phys. Lett. 245, 215 (1995). DOI: [10.1016/0009-2614(95)01010-7](https://doi.org/10.1016/0009-2614(95)01010-7)
- N. Forsberg, P.-Å. Malmqvist, "Multiconfiguration perturbation theory with imaginary level shift", Chem. Phys. Lett. 274, 196 (1997). DOI: [10.1016/S0009-2614(97)00669-6](https://doi.org/10.1016/S0009-2614(97)00669-6)
- G. Ghigo, B. O. Roos, P.-Å. Malmqvist, "A modified definition of the zeroth-order Hamiltonian in multiconfigurational perturbation theory (CASPT2)", Chem. Phys. Lett. 396, 142 (2004). DOI: [10.1016/j.cplett.2004.08.032](https://doi.org/10.1016/j.cplett.2004.08.032)
- J. P. Zobel, J. J. Nogueira, L. González, "The IPEA dilemma in CASPT2", Chem. Sci. 8, 1482 (2017). DOI: [10.1039/c6sc03759c](https://doi.org/10.1039/c6sc03759c)
- R. Sarkar, P.-F. Loos, M. Boggio-Pasqua, D. Jacquemin, "Assessing the Performances of CASPT2 and NEVPT2 for Vertical Excitation Energies", J. Chem. Theory Comput. 18, 2418 (2022). DOI: [10.1021/acs.jctc.1c01197](https://doi.org/10.1021/acs.jctc.1c01197)
- K. G. Dyall, "The choice of a zeroth-order Hamiltonian for second-order perturbation theory with a complete active space self-consistent-field reference function", J. Chem. Phys. 102, 4909 (1995). DOI: [10.1063/1.469539](https://doi.org/10.1063/1.469539)
- C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J.-P. Malrieu, "Introduction of n-electron valence states for multireference perturbation theory", J. Chem. Phys. 114, 10252 (2001). DOI: [10.1063/1.1361246](https://doi.org/10.1063/1.1361246)
- C. Angeli, R. Cimiraglia, J.-P. Malrieu, "n-electron valence state perturbation theory: A spinless formulation and an efficient implementation of the strongly contracted and of the partially contracted variants", J. Chem. Phys. 117, 9138 (2002). DOI: [10.1063/1.1515317](https://doi.org/10.1063/1.1515317)
- J. Finley, P.-Å. Malmqvist, B. O. Roos, L. Serrano-Andrés, "The multi-state CASPT2 method", Chem. Phys. Lett. 288, 299 (1998). DOI: [10.1016/S0009-2614(98)00252-8](https://doi.org/10.1016/S0009-2614(98)00252-8)
- A. A. Granovsky, "Extended multi-configuration quasi-degenerate perturbation theory: The new approach to multi-state multi-reference perturbation theory", J. Chem. Phys. 134, 214113 (2011). DOI: [10.1063/1.3596699](https://doi.org/10.1063/1.3596699)
- T. Shiozaki, W. Győrffy, P. Celani, H.-J. Werner, "Communication: Extended multi-state complete active space second-order perturbation theory: Energy and nuclear gradients", J. Chem. Phys. 135, 081106 (2011). DOI: [10.1063/1.3633329](https://doi.org/10.1063/1.3633329)
- C. Angeli, S. Borini, M. Cestari, R. Cimiraglia, "A quasidegenerate formulation of the second order n-electron valence state perturbation theory approach", J. Chem. Phys. 121, 4043 (2004). DOI: [10.1063/1.1778711](https://doi.org/10.1063/1.1778711)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- PySCF の多参照摂動論のモジュール（`pyscf.mrpt`）：https://pyscf.org/user/mrpt.html
