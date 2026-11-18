# 第 5 回　Gaussian（2）：振動解析、熱化学、溶媒効果

!!! abstract "この回の道筋"
    Gaussian の `Freq` の出力から、振動数、零点エネルギー、エンタルピー、自由エネルギーを読み取る方法を学びます。
    次に、自由エネルギーを比べるときに必要な標準状態の補正と、低い振動数のモードの扱いを考えます。
    最後に、`SCRF` キーワードで溶媒の効果を取り入れる方法（PCM、SMD）を学びます。

    | 節 | 内容 |
    |---|---|
    | 1 | 振動数の出力を読む |
    | 2 | 熱化学の出力を読む |
    | 3 | 標準状態：1 atm と 1 M |
    | 4 | 低い振動数のモードと準調和近似 |
    | 5 | 溶媒の効果：SCRF、PCM、SMD |

    **前に読んでおく回**：計算パッケージの使い方 [第 1 回](01-workflow.md)、[第 4 回](04-gaussian-basics.md)

!!! tip "物理的なイメージ"
    振動解析が求めているのは、谷底のまわりでのエネルギーの曲がり具合（ばね定数の行列、ヘッセ行列）と、そこから決まる原子の振動数です。各振動を調和振動子とみなすと、温度が 0 でも振動数 ν ごとに ½hν のエネルギー（零点エネルギー）が残ります。温度を上げると、振動、回転、並進にエネルギーが分配され、その分配の数え方（分配関数）から、エンタルピー、エントロピー、自由エネルギーが求まります。振動数の低いモードは、調和振動子の近似が悪いのに、エントロピーには大きく効くので注意が必要です。

## 1. 振動数の出力を読む

この回でも、第 4 回の水分子の `HF/STO-3G Opt Freq` の計算を例にします。
`Freq` の部分の出力を、上から順に読みます。

### 1.1 Low frequencies の行

振動数の表の前に、次の 2 行が出ます。

```text
 Low frequencies ---  -51.8680  -51.7821  -51.6444   -0.0019   -0.0007   -0.0004
 Low frequencies --- 2169.8344 4141.7948 4392.8972
```

1 行目の 6 個は、分子全体の並進（3 個）と回転（3 個）にあたる値です。
直線分子では回転が 2 個なので、並進と回転にあたるのは 5 個です。
並進と回転では分子の形が変わらないので、これらは理論上は 0 になるはずです。
実際には、構造最適化の収束の甘さや数値計算の誤差のため、0 から少しずれます。
この例では、絶対値で 52 cm⁻¹ 以下でした。
公式の文書は、これらの値が 0 から大きくずれているときは、収束が不十分である可能性があり、最適化の条件を厳しくするか、DFT なら積分のグリッドを細かくするよう勧めています。

2 行目は、並進と回転を除いた本当の振動の振動数を、低いほうから並べたものです。
水では、振動は $3N - 6 = 3$ 個なので、3 個すべてが出ています。

### 1.2 振動数の表

その後ろに、振動ごとの詳しい表が出ます。
主な行は次のとおりです。

| 行 | 意味 |
|---|---|
| `Frequencies --` | 振動数（cm⁻¹） |
| `Red. masses --` | その振動の換算質量 |
| `IR Inten --` | 赤外吸収の強度 |
| 原子ごとの 3 列の数字 | その振動で各原子が動く向き（基準振動の形） |

この例の振動数は、2169.8、4141.8、4392.9 cm⁻¹ でした。
それぞれ、変角振動、対称伸縮振動、逆対称伸縮振動にあたります。
実験の値（変角はおよそ 1600 cm⁻¹）より大きいのは、電子相関を含まない HF 法と、小さすぎる基底関数 STO-3G を使ったうえ、調和近似（1.4 節）で計算したからです。
ここでは出力の読み方だけに注目します。

振動の形は、数字の表より、可視化ソフト（第 2 回）でアニメーションにして見るほうが分かりやすくなります。

### 1.3 虚振動

第 1 回で見たとおり、虚振動は負の振動数として出力されます。
虚振動があると、表の前に次の形の行が出ます。

```text
 ******    2 imaginary frequencies (negative Signs) ******
```

これは、O–H の長さを最適化した値（0.989 Å）のまま、わざと直線にした水分子で振動解析をしたときの出力です。
変角振動にあたる 2 個の縮重した振動数が、どちらも約 $-2270$ cm⁻¹ になりました。
極小を求めたつもりで虚振動が出たときは、その振動の向きに構造を少しずらして、最適化をやり直します。
遷移状態では、虚振動が 1 個で、その形が目的の反応の向きになっていることを確かめます。

### 1.4 振動数のスケーリング

Gaussian が出す振動数は、調和近似（エネルギーを 2 次の式で近似したもの）で計算した値です。
調和近似の振動数は、実験の振動数より系統的に大きくなる傾向があります。
そこで、計算レベルごとに決めた係数（スケーリング因子）を掛けて、実験と比べることがあります（Scott, Radom, 1996）。
Gaussian では、`Scale` キーワードで係数を指定できます。
係数は計算レベルごとに違うので、使う計算レベルについて決められた値を文献で確かめます。

**ここまでで分かったこと**：`Low frequencies` の最初の 6 個は並進と回転で 0 に近いはずの値、`Frequencies` が本当の振動数です。虚振動は負の値で表されます。

## 2. 熱化学の出力を読む

### 2.1 何を計算しているか

振動数の表の後ろに、熱化学（thermochemistry）の出力が続きます。
Gaussian は、分子を理想気体として扱い、統計力学の式で温度による補正を計算します。
このとき、並進、回転、振動、電子の運動を分けて考えます。
振動はすべて調和振動子として、回転は剛体の回転として扱います。
この近似を、剛体回転子・調和振動子（RRHO）近似と呼びます。

公式の文書によると、既定の条件は 298.15 K、1 atm です。
温度は `Temperature`、圧力は `Pressure` キーワードで変えられます。
たとえば、ルートセクションに `Temperature=373.15` を加えると、373.15 K での値が出ます。

### 2.2 出力の行

出力のうち、よく使う 8 行は次のとおりです（単位はハートリー）。

```text
 Zero-point correction=                           0.024387 (Hartree/Particle)
 Thermal correction to Energy=                    0.027220
 Thermal correction to Enthalpy=                  0.028164
 Thermal correction to Gibbs Free Energy=         0.006649
 Sum of electronic and zero-point Energies=            -74.941514
 Sum of electronic and thermal Energies=               -74.938682
 Sum of electronic and thermal Enthalpies=             -74.937737
 Sum of electronic and thermal Free Energies=          -74.959252
```

上の 4 行が補正の量、下の 4 行が、電子エネルギー $E_\mathrm{el}$（`SCF Done` の値）にそれぞれの補正を足したものです。
この例では $E_\mathrm{el} = -74.965901$ $E_\mathrm{h}$ なので、たとえば $-74.965901 + 0.006649 = -74.959252$ です。

補正の 4 行の意味を、上から順に説明します。

**零点エネルギー（ZPE）**：量子力学では、振動の最も低い状態でも、エネルギーは谷の底より高くなります。
振動数 $\nu_i$ の振動 1 個あたり $h\nu_i/2$ だけ高く、これを全部の振動について足したものです。

$$
E_\mathrm{ZPE} = \sum_i \frac{1}{2} h \nu_i
$$

この例で確かめると、$(2169.83 + 4141.79 + 4392.90)/2 = 5352.26$ cm⁻¹ です。
1 $E_\mathrm{h}$ = 219474.6 cm⁻¹ で換算すると 0.024387 $E_\mathrm{h}$ となり、出力と一致します。

**内部エネルギーの補正**（`Thermal correction to Energy`）：零点エネルギーに、温度 $T$ で並進、回転、振動が励起されることによるエネルギーを加えたものです。

**エンタルピーの補正**：理想気体では、1 分子あたりのエンタルピーは $H = U + k_\mathrm{B} T$ です（$k_\mathrm{B}$ はボルツマン定数）。
この例でも、$0.028164 - 0.027220 = 0.000944$ $E_\mathrm{h}$ で、298.15 K の $k_\mathrm{B} T$（0.000944 $E_\mathrm{h}$）に一致します。

**ギブズ自由エネルギーの補正**：エンタルピーの補正から、温度とエントロピーの積を引いたものです。

$$
G_\mathrm{corr} = H_\mathrm{corr} - T S
$$

エントロピー $S$ の内訳（並進、回転、振動、電子）は、そのすぐ後ろの表に、cal/(mol K) の単位で出ます。

### 2.3 自由エネルギーを比べるときの注意

反応の自由エネルギーは、各分子の `Sum of electronic and thermal Free Energies` の差から求めます。
第 1 回で見たとおり、すべての分子を同じ計算レベルで計算します。
よく使われるのは、構造と振動数を軽いレベルで計算し、電子エネルギーだけを重いレベルで計算し直す方法です。
このときは、次のように足し合わせます。

$$
G = E_\mathrm{el}(\text{重いレベル}) + G_\mathrm{corr}(\text{軽いレベル})
$$

遷移状態の熱化学では、虚振動は計算から除かれます。
出力にも、その旨の行（`imaginary frequencies ignored`）が出ます。

**ここまでで分かったこと**：熱化学の出力は、電子エネルギーに足す補正（ZPE、内部エネルギー、エンタルピー、自由エネルギー）と、その和です。既定は理想気体、298.15 K、1 atm です。

## 3. 標準状態：1 atm と 1 M

2 節の自由エネルギーは、1 atm の理想気体の値です。
一方、溶液中の反応では、各分子の濃度が 1 mol/L（1 M）の状態を標準にすることがふつうです。
二つの標準状態では、同じ分子でも並進のエントロピーが違うので、自由エネルギーが違います。

導出は付録 A に回しますが、1 atm から 1 M に移すときの自由エネルギーの変化は次のとおりです。

$$
\Delta G^{\circ}_{\mathrm{1\,atm} \to \mathrm{1\,M}} = RT \ln\frac{c^{\circ} RT}{p^{\circ}} = RT \ln 24.46 = 1.89\ \text{kcal/mol} \quad (T = 298.15\ \text{K})
$$

ここで $p^{\circ}$ = 1 atm、$c^{\circ}$ = 1 mol/L、$R$ は気体定数です。
24.46 は、298.15 K、1 atm の理想気体 1 mol の体積（L）です。
この補正は、分子 1 個ごとに同じ量だけ自由エネルギーを上げます（Bursch ら, 2022）。

そのため、反応の前後で分子の数が変わらなければ、補正は打ち消し合います。
分子の数が変わる反応では、結果が変わります。
たとえば A + B → C では、反応物が 2 個、生成物が 1 個なので、反応の自由エネルギーは 1 M の標準で $1.89$ kcal/mol だけ低くなります。
二つの分子が結合する反応や、配位子が外れる反応では、この補正を忘れないようにします。

**ここまでで分かったこと**：Gaussian の自由エネルギーは 1 atm の標準です。溶液の 1 M の標準にするには、分子 1 個ごとに 298.15 K で 1.89 kcal/mol を加えます。

## 4. 低い振動数のモードと準調和近似

### 4.1 何が問題か

大きく柔らかい分子には、振動数が 100 cm⁻¹ より低いモードがたくさんあります。
メチル基の回転や、分子全体がゆっくり曲がる運動などです。
このようなモードは、調和振動子というより、ほとんど自由に回る回転に近い運動です。

調和振動子のエントロピーは、振動数が 0 に近づくと際限なく大きくなります。
振動数 $\nu$ の調和振動子 1 個のエントロピーは、$x = h\nu/(k_\mathrm{B}T)$ として次の式で表されます。

$$
S_\mathrm{vib} = R \left[ \frac{x}{e^{x} - 1} - \ln\left(1 - e^{-x}\right) \right]
$$

$x \to 0$ では、右辺の第 2 項が $-\ln x$ のようにふるまうので、$S_\mathrm{vib}$ は発散します。
この式で、298.15 K での 1 モードあたりの値を計算してみます。

```python
import numpy as np

h = 6.62607015e-34      # プランク定数（J s）
c = 2.99792458e10       # 光速（cm/s）
kB = 1.380649e-23       # ボルツマン定数（J/K）
R = 1.987204e-3         # 気体定数（kcal/(mol K)）
T = 298.15              # 温度（K）

def s_harmonic(nu):
    """振動数 nu（cm^-1）の調和振動子 1 個のエントロピー（kcal/(mol K)）"""
    x = h * c * nu / (kB * T)
    return R * (x / np.expm1(x) - np.log(-np.expm1(-x)))

for nu in [5, 10, 20, 50, 100, 200, 500, 1000]:
    s = s_harmonic(nu)
    print(f"{nu:5d} cm^-1:  S = {s*1000:6.2f} cal/(mol K),  T*S = {T*s:5.2f} kcal/mol")
```

出力は次のとおりです。

```text
    5 cm^-1:  S =   9.39 cal/(mol K),  T*S =  2.80 kcal/mol
   10 cm^-1:  S =   8.01 cal/(mol K),  T*S =  2.39 kcal/mol
   20 cm^-1:  S =   6.63 cal/(mol K),  T*S =  1.98 kcal/mol
   50 cm^-1:  S =   4.82 cal/(mol K),  T*S =  1.44 kcal/mol
  100 cm^-1:  S =   3.45 cal/(mol K),  T*S =  1.03 kcal/mol
  200 cm^-1:  S =   2.13 cal/(mol K),  T*S =  0.64 kcal/mol
  500 cm^-1:  S =   0.66 cal/(mol K),  T*S =  0.20 kcal/mol
 1000 cm^-1:  S =   0.09 cal/(mol K),  T*S =  0.03 kcal/mol
```

500 cm⁻¹ を超えるモードの寄与は小さいのに対し、低いモードは 1 個で 1〜3 kcal/mol も $-TS$ を変えます。
しかも、低い振動数の値そのものが、最適化の収束やグリッドでわずかに変わるだけで、大きく動きます。
低いモードが多い分子では、自由エネルギーの誤差がたまりやすいということです。

### 4.2 準調和近似の考え方

そこで、低いモードのエントロピーだけを、別の扱いに置き換える方法が提案されています。
まとめて準調和近似（quasi-harmonic approximation）と呼ばれ、よく使われるのは次の二つです。

- **振動数を引き上げる方法**：決めた値（100 cm⁻¹ がよく使われます）より低い振動数を、その値に置き換えてエントロピーを計算します。Cramer、Truhlar らの準調和近似（Ribeiro ら, 2011）に由来する方法として知られています。
- **自由回転子と混ぜる方法**：低いモードのエントロピーを、調和振動子と自由回転子の値の重みつきの平均で表します（Grimme, 2012）。重みは、振動数 $\nu$ が基準の値 $\nu_0$ より十分高ければ調和振動子に、十分低ければ自由回転子に近づくように決めます。Q-Chem の文書と GoodVibes の実装では、重みの関数は $w(\nu) = 1/[1 + (\nu_0/\nu)^4]$ で、$\nu_0 = 100$ cm⁻¹ が既定です。

どちらの方法も、主に低いモードのエントロピーを直すためのものです。
Gaussian の `Freq` は、既定では調和振動子の値をそのまま出力します。
準調和近似を使うには、Gaussian の出力を読んで補正を計算するプログラム（たとえば GoodVibes, Luchini ら, 2020）を使うのが手軽です。
論文では、どの方法と基準の振動数を使ったかを書きます。

**ここまでで分かったこと**：調和振動子のエントロピーは振動数が 0 に近いと大きくなりすぎるので、低いモードは準調和近似で補正することがあります。使った方法を明記します。

## 5. 溶媒の効果：SCRF、PCM、SMD

### 5.1 連続誘電体モデル

溶液中の分子を計算するとき、溶媒の分子を一つずつ入れると計算が重くなります。
そこで、溶媒を、誘電率をもつ連続した媒質で置き換える方法がよく使われます。
分子を、原子の半径をもとに作った空洞（cavity）に入れ、そのまわりを媒質が取り囲むと考えます。
分子の電荷の分布が媒質を分極させ、その分極が逆に分子の電子状態を変えます。
この効果を、SCF の中で自己無撞着に取り込みます。

このような方法を連続誘電体モデルと呼び、代表が分極連続体モデル（PCM）です（Tomasi ら, 2005）。
Gaussian では、`SCRF` キーワードで指定します。

### 5.2 PCM と SMD の指定

公式の文書によると、`SCRF` の既定のモデルは IEFPCM（積分方程式の形で書いた PCM）で、既定の溶媒は水です。
溶媒は `Solvent=` で指定します。
次の 2 行は、水分子で実際に動かして文法を確かめたものです。

```text
# HF/STO-3G SCRF=(PCM,Solvent=Water)
# HF/STO-3G SCRF=(SMD,Solvent=Water)
```

SMD は、PCM の静電的な相互作用に、空洞の大きさや溶媒との分散力などによる静電的でない寄与を加えたモデルです（Marenich ら, 2009）。
静電的でない寄与の大きさは、出力の `SMD-CDS (non-electrostatic) energy` の行で確かめられます。
公式の文書は、溶媒和自由エネルギー（気相から溶液に移すときの自由エネルギーの変化）を求めるには SMD を勧めています。

### 5.3 溶液中の自由エネルギーを組み立てる

公式の文書によると、溶媒和自由エネルギー $\Delta G_\mathrm{solv}$ は、気相の計算と `SCRF=SMD` の計算をして、エネルギーの差をとって求めます。
これを使うと、溶液中の自由エネルギーは次のように組み立てられます。

$$
G^{\circ}(\text{溶液, 1 M}) = G^{\circ}(\text{気相, 1 atm}) + \Delta G_\mathrm{solv} + 1.89\ \text{kcal/mol}
$$

$G^{\circ}(\text{気相, 1 atm})$ は、気相の `Freq` の計算で得た自由エネルギーです。
最後の項は、3 節の標準状態の補正です。
SMD のパラメータは、気相も溶液も 1 M にそろえた標準状態の実験値（Minnesota Solvation Database）に合わせて決められています。
そのため、SMD の $\Delta G_\mathrm{solv}$ はこの標準状態の値で、気相の 1 atm の値と組み合わせるときに、この補正が必要になります（Bursch ら, 2022）。

構造最適化や振動解析を、初めから `SCRF` をつけて溶液中で行うこともできます。
その場合も、`Freq` の熱化学は理想気体として 1 atm で計算されるので、1 M の標準にするには同じ補正を加えます（Bursch ら, 2022）。

### 5.4 連続誘電体モデルの限界

連続誘電体モデルには、溶媒の分子が一つずつの形で現れません。
そのため、溶質と溶媒の間の水素結合や、溶媒が金属に配位する効果のような、特定の溶媒分子との相互作用は表せません。
このような効果が大切な系では、溶媒の分子を何個か露わに加え、その外側を連続誘電体で囲む方法が使われます。

**ここまでで分かったこと**：溶媒の効果は `SCRF=(SMD,Solvent=…)` などで取り入れます。溶液中の自由エネルギーは、気相の自由エネルギー、溶媒和自由エネルギー、標準状態の補正を足して組み立てます。

## まとめ

- `Freq` の出力の `Low frequencies` の最初の 6 個は並進と回転で 0 に近いはず、`Frequencies` が本当の振動数で、虚振動は負の値で出ます。
- 熱化学の出力は、電子エネルギーに ZPE、内部エネルギー、エンタルピー、自由エネルギーの補正を足したもので、既定は理想気体、298.15 K、1 atm です。
- 溶液中の反応を 1 M の標準で比べるときは、分子 1 個ごとに 1.89 kcal/mol（298.15 K）を加えます。分子の数が変わる反応で効きます。
- 低いモードの多い分子では、調和振動子のエントロピーが大きくなりすぎるので、準調和近似で補正することがあります。
- 溶媒の効果は `SCRF` で取り入れ、溶媒和自由エネルギーには SMD が勧められています。

## 付録 A　標準状態の補正の導出

理想気体 1 mol の化学ポテンシャルは、圧力 $p$ のとき次の式で表されます。

$$
\mu(p) = \mu^{\circ}(p^{\circ}) + RT \ln\frac{p}{p^{\circ}}
$$

標準状態 $p^{\circ}$ = 1 atm から、濃度 $c$ の状態に移ることを考えます。
理想気体では $p = cRT$ なので、濃度 $c^{\circ}$ = 1 mol/L の状態の圧力は $p_1 = c^{\circ} R T$ です。
これを上の式に入れると、次のようになります。

$$
\mu(c^{\circ}) - \mu^{\circ}(p^{\circ}) = RT \ln\frac{c^{\circ} R T}{p^{\circ}}
$$

298.15 K では、$RT/p^{\circ}$ = 0.0820574 L atm/(mol K) × 298.15 K ÷ 1 atm = 24.465 L/mol です。
これに $c^{\circ}$ = 1 mol/L を掛けた 24.465 は単位のない数で、1 M の気体の圧力が 24.465 atm であることを表します。
したがって、変化は次の値になります。

$$
RT \ln 24.465 = 1.987204 \times 10^{-3}\ \text{kcal/(mol K)} \times 298.15\ \text{K} \times 3.1973 = 1.894\ \text{kcal/mol}
$$

濃度を上げると、分子が動ける体積が狭くなり、並進のエントロピーが下がります。
そのため、自由エネルギーは上がる向き（正の値）になります。

## 参考文献

- A. P. Scott, L. Radom, "Harmonic Vibrational Frequencies: An Evaluation of Hartree−Fock, Møller−Plesset, Quadratic Configuration Interaction, Density Functional Theory, and Semiempirical Scale Factors", J. Phys. Chem. 100, 16502 (1996). DOI: [10.1021/jp960976r](https://doi.org/10.1021/jp960976r)
- C. J. Cramer, "Essentials of Computational Chemistry" (Wiley, 2004). ISBN: 9780470091821
- M. Bursch, J.-M. Mewes, A. Hansen, S. Grimme, "Best-Practice DFT Protocols for Basic Molecular Computational Chemistry", Angew. Chem. Int. Ed. 61, e202205735 (2022). DOI: [10.1002/anie.202205735](https://doi.org/10.1002/anie.202205735)
- R. F. Ribeiro, A. V. Marenich, C. J. Cramer, D. G. Truhlar, "Use of Solution-Phase Vibrational Frequencies in Continuum Models for the Free Energy of Solvation", J. Phys. Chem. B 115, 14556 (2011). DOI: [10.1021/jp205508z](https://doi.org/10.1021/jp205508z)
- S. Grimme, "Supramolecular Binding Thermodynamics by Dispersion-Corrected Density Functional Theory", Chem. Eur. J. 18, 9955 (2012). DOI: [10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)
- G. Luchini, J. V. Alegre-Requena, I. Funes-Ardoiz, R. S. Paton, "GoodVibes: automated thermochemistry for heterogeneous computational chemistry data", F1000Research 9, 291 (2020). DOI: [10.12688/f1000research.22758.1](https://doi.org/10.12688/f1000research.22758.1)
- J. Tomasi, B. Mennucci, R. Cammi, "Quantum Mechanical Continuum Solvation Models", Chem. Rev. 105, 2999 (2005). DOI: [10.1021/cr9904009](https://doi.org/10.1021/cr9904009)
- A. V. Marenich, C. J. Cramer, D. G. Truhlar, "Universal Solvation Model Based on Solute Electron Density and on a Continuum Model of the Solvent Defined by the Bulk Dielectric Constant and Atomic Surface Tensions", J. Phys. Chem. B 113, 6378 (2009). DOI: [10.1021/jp810292n](https://doi.org/10.1021/jp810292n)

## 関連リンク

- Gaussian 16, Freq キーワード：<https://gaussian.com/freq/>
- Gaussian, Thermochemistry in Gaussian（熱化学の計算方法の解説）：<https://gaussian.com/thermo/>
- Gaussian 16, SCRF キーワード：<https://gaussian.com/scrf/>
- Q-Chem マニュアル, 低い振動数のモードの扱い（準 RRHO）：<https://manual.q-chem.com/latest/sect_qrrho.html>
- Minnesota Solvation Database の説明書（標準状態の換算）：<https://comp.chem.umn.edu/mnsol/mnsol_Manual_v2.0.1-2020.5.20.pdf>
