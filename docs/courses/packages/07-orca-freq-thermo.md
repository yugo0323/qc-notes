# 第 7 回　ORCA（2）：振動解析と熱化学

!!! abstract "この回の道筋"
    ORCA で振動解析をする二つの方法（`Freq` と `NumFreq`）と、振動数の出力の読み方を学びます。
    次に、熱化学の出力を上から読み、ORCA が既定で使う準 RRHO 近似がエントロピーをどう変えるかを、実際の計算で確かめます。
    最後に、溶媒の効果を CPCM と SMD で取り入れます。入力例はすべて ORCA 6.1.1 で動かしたものです。

    | 節 | 内容 |
    |---|---|
    | 1 | Freq と NumFreq |
    | 2 | 振動数の出力を読む |
    | 3 | 熱化学の出力を読む |
    | 4 | 準 RRHO 近似：ORCA の既定 |
    | 5 | 溶媒の効果：CPCM と SMD |

    **前に読んでおく回**：計算パッケージの使い方 [第 5 回](05-gaussian-freq-thermo.md)、[第 6 回](06-orca-basics.md)

!!! tip "物理的なイメージ"
    第 5 回と同じく、振動解析は、エネルギーの曲がり具合から振動数を求め、そこから熱力学量を数えます。`Freq` は曲がり具合を解析的な微分で求め、`NumFreq` は原子を少しずつずらしたときの力の変化から数値的に求めます。準 RRHO 近似は、調和振動子の近似が悪い振動数の低いモードを、自由な回転に近いものとして扱い直し、エントロピーの過大評価を抑える工夫です。

## 1. Freq と NumFreq

振動解析では、エネルギーを原子の座標で 2 回微分した行列（ヘシアン）を求め、それを対角化して振動数を得ます（第 1 回）。
ORCA には、ヘシアンの求め方が二つあります。

| キーワード | ヘシアンの求め方 | 使える方法 |
|---|---|---|
| `Freq` | 解析的に微分する | HF と DFT（ORCA 6.1 では二重混成汎関数と RI-JK を除く） |
| `NumFreq` | 原子を少しずつ動かして勾配を計算し、その差から求める | 勾配が計算できる方法すべて |

`Freq` のほうが速く、精度も高いので、HF や DFT では `Freq` を使います。
`NumFreq` は、解析的なヘシアンがない方法で振動数が必要なときに使います。

`NumFreq` では、各原子を $x, y, z$ の正と負の向きに少しずつ動かして（既定は中心差分）、そのたびに勾配を計算します。
水分子を HF/STO-3G で計算すると、出力に次の行が出ました。

```text
Central differences            ... used
Number of displacements        ... 18 - 6
```

座標は $3 \times 3 = 9$ 個あり、それぞれ正負に動かすので 18 回です。
公式の文書によると、並進に対する不変性を使って 6 回分を省いています。
1 コアでの所要時間は、`Freq` が 1 秒ほど、`NumFreq` が 12 秒ほどでした。
得られた振動数は、`Freq` で 2169.82、4139.63、4390.66 cm⁻¹、`NumFreq` で 2169.59、4139.76、4390.70 cm⁻¹ で、差は 1 cm⁻¹ 以下でした。

`NumFreq` は、途中で止まっても、`%freq Restart true end` を書いて同じ場所で実行し直すと、続きから計算できます。
このとき、計算方法、基底関数、構造を変えてはいけません。

第 4 回の Gaussian と同じく、最適化と振動解析は一つの入力にまとめられます。

```text
! HF STO-3G Opt Freq
```

この場合、振動解析は最適化した構造で行われます。

**ここまでで分かったこと**：HF と DFT では解析的な `Freq` を使います。`NumFreq` は勾配の差からヘシアンを作る方法で、どの方法にも使えますが時間がかかります。

## 2. 振動数の出力を読む

### 2.1 VIBRATIONAL FREQUENCIES の表

振動数は、`VIBRATIONAL FREQUENCIES` の見出しの下に出ます。
水分子の HF/STO-3G の例では、次のとおりでした。

```text
     0:       0.00 cm**-1
     1:       0.00 cm**-1
     2:       0.00 cm**-1
     3:       0.00 cm**-1
     4:       0.00 cm**-1
     5:       0.00 cm**-1
     6:    2169.82 cm**-1
     7:    4139.63 cm**-1
     8:    4390.66 cm**-1
```

番号は 0 から数えます。
最初の 6 個は並進と回転にあたり、ちょうど 0.00 と出ます。
第 5 回の Gaussian では、この 6 個は 0 から少しずれた値として出ました。
ORCA は、ヘシアンから並進と回転の成分を取り除いてから対角化するので、ちょうど 0 になります。
そのため公式の文書が注意しているとおり、ORCA ではこの 6 個を見ても、ヘシアンの数値の精度は分かりません。
精度を確かめたいときは、最適化の収束条件や DFT の積分の格子を変えて計算し、振動数がどれだけ動くかを比べます。

### 2.2 虚振動

虚振動は、負の値に印をつけて出力されます。
たとえば、次の行は、水の二量体を B3LYP-D3(BJ)/def2-SVP で、既定の条件で最適化したあとの振動解析で実際に出たものです。

```text
     6:     -20.05 cm**-1  ***imaginary mode***
```

第 5 回と同じく、極小を求めたつもりで虚振動が出たら、その振動の向きに構造を少しずらして最適化をやり直します。

### 2.3 振動数のスケーリング

調和近似の振動数に係数を掛けて実験と比べるときは（第 5 回）、`%freq ScalFreq 0.96 end` のように係数を指定します。
出力の表の上には `Scaling factor for frequencies =  1.000000000  (already applied!)` の行が出て、どの係数が掛かっているかを確かめられます。
係数を掛けた振動数は、熱化学の計算にもそのまま使われます。
実際に、最適化した水分子の構造（`.xyz` ファイル）で `ScalFreq 0.90` とすると、振動数は 1952.95 cm⁻¹（係数なしの 2169.94 cm⁻¹ の 0.90 倍）などになりました。
零点エネルギーも、0.024378 から 0.021940 $E_\mathrm{h}$ に下がりました。

**ここまでで分かったこと**：ORCA の振動数の表では、並進と回転の 6 個がちょうど 0.00 と出ます。虚振動は負の値と `***imaginary mode***` の印で表されます。

## 3. 熱化学の出力を読む

### 3.1 何を計算しているか

振動数の表のあとに、`THERMOCHEMISTRY AT 298.15K` の見出しで熱化学の結果が出ます。
第 5 回と同じく、分子を理想気体とみなし、並進、回転、振動、電子の寄与を統計力学の式で計算します。
見出しのすぐ下に、条件が出ます。

```text
Temperature         ...   298.15 K
Pressure            ...     1.00 atm
Total Mass          ...    18.02 AMU
Quasi RRHO          ...     True
Cut-Off Frequency   ...     1.00 cm^-1
```

既定は 298.15 K、1 atm です。
`Quasi RRHO ... True` は、振動のエントロピーに準 RRHO 近似を使っていることを表します（4 節）。

### 3.2 内部エネルギー、エンタルピー、自由エネルギー

出力は、内部エネルギー、エンタルピー、エントロピー、ギブズ自由エネルギーの順に、ブロックに分かれて出ます。
水分子の HF/STO-3G の例から、主な行を抜き出します（単位はハートリー）。

```text
Electronic energy                ...    -74.96590119 Eh
Zero point energy                ...      0.02437664 Eh      15.30 kcal/mol
Thermal vibrational correction   ...      0.00000028 Eh       0.00 kcal/mol
Thermal rotational correction    ...      0.00141627 Eh       0.89 kcal/mol
Thermal translational correction ...      0.00141627 Eh       0.89 kcal/mol
Total thermal energy                    -74.93869173 Eh
...
Thermal Enthalpy correction       ...      0.00094421 Eh       0.59 kcal/mol
Total Enthalpy                    ...    -74.93774752 Eh
...
Electronic entropy                ...      0.00000000 Eh      0.00 kcal/mol
Vibrational entropy               ...      0.00000031 Eh      0.00 kcal/mol
Rotational entropy                ...      0.00507112 Eh      3.18 kcal/mol
Translational entropy             ...      0.01644380 Eh     10.32 kcal/mol
Final entropy term                ...      0.02151523 Eh     13.50 kcal/mol
...
Final Gibbs free energy         ...    -74.95926275 Eh
G-E(el)                           ...      0.00663844 Eh      4.17 kcal/mol
```

上から順に読みます。

- **Electronic energy**：電子エネルギー $E_\mathrm{el}$ です。最後の一点計算の `FINAL SINGLE POINT ENERGY` と同じ値です。
- **Zero point energy** から **Thermal translational correction**：零点エネルギーと、温度による振動、回転、並進のエネルギーです。これらを $E_\mathrm{el}$ に足したものが `Total thermal energy`（内部エネルギー $U$）です。
- **Thermal Enthalpy correction**：理想気体の $H = U + k_\mathrm{B}T$ の $k_\mathrm{B}T$ です（第 5 回）。`Total Enthalpy` が $H$ です。
- **エントロピーの 4 行**：$T S$ の値で、エネルギーの単位で出ます。単位がエントロピー（cal/(mol K)）ではなく、温度を掛けたエネルギーである点に注意します。
- **Final Gibbs free energy**：$G = H - TS$ です。
- **G-E(el)**：$G - E_\mathrm{el}$ です。Gaussian の `Thermal correction to Gibbs Free Energy` にあたります。

第 5 回で見た「軽いレベルの熱補正と重いレベルの電子エネルギーを足す」使い方では、この `G-E(el)` の値を使います。
自由エネルギーは 1 atm の標準状態の値なので、溶液の 1 M の標準にそろえるときは、第 5 回の 1.89 kcal/mol（298.15 K）の補正を加えます。

### 3.3 回転の対称数と温度の変更

エントロピーのブロックには、次の行も出ます。

```text
Point Group:  C2v, Symmetry Number:   2
```

回転のエントロピーは、分子の対称数（回転させて重なる向きの数）で変わります。
ORCA は点群から対称数を決め、さらに対称数を 1 から 12 にした場合の回転のエントロピーを表にして出力します。
構造がわずかに歪んで、点群が正しく判定されなかったときは、この表から正しい対称数の値を選べます。

温度を変えるときは、`%freq` ブロックの `Temp` に書きます。
カンマで区切ると、複数の温度で計算できます。

```text
%freq
  Temp 298.15, 373.15
end
```

この入力では、`THERMOCHEMISTRY AT 298.15K` と `THERMOCHEMISTRY AT 373.15K` の二つのブロックが出ました。

**ここまでで分かったこと**：熱化学の出力は U、H、TS、G の順に出ます。エントロピーは $TS$ の形で出ること、`G-E(el)` が自由エネルギーの補正であることに注意します。

## 4. 準 RRHO 近似：ORCA の既定

### 4.1 ORCA は既定で準 RRHO を使う

第 5 回で見たとおり、調和振動子のエントロピーは、振動数が 0 に近いと大きくなりすぎます。
Grimme（2012）は、振動数が基準の値 $\omega_0$ より十分低いモードを自由回転子のように、十分高いモードを調和振動子のように扱い、その間を滑らかにつなぐ方法を提案しました。
これを準 RRHO（quasi-RRHO）近似と呼びます。

ORCA 6.1 の文書によると、ORCA の振動のエントロピーは、**既定で**この準 RRHO で計算されます。
基準の振動数 $\omega_0$ の既定は、原論文と同じ 100 cm⁻¹ で、`%freq` の `QRRHORefFreq` で変えられます。
出力にも、次の行で基準の振動数が示されます。

```text
Vibrational entropy computed according to the QRRHO of S. Grimme
Chem.Eur.J. 2012 18 9955 using a reference frequency of 100.0 cm-1
```

ORCA の文書は、原論文との違いも一つ挙げています。
自由回転子の式に入れる平均の慣性モーメントを、原論文は分子によらない一定の値にしましたが、ORCA は計算している分子の慣性モーメントの平均を使います。
文書によると、これによる自由エネルギーの差は、多くの場合 0.1 kcal/mol 以内です。

### 4.2 版による違い

手元にある ORCA の三つの版で、同じ水分子の振動解析をして、出力を比べました。

| 版 | 準 RRHO の既定 | 出力の表示 |
|---|---|---|
| 4.2.1 | 使う | `Vibrational entropy computed according to the QRRHO of S. Grimme` の行 |
| 5.0.4 | 使う | 4.2.1 と同じ |
| 6.1.1 | 使う | 上の行に加えて、`Quasi RRHO ... True` と基準の振動数を表示 |

三つの版のどれでも、`%freq QuasiRRHO false end` で準 RRHO を止められました。
`QRRHORefFreq` は、ORCA 6.0 と 6.1 の文書に記載があります。
古い版の結果と比べるときは、出力のこれらの行で、どの扱いだったかを確かめます。

ここで大切なのは、**Gaussian の `Freq` は既定で調和振動子の値を出し、ORCA は既定で準 RRHO の値を出す**という違いです（第 5 回）。
低い振動数のモードが多い分子では、同じ構造、同じ振動数でも、二つのプログラムの自由エネルギーが違ってきます。
論文には、どちらの扱いを使ったかを書きます。

### 4.3 低いモードのある分子で比べる

準 RRHO の効果を、メチル基がほぼ自由に回るトルエンで確かめます。
HF/STO-3G で最適化と振動解析をすると、最も低い振動数（メチル基の回転）は 6.35 cm⁻¹ でした。
最適化した構造（`.xyz` ファイル）で振動解析だけをやり直すと、このモードは 30.50 cm⁻¹ になりました。
二つの計算の電子エネルギーは $10^{-8}$ $E_\mathrm{h}$ しか違いません。
二つの計算では、SCF の収束条件が違いました。
出力の `Energy Change TolE` の行を見ると、最適化と続けた計算では、最適化のために厳しくした $10^{-8}$ $E_\mathrm{h}$ が使われていました（第 6 回 3 節）。
振動解析だけの計算では、既定の $10^{-6}$ $E_\mathrm{h}$ のままでした。
実際、振動解析だけの計算に `VeryTightSCF` を加えると、このモードは 12.92 cm⁻¹ になりました。
このことから、食い違いの大きな原因は SCF の収束の違いだと推定できます。
ただし、12.92 cm⁻¹ は 6.35 cm⁻¹ とも一致しないので、原因をすべて切り分けたわけではありません。
いずれにしても、このように平らなモードの振動数は、計算条件のわずかな違いで大きく動きます。

それぞれの計算について、`QuasiRRHO false`（調和振動子）と既定（準 RRHO）で、振動のエントロピーの項 $T S_\mathrm{vib}$ を比べました（298.15 K）。

| 最も低い振動数 | 調和振動子 | 準 RRHO |
|---|---|---|
| 6.35 cm⁻¹ | 4.40 kcal/mol | 3.41 kcal/mol |
| 30.50 cm⁻¹ | 3.48 kcal/mol | 2.96 kcal/mol |
| 二つの差 | 0.92 kcal/mol | 0.45 kcal/mol |

準 RRHO では、$T S_\mathrm{vib}$ が小さくなり、低いモードの振動数のぶれによる差も半分ほどに減りました。
ただし、差が消えるわけではありません。
低いモードの多い分子の自由エネルギーには、このくらいの不確かさがあると考えておきます。

**ここまでで分かったこと**：ORCA は 4.2.1、5.0.4、6.1.1 のどれでも既定で準 RRHO（$\omega_0$ = 100 cm⁻¹）を使います。Gaussian の既定とは違うので、使った扱いを記録します。

## 5. 溶媒の効果：CPCM と SMD

### 5.1 指定のしかた

第 5 回の連続誘電体モデルは、ORCA では「!」行に書きます。
次の 2 行は、水分子で実際に動かして確かめたものです。

```text
! B3LYP D3BJ def2-SVP CPCM(Water)
! B3LYP D3BJ def2-SVP SMD(Water)
```

`CPCM(溶媒名)` は導体様分極連続体モデル（C-PCM）です（Barone, Cossi, 1998）。
導体様とは、まず溶媒を導体（誘電率が無限大）とみなして表面の電荷を求め、それを実際の誘電率に合わせて縮める考え方です。
`SMD(溶媒名)` は、C-PCM に、空洞の形成や分散力などの静電的でない寄与を加えたモデルです（Marenich ら, 2009）。
ORCA の C-PCM の実装は、Garcia-Ratés と Neese（2020）に書かれています。

### 5.2 出力を読む

`CPCM(Water)` の出力では、まず溶媒の設定が出ます。

```text
CPCM parameters:
  Epsilon                                         ...      80.1510
...
Solvent:                                          ... WATER
```

SCF が収束したあと、溶媒との静電的な相互作用のエネルギーが `CPCM Dielectric` の行に出ます。
この例では $-0.010658$ $E_\mathrm{h}$ でした。
`Free-energy (cav+disp)` の行には `This term is not implemented in the current solvation scheme` と出ました。
C-PCM だけでは、静電的でない寄与が入っていないということです。

`SMD(Water)` の出力では、静電的でない寄与が加わります。

```text
SMD CDS free energy correction energy :                 1.45916     Kcal/mol
```

CDS は、空洞（cavity）、分散（dispersion）、溶媒の構造（solvent structure）の頭文字です。
最後の `FINAL SINGLE POINT ENERGY` には、この寄与も含まれています。

### 5.3 溶媒和自由エネルギーを求める

溶媒和自由エネルギー $\Delta G_\mathrm{solv}$ は、同じ構造で気相と SMD の一点計算をして、`FINAL SINGLE POINT ENERGY` の差をとって求めます。
B3LYP-D3(BJ)/def2-SVP の水分子では、次のとおりでした（構造は、気相の HF/STO-3G で最適化したもの）。

| 計算 | FINAL SINGLE POINT ENERGY |
|---|---|
| 気相 | $-76.320892$ $E_\mathrm{h}$ |
| SMD(Water) | $-76.332865$ $E_\mathrm{h}$ |
| 差 $\Delta G_\mathrm{solv}$ | $-0.011973$ $E_\mathrm{h}$ = $-7.51$ kcal/mol |

ORCA 6.1 の文書も注意しているとおり、SMD は気相と溶液の両方で濃度を 1 mol/L にそろえた標準状態で作られています。
したがって、第 5 回と同じく、溶液中の 1 M の自由エネルギーは次のように組み立てます。

$$
G^{\circ}(\text{溶液, 1 M}) = G^{\circ}(\text{気相, 1 atm}) + \Delta G_\mathrm{solv} + 1.89\ \text{kcal/mol}
$$

$G^{\circ}(\text{気相, 1 atm})$ は、3 節の `Final Gibbs free energy` です。
この式の考え方はプログラムによらないので、Gaussian で求めた場合と同じです。

**ここまでで分かったこと**：溶媒は「!」行の `CPCM(溶媒名)` や `SMD(溶媒名)` で指定します。C-PCM には静電的でない寄与が入らないので、溶媒和自由エネルギーには SMD を使います。

## まとめ

- HF と DFT の振動解析は解析的な `Freq` で、それ以外の方法は数値微分の `NumFreq` で行います。
- ORCA の振動数の表では並進と回転の 6 個がちょうど 0.00 と出るので、これで精度は判断できません。虚振動は `***imaginary mode***` の印で分かります。
- 熱化学の出力は U、H、TS、G の順で、エントロピーは $TS$ の形で出ます。`G-E(el)` が自由エネルギーの補正です。
- ORCA は 4.2.1、5.0.4、6.1.1 のどれでも既定で準 RRHO を使います。Gaussian の既定（調和振動子）と違うので、使った扱いを記録します。
- 溶媒は `CPCM(溶媒名)`、`SMD(溶媒名)` で指定し、溶液中の 1 M の自由エネルギーは、気相の自由エネルギー、溶媒和自由エネルギー、1.89 kcal/mol の和で組み立てます。

## 参考文献

- S. Grimme, "Supramolecular Binding Thermodynamics by Dispersion-Corrected Density Functional Theory", Chem. Eur. J. 18, 9955 (2012). DOI: [10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)
- V. Barone, M. Cossi, "Quantum Calculation of Molecular Energies and Energy Gradients in Solution by a Conductor Solvent Model", J. Phys. Chem. A 102, 1995 (1998). DOI: [10.1021/jp9716997](https://doi.org/10.1021/jp9716997)
- A. V. Marenich, C. J. Cramer, D. G. Truhlar, "Universal Solvation Model Based on Solute Electron Density and on a Continuum Model of the Solvent Defined by the Bulk Dielectric Constant and Atomic Surface Tensions", J. Phys. Chem. B 113, 6378 (2009). DOI: [10.1021/jp810292n](https://doi.org/10.1021/jp810292n)
- M. Garcia-Ratés, F. Neese, "Effect of the Solute Cavity on the Solvation Energy and its Derivatives within the Framework of the Gaussian Charge Scheme", J. Comput. Chem. 41, 922 (2020). DOI: [10.1002/jcc.26139](https://doi.org/10.1002/jcc.26139)
- F. Neese, "Software Update: The ORCA Program System—Version 6.0", WIREs Comput. Mol. Sci. 15, e70019 (2025). DOI: [10.1002/wcms.70019](https://doi.org/10.1002/wcms.70019)

## 関連リンク

- ORCA 6.1 マニュアル, 振動数の計算：<https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html>
- ORCA 6.1 マニュアル, 熱化学（準 RRHO）：<https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html>
- ORCA 6.1 マニュアル, 連続誘電体モデル（CPCM、SMD）：<https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html>
