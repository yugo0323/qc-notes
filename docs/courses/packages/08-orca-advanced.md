# 第 8 回　ORCA（3）：DLPNO-CCSD(T) と CASSCF/NEVPT2

!!! abstract "この回の道筋"
    ORCA でよく使われる二つの発展的な計算、DLPNO-CCSD(T) と CASSCF/NEVPT2 の入力と出力を学びます。
    前半では、DLPNO-CCSD(T) の近似の考え方を一言でつかんだうえで、補助基底、PNO の閾値、開殻の扱いを確かめます。
    後半では、CASSCF の活性空間の指定と、軌道の選び方を誤ると何が起きるかを、N₂ の実際の計算で見ます。入力例はすべて ORCA 6.1.1 で動かしたものです。

    | 節 | 内容 |
    |---|---|
    | 1 | DLPNO-CCSD(T) の考え方 |
    | 2 | DLPNO-CCSD(T) の入力と出力 |
    | 3 | PNO の閾値と三重励起の扱い |
    | 4 | 開殻の DLPNO-CCSD(T) |
    | 5 | CASSCF の入力：活性空間の指定 |
    | 6 | 活性軌道の選び方 |
    | 7 | NEVPT2 |

    **前に読んでおく回**：計算パッケージの使い方 [第 6 回](06-orca-basics.md)

!!! tip "物理的なイメージ"
    DLPNO-CCSD(T) は、「電子の避け合いは近くでしか大きく効かない」という性質を使って、CCSD(T) の計算を大幅に減らす方法です（量子化学の理論 第 11 回）。CASSCF/NEVPT2 は、電子の入り方が一つに決まらない軌道の組だけをすべての入り方で扱い（CASSCF）、そのうえで、その外への電子の跳び出しを摂動で足す方法です（第 12・13 回）。どちらも、何を近似して捨てているか（PNO の閾値、活性空間の選び方）を理解して指定することが、結果の信頼性を決めます。

## 1. DLPNO-CCSD(T) の考え方

CCSD(T) は、単一の Hartree–Fock の配置がよい出発点になる分子で、エネルギーを高い精度で求められる方法です。
ただし、計算時間は分子の大きさ $N$ のおよそ 7 乗で増えるので、ふつうの（正準の）CCSD(T) は数十原子が限界です。

DLPNO-CCSD(T)（domain-based local pair natural orbital）は、電子相関が近くの電子どうしで主に生じることを使って、この計算を軽くする方法です（Riplinger, Neese, 2013; Riplinger ら, 2013）。
考え方は、次の三つの段階にまとめられます。

1. 占有軌道を、空間的に局在した軌道に変換します。
2. 占有軌道の組（電子対）ごとに、まず安価な MP2 型の見積もりで相関エネルギーの大きさを調べます。小さい電子対は、粗い扱いで済ませます。
3. 残った電子対ごとに、その電子対の相関を表すのに適した仮想軌道（pair natural orbital、PNO）を作り、占有数の大きいものだけを使います。

つまり、「相関に効く電子対」と「その電子対に効く仮想軌道」だけを残すことで、計算量を減らします。
ORCA の文書によると、既定の設定で、正準の CCSD(T) の相関エネルギーの 99.8〜99.9% 程度が得られます。

**ここまでで分かったこと**：DLPNO-CCSD(T) は、局在化した占有軌道の電子対ごとに、相関に効く仮想軌道（PNO）だけを残して CCSD(T) を近似します。

## 2. DLPNO-CCSD(T) の入力と出力

### 2.1 入力

水分子の DLPNO-CCSD(T)/cc-pVDZ の入力は次のとおりです。

```text
! DLPNO-CCSD(T) cc-pVDZ cc-pVDZ/C TightSCF
%maxcore 1000
* xyzfile 0 1 h2o_min.xyz
```

`cc-pVDZ/C` は、相関の計算に使う RI の補助基底です（第 6 回 5.3 節）。
DLPNO の方法は RI 近似と一体になっているので、`/C` の補助基底は**必ず**指定します。
試しに `cc-pVDZ/C` を外して動かすと、次のメッセージを出して止まりました。

```text
WARNING: For DLPNO calculations we need at least one auxiliary basis set
```

def2 系の基底関数なら、`def2-TZVP def2-TZVP/C` のように、同じ名前に `/C` をつけたものを使います。
参照の HF の計算も速くしたいときは、`def2/J RIJCOSX` を加えます。

`TightSCF` は、HF の収束条件を厳しくするキーワードです。
相関の計算は HF の軌道の上に組み立てるので、ORCA の文書の例でも DLPNO の計算には `TightSCF` が付いています。

公式の文書によると、`%maxcore`（第 6 回 3 節）を大きくすると、最初の積分の変換が速くなります。必須ではありませんが、相関の計算ではメモリーを多めに与えます。

### 2.2 出力

計算の最後に、次の行が出ます。

```text
Triples Correction (T)                     ...     -0.003118474
E(CCSD)                                    ...    -76.237207643
E(CCSD(T))                                 ...    -76.240326117
```

`E(CCSD(T))` が全エネルギーで、`FINAL SINGLE POINT ENERGY` にも同じ値が出ます。
その少し前には、相関エネルギーの内訳が出ます。

```text
E(CORR)(strong-pairs)                      ...     -0.214076901
E(CORR)(weak-pairs)                        ...     -0.000007793
```

`weak-pairs` が、1 節の段階 2 で粗い扱いにした電子対の寄与です。

同じ構造で正準の CCSD(T)（`! CCSD(T) cc-pVDZ TightSCF`）も計算しました。
全エネルギーは $-76.240343$ $E_\mathrm{h}$ で、DLPNO との差は $1.7 \times 10^{-5}$ $E_\mathrm{h}$（0.011 kcal/mol）でした。
水のような小さな分子では DLPNO の利点はありませんが、近似による差が小さいことは確かめられます。

出力には `T1 diagnostic` の行も出ます（この例では 0.0065）。
これは、CCSD の一電子励起の振幅の大きさから、単一の HF 配置がよい出発点かどうかを見る目安です。
値が大きい分子は、5 節以降の多配置の方法を検討する手がかりになります。

どちらの計算でも、出力に `Frozen core treatment ... chemical core (2 el)` と出ました。
ORCA は既定で、O の 1s のような内殻の電子を相関の計算から外します（凍結内殻近似）。
内殻も相関させたいときは `NoFrozenCore` を書きます。

**ここまでで分かったこと**：DLPNO-CCSD(T) には `/C` の補助基底が必須です。全エネルギーは `E(CCSD(T))` の行で読み、既定では内殻の電子は相関させません。

## 3. PNO の閾値と三重励起の扱い

### 3.1 LoosePNO、NormalPNO、TightPNO

DLPNO の精度は、主に次の閾値で決まります。

| 閾値 | 意味 |
|---|---|
| `TCutPairs` | 段階 2 で、この値より相関エネルギーの見積もりが小さい電子対を粗く扱う |
| `TCutPNO` | 段階 3 で、占有数がこの値より大きい PNO だけを残す |
| `TCutMKN` | 各電子対の補助基底のフィッティングに使う範囲を決める技術的な値 |

ORCA は、これらをまとめて設定する三つのキーワードを用意しています（Liakos ら, 2015）。
ORCA 6.1 の文書の表から、主な値を抜き出します。

| キーワード | `TCutPairs` | `TCutPNO` | `TCutMKN` |
|---|---|---|---|
| `LoosePNO` | $10^{-3}$ | $1.00 \times 10^{-6}$ | $10^{-3}$ |
| `NormalPNO`（既定） | $10^{-4}$ | $3.33 \times 10^{-7}$ | $10^{-3}$ |
| `TightPNO` | $10^{-5}$ | $1.00 \times 10^{-7}$ | $10^{-3}$ |

2 節の計算の出力にも `TCutPNO: 3.330e-07` と `TCutPairs: 1.000e-04` が出ていて、既定が `NormalPNO` であることを確かめられます。
公式の文書は、閾値を一つずつ変えるより、これらのキーワードで精度を選ぶことを勧めています。
文書では、弱い相互作用や配座のエネルギー差のように、小さなエネルギー差を比べるときに `TightPNO` が勧められています。

```text
! DLPNO-CCSD(T) def2-TZVP def2-TZVP/C TightPNO TightSCF
```

比べるエネルギーは、すべて同じ閾値で計算します。
閾値が違うと、近似による誤差が打ち消し合わなくなるからです。

### 3.2 三重励起の (T) と (T1)

ORCA の文書によると、キーワード `DLPNO-CCSD(T)` の三重励起の補正は、計算を簡単にするための近似（半正準近似）を使ったものです。
多くの分子ではよい近似ですが、HOMO と LUMO のエネルギーの差が小さい分子などでは、まれに大きな誤差が出ます。
そこで、三重励起の振幅を繰り返し計算で求める改良版が、`DLPNO-CCSD(T1)` のキーワードで使えます（Guo ら, 2018）。
この論文では、半正準近似のほうを (T0)、改良版を (T) と呼んでいて、ORCA のキーワードとは名前が違うことに注意します。

また、ORCA 6.1 の文書によると、ORCA 6.0 からは、第一遷移金属の 3s、3p 軌道（準内殻）を含む電子対に、既定でより厳しい `TCutPNO` が使われます。
そのため、遷移金属を含む分子では、ORCA 5 と 6 で DLPNO のエネルギーが少し変わることがあります。

**ここまでで分かったこと**：精度は `LoosePNO`、`NormalPNO`（既定）、`TightPNO` で選びます。三重励起は、キーワード `DLPNO-CCSD(T)` が半正準近似、`DLPNO-CCSD(T1)` が改良版です。

## 4. 開殻の DLPNO-CCSD(T)

開殻の分子（二重項、高スピンの三重項など）にも、DLPNO-CCSD(T) を使えます（Saitow ら, 2017）。
参照の軌道は、`UHF` か `ROHF` で指定します。
次の入力は、OH ラジカル（二重項）の例です。

```text
! UHF DLPNO-CCSD(T) cc-pVDZ cc-pVDZ/C TightSCF
%maxcore 1000
* xyz 0 2
O  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.970000
*
```

UHF を参照にした場合、ORCA は UHF の軌道から準制限軌道（QRO、quasi-restricted orbitals）を作り、それを使って相関の計算をします。
QRO は、α と β の電子が同じ空間軌道を使う形に UHF の軌道を作り直したものです。
出力にも `Localizing the valence QROs` の行が出て、この手順が行われたことが分かります。
この例の全エネルギーは $-75.559257$ $E_\mathrm{h}$ でした。

ORCA の文書は、開殻の DLPNO-CCSD を、閉殻と**高スピン**の開殻のための方法としています。
開殻一重項（二つの軌道に 1 個ずつ入った電子が、一重項に組んだ状態）や、反強磁性的に結合した金属クラスターのように、スピンの向きが場所によって分かれた状態は対象外です。
文書は、SCF が対称性の破れた解（broken-symmetry 解）に収束した場合、DLPNO-CCSD の結果が物理的に意味をもつとは限らないと注意しています。
このような系には、DLPNO-NEVPT2（7 節）などのほうが適しているとされています。

**ここまでで分かったこと**：開殻では `UHF` か `ROHF` の参照を指定し、ORCA は QRO を使って相関を計算します。開殻一重項や反強磁性的な状態は対象外です。

## 5. CASSCF の入力：活性空間の指定

### 5.1 活性空間とは

結合が切れかけた分子や、ジラジカル、多くの遷移金属錯体では、一つの電子配置だけでは波動関数をよく表せません。
このような場合には、重要な軌道の中で電子の入り方をすべて考え、それらの配置を重ね合わせた波動関数を使います。

CASSCF（complete active space SCF）では、軌道を次の三つに分けます（Roos ら, 1980）。

| 軌道 | 電子の入り方 |
|---|---|
| 内側の軌道（inactive） | 常に 2 個の電子が入る |
| 活性軌道（active） | 指定した数の電子を、あらゆる入り方で入れる |
| 外側の軌道（external） | 常に空 |

活性軌道の中の電子配置をすべて重ね合わせ（完全な CI）、その係数と軌道の形を同時に最適化します。
活性空間は「CAS(電子数, 軌道数)」と書きます。

### 5.2 %casscf ブロック

ORCA では、`%casscf` ブロックに活性空間を書きます。
N₂ の三重結合の 6 個の電子を、結合性の 3 個の軌道（σ、π が 2 個）と反結合性の 3 個の軌道（π* が 2 個、σ*）に入れる CAS(6,6) の例です。

```text
! def2-SVP
%maxcore 1000
%casscf
  nel    6
  norb   6
  mult   1
  nroots 1
end
* xyz 0 1
N  0.0  0.0  0.0
N  0.0  0.0  1.0977
*
```

| キーワード | 意味 |
|---|---|
| `nel` | 活性空間の電子数 |
| `norb` | 活性軌道の数 |
| `mult` | 求める状態のスピン多重度 |
| `nroots` | 求める状態の数（基底状態だけなら 1） |

`mult` と `nroots` には、カンマで区切って複数の値を書けます。
たとえば `mult 1,3` と `nroots 4,2` で、一重項を 4 個、三重項を 2 個求めます。
複数の状態を求めると、それらの平均のエネルギーで軌道を最適化する状態平均 CASSCF（SA-CASSCF）になります。
既定では、各状態に等しい重みがつきます。

### 5.3 出力

計算が収束すると、`Final CASSCF energy` の行にエネルギーが出ます。
状態の中身は、次のような表で出ます（6 節の、正しい活性空間の計算から）。

```text
ROOT   0:  E=    -108.9888867002 Eh
      0.93462 [     0]: 222000
      0.01329 [    14]: 211110
      0.01174 [     4]: 220200
```

`222000` は、6 個の活性軌道に入る電子の数を、軌道の順に並べたものです。
左の数は、その配置の重み（CI 係数の 2 乗）です。
この例では、HF の配置（結合性の 3 個の軌道に 2 個ずつ）が 93% を占め、残りを反結合性の軌道に電子が移った配置が分け合っています。

活性軌道の自然軌道の占有数は、`ORBITAL ENERGIES` の表の `OCC` の列に出ます。
この例では、1.9826、1.9408、1.9408、0.0591、0.0591、0.0176 でした。
2 や 0 からずれた値が、静的な相関の大きさを表します。

**ここまでで分かったこと**：`%casscf` の `nel`、`norb`、`mult`、`nroots` で活性空間と求める状態を指定します。出力の配置の表と自然軌道の占有数で、波動関数の中身を確かめます。

## 6. 活性軌道の選び方

### 6.1 ORCA は軌道の番号の順に活性軌道を選ぶ

`%casscf` で決めるのは、活性軌道の**数**だけです。
**どの**軌道を活性にするかは、ORCA が軌道の番号の順に機械的に決めます。
5.2 節の N₂ の例では、出力に次のように出ました。

```text
Determined orbital ranges:
   Internal       0 -    3 (   4 orbitals)
   Active         4 -    9 (   6 orbitals)
   External      10 -   27 (  18 orbitals)
```

軌道の番号は 0 から数えます。
N₂ の電子は 14 個なので、内側に 4 個の軌道（8 電子）、その次の 6 個（4〜9 番）が活性軌道になります。
公式の文書は、望む軌道が活性になる位置に並んでいることを確かめるのは利用者の責任だと、強く注意しています。

### 6.2 実際に確かめる

そこで、N₂ の RHF/def2-SVP の軌道を、`%output Print[P_MOs] 1 end` で係数を出力して調べました。
4〜8 番は、σ 結合、π 結合 2 個、π* 2 個でした。
ところが 9 番は σ* ではなく、二つの N の s 関数が同じ符号で組み合わさった σ 型の軌道（結合の中点に対して対称な軌道）でした。
σ* にあたる、結合の中点に対して反対称な軌道は 10 番にありました。

このまま RHF の軌道を `MORead` で読み込み、CAS(6,6) を計算すると、9 回の繰り返しで収束し、警告もなく終わりました。
しかし、エネルギーは $-108.950712$ $E_\mathrm{h}$ でした。
次に、9 番と 10 番を入れ替えてから計算すると、エネルギーは $-108.988887$ $E_\mathrm{h}$ になりました。
入れ替えないと、38 m$E_\mathrm{h}$（24 kcal/mol）も高い、意図と違う活性空間の解に収束していたことになります。

入れ替えには、`%scf` ブロックの `rotate` を使います。

```text
! def2-SVP MORead
%moinp "n2_rhf.gbw"
%maxcore 1000
%scf
  rotate {9, 10, 90} end
end
%casscf
  nel    6
  norb   6
  mult   1
  nroots 1
end
* xyz 0 1
N  0.0  0.0  0.0
N  0.0  0.0  1.0977
*
```

`MORead` と `%moinp` で、前の RHF の計算の軌道（`.gbw` ファイル）を読みます。
`rotate {9, 10, 90}` は、9 番と 10 番の軌道を 90° 回転させる、つまり入れ替える指定です。
出力には `Orbital energies were swapped.` の行が出て、入れ替わったことが分かります。

間違った解を見分ける手がかりは、自然軌道の占有数でした。
入れ替えなかった計算では、6 番目の活性軌道の占有数が 0.0019 で、ほとんど空でした。
正しい計算では 0.0176 です。
占有数が 2 や 0 にごく近い活性軌道は、活性空間に入れる意味がほとんどない軌道です。
こうした軌道があったら、意図しない軌道が入っていないかを疑います。

なお、初期軌道を指定せず、ORCA の既定の推測（モデル・ポテンシャル）から始めた計算は、正しいエネルギーに収束しました。
しかし、繰り返しは 40 回かかり、途中で `Active space composition changed by more than 50%` という警告が出ました。
どの初期軌道でもうまくいくとは限らないので、活性軌道は必ず自分で確かめます。

### 6.3 軌道を選ぶための道具

軌道を確かめるには、軌道を可視化ソフトで描くのが確実です。
ほかに、ORCA の文書では次の方法が紹介されています。

- **自然軌道**：MP2 などの簡単な相関計算の自然軌道を求め、占有数が 2 や 0 から離れた軌道を活性の候補にします。
- **AVAS**：π 系の p 軌道や金属の d 軌道など、活性にしたい原子軌道を指定すると、それとよく重なる分子軌道を自動で選びます（Sayfutyarova ら, 2017）。ORCA では `! AVAS(Valence-D)` のようなキーワードで使えます。AVAS は `%casscf` の `nel` と `norb` を上書きする点に注意します。

**ここまでで分かったこと**：ORCA は軌道の番号の順に活性軌道を選ぶので、望む軌道が正しい位置にあるかを必ず確かめます。間違っていても警告なしに収束することがあり、自然軌道の占有数が手がかりになります。

## 7. NEVPT2

CASSCF は、活性空間の中の相関（静的相関）しか取り込みません。
活性空間の外の軌道との相関（動的相関）を加える方法の一つが、NEVPT2（n-electron valence state perturbation theory、2 次）です（Angeli ら, 2001）。

ORCA 6.1 では、CASSCF の入力に、次のどれかのキーワードを加えるだけで NEVPT2 を計算できます。

| キーワード | 方法 |
|---|---|
| `SC-NEVPT2` | 強縮約（strongly contracted）NEVPT2（Angeli, Cimiraglia, Malrieu, 2001） |
| `FIC-NEVPT2` | 完全内部縮約（fully internally contracted）NEVPT2。原論文の部分縮約 NEVPT2 にあたる |
| `DLPNO-NEVPT2` | FIC-NEVPT2 に DLPNO 近似を使ったもの。大きな分子向け |

公式の文書は、まず CASSCF を収束させて軌道を確かめ、そのあとで軌道を読み込んで NEVPT2 を計算する、二段階の進め方を勧めています。
次の入力は、6 節の正しい CASSCF の軌道（`n2_cas.gbw`）を読み込んで、FIC-NEVPT2 を計算するものです。

```text
! FIC-NEVPT2 def2-SVP MORead
%moinp "n2_cas.gbw"
%maxcore 1000
%casscf
  nel    6
  norb   6
  mult   1
  nroots 1
end
* xyz 0 1
N  0.0  0.0  0.0
N  0.0  0.0  1.0977
*
```

出力の最後に、次の結果が出ました。

```text
 	 Total Energy Correction : dE = -0.15460093254147
 	 Reference  Energy       : E0 = -108.98888670021809
 	 Total Energy (E0+dE)    : E  = -109.14348763275956
```

`Reference Energy` が CASSCF のエネルギー、`Total Energy Correction` が NEVPT2 の補正、その和が NEVPT2 の全エネルギーです。
補正の上には、励起のしかたごとの内訳（`Class V0_ijab` など）が出ます。
NEVPT2 でも、既定で内殻の電子は相関させません（凍結内殻近似）。

NEVPT2 の補正は CASSCF の波動関数の上に計算するので、6 節のように活性空間を誤ると、NEVPT2 の結果も意味を失います。
活性空間を確かめることが、多配置の計算で最も大切な手順です。

**ここまでで分かったこと**：NEVPT2 は、CASSCF の入力に `FIC-NEVPT2` などを加えるだけで計算できます。CASSCF を先に収束させて軌道を確かめる二段階の進め方が勧められています。

## まとめ

- DLPNO-CCSD(T) は、局在化した電子対ごとに相関に効く PNO だけを残す近似で、`/C` の補助基底が必須です。
- 精度は `LoosePNO`、`NormalPNO`（既定）、`TightPNO` で選びます。キーワード `DLPNO-CCSD(T)` の三重励起は半正準近似で、改良版は `DLPNO-CCSD(T1)` です。
- 開殻では `UHF` か `ROHF` の参照から QRO を作って計算します。開殻一重項や反強磁性的な状態は対象外です。
- CASSCF は `%casscf` の `nel`、`norb`、`mult`、`nroots` で指定します。ORCA は軌道の番号の順に活性軌道を選ぶので、`rotate` などで並びを直し、占有数で確かめます。
- NEVPT2 は `SC-NEVPT2`、`FIC-NEVPT2`、`DLPNO-NEVPT2` のキーワードで、確かめた CASSCF の軌道の上に計算します。

## 参考文献

- C. Riplinger, F. Neese, "An efficient and near linear scaling pair natural orbital based local coupled cluster method", J. Chem. Phys. 138, 034106 (2013). DOI: [10.1063/1.4773581](https://doi.org/10.1063/1.4773581)
- C. Riplinger, B. Sandhoefer, A. Hansen, F. Neese, "Natural triple excitations in local coupled cluster calculations with pair natural orbitals", J. Chem. Phys. 139, 134101 (2013). DOI: [10.1063/1.4821834](https://doi.org/10.1063/1.4821834)
- D. G. Liakos, M. Sparta, M. K. Kesharwani, J. M. L. Martin, F. Neese, "Exploring the Accuracy Limits of Local Pair Natural Orbital Coupled-Cluster Theory", J. Chem. Theory Comput. 11, 1525 (2015). DOI: [10.1021/ct501129s](https://doi.org/10.1021/ct501129s)
- Y. Guo, C. Riplinger, U. Becker, D. G. Liakos, Y. Minenkov, L. Cavallo, F. Neese, "Communication: An improved linear scaling perturbative triples correction for the domain based local pair-natural orbital based singles and doubles coupled cluster method [DLPNO-CCSD(T)]", J. Chem. Phys. 148, 011101 (2018). DOI: [10.1063/1.5011798](https://doi.org/10.1063/1.5011798)
- M. Saitow, U. Becker, C. Riplinger, E. F. Valeev, F. Neese, "A new near-linear scaling, efficient and accurate, open-shell domain-based local pair natural orbital coupled cluster singles and doubles theory", J. Chem. Phys. 146, 164105 (2017). DOI: [10.1063/1.4981521](https://doi.org/10.1063/1.4981521)
- B. O. Roos, P. R. Taylor, P. E. M. Siegbahn, "A complete active space SCF method (CASSCF) using a density matrix formulated super-CI approach", Chem. Phys. 48, 157 (1980). DOI: [10.1016/0301-0104(80)80045-0](https://doi.org/10.1016/0301-0104(80)80045-0)
- E. R. Sayfutyarova, Q. Sun, G. K.-L. Chan, G. Knizia, "Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals", J. Chem. Theory Comput. 13, 4063 (2017). DOI: [10.1021/acs.jctc.7b00128](https://doi.org/10.1021/acs.jctc.7b00128)
- C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J.-P. Malrieu, "Introduction of n-electron valence states for multireference perturbation theory", J. Chem. Phys. 114, 10252 (2001). DOI: [10.1063/1.1361246](https://doi.org/10.1063/1.1361246)
- C. Angeli, R. Cimiraglia, J.-P. Malrieu, "N-electron valence state perturbation theory: a fast implementation of the strongly contracted variant", Chem. Phys. Lett. 350, 297 (2001). DOI: [10.1016/S0009-2614(01)01303-3](https://doi.org/10.1016/S0009-2614(01)01303-3)

## 関連リンク

- ORCA 6.1 マニュアル, 結合クラスター法と DLPNO：<https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html>
- ORCA 6.1 マニュアル, CASSCF：<https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/CASSCF.html>
- ORCA 6.1 マニュアル, NEVPT2：<https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/NEVPT2.html>
- ORCA 6.1 マニュアル, 初期軌道と rotate：<https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html>
