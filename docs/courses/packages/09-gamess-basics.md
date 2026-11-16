# 第 9 回　GAMESS（1）：入力ファイルの書き方

!!! abstract "この回の道筋"
    GAMESS の入力ファイルは、`$CONTRL` や `$DATA` などの「グループ」を並べて作ります。
    この回では、まず書式の決まり（2 列目の `$`、`$END`）を学び、次に主なグループの役割と書き方を一つずつ見ます。
    最後に、HF、UHF、ROHF、DFT の指定と、構造最適化の入力と出力の読み方をまとめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 入力ファイルの全体の形 |
    | 2 | 書式の決まり |
    | 3 | $CONTRL：計算の種類を決める |
    | 4 | $SYSTEM と $BASIS |
    | 5 | $DATA：分子の構造と対称性 |
    | 6 | $SCF と DFT の指定 |
    | 7 | 構造最適化：RUNTYP=OPTIMIZE と $STATPT |

    **前に読んでおく回**：計算パッケージの使い方 [第 3 回](03-common-concepts.md)、[第 6 回](06-orca-basics.md)

## 1. 入力ファイルの全体の形

GAMESS は、アイオワ州立大学の Gordon らのグループを中心に開発されている量子化学計算プログラムです（Schmidt ら, 1993; Barca ら, 2020）。
入力は、拡張子 `.inp` のテキストファイルに書きます。
まず、水分子を RHF/6-31G(d) で構造最適化する入力の全体を示します。

```text
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE ICHARG=0 MULT=1 $END
 $SYSTEM MWORDS=50 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $STATPT OPTTOL=0.0001 NSTEP=50 $END
 $DATA
H2O RHF/6-31G(d) optimization
C1
O   8.0   0.000000   0.000000   0.117300
H   1.0   0.000000   0.757200  -0.469200
H   1.0   0.000000  -0.757200  -0.469200
 $END
```

`$` で始まる名前のまとまりを**グループ**と呼びます。
各グループの役割は次のとおりです。

| グループ | 役割 |
|---|---|
| `$CONTRL` | 波動関数の種類、計算の種類、電荷、多重度など、計算の中心になる設定 |
| `$SYSTEM` | メモリーなど、計算機に関する設定 |
| `$BASIS` | 基底関数 |
| `$STATPT` | 構造最適化の設定 |
| `$DATA` | 題名、対称性、原子の種類と座標 |

グループはどの順に並べてもかまいません。
省略したグループの設定には、既定の値が使われます。
`$DATA` だけは必ず書きます。

実行には、GAMESS をインストールしたときに用意される `rungms` というスクリプトを使います。
引数の与え方は計算機ごとの設定で変わることがあるので、管理者の説明に従います。

**ここまでで分かったこと**：GAMESS の入力は、`$名前` で始まるグループを並べたものです。`$DATA` は必須で、ほかのグループは省略すると既定の値になります。

## 2. 書式の決まり

GAMESS の入力には、ほかのプログラムにない厳しい決まりがあります。
INPUT.DOC と INTRO.DOC（GAMESS の公式の文書）から、要点をまとめます。

1. **`$` は必ず 2 列目に置きます。** 1 列目は空白にします。1 列目や 3 列目に `$` があると、そのグループは読まれません。
2. **グループは `$END` で終わります。** 1 節の例のように、`$END` は同じ行の終わりに書いても、別の行に書いてもかまいません。
3. **キーワードは「名前=値」の形で、空白かカンマで区切ります。** キーワードを省略形で書くことはできません。
4. **大文字と小文字は区別しません。**
5. **論理値は `.TRUE.` か `.FALSE.`（`.T.`、`.F.`）で書きます。** 前後のピリオドが必要です。
6. **グループの外の行は注釈として無視されます。** グループの中では、`!` で注釈が始まり、次の `!` で終わります。

決まり 1 を使うと、行の先頭に文字を加えるだけで、グループを一時的に無効にできます。

```text
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE $END
--$CONTRL EXETYP=CHECK $END
```

2 行目は 2 列目が `$` ではないので、読まれません。
ちなみに `EXETYP=CHECK` は、エネルギーを計算せずに入力とメモリーの量だけを確かめる設定です。
大きな計算の前に、入力の誤りを見つけるのに使えます。

`$DATA` だけは、キーワードを使わず、行の位置で意味が決まる形式です（5 節）。

**ここまでで分かったこと**：`$` は 2 列目、グループは `$END` で終わります。グループの中は「名前=値」を空白で区切って並べます。

## 3. $CONTRL：計算の種類を決める

`$CONTRL` は、計算の中心になる設定を書くグループです。
綴りは CONTROL ではなく、O が一つの `CONTRL` です（GAMESS のグループ名とキーワードは 6 文字以内です）。
よく使うキーワードは次のとおりです。

| キーワード | 意味 | 既定 |
|---|---|---|
| `SCFTYP` | 波動関数の種類：`RHF`、`UHF`、`ROHF`、`MCSCF` など | `RHF` |
| `RUNTYP` | 計算の種類：`ENERGY`（一点）、`GRADIENT`、`HESSIAN`（振動解析）、`OPTIMIZE`、`SADPOINT`（遷移状態）など | `ENERGY` |
| `ICHARG` | 分子全体の電荷 | 0 |
| `MULT` | スピン多重度 | 1 |
| `DFTTYP` | DFT の汎関数（6 節） | `NONE`（DFT を使わない） |
| `MAXIT` | SCF の繰り返しの上限 | 30 |
| `ISPHER` | d 関数などを球面調和関数で使うなら 1 | −1（直交座標の関数） |

第 3 回で学んだ RHF、UHF、ROHF の区別は、`SCFTYP` で指定します。
INPUT.DOC では、`SCFTYP` の既定は単に `RHF` とされていて、Gaussian や ORCA のように多重度に応じて UHF に切り替わるという記述はありません。
開殻の分子では、`SCFTYP` を必ず書きます。
たとえば、OH ラジカルを UHF と ROHF で計算するときの `$CONTRL` は、次のようになります。

```text
 $CONTRL SCFTYP=UHF  MULT=2 $END
 $CONTRL SCFTYP=ROHF MULT=2 $END
```

SCF の繰り返しの上限 `MAXIT` は、`$SCF` ではなく `$CONTRL` に書く点にも注意します。
既定の 30 回は、ほかのプログラムより少なめです。

**ここまでで分かったこと**：`$CONTRL` の `SCFTYP` で RHF、UHF、ROHF を、`RUNTYP` で計算の種類を選びます。既定は RHF なので、開殻では `SCFTYP` を必ず書きます。

## 4. $SYSTEM と $BASIS

### 4.1 $SYSTEM：メモリー

`$SYSTEM` の `MWORDS` は、1 コアが使ってよいメモリーの量です。
単位は「100 万語」で、1 語は 64 ビット（8 バイト）です。
そのため、`MWORDS=50` は $50 \times 10^6 \times 8$ バイト = 400 MB にあたります。
既定は `MWORDS=1`（8 MB）なので、小さな分子でも足りなくなることがあり、たいてい明示します。

MP2 や結合クラスター法などを並列で計算するときは、`MEMDDI` で、全体で共有するメモリーの量も指定します。
この場合、1 コアあたりに必要な量は、おおよそ `MEMDDI`/コア数 + `MWORDS` です。

### 4.2 $BASIS：基底関数

`$BASIS` では、基底関数を系列の名前と追加の関数の数で組み立てます。

| 書き方 | 基底関数 |
|---|---|
| `GBASIS=STO NGAUSS=3` | STO-3G |
| `GBASIS=N31 NGAUSS=6` | 6-31G |
| `GBASIS=N31 NGAUSS=6 NDFUNC=1` | 6-31G(d) |
| `GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1` | 6-31G(d,p) |
| `GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.TRUE.` | 6-31+G(d) |
| `GBASIS=CCD` | cc-pVDZ（`CCT` で cc-pVTZ） |
| `GBASIS=ACCD` | aug-cc-pVDZ |

`NDFUNC` は重い原子（この場合は Li 以降）の d 関数の数、`NPFUNC` は水素とヘリウムの p 関数の数、`DIFFSP` は重い原子の分散関数です。

cc-pVDZ などの Dunning 系の基底関数は、球面調和関数（d 関数 5 個）で使うように作られています。
INPUT.DOC によると、これらは `$CONTRL` に `ISPHER=1` を書かないと計算が始まりません。
一方、6-31G(d) は、既定の `ISPHER=-1`（d 関数 6 個）のままで使うのがふつうです。
第 3 回で見たとおり、d 関数の数でエネルギーが少し変わるので、ほかのプログラムと比べるときはそろえます。

**ここまでで分かったこと**：`MWORDS` は 1 コアあたりのメモリーで、単位は 800 万バイトです。基底関数は `GBASIS` と追加の関数の数で組み立て、Dunning 系には `ISPHER=1` が必要です。

## 5. $DATA：分子の構造と対称性

### 5.1 C1 で書く

`$DATA` は、行の位置で意味が決まるグループです。
1 節の例をもう一度示します。

```text
 $DATA
H2O RHF/6-31G(d) optimization
C1
O   8.0   0.000000   0.000000   0.117300
H   1.0   0.000000   0.757200   -0.469200
H   1.0   0.000000  -0.757200   -0.469200
 $END
```

| 行 | 内容 |
|---|---|
| 1 行目 | `$DATA`（2 列目から） |
| 2 行目 | 題名（1 行） |
| 3 行目 | 点群の記号 |
| 4 行目以降 | 原子の名前、核電荷、$x, y, z$（Å） |
| 最後 | `$END` |

原子の名前は出力に表示されるだけで、原子の種類は**核電荷**の数（O なら 8.0）で決まります。
座標の単位は Å が既定です（`$CONTRL` の `UNITS=BOHR` で原子単位に変えられます）。

点群を `C1`（対称性なし）にすると、すべての原子を書きます。
C1 のときは、点群の行のすぐ次から原子を書き、**空行を入れません**。

### 5.2 対称性を使って書く

GAMESS は分子の対称性を使って計算を速くします。
対称性を使うときは、点群の行に記号と主軸の次数を書き、**対称性で区別される原子（対称独立な原子）だけ**を書きます。
残りの原子は、GAMESS が対称操作で作ります。
水分子（C₂ᵥ）なら、次のようになります。

```text
 $DATA
H2O RHF/6-31G(d) optimization, C2v
Cnv 2

O   8.0   0.000000   0.000000   0.117300
H   1.0   0.757200   0.000000  -0.469200
 $END
```

`Cnv 2` は「$n = 2$ の Cnv」、つまり C₂ᵥ を表します。
点群の行の次の空行は、分子を置く座標系として標準のもの（INPUT.DOC が master frame と呼ぶもの）を使う、という指定です。
C1 以外では、この空行が必要です。
標準の座標系では、主軸が $z$ 軸、C₂ᵥ の鏡映面の一つが $xz$ 面です。
この例では、水分子が $xz$ 面にくるように、H の座標を $x$ 方向に書きました。
GAMESS が対称操作でどの原子を作ったかは、`EXETYP=CHECK` の計算の出力で、エネルギーを計算する前に確かめられます。

対称性の指定は間違えやすいので、慣れないうちは `C1` で書くのが安全です。
`C1` ならすべての原子を書くので、xyz ファイルの座標に核電荷の列を加えるだけで使えます。
ただし、対称性を使わない分、計算は遅くなります。

**ここまでで分かったこと**：`$DATA` は「題名、点群、原子」の順に書きます。C1 では空行を入れずにすべての原子を書き、C1 以外では空行のあとに対称独立な原子だけを書きます。

## 6. $SCF と DFT の指定

### 6.1 $SCF

`$SCF` は、SCF の収束を細かく制御するグループです。
よく使うキーワードは次のとおりです。

| キーワード | 意味 |
|---|---|
| `DIRSCF=.TRUE.` | 積分をディスクに保存せず、毎回計算し直す（直接 SCF） |
| `DIIS`、`SOSCF` | 収束を速める方法（どちらか一方） |
| `DAMP=.TRUE.` | ダンピング（第 3 回） |
| `SHIFT=.TRUE.` | レベルシフト（第 3 回） |
| `CONV` | 密度行列の変化の収束の閾値 |

INPUT.DOC によると、収束を速める方法の既定は、RHF、UHF、ROHF では 2 次の収束法 SOSCF、DFT では DIIS です。
`CONV` の既定は $10^{-5}$ で、MP2 や TDDFT などの計算を続けるときは $10^{-6}$ になります。
既定では積分をディスクに保存するので、大きな分子では `DIRSCF=.TRUE.` を使います。

### 6.2 DFT の指定

GAMESS では、DFT は「SCF の種類」ではなく、`SCFTYP` に `DFTTYP` を加えて指定します。

```text
 $CONTRL SCFTYP=RHF DFTTYP=B3LYP $END
```

これで、閉殻の B3LYP（RKS）になります。
`SCFTYP=UHF` や `SCFTYP=ROHF` と組み合わせれば、非制限や制限開殻の DFT になります。

第 3 回で、B3LYP の中身がプログラムによって少し違うと述べました。
GAMESS の `DFTTYP=B3LYP` は、局所密度近似の相関の部分に VWN の式 5（VWN5）を使います。
INPUT.DOC によると、VWN の別の式を使うプログラムと結果を合わせたいときは、`DFTTYP=B3LYPV1R` を使います。

分散補正は、`$DFT` グループで加えます。

```text
 $DFT DC=.TRUE. IDCVER=4 $END
```

`DC=.TRUE.` で Grimme の分散補正を加え、`IDCVER=4` で D3(BJ) を選びます（`IDCVER=3` が D3）。
INPUT.DOC によると、`IDCVER` を書けば `DC=.TRUE.` も自動で設定されます。

**ここまでで分かったこと**：SCF の細かな制御は `$SCF` で行います。DFT は `SCFTYP` と `DFTTYP` の組み合わせで指定し、分散補正は `$DFT` の `DC` と `IDCVER` で加えます。

## 7. 構造最適化：RUNTYP=OPTIMIZE と $STATPT

### 7.1 入力

構造最適化は、`$CONTRL` に `RUNTYP=OPTIMIZE` を書いて行います。
細かな設定は `$STATPT` に書きます。

| キーワード | 意味 | 既定 |
|---|---|---|
| `OPTTOL` | 勾配の収束の閾値（hartree/bohr） | 0.0001 |
| `NSTEP` | 最適化の繰り返しの上限 | 50（遷移状態の探索では 20） |
| `METHOD` | 最適化の方法 | `QA` |
| `HSSEND=.TRUE.` | 最適化が終わったあとに、ヘシアンを計算して振動解析をする | `.FALSE.` |

INPUT.DOC によると、勾配の最大の成分が `OPTTOL` より小さく、かつ勾配の 2 乗平均の平方根が `OPTTOL` の 1/3 より小さくなると、収束と判定されます。
Gaussian（第 4 回）や ORCA（第 6 回）と違い、判定に使うのは勾配だけです。

`HSSEND=.TRUE.` を書くと、Gaussian の `Opt Freq` のように、最適化した構造で振動解析まで続けて行えます。
振動解析だけをするときは、`RUNTYP=HESSIAN` を使います。

B3LYP-D3(BJ)/6-31G(d) で最適化し、続けて振動解析をする入力は、たとえば次のようになります。

```text
 $CONTRL SCFTYP=RHF DFTTYP=B3LYP RUNTYP=OPTIMIZE ICHARG=0 MULT=1 $END
 $SYSTEM MWORDS=50 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $DFT    DC=.TRUE. IDCVER=4 $END
 $STATPT OPTTOL=0.0001 NSTEP=50 HSSEND=.TRUE. $END
 $DATA
H2O B3LYP-D3(BJ)/6-31G(d) optimization and frequencies
C1
O   8.0   0.000000   0.000000   0.117300
H   1.0   0.000000   0.757200  -0.469200
H   1.0   0.000000  -0.757200  -0.469200
 $END
```

### 7.2 出力を読む

GAMESS の出力では、最適化の各段階に `NSERCH` という番号がつきます（0 から数えます）。
例として、アイオワ州立大学の GAMESS の講習資料に載っている、アンモニアの最適化の出力から抜き出します。
これは 1998 年の版の GAMESS の出力なので、現在の版とは行の形が違うところがあるかもしれません。ここでは、どの語を探せばよいかの例として読んでください。

```text
           NSERCH=   0           ENERGY=       -55.3983197
           MAXIMUM GRADIENT =    .1098984     RMS GRADIENT = .0550384
...
             NSERCH=   5    ENERGY=     -55.4554198
             MAXIMUM GRADIENT =   .0000053    RMS GRADIENT =   .0000028
1        ***** EQUILIBRIUM GEOMETRY LOCATED *****
```

各段階のエネルギーと、勾配の最大の成分（`MAXIMUM GRADIENT`）、2 乗平均の平方根（`RMS GRADIENT`）が出ます。
7.1 節の判定の条件を満たすと、`EQUILIBRIUM GEOMETRY LOCATED` の行が出て、そのあとに最適化した構造の座標が続きます。
この例では、`NSERCH=5` で、勾配の最大の成分が 0.0000053 と `OPTTOL` の 0.0001 より十分小さくなっています。

出力の最後に `EXECUTION OF GAMESS TERMINATED NORMALLY` があれば、計算は正常に終わっています。
行の形は版によって違うことがあるので、`NSERCH`、`EQUILIBRIUM GEOMETRY LOCATED`、`TERMINATED NORMALLY` などの語で検索して読みます。

`NSTEP` の上限に達して止まった場合は、最後の構造から、入力を作り直して続けます。

**ここまでで分かったこと**：最適化は `RUNTYP=OPTIMIZE` と `$STATPT` で行い、勾配だけで収束を判定します。出力では `NSERCH` ごとの勾配と、`EQUILIBRIUM GEOMETRY LOCATED` の行を確かめます。

## まとめ

- GAMESS の入力は、2 列目の `$` で始まり `$END` で終わるグループを並べたものです。`$DATA` は必須です。
- `$CONTRL` の `SCFTYP` で RHF、UHF、ROHF を、`RUNTYP` で計算の種類を選びます。既定は RHF なので、開殻では `SCFTYP` を必ず書きます。
- `$SYSTEM` の `MWORDS` は 1 コアあたりのメモリー（100 万語単位）、基底関数は `$BASIS` の `GBASIS` などで組み立てます。Dunning 系には `ISPHER=1` が必要です。
- `$DATA` は「題名、点群、原子」の順に書き、C1 では空行なしですべての原子を、それ以外では空行のあとに対称独立な原子だけを書きます。
- DFT は `SCFTYP` と `DFTTYP` の組み合わせで、最適化は `RUNTYP=OPTIMIZE` と `$STATPT` で指定し、`EQUILIBRIUM GEOMETRY LOCATED` で収束を確かめます。

## 参考文献

- M. W. Schmidt, K. K. Baldridge, J. A. Boatz, S. T. Elbert, M. S. Gordon, J. H. Jensen, S. Koseki, N. Matsunaga, K. A. Nguyen, S. Su, T. L. Windus, M. Dupuis, J. A. Montgomery, "General atomic and molecular electronic structure system", J. Comput. Chem. 14, 1347 (1993). DOI: [10.1002/jcc.540141112](https://doi.org/10.1002/jcc.540141112)
- G. M. J. Barca ら, "Recent developments in the general atomic and molecular electronic structure system", J. Chem. Phys. 152, 154102 (2020). DOI: [10.1063/5.0005188](https://doi.org/10.1063/5.0005188)

## 関連リンク

- GAMESS の文書の一覧：<https://www.msg.chem.iastate.edu/gamess/documentation.html>
- GAMESS, INPUT.DOC（入力の説明）：<https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-input.txt>
- GAMESS, INTRO.DOC（入力の書式の決まり）：<https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-intro.txt>
- アイオワ州立大学 GAMESS の講習資料（アンモニアの最適化の出力など）：<https://www.msg.chem.iastate.edu/tutorials/tutorials.html>
