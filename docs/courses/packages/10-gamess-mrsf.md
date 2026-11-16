# 第 10 回　GAMESS（2）：MRSF-TDDFT で励起状態を計算する

!!! abstract "この回の道筋"
    MRSF-TDDFT は、三重項を出発点にして電子のスピンを反転させることで、基底状態と励起状態を同じ扱いで求める方法です。
    この回では、ふつうの TDDFT とスピン反転 TDDFT の弱点から始めて、MRSF-TDDFT の考え方を一言でつかみます。
    次に、GAMESS での入力（`TDDFT=MRSF`、ROHF の三重項の参照、BHHLYP）を INPUT.DOC に沿って組み立て、出力の読み方を学びます。

    | 節 | 内容 |
    |---|---|
    | 1 | ふつうの TDDFT の弱点 |
    | 2 | スピン反転 TDDFT（SF-TDDFT） |
    | 3 | MRSF-TDDFT の考え方 |
    | 4 | GAMESS での入力 |
    | 5 | 出力を読む |
    | 6 | 構造最適化と注意点 |

    **前に読んでおく回**：計算パッケージの使い方 [第 3 回](03-common-concepts.md)、[第 9 回](09-gamess-basics.md)

## 1. ふつうの TDDFT の弱点

TDDFT（時間依存密度汎関数理論）の線形応答は、励起エネルギーを計算する最もよく使われる方法です。
ふつうの TDDFT では、閉殻の基底状態（一つの Kohn–Sham 行列式）を参照にし、そこから電子を 1 個ずつ励起した配置を重ね合わせて、励起状態を表します。

この方法は、基底状態が一つの配置でよく表せる分子では、うまく働きます。
しかし、次のような場合には苦手です。

- 結合が切れかけた構造や、ジラジカルのように、基底状態そのものが複数の配置の重ね合わせになる場合。
- 基底状態と第一励起状態のエネルギーが近づき、交わる構造（円錐交差）の付近。

どちらも原因は同じで、基底状態だけが参照として特別扱いされ、励起状態と同じ扱いで表されないからです。
基底状態に重要な配置（たとえば二電子励起の配置）は、ふつうの TDDFT の応答には入りません。

**ここまでで分かったこと**：ふつうの TDDFT は閉殻の基底状態を特別扱いするので、基底状態が多配置的な場合や、基底状態と励起状態が交わる付近を苦手とします。

## 2. スピン反転 TDDFT（SF-TDDFT）

この弱点を避けるために、Shao、Head-Gordon、Krylov（2003）は、参照を高スピンの三重項に変えたうえで、電子のスピンを反転させる励起を考える TDDFT を作りました。
これがスピン反転（spin-flip、SF）TDDFT です。

考え方を、二つの不対電子が入る軌道を $\phi_1$、$\phi_2$ として説明します。
三重項の $M_S = +1$ の成分では、$\phi_1$ と $\phi_2$ に α 電子が 1 個ずつ入っています。
ここから、α 電子の 1 個を β 電子に変えながら移す（スピンを反転させる）と、$M_S = 0$ の配置ができます。

- $\phi_2$ の α 電子を $\phi_1$ の β に移すと、$\phi_1$ に 2 個入った閉殻の配置。
- $\phi_1$ の α 電子を $\phi_2$ の β に移すと、$\phi_2$ に 2 個入った配置。
- $\phi_1$ と $\phi_2$ に 1 個ずつ入った開殻の配置。

基底の一重項に必要な閉殻の配置と、励起状態に必要な配置が、どれも「三重項から 1 回スピンを反転した配置」として同じ扱いで得られます。
そのため、基底状態も応答の状態の一つとして求められ、基底状態と励起状態の交わりも表せます。

ただし、SF-TDDFT には弱点があります。
$M_S = +1$ の成分からの反転（α → β）だけでは、純粋なスピン状態を作るのに必要な配置の一部が足りません。
その結果、得られる状態に一重項と三重項が混ざります（スピン汚染、第 3 回）。
GAMESS の INPUT.DOC も、SF-TDDFT の状態はスピン汚染していて、どの多重度の状態かは出力を見て判断する必要があると注意しています。

**ここまでで分かったこと**：SF-TDDFT は、三重項の $M_S=+1$ の成分から 1 回スピンを反転させて、基底状態と励起状態を同じ扱いで求めます。ただし、反転の向きが一つだけなので、スピン汚染が残ります。

## 3. MRSF-TDDFT の考え方

MRSF-TDDFT（mixed-reference spin-flip TDDFT）は、Lee、Filatov、Lee、Choi（2018）が提案した、SF-TDDFT の改良版です。
考え方を一言で言うと、次のとおりです。

> 三重項の $M_S=+1$ と $M_S=-1$ の成分を混ぜた参照（混合参照）から、スピンを反転させる。

$M_S = +1$ の成分からは α → β の反転で、$M_S = -1$ の成分からは β → α の反転で、$M_S = 0$ の配置を作ります。
二つの向きの反転を合わせることで、SF-TDDFT で足りなかった配置の多くが補われます。
原論文では、この混合参照を、二つの成分の密度行列を組み合わせた仮想的な参照の密度行列として定式化しています。
一重項と三重項を混ぜるのではなく、**同じ三重項の二つの成分**を混ぜる点に注意します。
原論文の要旨によると、これによって SF-TDDFT の応答状態のスピン汚染はほぼ取り除かれ、状態の見分けが簡単になります。

INPUT.DOC によると、GAMESS の MRSF-TDDFT では、求める状態のスピンを指定できます。
$M_S = 0$ の一重項と三重項、$M_S = 2$ の五重項を、それぞれ別々に求めます。
一重項と三重項がはっきり分かれるので、一重項どうしの円錐交差の探索も簡単になります（Lee ら, 2019）。

**ここまでで分かったこと**：MRSF-TDDFT は、三重項の $M_S=\pm 1$ の成分を混ぜた参照から両方の向きにスピンを反転させ、SF-TDDFT のスピン汚染の大部分を取り除きます。求める状態のスピンを指定できます。

## 4. GAMESS での入力

### 4.1 入力の例

次は、エチレンの基底状態と低い一重項の励起状態を、MRSF-TDDFT（BHHLYP/6-31G(d)）で計算する入力の例です。
キーワードはすべて、GAMESS の INPUT.DOC（2024 年 7 月 15 日版）の記述と照らし合わせたものです。
GAMESS を手元で動かせなかったため、この入力は実行していません。

```text
 $CONTRL SCFTYP=ROHF MULT=3 DFTTYP=BHHLYP TDDFT=MRSF RUNTYP=ENERGY $END
 $SYSTEM MWORDS=100 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $TDDFT  NSTATE=3 IROOT=1 MULT=1 $END
 $DATA
Ethylene MRSF-TDDFT BHHLYP/6-31G(d)
C1
C   6.0   0.000000   0.000000   0.666000
C   6.0   0.000000   0.000000  -0.666000
H   1.0   0.000000   0.923000   1.239000
H   1.0   0.000000  -0.923000   1.239000
H   1.0   0.000000   0.923000  -1.239000
H   1.0   0.000000  -0.923000  -1.239000
 $END
```

`$SYSTEM`、`$BASIS`、`$DATA` の書き方は第 9 回と同じです。
MRSF-TDDFT に特有の部分を、次に説明します。

### 4.2 $CONTRL：参照の三重項

| キーワード | 意味 |
|---|---|
| `SCFTYP=ROHF` | 参照の軌道を ROHF（DFT なら ROKS）で求める |
| `MULT=3` | 参照は三重項 |
| `DFTTYP=BHHLYP` | 汎関数 |
| `TDDFT=MRSF` | MRSF-TDDFT の計算をする |

`$CONTRL` の `MULT=3` は、**参照**の多重度です。
求めたい状態が一重項でも、ここは 3 にします。
INPUT.DOC によると、現在の MRSF-TDDFT が対応している参照は、三重項（`MULT=3`）だけです。

参照の軌道は UHF でも ROHF でも求められますが、INPUT.DOC によると、エネルギーの勾配（6 節）が計算できるのは ROHF の参照だけです。
そのため、ここでは `SCFTYP=ROHF` にしました。

### 4.3 汎関数：BHHLYP

`DFTTYP=BHHLYP` は、HF の交換と Becke（B88）の交換を半分ずつ混ぜ、相関に LYP を使うハイブリッド汎関数です。
INPUT.DOC の `$LIBXC` の汎関数の一覧には、同じ名前の汎関数が「0.5HF + 0.5Becke88 + LYP」と書かれています。
HF の交換を半分取り入れる考え方は、Becke（1993）の「half-and-half」の汎関数にさかのぼります（Becke の原論文の式そのものは、B88 や LYP を使う BHHLYP とは違います）。

スピン反転の方法でこのような汎関数がよく使われるのには、理由があります。
GAMESS の SF-TDDFT と MRSF-TDDFT は、「collinear 近似」と呼ばれる形で作られています。
INPUT.DOC によると、この近似では、スピンを反転させる励起のエネルギーに大きく効くのは HF の交換の項だけです。
そのため、HF の交換を多く含むハイブリッド汎関数が使われます。
SF-TDDFT の原論文でも、HF の交換が 50% のハイブリッド汎関数で、結果が大きく改善したと報告されています（Shao ら, 2003）。

### 4.4 $TDDFT：求める状態

| キーワード | 意味 | 既定 |
|---|---|---|
| `NSTATE` | 求める状態の数（参照は数えない） | 1 |
| `IROOT` | 性質の計算や構造最適化に使う状態の番号 | 1 |
| `MULT` | 求める状態の多重度：1、3、5 | 1 |

`$TDDFT` の `MULT` は、**求める状態**の多重度です。
`$CONTRL` の `MULT`（参照の多重度）と同じ名前なので、取り違えないようにします。
一重項の励起状態を求めるなら 1、三重項なら 3 です。

INPUT.DOC によると、MRSF と SF の計算は Tamm–Dancoff 近似（励起の振幅だけを使い、脱励起を無視する近似）でだけ作られているので、`TAMMD` は自動で `.TRUE.` になります。

**ここまでで分かったこと**：`$CONTRL` に `SCFTYP=ROHF MULT=3 TDDFT=MRSF` と汎関数を、`$TDDFT` に `NSTATE`、`IROOT` と、求める状態の `MULT` を書きます。二つの `MULT` の意味の違いに注意します。

## 5. 出力を読む

### 5.1 状態の番号のつけ方

MRSF-TDDFT の出力を読むときに最も大切なのは、状態の番号の意味です。

- **状態 0** は参照、つまり ROHF の三重項です。
- **状態 1、2、…** が応答の状態です。

INPUT.DOC は、SF と MRSF では、**最初の応答の状態（状態 1）が基底状態にあたる**と強調しています。
`MULT=1` で求めた一重項には、基底の一重項も含まれるからです。
つまり、状態 1 が S₀、状態 2 が S₁、状態 3 が S₂ です（構造によっては順番が変わることもあるので、状態の中身も確かめます）。

### 5.2 エネルギーの表

この節で示す出力の見出しと数値は、手元で GAMESS を動かして得たものではなく、開発グループの解説（関連リンクの MQCP の文書）に載っている出力の例から取ったものです。
この例は、4.1 節のエチレンではなく、Be 原子を MRSF-TDDFT（BHHLYP/6-31G）で計算したものです。
この例では、計算の最後に次の見出しの表が出ます。

```text
SUMMARY OF MRSF-DFT RESULTS
```

表には、状態ごとに、全エネルギー（ハートリー）、励起エネルギー（eV）、$\langle \hat S^2 \rangle$、遷移双極子、振動子強度が並びます。
ここで、励起エネルギーの列は、**参照（状態 0、三重項）からの**エネルギーです。
基底の一重項は参照の三重項より低いことが多いので、状態 1 の励起エネルギーは負の値になることがあります。
解説の Be 原子の例でも、状態 1 の値は $-2.296$ eV でした。

S₀ からの鉛直励起エネルギーを求めるには、状態 1 を基準にします。
たとえば S₁ の鉛直励起エネルギーは、次の差です。

$$
\Delta E(\mathrm{S}_0 \to \mathrm{S}_1) = E(\text{状態 2}) - E(\text{状態 1})
$$

全エネルギーの列（ハートリー）で差をとり、1 $E_\mathrm{h}$ = 27.211 eV で換算します。
解説の例では、この表のすぐ下に、状態どうしの励起エネルギーの表（`1  ->  2` などの行）も出ていて、状態 1 からの励起エネルギーを直接読めます。
同じ解説は、参照の三重項は変分法で、応答の状態は線形応答で求めた別の種類の量なので、参照の三重項のエネルギーを解析に使わないよう勧めています。
三重項の状態が必要なら、`$TDDFT` の `MULT=3` で、三重項の応答の状態として求めます。

### 5.3 IROOT の状態

同じ解説の例（Be 原子）では、`IROOT` で選んだ状態が、次の形の行で示されます。

```text
SELECTING EXCITED STATE IROOT=   1 AT E=      -14.6509883896
```

構造最適化（6 節）や性質の計算は、この状態について行われます。
4.1 節の入力では `IROOT=1` なので、基底状態 S₀ が選ばれます。

出力の細かな形は GAMESS の版によって変わることがあるので、`MRSF` や `IROOT` などの語で検索して読みます。

**ここまでで分かったこと**：状態 0 は参照の三重項、状態 1 が基底状態 S₀ です。励起エネルギーの列は参照からの値なので、S₀ からの励起エネルギーは状態 1 との差で求めます。

## 6. 構造最適化と注意点

### 6.1 励起状態の構造最適化

INPUT.DOC によると、MRSF-TDDFT では、ROHF の参照について、エネルギーの勾配が解析的に計算できます。
そのため、第 9 回と同じ `RUNTYP=OPTIMIZE` で、励起状態の構造最適化ができます。
最適化する状態は、`$TDDFT` の `IROOT` で選びます。
たとえば S₁ を最適化するなら、状態 2 にあたるので `IROOT=2` にします。

```text
 $CONTRL SCFTYP=ROHF MULT=3 DFTTYP=BHHLYP TDDFT=MRSF RUNTYP=OPTIMIZE $END
 $TDDFT  NSTATE=3 IROOT=2 MULT=1 $END
```

INPUT.DOC には、ほかに、二つや三つの状態の円錐交差の探索（`$CONICL`）や、非断熱の分子動力学（`$MD`）も、MRSF-TDDFT で使えると書かれています。

### 6.2 参照の三重項に注意する

MRSF-TDDFT の結果は、参照の ROHF の三重項の軌道の上に組み立てられます。
INPUT.DOC の `$SCF` の説明によると、エネルギーの近い三重項が複数あると、SCF が構造ごとに違う三重項に収束し、ポテンシャルエネルギー面が不連続になることがあります。
このようなときは、`$SCF` の `MOM=.TRUE.` で、読み込んだ初期軌道と同じ性格の参照を保つようにできます。
最適化の途中でエネルギーが急に跳んだら、まず参照の三重項が変わっていないかを確かめます。

**ここまでで分かったこと**：MRSF-TDDFT では、`RUNTYP=OPTIMIZE` と `IROOT` で励起状態の構造最適化ができます。参照の三重項が構造の途中で入れ替わらないかに注意します。

## まとめ

- ふつうの TDDFT は閉殻の基底状態を特別扱いするので、基底状態が多配置的な場合や、基底状態と励起状態の交わりの付近を苦手とします。
- SF-TDDFT は三重項の $M_S=+1$ の成分からスピンを反転させて基底状態と励起状態を同じ扱いで求めますが、スピン汚染が残ります。
- MRSF-TDDFT は、$M_S=+1$ と $M_S=-1$ の成分を混ぜた参照から両方の向きに反転させ、スピン汚染の大部分を取り除きます。
- GAMESS では `$CONTRL` に `SCFTYP=ROHF MULT=3 TDDFT=MRSF` と `DFTTYP=BHHLYP` などを、`$TDDFT` に `NSTATE`、`IROOT`、求める状態の `MULT` を書きます。
- 状態 0 は参照の三重項、状態 1 が基底状態 S₀ です。励起エネルギーは状態 1 との差で求め、構造最適化は `IROOT` で状態を選びます。

## 参考文献

- Y. Shao, M. Head-Gordon, A. I. Krylov, "The spin–flip approach within time-dependent density functional theory: Theory and applications to diradicals", J. Chem. Phys. 118, 4807 (2003). DOI: [10.1063/1.1545679](https://doi.org/10.1063/1.1545679)
- S. Lee, M. Filatov, S. Lee, C. H. Choi, "Eliminating spin-contamination of spin-flip time dependent density functional theory within linear response formalism by the use of zeroth-order mixed-reference (MR) reduced density matrix", J. Chem. Phys. 149, 104101 (2018). DOI: [10.1063/1.5044202](https://doi.org/10.1063/1.5044202)
- S. Lee, E. E. Kim, H. Nakata, S. Lee, C. H. Choi, "Efficient implementations of analytic energy gradient for mixed-reference spin-flip time-dependent density functional theory (MRSF-TDDFT)", J. Chem. Phys. 150, 184111 (2019). DOI: [10.1063/1.5086895](https://doi.org/10.1063/1.5086895)
- A. D. Becke, "A new mixing of Hartree–Fock and local density-functional theories", J. Chem. Phys. 98, 1372 (1993). DOI: [10.1063/1.464304](https://doi.org/10.1063/1.464304)
- G. M. J. Barca ら, "Recent developments in the general atomic and molecular electronic structure system", J. Chem. Phys. 152, 154102 (2020). DOI: [10.1063/5.0005188](https://doi.org/10.1063/5.0005188)

## 関連リンク

- GAMESS, INPUT.DOC（`TDDFT=MRSF` と `$TDDFT` の説明）：<https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-input.txt>
- MQCP の文書, MRSF-TDDFT の解説（GAMESS の入力と出力の例）：<https://mqcp.readthedocs.io/en/latest/mrsf.html>
