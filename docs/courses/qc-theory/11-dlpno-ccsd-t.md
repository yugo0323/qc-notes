# 第 11 回　局所相関法：DLPNO-CCSD(T)

!!! abstract "この回の道筋"
    CCSD(T) は正確ですが、計算量が系の大きさの 7 乗で増えるため、大きな分子には使えません。この回では、「電子相関は短い距離でしか効かない」という性質を使って計算を大幅に減らす局所相関法を学びます。その代表である DLPNO-CCSD(T) について、ペア自然軌道（PNO）の考え方、精度を決める閾値（LoosePNO / NormalPNO / TightPNO）、そして使うときの注意点（PNO の打ち切り誤差とその外挿）を、原論文とベンチマークの論文の数値で確かめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 局所相関の考え方：相関は短距離 |
    | 2 | 電子対に分ける：ペア相関エネルギー |
    | 3 | ペア自然軌道（PNO） |
    | 4 | DLPNO-CCSD(T) の組み立て |
    | 5 | 閾値：LoosePNO / NormalPNO / TightPNO |
    | 6 | 精度と注意点 |

    **前に読んでおく回**：[量子化学の理論 第 9 回](09-mp2.md)、[第 10 回](10-coupled-cluster.md)

!!! tip "物理的なイメージ"
    局所相関法が利用しているのは、「電子が避け合う動きは、近くの電子どうしでしか大きく効かない」という性質です。局在化した軌道で見ると、離れた電子対の相関エネルギーは、距離の 6 乗に反比例して弱まります（分散力と同じ形）。そのため、ほとんどの対は安く見積もるだけで十分です。大事な対についても、その対の避け合いに必要な空軌道の組み合わせはわずかで、それを対ごとに選び直したものがペア自然軌道（PNO）です。閾値は、どこまで小さい寄与を捨てるかを決めるものさしです。

この回では、閉殻の分子を考えます。$i, j$ は占有の空間軌道、$a, b$ は仮想の空間軌道です。

## 1. 局所相関の考え方：相関は短距離

第 10 回で見たように、CCSD(T) の計算量は $N^7$ に比例します（$N$ は系の大きさ）。分子が 10 倍になると、計算時間は $10^7$ 倍です。これは手法そのものの性質ではなく、HF の**正準軌道**（Fock 行列を対角化して得た軌道）を使って式を書いたことから生じています。

正準軌道は、たいてい分子全体に広がっています。そのため、どの軌道の組を選んでも、2 電子積分や振幅が 0 になりません。ところが、物理的には電子相関は短い距離の現象です。第 7 回で見たように、動的相関は「近づいた電子どうしが互いを避ける」動きでした。

離れた二つの電子対の相関がどう減るかは、次のように考えると分かります。原子 A の近くの電子対と、遠く離れた原子 B の近くの電子対があるとします。二つの電子対の電子が同時にゆらぐと、それぞれに瞬間的な双極子が生じます。双極子どうしの相互作用は距離 $R$ の 3 乗に反比例し、エネルギーへの寄与は 2 次の摂動なのでその 2 乗、つまり $R^{-6}$ に比例します。これは、ロンドンの分散力と同じ形です。距離が 2 倍になると、寄与は 64 分の 1 になります。

そこで、軌道を空間的に局在させ、遠く離れた軌道どうしの相関を省いたり、粗く扱ったりすれば、計算量を減らせるはずです。この考え方を**局所相関**（local correlation）と呼び、Pulay が 1980 年代に提案しました [Pulay 1983, Saebø 1993]。

**ここまでで分かったこと**：CCSD(T) の高い計算量は、広がった正準軌道を使うことから来ています。相関は短距離なので、局在化した軌道を使えば、遠く離れた部分の計算を省けます。

## 2. 電子対に分ける：ペア相関エネルギー

占有軌道は、互いに混ぜ合わせても（ユニタリー変換しても）、HF の決定子は変わりません。そこで、占有軌道を空間的に局在した軌道（**局在化軌道**）に変換します。局在化軌道は、化学で言う「C–H 結合」「O の孤立電子対」などに近い形になります。

CC 法の相関エネルギーは、占有軌道の組 $(i, j)$ ごとの和に分けられます。

$$
E_\mathrm{corr} = \sum_{i \le j} \varepsilon_{ij}
$$

$\varepsilon_{ij}$ を**ペア相関エネルギー**と呼びます。第 9 回の MP2 の式でも、$i, j$ についての和の中身を取り出せば、そのまま $\varepsilon_{ij}$ になります。局在化軌道を使うと、$i$ と $j$ が遠く離れた組では、1 節の理由で $\varepsilon_{ij}$ が非常に小さくなります。

電子対の数は占有軌道の数の 2 乗に比例します。しかし、遠い組の寄与を無視すれば、意味のある組の数は系の大きさに比例するだけになります。1 個の軌道のまわりにある軌道の数は、系がどれだけ大きくても、ほぼ一定だからです。

**ここまでで分かったこと**：局在化軌道を使うと、相関エネルギーを電子対ごとの寄与に分けられ、遠い電子対の寄与は小さくなります。

## 3. ペア自然軌道（PNO）

電子対の数を減らすだけでは、まだ不十分です。電子対 $(i, j)$ ごとに、仮想軌道 $a, b$ のすべての組についての振幅 $t_{ij}^{ab}$ が必要だからです。仮想軌道の数 $v$ は系の大きさに比例するので、1 組あたりの振幅の数は $v^2$ 個になります。

そこで、**電子対ごとに、その電子対の相関に必要な仮想軌道だけを選び直します**。手順は次のとおりです。

1. まず、安い近似（MP2 の振幅）で $t_{ij}^{ab}$ を見積もります。
2. その振幅から、電子対 $(i, j)$ の励起先の電子の分布を表す密度行列（ペア密度行列）をつくります。
3. これを対角化します。固有ベクトルが新しい仮想軌道、固有値がその軌道の**占有数**（その軌道がどれだけ使われているか）です。
4. 占有数が閾値 $T_\mathrm{CutPNO}$ より大きい軌道だけを残します。

こうして得た軌道を、**ペア自然軌道**（pair natural orbital, PNO）と呼びます。占有数の大きな少数の PNO が、その電子対の相関のほとんどを担います。そのため、少ない数の PNO で、元の仮想軌道すべてを使った場合とほぼ同じ相関エネルギーが得られます。

PNO の考え方は古く、1960 年代に提案されました [Edmiston 1966]。その後しばらく使われませんでしたが、Neese らが、効率よく実装できることを示して復活させました [Neese 2009]。この最初の方法では、既定の閾値で正準な計算の相関エネルギーの約 99.8% が得られたと報告されています。

**ここまでで分かったこと**：PNO は、電子対ごとに選び直した、その電子対専用の仮想軌道です。占有数が $T_\mathrm{CutPNO}$ より大きいものだけを残すことで、振幅の数を大きく減らせます。

## 4. DLPNO-CCSD(T) の組み立て

PNO を使った CCSD を、さらに大きな分子で使えるようにしたのが、Riplinger と Neese の **DLPNO-CCSD**（domain-based local PNO CCSD）です [Riplinger 2013a]。その主な工夫を、計算の流れに沿って並べます。

1. **電子対のふるい分け**：局在化軌道のすべての組について、ペア相関エネルギーを安い方法で見積もります。値が閾値 $T_\mathrm{CutPairs}$ より小さい組は「弱い電子対」として、CCSD の計算から外します。外した組の寄与は、MP2 の程度の近似で見積もって足します。
2. **ドメイン**：各 PNO を、その電子対の近くにある原子の**射影原子軌道**（PAO）だけで展開します。PAO は、各原子の基底関数から占有軌道の成分を取り除き、仮想軌道の空間に射影したものです。Pulay が局所相関のために導入しました [Pulay 1983]。この原子の範囲を**ドメイン**と呼びます。原論文では、ドメインによる相関エネルギーの損失は 0.05% 未満で、ドメインは平均 15〜20 原子（ときに 30 原子まで）でした [Riplinger 2013a]。
3. **CCSD**：残った「強い電子対」について、それぞれの PNO の空間で CCSD の方程式を解きます。CCSD の方程式の項は省きません。
4. **(T)**：3 電子励起の補正は、3 個の電子対の PNO の空間を合わせた「3 電子励起の自然軌道」（TNO）を使って計算します [Riplinger 2013b]。

原論文の冒頭では、PNO を使うと、大きな分子では振幅の数が、同じ系・同じ基底関数の正準 CCSD に比べて $10^5$〜$10^6$ 分の 1 に減ると述べられています [Riplinger 2013a]。計算量は、系の大きさにほぼ比例するところまで下がります。(T) の計算は全体の時間を 30〜50% 増やす程度で、タンパク質クランビン（644 原子、6400 以上の基底関数）全体の CCSD(T) 計算が示されました [Riplinger 2013b]。

### (T0) と (T)

最初の DLPNO-CCSD(T) では、(T) の計算に**半正準近似**を使っていました。後に、この近似を外した反復的な (T) が開発され、元の方法は区別のため DLPNO-CCSD(T0) と呼ばれるようになりました [Guo 2018]。Guo らによると、(T0) はほとんどの場合、相対エネルギーでは十分に正確です。しかし、まれに、特に HOMO–LUMO のエネルギー差が小さい系では (T0) が破綻し、絶対エネルギーの誤差も (T) の補正が大きいほど大きくなります。反復的な (T) の計算時間は、(T0) のおよそ 2 倍です。

ORCA では、キーワード `DLPNO-CCSD(T)` が (T0) を、`DLPNO-CCSD(T1)` が反復的な (T) を意味します（ORCA 6.1 で、反復の (T) の計算が出力に現れるかどうかで確かめました）。

**ここまでで分かったこと**：DLPNO-CCSD(T) は、電子対のふるい分け、ドメイン、PNO、TNO を組み合わせて、CCSD(T) の計算量をほぼ線形に下げた方法です。(T) の扱い方には (T0) と反復的な (T) の 2 種類があります。

## 5. 閾値：LoosePNO / NormalPNO / TightPNO

DLPNO-CCSD(T) の精度と速さは、主に次の三つの閾値で決まります。

- $T_\mathrm{CutPNO}$：PNO を残すかどうかを決める、占有数の閾値です。小さいほど PNO が増えます。
- $T_\mathrm{CutPairs}$：電子対を CCSD で扱うかどうかを決める、ペア相関エネルギーの見積もりの閾値です。
- $T_\mathrm{CutDO}$：ドメインに原子を入れるかどうかを決める閾値です。

利用者が個々の値を選ばなくて済むように、三つの組み合わせ LoosePNO、NormalPNO、TightPNO が定義されています [Liakos 2015]。次の表は、ORCA 6.1 で H₂O の小さな計算を実行し、出力の「DLPNO SETTINGS」の欄に表示された値を写したものです。出力には、前段の粗い MP2 や (T) のために、これを 2 倍、4 倍にゆるめた $T_\mathrm{CutDO}$ も表示されるので、読み違えないようにしてください。版によって値が変わることがあるので、使う版の出力で確かめてください。

| 設定 | $T_\mathrm{CutPNO}$ | $T_\mathrm{CutPairs}$ | $T_\mathrm{CutDO}$ |
|---|---|---|---|
| LoosePNO | $1 \times 10^{-6}$ | $10^{-3}$ | 0.02 |
| NormalPNO | $3.33 \times 10^{-7}$ | $10^{-4}$ | 0.01 |
| TightPNO | $1 \times 10^{-7}$ | $10^{-5}$ | 0.005 |

ORCA の入力は、たとえば次のようになります。手法、閾値の設定、基底関数、RI 近似の補助基底関数を 1 行目に並べます。

```
! DLPNO-CCSD(T) NormalPNO cc-pVDZ cc-pVDZ/C
* xyz 0 1
O 0.000 0.000 0.000
H 0.000 0.759 0.588
H 0.000 -0.759 0.588
*
```

`NormalPNO` を `TightPNO` に替えれば、閾値の組がまとめて替わります。

**ここまでで分かったこと**：DLPNO-CCSD(T) の精度は主に $T_\mathrm{CutPNO}$、$T_\mathrm{CutPairs}$、$T_\mathrm{CutDO}$ で決まり、LoosePNO、NormalPNO、TightPNO の三つの組が用意されています。

## 6. 精度と注意点

DLPNO-CCSD(T) の誤差は、「同じ基底関数の正準 CCSD(T) からのずれ」として測ります。基底関数の不完全さの誤差や、CCSD(T) そのものの誤差は、これとは別にあります。

### ベンチマークの結果

- Liakos らは、既定の閾値（NormalPNO）で、すべてのテストセットについて 1 kcal/mol 以内の精度が得られることを示しました。さらに閾値を厳しくすると、正準な結果からの平均絶対誤差と標準偏差が 0.25 kcal/mol 未満になりました [Liakos 2015]。
- Altun らの論文の序論によると、TightPNO（$T_\mathrm{CutPNO} = 10^{-7}$）では、正準 CCSD(T) の相関エネルギーのおよそ 99.9% が得られます。また、GMTKN55 と呼ばれる 1505 個の反応からなるベンチマークでは、TightPNO の平均絶対誤差は、55 個の部分集合のうち 39 個で 0.2 kcal/mol 未満、14 個で 0.2〜0.4 kcal/mol、2 個で 0.4 kcal/mol を超えました [Altun 2020]。
- 同じ論文によると、LoosePNO（$T_\mathrm{CutPNO} = 10^{-6}$）は予備的な計算にだけ勧められる設定です。$T_\mathrm{CutPNO} = 10^{-6}$ での平均絶対誤差は、多くの部分集合で 1.5〜3 kcal/mol で、人工的な分子からなる難しい集合では最大 6 kcal/mol に達しました [Altun 2020]。

### 注意 1：PNO の打ち切り誤差

DLPNO の誤差の主な原因の一つは、PNO を打ち切ったことによる誤差（**PNO 打ち切り誤差**）です。この誤差は、$T_\mathrm{CutPNO}$ を小さくすると、なめらかに小さくなっていきます [Riplinger 2013b, Altun 2020]。

この性質を使い、二つの $T_\mathrm{CutPNO}$ の結果から、PNO を打ち切らない極限（完全 PNO 空間の極限、CPS）を外挿する方法が提案されています [Altun 2020]。$T_\mathrm{CutPNO} = 10^{-X}$ と $10^{-Y}$（$Y = X + 1$）で求めた相関エネルギーを $E_X$、$E_Y$ とすると、

$$
E_\mathrm{CPS} = E_X + F\,(E_Y - E_X), \qquad F = 1.5
$$

です。$F$ は、基底関数の外挿の式と同じ形の式から出てくる係数で、ベンチマークに合うように 1.5 と決められました。$X = 6$、$Y = 7$ の組を CPS(6/7) と呼びます。たとえば、非共有結合の相互作用のベンチマーク S66 では、TightPNO の平均絶対誤差 0.20 kcal/mol が、CPS(6/7) では 0.11 kcal/mol に下がりました（aug-cc-pVDZ-DK 基底） [Altun 2020]。$T_\mathrm{CutPNO} = 10^{-6}$ の計算の時間は $10^{-7}$ の 30〜40% ほどなので、外挿の手間は TightPNO の計算より少し多い程度で済みます。

ただし、外挿で直るのは PNO 打ち切り誤差だけで、そのほかの近似（電子対のふるい分けやドメインなど）の誤差は残ります。また、Riplinger らは、(T) の打ち切り誤差を外挿すると絶対エネルギーは改善するが、相対エネルギーは改善しないと報告しています [Riplinger 2013b]。外挿の効果は、量の種類によって違うことに注意してください。

### 注意 2：遷移金属錯体

遷移金属錯体では、DLPNO の誤差が大きくなりやすいことが知られています。Altun らは、第一遷移系列の金属を含む反応のベンチマークで、二つの誤差の原因を調べました [Altun 2023]。一つは、3s・3p の内殻に近い電子（セミコア電子）の相関の扱いです。もう一つは、動的相関によって軌道が変化する効果を、DLPNO の出発点である局所 MP2 の見積もりが十分に捉えられないことです。この研究のベンチマーク（MOBH35、def2-SVP 基底、正準 CCSD(T) との差）で、反応エネルギーと障壁の平均絶対誤差は、従来の設定で $T_\mathrm{CutPNO} = 10^{-6}$ のとき 1.08 kcal/mol、$10^{-7}$ のとき 0.71 kcal/mol でした。改良した設定と補正を使うと、これが小さくなります。また、T1 診断（第 10 回）が 0.02 を超える系では、PNO 空間の大きさに対する収束が遅いことも報告されています。

### 注意 3：CCSD(T) の限界はそのまま残る

DLPNO-CCSD(T) は、正準 CCSD(T) を**近似する**方法です。第 10 回で見たように、多参照性が強い系では CCSD(T) そのものが破綻します。DLPNO にしても、この問題は直りません。T1・D1 診断などで、まず CCSD(T) を使ってよい系かどうかを確かめてください。

### 使い方のまとめ

ここまでの文献の結果をまとめると、次のような使い方が考えられます。

- LoosePNO は予備的な計算だけに使い、通常の計算には NormalPNO を使います。
- 非共有結合の相互作用など、小さなエネルギー差を議論するときは TightPNO を使うか、CPS(6/7) で外挿します。
- 遷移金属錯体では、$T_\mathrm{CutPNO}$ を変えて結果の収束を確かめます。
- どの場合も、DLPNO の誤差とは別に、基底関数の誤差と CCSD(T) そのものの誤差があることを忘れないでください。

**ここまでで分かったこと**：DLPNO-CCSD(T) は、NormalPNO で 1 kcal/mol 程度、TightPNO ではさらに高い精度で正準 CCSD(T) を再現します。主な誤差である PNO 打ち切り誤差は、二つの閾値の結果から外挿して減らせます。遷移金属錯体では誤差が大きくなりやすく、多参照性による CCSD(T) の破綻は DLPNO でも直りません。

## まとめ

- 電子相関は短距離（離れた電子対の寄与は $R^{-6}$ で減る）なので、局在化軌道を使えば遠い部分の計算を省けます。
- PNO は電子対ごとに選び直した仮想軌道で、占有数が $T_\mathrm{CutPNO}$ より大きいものだけを残します。
- DLPNO-CCSD(T) は、電子対のふるい分け、ドメイン、PNO、TNO を組み合わせ、CCSD(T) の計算量をほぼ線形に下げます。
- 閾値の組は LoosePNO、NormalPNO、TightPNO で、主な誤差である PNO 打ち切り誤差は CPS 外挿で減らせます。
- 遷移金属錯体では誤差の確認が必要で、多参照性による CCSD(T) の破綻は DLPNO でも直りません。

## 参考文献

- P. Pulay, "Localizability of dynamic electron correlation", Chem. Phys. Lett. 100, 151 (1983). DOI: [10.1016/0009-2614(83)80703-9](https://doi.org/10.1016/0009-2614(83)80703-9)
- S. Saebø, P. Pulay, "Local Treatment of Electron Correlation", Annu. Rev. Phys. Chem. 44, 213 (1993). DOI: [10.1146/annurev.pc.44.100193.001241](https://doi.org/10.1146/annurev.pc.44.100193.001241)
- C. Edmiston, M. Krauss, "Pseudonatural Orbitals as a Basis for the Superposition of Configurations. I. He2+", J. Chem. Phys. 45, 1833 (1966). DOI: [10.1063/1.1727841](https://doi.org/10.1063/1.1727841)
- F. Neese, F. Wennmohs, A. Hansen, "Efficient and accurate local approximations to coupled-electron pair approaches: An attempt to revive the pair natural orbital method", J. Chem. Phys. 130, 114108 (2009). DOI: [10.1063/1.3086717](https://doi.org/10.1063/1.3086717)
- C. Riplinger, F. Neese, "An efficient and near linear scaling pair natural orbital based local coupled cluster method", J. Chem. Phys. 138, 034106 (2013). DOI: [10.1063/1.4773581](https://doi.org/10.1063/1.4773581)
- C. Riplinger, B. Sandhoefer, A. Hansen, F. Neese, "Natural triple excitations in local coupled cluster calculations with pair natural orbitals", J. Chem. Phys. 139, 134101 (2013). DOI: [10.1063/1.4821834](https://doi.org/10.1063/1.4821834)
- Y. Guo, C. Riplinger, U. Becker, D. G. Liakos, Y. Minenkov, L. Cavallo, F. Neese, "Communication: An improved linear scaling perturbative triples correction for the domain based local pair-natural orbital based singles and doubles coupled cluster method [DLPNO-CCSD(T)]", J. Chem. Phys. 148, 011101 (2018). DOI: [10.1063/1.5011798](https://doi.org/10.1063/1.5011798)
- D. G. Liakos, M. Sparta, M. K. Kesharwani, J. M. L. Martin, F. Neese, "Exploring the Accuracy Limits of Local Pair Natural Orbital Coupled-Cluster Theory", J. Chem. Theory Comput. 11, 1525 (2015). DOI: [10.1021/ct501129s](https://doi.org/10.1021/ct501129s)
- A. Altun, F. Neese, G. Bistoni, "Extrapolation to the Limit of a Complete Pair Natural Orbital Space in Local Coupled-Cluster Calculations", J. Chem. Theory Comput. 16, 6142 (2020). DOI: [10.1021/acs.jctc.0c00344](https://doi.org/10.1021/acs.jctc.0c00344)
- A. Altun, C. Riplinger, F. Neese, G. Bistoni, "Exploring the Accuracy Limits of PNO-Based Local Coupled-Cluster Calculations for Transition-Metal Complexes", J. Chem. Theory Comput. 19, 2039 (2023). DOI: [10.1021/acs.jctc.3c00087](https://doi.org/10.1021/acs.jctc.3c00087)
