# 第 10 回　多参照性の診断と多参照法

!!! abstract "この回の道筋"
    単一参照の CCSD(T) や DFT が使えるかどうかは、静的相関（多参照性）の強さで決まります。この回では、まず多参照性の診断の定義を確かめ、それらがどこまで当てになるかを、文献と自分の計算（Cr$_2$ と TiCl$_4$）で調べます。次に、診断で多参照性が強いと分かったときに使う CASSCF の活性空間の選び方、CASPT2・NEVPT2 を遷移金属に使うときの注意、大きな活性空間を扱う DMRG の利用を順に見ます。

    | 節 | 内容 |
    |---|---|
    | 1 | 多参照性の診断：定義 |
    | 2 | 診断はどこまで当てになるか |
    | 3 | 計算例：Cr$_2$ と TiCl$_4$ |
    | 4 | 遷移金属の CASSCF の活性空間 |
    | 5 | CASPT2・NEVPT2 を使うときの注意 |
    | 6 | DMRG の利用 |

    **前に読んでおく回**：[第 4 回 遷移金属の電子相関](04-electron-correlation-tm.md)、[量子化学の理論 第 10 回](../qc-theory/10-coupled-cluster.md)、[第 12 回](../qc-theory/12-casscf.md)、[第 13 回](../qc-theory/13-caspt2-nevpt2.md)、[第 14 回](../qc-theory/14-dmrg-sci.md)

## 1. 多参照性の診断：定義

第 4 回で見たように、遷移金属の錯体では静的相関が強くなりやすく、その強さは系によって大きく違います。そこで、計算の前に「この系に単一参照の方法を使ってよいか」を見積もる量が、いくつも提案されてきました。これを**多参照性の診断**（multireference diagnostic）と呼びます。主なものを、もとにする量で分けて定義します。

### 1.1 結合クラスター法の振幅から

CCSD の 1 電子励起の振幅 $t_i^a$（$i$ は占有軌道、$a$ は仮想軌道）は、静的相関が強いと大きくなります。HF の軌道がよい出発点でないとき、それを補うように軌道を回し直す役目を、1 電子励起が担うからです。

- **T1 診断**：$t_i^a$ を並べたベクトルの長さを、相関させた電子の数 $N$ の平方根で割った量 $\mathcal T_1 = \sqrt{\sum_{ia}(t_i^a)^2/N}$ です（Lee, Taylor 1989、関連リンク）。定義と、主要族の分子での目安 0.02 は、[量子化学の理論 第 10 回](../qc-theory/10-coupled-cluster.md)で扱いました。
- **D1 診断**：$t_i^a$ を行列と見たときの最大の特異値です（Janssen, Nielsen 1998）。T1 は全体の平均なので、大きな分子の一部だけで振幅が大きいと薄まります。D1 は最も悪い部分の値を拾います。
- **振幅の最大値** $\max|t_i^a|$：1 つの励起だけが大きい場合を拾います。

### 1.2 エネルギーから

- **%TAE[(T)]**：全原子化エネルギー（TAE。分子をばらばらの原子にするのに要るエネルギー）のうち、CCSD(T) の三重励起 (T) が占める割合です。$100 \times [\mathrm{TAE}(\mathrm{CCSD(T)}) - \mathrm{TAE}(\mathrm{CCSD})]/\mathrm{TAE}(\mathrm{CCSD(T)})$ で定義します。この量は、Karton らが W4 理論の論文で提案したもので、CCSD(T) を超える高次の励起の寄与の大きさのよい目安になることが、繰り返し確かめられています（Karton ら 2006; Semidalas, Martin 2022）。ただし、原子化エネルギーを使うので、大きな分子には使いにくい量です（Semidalas, Martin 2022）。

### 1.3 自然軌道の占有数から

自然軌道とその占有数は、[量子化学の理論 第 12 回](../qc-theory/12-casscf.md)で扱いました。閉殻の分子で、静的相関がなければ、占有数は 2 か 0 に近くなります。2 や 0 から大きくずれた軌道があれば、静的相関が強いことを示します。

- **CASSCF の主配置の重み** $C_0^2$：CASSCF の波動関数で、最も重い配置の係数の 2 乗です。1 に近ければ単一参照に近いことを示します（Jiang ら 2012）。
- **M 診断**：Tishchenko らが提案した、自然軌道の占有数を使う診断です（Tishchenko ら 2008）。閉殻の分子では、HOMO と LUMO にあたる自然軌道の占有数 $n_\mathrm{HOMO}$、$n_\mathrm{LUMO}$ を使って $(2 - n_\mathrm{HOMO} + n_\mathrm{LUMO})/2$ になります（Semidalas, Martin 2022）。

### 1.4 DFT と HF の計算から

相関の計算をしなくても見積もれる診断もあります。

- **FOD（分数占有数で重みづけた密度）**：Grimme と Hansen は、電子の温度を高く設定した有限温度の DFT（たとえば TPSS で 5000 K）で軌道の占有数を分数にし、その分数の占有数で重みづけた電子密度 $\rho^\mathrm{FOD}$ を作りました（Grimme, Hansen 2015）。この密度は、強く相関した「熱い」電子がどこにあるかを示します。空間で積分した数 $N_\mathrm{FOD}$ は、静的相関の大きさの目安になります。
- **スピン対称性の破れ**：閉殻の RHF の解が、スピン対称性を破った UHF の解へ不安定になるかどうかです（[量子化学の理論 第 7 回](../qc-theory/07-electron-correlation.md)）。結合を伸ばした H$_2$ のような、近縮退による静的相関の典型的な兆候です（Semidalas, Martin 2022）。
- **厳密交換の割合への感度**：DFT の結果が、厳密交換の割合を変えたときにどれだけ変わるかです。第 5 回で見たスピン状態の傾きも、この種の量です。Semidalas と Martin は、原子化エネルギーの、HF 交換の割合に対する傾きを使う診断を調べています（Semidalas, Martin 2022）。

**ここまでで分かったこと**：多参照性の診断には、CCSD の振幅（T1、D1、振幅の最大値）、エネルギー（%TAE[(T)]）、自然軌道の占有数（$C_0^2$、M 診断）、DFT と HF の計算（FOD、スピン対称性の破れ）にもとづくものがあります。

## 2. 診断はどこまで当てになるか

### 2.1 目安の値は系の種類で違う

主要族の分子の T1 の目安は 0.02 でした。Jiang らは、3d 金属を含む多くの分子を調べ、遷移金属の化合物では T1 > 0.05、D1 > 0.15、|%TAE| > 10 を、静的相関が強いことの目安として提案しました（Jiang ら 2012）。主要族の目安を遷移金属にそのまま当てはめると、多くの錯体が「多参照」と判定されてしまいます。

ただし、Jiang ら自身が、目安を下回っても単一参照とは限らないことを指摘しています。静的相関が強いことで知られる Cr$_2$ や Mn$_2$ で、T1 < 0.05、D1 < 0.15 となる例があったからです。

### 2.2 診断どうしが食い違う

Duan らは、小さな有機分子 3165 個（平衡構造とひずんだ構造）について、15 個の診断を計算しました（Duan ら 2020）。その結果、診断どうしで多参照性の判定が食い違い、2 つずつの診断の間の直線的な相関も弱いものでした。また、相関エネルギーの回復率（%$E_\mathrm{corr}$）を予測する力は、DFT にもとづく安い診断のどれも、波動関数にもとづく最良の診断に遠く及びませんでした。

Semidalas と Martin は、MOBH35（第 7 回）の多くの診断を統計的に分類し、診断が次のグループに分かれることを示しました（Semidalas, Martin 2022）。

1. 1 電子励起の大きさ（T1、D1 など）
2. HOMO–LUMO ギャップの小ささと、相関のエントロピー（自然軌道の占有数にもとづく量）
3. 相関エネルギーの熱化学的な大きさ（%TAE など）と、厳密交換の割合への感度

グループが違う診断は、静的相関の違う側面を見ています。1 つの診断だけで判断せず、違うグループから複数を見るのがよい、ということです。

### 2.3 診断と CCSD(T) の誤差の関係

診断の本当の目的は、CCSD(T) などの誤差を予測することです。その点で、次の 2 つの結果は対照的です。

- **相関しなかった例**：第 6 回の SSE17 では、CCSD(T) の誤差は、調べたどの診断（CCSD の振幅、(T) の寄与、CASSCF の主配置の重み、自然軌道の占有数、厳密交換の割合への感度）とも相関しませんでした。D1 が 0.15 を超える錯体が 6 割以上、0.20 を超える錯体が 4 個ありましたが、CCSD(T) は精度を保ちました（Radoń ら 2024）。
- **役に立った例**：第 7 回の MOBH35 の反応 9 では、生成物の T1 が 0.05、D1 が 0.52 と大きく、実際に局所相関法が破綻しました（Semidalas, Martin 2022）。また 3dTMV では、CCSD の波動関数と PBE0 の解のスピン対称性の破れが、CCSD(T) が使える範囲をよく区切りました（Neugebauer ら 2023）。

つまり、診断の値が大きいことは「注意して確かめるべき」という警告であり、CCSD(T) が間違っていることの証明ではありません。逆に、値が小さいことも、正しさの保証ではありません。

**ここまでで分かったこと**：診断の目安の値は系の種類で違い、診断どうしも食い違います。診断は警告として使い、違う種類の診断を複数見ます。値の大小だけで、CCSD(T) の正否は決まりません。

## 3. 計算例：Cr$_2$ と TiCl$_4$

診断が出発点の選び方にどれだけ左右されるかを、PySCF で確かめます。比べるのは、静的相関が強いことで知られる Cr$_2$（Cr–Cr 1.68 Å）と、d 電子をもたない閉殻の TiCl$_4$（Ti–Cl 2.17 Å の正四面体）です。基底関数は def2-SVP で、Cr、Ti、Cl の 1s2s2p を凍結しました。スクリプトは付録 A に載せます。

```text
--- Cr2
  RHF：安定性解析でたどった解は、最初の解より 161.2 kcal/mol 低い
  CCSD（最初の RHF 解から）: T1 = 0.032, D1 = 0.090, max|t1| = 0.075
  CCSD（安定な RHF 解から）: T1 = 0.311, D1 = 1.119, max|t1| = 1.086
  UHF：E(UHF) - E(安定な RHF) = -223.9 kcal/mol, <S^2> = 4.86
  UHF 自然軌道で占有数が 0.02〜1.98 のもの: 1.65 1.49 1.49 1.43 1.17 1.17 0.83 0.83 0.57 0.51 0.51 0.35
--- TiCl4
  RHF：安定性解析でたどった解は、最初の解より 0.0 kcal/mol 低い
  CCSD（最初の RHF 解から）: T1 = 0.022, D1 = 0.072, max|t1| = 0.071
  CCSD（安定な RHF 解から）: T1 = 0.022, D1 = 0.072, max|t1| = 0.071
  UHF：E(UHF) - E(安定な RHF) = -0.0 kcal/mol, <S^2> = 0.00
  UHF 自然軌道で占有数が 0.02〜1.98 のもの: なし
```

結果から、次のことが分かります。

**TiCl$_4$** では、RHF の解は最初から安定で、UHF もスピン対称性を破りませんでした。T1 は 0.022 で、主要族の目安 0.02 をわずかに超えますが、遷移金属の目安 0.05 は大きく下回ります。どの診断も、単一参照でよいことを示しています。

**Cr$_2$** では、事情がまったく違います。

- **RHF の解が 1 つに決まらない**：最初に得た RHF の解は不安定で、安定性解析で軌道を回し直すと、161 kcal/mol も低い RHF の解が見つかりました。
- **診断の値が出発点で 10 倍変わる**：最初の RHF の解から出発した CCSD では、T1 = 0.032、D1 = 0.090 で、Jiang らの目安を下回りました。第 2.1 節で紹介した「Cr$_2$ でも目安を下回る例がある」という指摘と同じ状況です。ところが、安定な RHF の解から出発すると、T1 = 0.311、D1 = 1.119 と、目安を大きく超えました。1 電子励起の振幅は、出発点の軌道の回し直しに当たります。その最大値が 1 を超えるのは、軌道を大きく回し直さなければならないほど、出発点の決定子が目的の状態から離れていることを示します。このような場合、CCSD の結果そのものも信頼できません。
- **スピン対称性を大きく破る**：UHF は、安定な RHF より 224 kcal/mol も低く、$\langle S^2\rangle$ は一重項の 0 ではなく 4.86 でした。UHF の自然軌道には、占有数が 2 と 0 から大きく離れた軌道が 12 個ありました。

Cr$_2$ のような系では、T1 や D1 の値は、どの RHF の解から出発したかで大きく変わり、それだけでは判断できません。スピン対称性の破れと自然軌道の占有数は、どちらの RHF の解から見ても、強い静的相関をはっきり示しています。遷移金属の計算では、まず SCF の解の安定性を確かめ、そのうえで違う種類の診断を見ることが大切です。

**ここまでで分かったこと**：TiCl$_4$ ではどの診断も単一参照を示しました。Cr$_2$ では、T1・D1 が出発点の RHF の解で 10 倍変わり、目安を下回ることもありましたが、スピン対称性の破れと自然軌道の占有数は強い静的相関を示しました。

## 4. 遷移金属の CASSCF の活性空間

診断で静的相関が強いと分かったら、CASSCF などの多参照の方法を使います。CASSCF の結果は、活性空間の選び方で大きく変わります（[量子化学の理論 第 12 回](../qc-theory/12-casscf.md)）。ここでは、遷移金属に特有の選び方を見ます。

### 4.1 単核の錯体の標準的な選び方

Radoń らは、SSE17 の CASPT2 の計算で、Pierloot の指針に基づいて、次の軌道を活性空間に入れました（Radoń ら 2024）。

1. **金属の 5 個の 3d 軌道**。
2. **配位子の軌道 1〜2 個**：主に二重に占有され、3d 軌道と大きく重なって、金属–配位子の共有結合をつくる軌道です。
3. **4d の性格をもつ空軌道（最大 5 個）**：第 4 回の二重殻効果（3d′ 軌道）を表すためです。π 逆供与を表す役も兼ねることがあります。

正八面体の錯体では、これで 10〜12 個の活性軌道になり、金属–配位子の共有結合と二重殻効果を扱える標準的な選び方とされています。ただし、低いスピンの状態で、対応する 3d 軌道がほとんど空のときは、4d の軌道を 3 個に減らしています。入れたままにすると、軌道が意図しない形に回ってしまう（制御できない軌道の回転）ためです。また、有機金属の Ni の錯体では、14 個の活性軌道を使っています。

このように、「何個の 3d′ 軌道を入れるか」は、状態ごとに確かめる必要があります。

### 4.2 反応に沿って活性空間をそろえる

反応エネルギーや障壁を多参照の方法で計算するときは、反応物、遷移状態、生成物で、同じ性格の活性空間を使う必要があります。活性空間の中身が途中で入れ替わると、エネルギーの差に大きな誤差が入ります。Hehn らは、Fe(III) を含む反応では、反応全体で一貫した活性空間を選べず、NEVPT2 の利点を生かせなかったと報告しています（Hehn ら 2024）。

### 4.3 選ぶための道具

活性空間の選び方は、経験と専門知識を要し、自動で選ぶ方法は長く難しいとされてきました（Veryazov ら 2011; Stein, Reiher 2016）。現在は、次のような道具があります。

- **AVAS**：金属の 3d 軌道など、指定した原子軌道に射影して活性軌道を選びます（[量子化学の理論 第 12 回](../qc-theory/12-casscf.md)）。
- **UHF の自然軌道の占有数**：第 3 節の Cr$_2$ のように、占有数が 2 と 0 から離れた軌道を選びます（同上）。
- **FOD**：$\rho^\mathrm{FOD}$ に寄与する軌道を、活性空間の候補にします（Bauer ら 2017）。
- **autoCAS**：DMRG で 100 個程度までの軌道を含む大きな空間を計算し、軌道どうしのもつれの大きさから活性軌道を選びます（Stein, Reiher 2016）。

どの道具を使っても、選ばれた軌道の形と占有数を目で確かめ、活性空間を変えて結果が変わらないかを確かめることが必要です。

**ここまでで分かったこと**：単核の錯体では、3d 軌道 5 個、共有結合をつくる配位子の軌道 1〜2 個、二重殻効果のための 4d の軌道を入れるのが標準的です。反応に沿って活性空間の性格をそろえ、自動の道具で選んだ結果も目で確かめます。

## 5. CASPT2・NEVPT2 を使うときの注意

CASSCF は動的相関をほとんど含まないので、その上に CASPT2 や NEVPT2 で動的相関を加えます（[量子化学の理論 第 13 回](../qc-theory/13-caspt2-nevpt2.md)）。遷移金属に使うときの注意を、第 4〜6 回の結果とあわせてまとめます。

**IPEA シフトの値を明記する**：CASPT2 には、開殻の軌道の扱いを補正する経験的な IPEA シフトがあり、標準の値は 0.25 $E_\mathrm{h}$ です（Ghigo ら 2004）。スピン状態の問題では、この値によって結果が数 kcal/mol 変わり、標準の値の CASPT2 は高スピンを安定化しすぎる傾向があります（第 5 回、第 6 回）。値を大きくすると結果がよくなる場合もありますが、その一部は誤差の打ち消しです（第 5 回）。比べたい文献と値をそろえて使います。

**NEVPT2 にはこの経験的なパラメータがない**：NEVPT2 は、活性電子どうしの 2 電子相互作用を 0 次に含む Dyall のハミルトニアンを使い、侵入者状態が起こらないように作られています（Angeli ら 2001、[量子化学の理論 第 13 回](../qc-theory/13-caspt2-nevpt2.md)）。CASPT2 の IPEA シフトにあたる調整の項はありません。ただし、パラメータがないことは、精度が高いことを意味しません。実験から作った 4 つの Fe 錯体の参照値との比較では、NEVPT2 は最大 7 kcal/mol の誤差を出し、CASPT2 より悪い結果でした（Radoń 2019）。

**3s3p の相関を確かめる**：CASPT2 の高スピンへの偏りの一因は、3s3p の準内殻相関の記述の不正確さでした（第 4 回）。この部分を CCSD(T) で置き換える CASPT2/CC が、偏りを一部取り除きます（Phung ら 2018）。

**さらに上の方法との比較**：CASPT2 の上に多参照 CI の補正を加える CASPT2+δMRCI（Reimann, Kaupp 2023）、CASPT2 と CASPT3 の平均をとる CASPT2.5（Zhang, Truhlar 2020）などの改良も提案されています。ただし、SSE17 では、どの多参照の方法も、単一参照の CCSD(T) を一貫して上回ることはありませんでした（Radoń ら 2024、第 6 回）。多参照 CI に大きさの補正（Davidson 補正など。[量子化学の理論 第 8 回](../qc-theory/08-configuration-interaction.md)）を加えた MRCI+Q では、補正の式の選び方で誤差が大きく変わりました（Radoń 2019; Radoń ら 2024）。

これらの結果から、単核の錯体で静的相関がそれほど強くない場合は、多参照の方法を使うことが、必ずしも精度の向上につながらないことが分かります。多参照の方法が本当に必要になるのは、次の節の多核の錯体のように、単一参照の方法ではうまく扱えないと考えられる系です。

**ここまでで分かったこと**：CASPT2 は IPEA シフトの値と 3s3p 相関の扱いに結果が左右され、NEVPT2 はパラメータがない代わりに、精度が高いとは限りません。単核の錯体では、多参照の方法が CCSD(T) より正確とは限りません。

## 6. DMRG の利用

金属を複数含む錯体では、活性空間に入れるべき軌道が多くなります。たとえば、鉄 4 個のクラスターでは、3d 軌道だけで 20 個になり、配位子の軌道も加えると、通常の CASSCF で扱える大きさ（20 軌道程度）を超えます。このような大きな活性空間を扱うのが、DMRG（密度行列繰り込み群）です。原理は[量子化学の理論 第 14 回](../qc-theory/14-dmrg-sci.md)で扱いました。ここでは、遷移金属への代表的な応用を紹介します。

- **光化学系 II の Mn$_4$CaO$_5$ クラスター**：Kurashige、Chan、Yanai は、DMRG で、このクラスターのほぼ厳密な多電子の波動関数を求めました（Kurashige ら 2013）。計算は、X 線結晶構造に対して提案されていた修正を支持し、構造のひずみに関係する低いエネルギーの複数のポテンシャル面を見いだしました。
- **鉄–硫黄クラスター**：Sharma らは、[2Fe-2S] と [4Fe-4S] のクラスターの低い電子状態を、模型を仮定せずに計算しました（Sharma ら 2014）。広く使われてきたハイゼンベルク–二重交換の模型は、状態の数を 1〜2 桁少なく見積もっており、その原因は Fe の d–d 励起を模型が含まないことにありました。

ただし、DMRG だけでは動的相関が足りません。Phung らは、鉄ポルフィリンのスピン状態などで、DMRG の上に CASPT2 を加えることが、遷移金属の錯体で定量的な結果を得るのに欠かせないことを示しました。計算量を減らすために 4 体の密度行列を近似する版（DMRG-cu(4)-CASPT2）は、小〜中程度の活性空間ではよい結果でしたが、30 軌道程度の大きな活性空間では精度が落ちました（Phung ら 2016）。第 6 回で紹介した Hehn らの研究では、大きな活性空間の DMRG-NEVPT2 が、CC や AFQMC と同じ結論を出しています（Hehn ら 2024）。

多核の錯体で、金属どうしの不対電子が反強磁性的に結合した状態は、単一参照の CCSD(T) ではうまく扱えないと考えられています（Radoń ら 2024）。このような系では、DMRG などの多参照の方法が主な選択肢になります。一方で、動的相関の扱い、活性空間の選び方、計算量のどれにも課題が残っており、単核の錯体ほど確立した手順はまだありません。

**ここまでで分かったこと**：DMRG は、多核の錯体のような大きな活性空間を扱え、Mn$_4$CaO$_5$ や鉄–硫黄クラスターに応用されてきました。定量的な結果には、DMRG の上に動的相関を加えることが必要です。

## まとめ

- 多参照性の診断には、CCSD の振幅（T1、D1）、エネルギー（%TAE[(T)]）、自然軌道の占有数、DFT と HF の計算（FOD、スピン対称性の破れ）にもとづくものがあり、それぞれ静的相関の違う側面を見ています。
- 診断の目安の値は系の種類で違い、診断どうしも食い違います。Cr$_2$ の計算では、T1・D1 が出発点の RHF の解で 10 倍変わりました。SCF の安定性を確かめ、違う種類の診断を複数見ます。
- 単核の錯体の活性空間は、3d 軌道、共有結合をつくる配位子の軌道、二重殻効果のための 4d の軌道が標準で、反応に沿って性格をそろえます。
- CASPT2 は IPEA シフトと 3s3p 相関の扱いに左右され、NEVPT2 も精度が高いとは限りません。単核の錯体では、多参照の方法が CCSD(T) より正確とは限りません。
- 多核の錯体には DMRG などの多参照の方法が必要で、その上に動的相関を加えます。

## 付録 A　Cr$_2$ と TiCl$_4$ の診断の計算のスクリプト

PySCF 2.15 で実行しました。1 コアで全体に約 1 分 20 秒かかりました。

```python
"""Cr2 と TiCl4 で、多参照性の診断の値を比べる（def2-SVP、Cr・Ti・Cl の 1s2s2p を凍結）。"""
import numpy as np
from scipy.linalg import eigh
from pyscf import gto, scf, cc
from pyscf.data.elements import chemcore          # 凍結する軌道の数を返す
from pyscf.cc.ccsd import get_t1_diagnostic, get_d1_diagnostic

def ccsd_diag(mf):
    mycc = cc.CCSD(mf, frozen=chemcore(mf.mol)).run()
    assert mycc.converged
    t1 = mycc.t1
    return get_t1_diagnostic(t1), get_d1_diagnostic(t1), np.abs(t1).max()

d = 2.17 / np.sqrt(3)   # Ti-Cl 2.17 Å の正四面体
systems = {
    "Cr2":   "Cr 0 0 0; Cr 0 0 1.68",
    "TiCl4": f"Ti 0 0 0; Cl {d} {d} {d}; Cl {-d} {-d} {d}; Cl {-d} {d} {-d}; Cl {d} {-d} {-d}",
}
for name, geom in systems.items():
    mol = gto.M(atom=geom, basis="def2-svp", verbose=0)
    print(f"--- {name}")
    # 1) RHF：最初に得た解と、安定性解析でたどった最も低い解
    rhf = scf.RHF(mol).run(conv_tol=1e-9)
    e_first = rhf.e_tot
    diag_first = ccsd_diag(rhf)
    while True:
        mo, stable = rhf.stability(internal=True, external=False, return_status=True)[::2]
        if stable:
            break
        rhf.kernel(rhf.make_rdm1(mo, rhf.mo_occ))
    print(f"  RHF：安定性解析でたどった解は、最初の解より {(e_first - rhf.e_tot)*627.5:.1f} kcal/mol 低い")
    for label, (T1, D1, mx) in [("最初の RHF 解", diag_first), ("安定な RHF 解", ccsd_diag(rhf))]:
        print(f"  CCSD（{label}から）: T1 = {T1:.3f}, D1 = {D1:.3f}, max|t1| = {mx:.3f}")
    # 2) スピン対称性を破った UHF を、安定性解析でたどって探す
    uhf = scf.UHF(mol).run(conv_tol=1e-9)
    while True:
        mo, stable = uhf.stability(return_status=True)[::2]
        if stable:
            break
        uhf.kernel(uhf.make_rdm1(mo, uhf.mo_occ))
    print(f"  UHF：E(UHF) - E(安定な RHF) = {(uhf.e_tot - rhf.e_tot)*627.5:.1f} kcal/mol, <S^2> = {uhf.spin_square()[0]:.2f}")
    # UHF の自然軌道の占有数：全密度行列 P を、重なり S のもとで対角化
    S = uhf.get_ovlp(); P = sum(uhf.make_rdm1())
    n = eigh(S @ P @ S, S, eigvals_only=True)[::-1]
    frac = [x for x in n if 0.02 < x < 1.98]
    print(f"  UHF 自然軌道で占有数が 0.02〜1.98 のもの: {' '.join(f'{x:.2f}' for x in frac) or 'なし'}")
```

PySCF の `get_t1_diagnostic` は、1 電子励起の振幅のノルムの 2 乗を、相関させた電子の数（占有軌道の数の 2 倍）で割った平方根を返します。`get_d1_diagnostic` は、$\mathbf t \mathbf t^\mathsf T$ と $\mathbf t^\mathsf T \mathbf t$ の最大固有値の平方根を返します。どちらも第 1.1 節の定義と同じです。`chemcore` は、各原子の凍結する軌道の数（ここでは Cr、Ti、Cl とも 1s2s2p の 5 個）の合計を返します。出力は本文の第 3 節に載せたとおりです。

## 参考文献

- C. L. Janssen, I. M. B. Nielsen, "New diagnostics for coupled-cluster and Møller–Plesset perturbation theory", Chem. Phys. Lett. 290, 423 (1998). DOI: [10.1016/S0009-2614(98)00504-1](https://doi.org/10.1016/S0009-2614(98)00504-1)
- C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J.-P. Malrieu, "Introduction of n-electron valence states for multireference perturbation theory", J. Chem. Phys. 114, 10252 (2001). DOI: [10.1063/1.1361246](https://doi.org/10.1063/1.1361246)
- G. Ghigo, B. O. Roos, P.-Å. Malmqvist, "A modified definition of the zeroth-order Hamiltonian in multiconfigurational perturbation theory (CASPT2)", Chem. Phys. Lett. 396, 142 (2004). DOI: [10.1016/j.cplett.2004.08.032](https://doi.org/10.1016/j.cplett.2004.08.032)
- A. Karton, E. Rabinovich, J. M. L. Martin, B. Ruscic, "W4 theory for computational thermochemistry: In pursuit of confident sub-kJ/mol predictions", J. Chem. Phys. 125, 144108 (2006). DOI: [10.1063/1.2348881](https://doi.org/10.1063/1.2348881)
- O. Tishchenko, J. Zheng, D. G. Truhlar, "Multireference Model Chemistries for Thermochemical Kinetics", J. Chem. Theory Comput. 4, 1208 (2008). DOI: [10.1021/ct800077r](https://doi.org/10.1021/ct800077r)
- V. Veryazov, P. Å. Malmqvist, B. O. Roos, "How to select active space for multiconfigurational quantum chemistry?", Int. J. Quantum Chem. 111, 3329 (2011). DOI: [10.1002/qua.23068](https://doi.org/10.1002/qua.23068)
- W. Jiang, N. J. DeYonker, A. K. Wilson, "Multireference Character for 3d Transition-Metal-Containing Molecules", J. Chem. Theory Comput. 8, 460 (2012). DOI: [10.1021/ct2006852](https://doi.org/10.1021/ct2006852)
- Y. Kurashige, G. K.-L. Chan, T. Yanai, "Entangled quantum electronic wavefunctions of the Mn4CaO5 cluster in photosystem II", Nat. Chem. 5, 660 (2013). DOI: [10.1038/nchem.1677](https://doi.org/10.1038/nchem.1677)
- S. Sharma, K. Sivalingam, F. Neese, G. K.-L. Chan, "Low-energy spectrum of iron–sulfur clusters directly from many-particle quantum mechanics", Nat. Chem. 6, 927 (2014). DOI: [10.1038/nchem.2041](https://doi.org/10.1038/nchem.2041)
- S. Grimme, A. Hansen, "A Practicable Real-Space Measure and Visualization of Static Electron-Correlation Effects", Angew. Chem. Int. Ed. 54, 12308 (2015). DOI: [10.1002/anie.201501887](https://doi.org/10.1002/anie.201501887)
- C. J. Stein, M. Reiher, "Automated Selection of Active Orbital Spaces", J. Chem. Theory Comput. 12, 1760 (2016). DOI: [10.1021/acs.jctc.6b00156](https://doi.org/10.1021/acs.jctc.6b00156)
- Q. M. Phung, S. Wouters, K. Pierloot, "Cumulant Approximated Second-Order Perturbation Theory Based on the Density Matrix Renormalization Group for Transition Metal Complexes: A Benchmark Study", J. Chem. Theory Comput. 12, 4352 (2016). DOI: [10.1021/acs.jctc.6b00714](https://doi.org/10.1021/acs.jctc.6b00714)
- C. A. Bauer, A. Hansen, S. Grimme, "The Fractional Occupation Number Weighted Density as a Versatile Analysis Tool for Molecules with a Complicated Electronic Structure", Chem. Eur. J. 23, 6150 (2017). DOI: [10.1002/chem.201604682](https://doi.org/10.1002/chem.201604682)
- Q. M. Phung, M. Feldt, J. N. Harvey, K. Pierloot, "Toward Highly Accurate Spin State Energetics in First-Row Transition Metal Complexes: A Combined CASPT2/CC Approach", J. Chem. Theory Comput. 14, 2446 (2018). DOI: [10.1021/acs.jctc.8b00057](https://doi.org/10.1021/acs.jctc.8b00057)
- M. Radoń, "Benchmarking quantum chemistry methods for spin-state energetics of iron complexes against quantitative experimental data", Phys. Chem. Chem. Phys. 21, 4854 (2019). DOI: [10.1039/C9CP00105K](https://doi.org/10.1039/C9CP00105K)
- C. Duan, F. Liu, A. Nandy, H. J. Kulik, "Data-Driven Approaches Can Overcome the Cost–Accuracy Trade-Off in Multireference Diagnostics", J. Chem. Theory Comput. 16, 4373 (2020). DOI: [10.1021/acs.jctc.0c00358](https://doi.org/10.1021/acs.jctc.0c00358)
- D. Zhang, D. G. Truhlar, "Spin Splitting Energy of Transition Metals: A New, More Affordable Wave Function Benchmark Method and Its Use to Test Density Functional Theory", J. Chem. Theory Comput. 16, 4416 (2020). DOI: [10.1021/acs.jctc.0c00518](https://doi.org/10.1021/acs.jctc.0c00518)
- E. Semidalas, J. M. L. Martin, "The MOBH35 Metal–Organic Barrier Heights Reconsidered: Performance of Local-Orbital Coupled Cluster Approaches in Different Static Correlation Regimes", J. Chem. Theory Comput. 18, 883 (2022). DOI: [10.1021/acs.jctc.1c01126](https://doi.org/10.1021/acs.jctc.1c01126)
- M. Reimann, M. Kaupp, "Spin-State Splittings in 3d Transition-Metal Complexes Revisited: Toward a Reliable Theory Benchmark", J. Chem. Theory Comput. 19, 97 (2023). DOI: [10.1021/acs.jctc.2c00925](https://doi.org/10.1021/acs.jctc.2c00925)
- H. Neugebauer, H. T. Vuong, J. L. Weber, R. A. Friesner, J. Shee, A. Hansen, "Toward Benchmark-Quality Ab Initio Predictions for 3d Transition Metal Electrocatalysts: A Comparison of CCSD(T) and ph-AFQMC", J. Chem. Theory Comput. 19, 6208 (2023). DOI: [10.1021/acs.jctc.3c00617](https://doi.org/10.1021/acs.jctc.3c00617)
- L. Hehn, P. Deglmann, M. Kühn, "Chelate Complexes of 3d Transition Metal Ions─A Challenge for Electronic-Structure Methods?", J. Chem. Theory Comput. 20, 4545 (2024). DOI: [10.1021/acs.jctc.3c01375](https://doi.org/10.1021/acs.jctc.3c01375)
- M. Radoń, G. Drabik, M. Hodorowicz, J. Szklarzewicz, "Performance of quantum chemistry methods for a benchmark set of spin-state energetics derived from experimental data of 17 transition metal complexes (SSE17)", Chem. Sci. 15, 20189 (2024). DOI: [10.1039/D4SC05471G](https://doi.org/10.1039/D4SC05471G)

## 関連リンク

- T. J. Lee, P. R. Taylor, "A diagnostic for determining the quality of single-reference electron correlation methods", Int. J. Quantum Chem. Symp. 23, 199 (1989)（NASA 技術報告サーバーの記録）：<https://ntrs.nasa.gov/citations/19900002238>
