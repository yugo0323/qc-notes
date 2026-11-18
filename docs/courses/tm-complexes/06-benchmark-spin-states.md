# 第 6 回　ベンチマーク（1）：スピン状態のエネルギー

!!! abstract "この回の道筋"
    計算法の良し悪しは、正しい「参照値」と比べて初めて分かります。この回では、スピン状態のエネルギー差の参照値を作る 2 つの道（実験から逆補正する道と、高精度の計算で求める道）を順に見て、それぞれの限界を確かめます。そのうえで、最近のベンチマークで各方法・汎関数がどの程度の誤差を出したかを、条件とともに整理し、現時点での実際的な勧めをまとめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 参照値を作る 2 つの道 |
    | 2 | 実験から作る参照値：逆補正の中身 |
    | 3 | 計算で作る参照値：高精度の方法とその限界 |
    | 4 | 最近のベンチマークの結果 |
    | 5 | 分かっていることと、まだ分からないこと |
    | 6 | 現時点での実際的な勧め |

    **前に読んでおく回**：[第 5 回 スピン状態のエネルギー差：なぜ計算が難しいか](05-spin-state-energetics.md)

この回でも、第 5 回と同じく、エネルギー差を「高いスピンの状態のエネルギー − 低いスピンの状態のエネルギー」で定義します。高スピンが基底状態なら負になります。単位は kcal/mol です。

## 1. 参照値を作る 2 つの道

ベンチマーク（benchmark）とは、多くの分子について、計算値を信頼できる参照値と比べ、方法ごとの誤差を調べる研究です。誤差は、平均絶対誤差（MAE。誤差の絶対値の平均）と最大誤差で表すのがふつうです。符号付きの平均誤差（MSE）は、方法がどちらかのスピン状態に偏っているかを示します。

ベンチマークの結論は、参照値の質を超えられません。スピン状態の参照値の作り方には、大きく 2 つの道があります（Radoń 2023）。

- **実験から作る**：スピンクロスオーバーのエンタルピーや、スピン禁制の吸収帯のエネルギーを測り、振動と環境の効果を計算で取り除きます（逆補正）。
- **計算で作る**：CCSD(T)、多参照の方法、量子モンテカルロ法など、最も精度が高いと考えられる方法で計算します。

第 5 回で見たように、高精度の方法どうしでも結果が食い違うことがあります。そのため、どちらの道にも長所と短所があります。次の 2 つの節で、順に中身を見ます。

**ここまでで分かったこと**：参照値は、実験から逆補正して作るか、高精度の計算で作ります。ベンチマークの結論は、参照値の質で決まります。

## 2. 実験から作る参照値：逆補正の中身

### 2.1 逆補正の式

Radoń らの SSE17 は、17 個の第一遷移系列の錯体（Fe(II)、Fe(III)、Co(II)、Co(III)、Mn(II)、Ni(II)）の実験値から作った参照値の集まりです（Radoń ら 2024）。参照値は、次の式で作ります。

$$
\Delta E_\mathrm{ref} = \Delta E_\mathrm{exptl} - \delta_\mathrm{vibr} - \delta_\mathrm{env} - \delta_\mathrm{subst}
$$

記号の意味は次のとおりです。

- $\Delta E_\mathrm{exptl}$：実験の生の値です。スピンクロスオーバーを起こす 9 個の錯体ではエンタルピー差 $\Delta H$、残りの 8 個では、最も低いスピン禁制の d–d 吸収帯の極大のエネルギーです。
- $\delta_\mathrm{vibr}$：振動の補正です。前者では零点振動と熱振動のエネルギーの差、後者では吸収帯の極大と垂直のエネルギー差のずれです。
- $\delta_\mathrm{env}$：溶媒や結晶の中にあることによる補正です。
- $\delta_\mathrm{subst}$：計算のモデルで置換基を省いたときの補正です（2 つの錯体だけ）。

こうして得る $\Delta E_\mathrm{ref}$ は、真空中の 1 分子の電子エネルギーの差です。前者からは断熱的な差（それぞれのスピン状態の最適構造どうしの差）、後者からは垂直の差（同じ構造での差）が得られます。

### 2.2 補正の大きさと、その求め方

補正の大きさの例を、SSE17 の表から抜き出します。

| 錯体 | 種類 | 環境 | $\Delta E_\mathrm{exptl}$ | $\delta_\mathrm{env}$ | $\delta_\mathrm{vibr}$ | $\Delta E_\mathrm{ref}$ |
|---|---|---|---|---|---|---|
| [Fe(HB(pz)$_3$)$_2$]（一重項–五重項） | 断熱 | CHCl$_3$ 溶液 | 5.7 | −0.2 | −1.0 | 6.9 |
| [Fe(tacn)$_2$]$^{2+}$（一重項–五重項） | 断熱 | 水など 3 種の溶媒 | 5.0〜5.7 | 1.8〜2.4 | −1.6〜−1.7 | 4.7 |
| [Fe(CN)$_6$]$^{4-}$（一重項–三重項） | 垂直 | 結晶 | 68.0 | −3.5 | −2.9 | 74.5 |
| [Fe(H$_2$O)$_6$]$^{2+}$（三重項–五重項） | 垂直 | 結晶 | −37.2 | 1.0 | −0.2 | −38.0 |

単位はすべて kcal/mol で、値は Radoń ら（2024）の Table 1 から写しました。[Fe(tacn)$_2$]$^{2+}$ は 3 つの溶媒の値を逆補正して平均したものです。

補正は、次のように求められています（Radoń ら 2024）。

- **振動の補正**：PBE0-D3(BJ)/def2-TZVP の調和振動数から求めます。垂直の差では、吸収帯の振動構造を Franck–Condon 近似でシミュレーションして求めます。SSE17 では、どの錯体でも符号が負でした。低いスピンから高いスピンに移ると、金属と配位子の結合が弱まり、振動数が下がるためです。
- **溶液の補正**：連続誘電体モデル COSMO に、実際の溶媒の誘電率を入れて求めます。水素結合が効く場合は、溶媒分子をあらわに加えます。
- **結晶の補正**：スピンクロスオーバーの錯体では、結晶全体の周期的な DFT 計算から求めます。吸収帯の錯体では、結晶構造から切り出したクラスターを CASPT2 で計算して求めます。

### 2.3 構造の選び方が補正を変える

垂直のエネルギー差は、計算に使う構造に敏感です。SSE17 では、計算に使う構造を、真空中ではなく COSMO の中で最適化した構造にしました。そのほうが結晶中の構造に近く、環境の補正が小さくなるからです（Radoń ら 2024）。

[Fe(CN)$_6$]$^{4-}$ がその例です。Fe–C 結合の長さは、真空中の最適化で 1.986 Å、結晶で 1.918 Å、COSMO の中で 1.912 Å でした。環境の補正は、COSMO の構造を基準にすると −3.5 kcal/mol ですが、真空中の構造を基準にすると 20 kcal/mol を超えます。つまり、SSE17 の垂直の参照値は、論文で配布された構造で計算したときにだけ使えます。構造を変えたら、補正を求め直す必要があります。

### 2.4 実験から作る参照値の限界

この道の限界は、次のとおりです。

- **不確かさ**：同じ錯体を複数の環境で測った場合、逆補正した値の平均からのずれは、最大で 0.7 kcal/mol でした。Radoń らは、補正と実験の誤差を合わせて、参照値の不確かさを 1〜3 kcal/mol と見積もっています。1 kcal/mol 未満の差を議論するのには使えません（Radoń ら 2024）。
- **補正そのものが計算に頼る**：振動と環境の補正は、DFT や CASPT2 で求めています。補正の大きさは数 kcal/mol に達し、前もって無視できるかどうかは分かりません。
- **使える錯体が限られる**：信頼できる実験値と結晶構造がそろった錯体は多くありません。SSE17 でも 17 個中 11 個が Fe の錯体です。
- **帰属の問題**：第 5 回で見た [Fe(H$_2$O)$_6$]$^{3+}$ のように、吸収帯の帰属が疑われる例があります（Reimann, Kaupp 2023）。SSE17 は、議論のあるこの錯体を除いています。

**ここまでで分かったこと**：実験から作る参照値は、振動・環境・置換基の補正を差し引いて得ます。不確かさは 1〜3 kcal/mol で、計算に使う構造と組にして使う必要があります。

## 3. 計算で作る参照値：高精度の方法とその限界

### 3.1 候補になる方法

計算で参照値を作るときに使われてきた方法を、簡単に紹介します。

- **CCSD(T)**：単一参照の結合クラスター法です（[量子化学の理論 第 10 回](../qc-theory/10-coupled-cluster.md)）。遷移金属では 3s3p 電子も相関させ、基底関数の極限（CBS）へ外挿する必要があります（第 4 回）。
- **CASPT2 とその改良版**：多参照の摂動論です。CASPT2/CC（3s3p の部分を CCSD(T) で置き換える、Phung ら 2018）、CASPT2+δMRCI（多参照 CI で補正する、Reimann, Kaupp 2023）、CASPT2.5（CASPT2 と CASPT3 の平均、Zhang, Truhlar 2020）などがあります。
- **拡散モンテカルロ法（DMC）**：乱数を使って、虚時間のシュレーディンガー方程式を解く方法です。波動関数の節（符号が変わる面）を試行関数で固定する「固定節近似」を使います。基底関数の不完全さの誤差がないことが利点です。
- **補助場量子モンテカルロ法（AFQMC）**：行列式の空間で乱数を使う方法です。符号の問題を避けるため、試行波動関数を使った「phaseless 近似」を入れます。

大きな錯体に使える局所相関の結合クラスター法や、計算量を抑えた多参照の方法の最近の進歩は、Feldt と Phung の総説にまとめられています（Feldt, Phung 2022）。

### 3.2 方法どうしの食い違い

これらの方法が同じ答えを出せば、参照値として安心して使えます。しかし、実際には大きく食い違うことがあります。

Mahajan らは、[Fe(H$_2$O)$_6$]$^{2+}$ の一重項と五重項のエネルギー差について、自分たちの AFQMC と文献の値を並べました（Mahajan ら 2024。2025 年に J. Chem. Theory Comput. に掲載。数値はプレプリント版の Table 5 から写しました）。次の表は、その一部です（$E_\mathrm{HS} - E_\mathrm{LS}$、kcal/mol）。

| 方法 | 基底関数など | 値 |
|---|---|---|
| AFQMC（試行関数 HF） | Fe: cc-pwCVTZ-DK、O・H: cc-pVDZ-DK | −44.2 |
| AFQMC（試行関数 CISD） | 同上 | −36.1 |
| 正準 CCSD(T) | 同上 | −37.7 |
| DLPNO-CCSD(T$_1$) | 同上 | −36.7 |
| CASPT2 | ANO-RCC triple-ζ | −46.6 |
| DMC（Song ら 2018 の値） | 基底関数の極限 | −41.0 |
| AFQMC（CISD）+ MP2 による CBS 補正 | 基底関数の極限の見積もり | −32 |
| DLPNO-CCSD(T$_1$)（Flöser ら 2020 の値） | QZ/5Z からの外挿 | −33.3 |

この表から、次のことが読み取れます。

- 同じ基底関数では、CISD を試行関数にした AFQMC と正準 CCSD(T) は 1.6 kcal/mol 以内で一致しました。一方、HF を試行関数にした AFQMC は、8 kcal/mol も離れています。AFQMC の値は、試行関数の選び方に左右されます。
- 基底関数の極限の見積もりどうし（AFQMC の −32 と、DLPNO-CCSD(T$_1$) の −33.3）はよく一致しましたが、DMC の値とは約 9 kcal/mol 違います。Mahajan らは、遷移金属で DMC を基準に使う難しさを示す例だとしています。
- 三重ゼータの基底関数の値（−36〜−37）と、極限の見積もり（−32〜−33）の間にも、数 kcal/mol の差があります。

基底関数の問題は、Flöser らも詳しく調べています。スピンクロスオーバーのモデル錯体では、5 重ゼータ（5Z）の基底関数まで使わないと極限に近づかない場合があり、TZ と QZ からの外挿では足りない場合がありました（Flöser ら 2020）。

食い違いは、この錯体だけではありません。Radoń らは、よく使われるモデル錯体 [Fe(NCH)$_6$]$^{2+}$ の一重項と五重項の差で、最良の DMC と CCSD(T) の値が 20 kcal/mol も違うことを指摘しています（Radoń ら 2024）。

### 3.3 計算で作る参照値の限界

まとめると、計算で作る参照値には、次の限界があります。

- **方法の近似**：DMC の固定節近似、AFQMC の phaseless 近似と試行関数、CASPT2 の IPEA シフト（第 5 回）、CCSD(T) の単一参照の性格など、それぞれに近似があります。どれが最も正しいかを、計算だけで決めるのは困難です。
- **基底関数の収束が遅い**：スピン状態の差は、基底関数を大きくしても、なかなか収束しません。
- **扱える大きさ**：正準 CCSD(T) や AFQMC は計算量が大きく、現実の配位子をもつ錯体には使いにくいです。SSE17 でも、正準 CCSD(T) が使えることが錯体を選ぶ条件の一つで、50 原子ほどが上限でした（Radoń ら 2024）。

一方で、方法の違う計算が一致すれば、その値はかなり信頼できます。Hehn らは、Fe(III) とニトリロ三酢酸の錯体について、CC、大きな活性空間の DMRG-NEVPT2、AFQMC という考え方の違う方法が、高スピン状態が大きく安定であることで一致したと報告しました。多くの DFT 汎関数と RPA は安定化を小さく見積もり、一部の汎関数は低スピンを基底状態と予測しました（Hehn ら 2024）。

**ここまでで分かったこと**：計算で作る参照値は、方法の近似と基底関数の収束の遅さのため、方法どうしで数〜20 kcal/mol 食い違うことがあります。考え方の違う方法が一致したときは、信頼度が上がります。

## 4. 最近のベンチマークの結果

### 4.1 SSE17：実験から作った参照値との比較

SSE17 のベンチマークは、実験から作った参照値を使い、波動関数法と DFT を同じ土俵で比べた点で、現時点で最も包括的な研究です（Radoń ら 2024）。計算の条件は次のとおりです。

- 構造：COSMO の中で最適化した PBE0-D3(BJ)/def2-TZVP の構造（すべての方法で共通）。
- 波動関数法：3s3p と価電子を相関させ、CBS の極限を見積もっています。局所相関の近似は使っていません。CASPT2 の IPEA シフトは 0.25 $E_\mathrm{h}$ です。
- DFT：def2-QZVPP 基底関数で、多くはスカラー相対論の補正を加えています。

結果を表にまとめます（MAE と最大誤差、kcal/mol）。

| 方法 | MAE | 最大誤差・偏り |
|---|---|---|
| CCSD(T)（HF 軌道） | 1.5 | 最大 −3.5、平均の符号付き誤差 −0.3 |
| CCSD(T)（PBE0 または PBE 軌道） | HF 軌道よりわずかに大きい | — |
| CASPT2/CC | 約 3 | 最大 −6 |
| CASPT2+δMRCI | 約 3 | 最大 +11 近く |
| CASPT2 | 4.1 | 最大 −7.3、平均の符号付き誤差 −3.3（高スピン寄り） |
| MRCI+Q（DSS 補正） | 約 4 | 最大 7〜9 |
| PWPB95-D3(BJ)（二重混成） | 2.4 | 6 以内 |
| B2PLYP-D3(BJ)（二重混成） | 2.8 | 6 以内 |
| LH14t-calPBE、LC-ωPBE、M06L、SSB-D など | 4 以内 | 約 8.5 以上 |
| B3LYP*-D3(BJ) | 5.1 | 14〜15 |
| TPSSh-D3(BJ) | 7.7 | 14〜15 |

この結果から、次のことが分かります。

- 単一参照の CCSD(T) が、調べた多参照の方法のどれよりも正確でした。誤差は、参照値の不確かさ（1〜3 kcal/mol）と同じ程度です。
- CCSD(T) の誤差は、多参照性の診断の値とは相関しませんでした。D1 診断が 0.15 を超える錯体が 6 割以上ありましたが、精度は落ちていませんでした。診断については[第 10 回](10-multireference-tm.md)で詳しく扱います。
- 二重混成汎関数のうち、PWPB95 と B2PLYP はよかったのですが、DSD-PBEB95 と DSD-PBEP86 は高スピンを安定化しすぎました。同じ段の汎関数でも、結果は大きく違います。
- これまでスピン状態に勧められてきた B3LYP* と TPSSh は、高スピンの錯体で低いスピンの状態を安定化しすぎ、大きな誤差を出しました。
- 分散補正を D3(BJ) から D4 に替えると、個々の錯体の値が最大 4.4 kcal/mol 変わりました。ただし、全体の成績はほぼ同じでした。

### 4.2 ほかのベンチマーク

SSE17 のほかにも、参照値の作り方が違うベンチマークがあります。主なものを、条件とともに挙げます。

| 研究 | 参照値 | 対象 | 主な結論 |
|---|---|---|---|
| Radoń 2019 | 実験から逆補正 | Fe 錯体 4 個 | CCSD(T)：最大誤差 2 未満、MAE 1 未満。CASPT2：高スピンを最大 5.5 安定化しすぎ。NEVPT2：最大 7。MRCI+Q は大きさの補正の式で MAE が 3 未満〜7 と変わる。すべてで釣り合う汎関数は少数（B2PLYP-D3、OPBE など） |
| Cirera ら 2018 | 実験のスピンクロスオーバーの挙動 | 第一遷移系列の錯体 | TPSSh（三重ゼータ基底）が、金属と酸化数によらず最もよい |
| Zhang, Truhlar 2020 | CASPT2.5（CBS 外挿） | FeF$_2$、FeCl$_2$、CoF$_2$ | 60 個の汎関数のうち、MN15-L、M06-SX、revM06 がよく一致。基底関数は def2-TZVP を推奨 |
| Neugebauer ら 2023（3dTMV） | ph-AFQMC | 28 個の 3d 金属の分子の垂直イオン化エネルギー（def2-SVP） | スピン対称性の破れで分けた「安全な領域」では、CCSD(T) の AFQMC からのずれは平均 2 kcal/mol 程度以下。その外では CCSD(T) が失敗しうる |

最後の 3dTMV はスピン状態ではなくイオン化エネルギーのベンチマークですが、「CCSD(T) をどこまで参照値に使えるか」を AFQMC で調べた研究として重要です。

表の Cirera らと SSE17 では、TPSSh の評価が逆になっています。前者は「スピンクロスオーバーが起こるかどうか、起こる温度の範囲を正しく予測できるか」を、後者は「逆補正した電子エネルギーの差を何 kcal/mol の誤差で再現するか」を見ています。物差しが違えば、汎関数の順位も変わります。ベンチマークの結論を使うときは、何を参照値にして、何を測ったのかを必ず確かめてください。

**ここまでで分かったこと**：SSE17 では、CCSD(T) が MAE 1.5 kcal/mol で最もよく、二重混成汎関数の一部が MAE 3 kcal/mol 未満でした。従来勧められてきた汎関数は大きな誤差を出しました。ただし、汎関数の評価は、参照値と物差しの選び方で変わります。

## 5. 分かっていることと、まだ分からないこと

ここまでの文献から、分かっていることと、まだ分からないことを分けて整理します。

**分かっていること**

- 単核の第一遷移系列の錯体では、3s3p を相関させ、基底関数の極限まで収束させた CCSD(T) が、実験から作った参照値を 1〜3 kcal/mol で再現します（Radoń 2019; Radoń ら 2024）。
- 標準の IPEA シフトの CASPT2 は、高スピンを安定化しすぎる傾向があります（Radoń 2019; Radoń ら 2024）。
- 近似の汎関数には、どの錯体でも正確な「万能の汎関数」はありません。厳密交換の割合の少ない混成汎関数（B3LYP*、TPSSh）も、SSE17 の範囲では万能ではありませんでした。

**まだ分からないこと**

- 多核の錯体：金属どうしの反強磁性的な結合を含むスピン状態は、単一参照の CCSD(T) ではうまく扱えないと考えられます。SSE17 の結論を多核の錯体に広げてはいけない、と著者自身が述べています（Radoń ら 2024）。
- 第二・第三遷移系列の金属、Mn・Co・Ni の錯体：参照値がまだ少なく、Fe に偏っています。
- 高精度の方法どうしの食い違い：DMC、AFQMC、CCSD(T) の値が、どの条件で一致し、どの条件で離れるかは、まだ研究が続いています（Mahajan ら 2024; Neugebauer ら 2023）。

## 6. 現時点での実際的な勧め

ここまでの結果から、スピン状態のエネルギーを計算するときの進め方をまとめます。これは、特定の方法を「正解」とするものではなく、文献の結果から言える範囲の勧めです。

1. **比べる量をそろえる**：断熱か垂直か、電子エネルギーかエンタルピーか、真空中か溶液中かを、まず決めます。実験と比べるなら、第 2 節の逆補正が必要です。
2. **小さなモデルで CCSD(T) を使う**：可能なら、配位子を簡単にしたモデルで、3s3p を相関させた CCSD(T) を計算し、基底関数の極限を見積もります。大きな系では DLPNO-CCSD(T) が使えますが、閾値と三重励起の扱いに注意が必要です（[第 9 回](09-dlpno-for-tm.md)）。Radoń らは、DFT で反応を調べるときに、簡単なモデルの CCSD(T) のスピン状態のエネルギーで補正する方法を、現実的な解決策の一つとして挙げています（Radoń ら 2024）。
3. **DFT は複数の汎関数で確かめる**：二重混成汎関数（PWPB95、B2PLYP など）を含め、厳密交換の割合の違う汎関数で計算し、スピン状態の順序が変わらないかを確かめます。厳密交換の割合を変えたときの傾き（第 5 回）も、結論の頑健さの目安になります。汎関数の選び方は[第 8 回](08-functional-choice.md)で扱います。
4. **多参照の方法を使うなら、偏りを知っておく**：CASPT2 を使うときは、IPEA シフトの値と活性空間を明記し、高スピン寄りの偏りがありうることを考えに入れます。
5. **多核の錯体は別に考える**：単核の錯体のベンチマークの結果を、そのまま当てはめないでください。

**ここまでで分かったこと**：比べる量をそろえ、可能なら CCSD(T) で確かめ、DFT は複数の汎関数で結論の頑健さを確かめるのが、現時点での実際的な進め方です。

## まとめ

- スピン状態の参照値は、実験から振動・環境の効果を逆補正して作るか、高精度の計算で作ります。前者の不確かさは 1〜3 kcal/mol で、計算に使う構造と組にして使います。
- 高精度の計算どうしでも、近似（固定節、phaseless、試行関数、IPEA シフト）と基底関数の収束の遅さのため、数〜20 kcal/mol 食い違うことがあります。
- SSE17 では、CCSD(T) が MAE 1.5 kcal/mol で最もよく、CASPT2 は高スピン寄りに偏り、二重混成汎関数の PWPB95 と B2PLYP が MAE 3 kcal/mol 未満でした。
- 汎関数の評価は、参照値と物差しの選び方で変わります。多核の錯体や第二・第三遷移系列は、まだ参照値が足りません。

## 参考文献

- M. Radoń, "Benchmarking quantum chemistry methods for spin-state energetics of iron complexes against quantitative experimental data", Phys. Chem. Chem. Phys. 21, 4854 (2019). DOI: [10.1039/C9CP00105K](https://doi.org/10.1039/C9CP00105K)
- M. Radoń, "Benchmarks for transition metal spin-state energetics: why and how to employ experimental reference data?", Phys. Chem. Chem. Phys. 25, 30800 (2023). DOI: [10.1039/D3CP03537A](https://doi.org/10.1039/D3CP03537A)
- M. Radoń, G. Drabik, M. Hodorowicz, J. Szklarzewicz, "Performance of quantum chemistry methods for a benchmark set of spin-state energetics derived from experimental data of 17 transition metal complexes (SSE17)", Chem. Sci. 15, 20189 (2024). DOI: [10.1039/D4SC05471G](https://doi.org/10.1039/D4SC05471G)
- Q. M. Phung, M. Feldt, J. N. Harvey, K. Pierloot, "Toward Highly Accurate Spin State Energetics in First-Row Transition Metal Complexes: A Combined CASPT2/CC Approach", J. Chem. Theory Comput. 14, 2446 (2018). DOI: [10.1021/acs.jctc.8b00057](https://doi.org/10.1021/acs.jctc.8b00057)
- S. Song, M.-C. Kim, E. Sim, A. Benali, O. Heinonen, K. Burke, "Benchmarks and Reliable DFT Results for Spin Gaps of Small Ligand Fe(II) Complexes", J. Chem. Theory Comput. 14, 2304 (2018). DOI: [10.1021/acs.jctc.7b01196](https://doi.org/10.1021/acs.jctc.7b01196)
- J. Cirera, M. Via-Nadal, E. Ruiz, "Benchmarking Density Functional Methods for Calculation of State Energies of First Row Spin-Crossover Molecules", Inorg. Chem. 57, 14097 (2018). DOI: [10.1021/acs.inorgchem.8b01821](https://doi.org/10.1021/acs.inorgchem.8b01821)
- B. M. Flöser, Y. Guo, C. Riplinger, F. Tuczek, F. Neese, "Detailed Pair Natural Orbital-Based Coupled Cluster Studies of Spin Crossover Energetics", J. Chem. Theory Comput. 16, 2224 (2020). DOI: [10.1021/acs.jctc.9b01109](https://doi.org/10.1021/acs.jctc.9b01109)
- D. Zhang, D. G. Truhlar, "Spin Splitting Energy of Transition Metals: A New, More Affordable Wave Function Benchmark Method and Its Use to Test Density Functional Theory", J. Chem. Theory Comput. 16, 4416 (2020). DOI: [10.1021/acs.jctc.0c00518](https://doi.org/10.1021/acs.jctc.0c00518)
- M. Feldt, Q. M. Phung, "Ab Initio Methods in First-Row Transition Metal Chemistry", Eur. J. Inorg. Chem. 2022, e202200014 (2022). DOI: [10.1002/ejic.202200014](https://doi.org/10.1002/ejic.202200014)
- M. Reimann, M. Kaupp, "Spin-State Splittings in 3d Transition-Metal Complexes Revisited: Toward a Reliable Theory Benchmark", J. Chem. Theory Comput. 19, 97 (2023). DOI: [10.1021/acs.jctc.2c00925](https://doi.org/10.1021/acs.jctc.2c00925)
- H. Neugebauer, H. T. Vuong, J. L. Weber, R. A. Friesner, J. Shee, A. Hansen, "Toward Benchmark-Quality Ab Initio Predictions for 3d Transition Metal Electrocatalysts: A Comparison of CCSD(T) and ph-AFQMC", J. Chem. Theory Comput. 19, 6208 (2023). DOI: [10.1021/acs.jctc.3c00617](https://doi.org/10.1021/acs.jctc.3c00617)
- L. Hehn, P. Deglmann, M. Kühn, "Chelate Complexes of 3d Transition Metal Ions─A Challenge for Electronic-Structure Methods?", J. Chem. Theory Comput. 20, 4545 (2024). DOI: [10.1021/acs.jctc.3c01375](https://doi.org/10.1021/acs.jctc.3c01375)
- A. Mahajan, J. H. Thorpe, J. S. Kurian, D. R. Reichman, D. A. Matthews, S. Sharma, "Beyond CCSD(T) accuracy at lower scaling with auxiliary field quantum Monte Carlo", arXiv:2410.02885 (2024).
