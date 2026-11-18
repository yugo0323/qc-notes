# 第 9 回　局所相関法 DLPNO-CCSD(T) を遷移金属に使うときの注意

!!! abstract "この回の道筋"
    DLPNO-CCSD(T) は、大きな遷移金属錯体に CCSD(T) 級の計算を届かせる方法として広く使われています。しかし第 7 回で見たように、遷移金属では正準 CCSD(T) から大きくずれることがあります。この回では、そのずれの原因を、PNO の閾値、三重励起の扱い、3s3p の準内殻相関、開殻の参照、多参照性、基底関数の順に一つずつ取り上げ、それぞれの対策と、正準 CCSD(T) と比べた研究の結果をまとめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 誤差を 3 つに分ける |
    | 2 | PNO の閾値の誤差と、その外挿 |
    | 3 | 三重励起：(T$_0$) と (T$_1$) |
    | 4 | 3s3p の準内殻相関と凍結内殻 |
    | 5 | 開殻の錯体：参照の行列式 |
    | 6 | T1・D1 などで多参照性を確かめる |
    | 7 | 基底関数の極限への外挿 |
    | 8 | 正準 CCSD(T) と比べた研究のまとめと、計算の進め方 |

    **前に読んでおく回**：[量子化学の理論 第 11 回 局所相関法：DLPNO-CCSD(T)](../qc-theory/11-dlpno-ccsd-t.md)、[第 4 回 遷移金属の電子相関](04-electron-correlation-tm.md)、[第 7 回 ベンチマーク（2）](07-benchmark-reactions.md)

## 1. 誤差を 3 つに分ける

DLPNO-CCSD(T) の考え方と、閾値の組（LoosePNO、NormalPNO、TightPNO）は、[量子化学の理論 第 11 回](../qc-theory/11-dlpno-ccsd-t.md)で扱いました。ここでは、それを遷移金属に使うときの注意に絞ります。

DLPNO-CCSD(T) で求めたエネルギーと、本当の値との差は、次の 3 つに分けて考えると整理しやすくなります。

1. **局所近似の誤差**：同じ基底関数の正準 CCSD(T) からのずれです。PNO の打ち切り、電子対のふるい分け、ドメイン（相関を考える原子の範囲）の近似から来ます。
2. **基底関数の誤差**：有限の基底関数を使うことによる、基底関数の極限（CBS）からのずれです。
3. **CCSD(T) そのものの誤差**：静的相関が強いと大きくなります。

この回で主に扱うのは 1 で、2 は第 7 節で扱います。3 は DLPNO にしても直らないので、第 6 節と[第 10 回](10-multireference-tm.md)で扱う診断で確かめます。遷移金属で難しいのは、静的相関が強い系ほど 1 も大きくなりやすく、1 と 3 が重なって現れることです（Semidalas, Martin 2022）。

**ここまでで分かったこと**：DLPNO-CCSD(T) の誤差は、局所近似、基底関数、CCSD(T) 自体の 3 つに分けて考えます。遷移金属では、局所近似の誤差と CCSD(T) 自体の誤差が重なりやすくなります。

## 2. PNO の閾値の誤差と、その外挿

### 2.1 反応エネルギーと障壁での大きさ

Altun らは、MOBH35（第 7 回）の反応エネルギーと障壁について、DLPNO-CCSD(T)/def2-SVP の結果を正準 CCSD(T)/def2-SVP と比べました（Altun ら 2023）。PNO を残すかどうかの閾値 $T_\mathrm{CutPNO} = 10^{-X}$ を変えたときの平均絶対誤差（MAE）は、次のとおりです（TightPNO の他の閾値と組み合わせた値、kcal/mol）。

| $T_\mathrm{CutPNO}$ | 従来の設定 | 準内殻の閾値を厳しくした設定 | さらに軌道緩和の補正を加えた設定 |
|---|---|---|---|
| $10^{-6}$ | 1.08 | 0.96 | 0.60 |
| $10^{-7}$ | 0.71 | 0.62 | 0.37 |
| $10^{-8}$ | 0.47 | 0.37 | 0.25 |
| $10^{-9}$ | 0.31 | 0.25 | 0.18 |
| CPS(6/7) | 0.56 | 0.49 | 0.35 |
| CPS(7/8) | 0.37 | 0.29 | 0.22 |

表の 3 列目と 4 列目の設定は、第 4 節と第 6 節で説明します。CPS($X$/$Y$) は、2 つの閾値の結果から、PNO を打ち切らない極限へ外挿した値です（外挿の式は[量子化学の理論 第 11 回](../qc-theory/11-dlpno-ccsd-t.md)）。

平均では、TightPNO の標準の $10^{-7}$ で 1 kcal/mol 未満です。しかし、平均が小さくても、特定の反応で大きくずれることがあります。第 7 回の MOBH35 の反応 9 では、DLPNO-CCSD(T) と DLPNO-CCSD(T$_1$) が大きく破綻しました（Semidalas, Martin 2022）。

### 2.2 スピン状態での大きさ

スピン状態のエネルギー差では、PNO の閾値の影響がさらに大きくなることがあります。理由は、スピン状態によって、正準の相関エネルギーのうち DLPNO で取り込める割合が違うことです。割合の違いが、そのままエネルギー差の誤差になります。

Drosou らは、鉄の錯体のスピン状態のエネルギー差を、CASPT2/CC（第 4 回）と正準 CCSD(T) の値と比べました（Drosou ら 2022）。主な結果は次のとおりです。

- ある Fe(III) の錯体では、取り込める相関エネルギーの割合がスピン状態ごとに違い、既定の NormalPNO の設定で、CCSD の部分だけで高スピンが 3.8 と 7.5 kcal/mol 有利になりすぎました。割合の違いは TightPNO でも見られました。
- 閾値による誤差のうち、電子対のふるい分けの誤差よりも、ドメインの誤差のほうがずっと大きく、誤差は主に $T_\mathrm{CutPNO}$ で決まっていました。
- TightPNO を基本にして、$T_\mathrm{CutPNO} = 10^{-6}$ と $10^{-7}$ の 2 点で CPS 外挿すると、スピン状態の間の取り込み割合の差がほぼなくなり、CCSD の部分の誤差は −0.4 と −0.8 kcal/mol に下がりました。
- 計算量を抑えるため、NormalPNO を基本にして $10^{-6}$ と $3.33 \times 10^{-7}$（NormalPNO の既定値）の 2 点で外挿する方法も提案されています。誤差の増え方は 1.5 kcal/mol 未満でした。

著者らは、スピン状態のエネルギーに DLPNO を使うときは、CPS 外挿を標準の手順とみなすべきだと述べています。

**ここまでで分かったこと**：PNO の閾値の誤差は、遷移金属の反応では平均 1 kcal/mol 程度ですが、スピン状態の差では、状態ごとの相関エネルギーの取り込み割合の違いから数 kcal/mol に達します。2 点の閾値からの CPS 外挿が有効です。

## 3. 三重励起：(T$_0$) と (T$_1$)

DLPNO-CCSD(T) の三重励起の補正には、2 つの版があります（Guo ら 2018）。

- **(T$_0$)**：初期の版です。占有軌道の Fock 行列の非対角要素を無視する「半正準近似」を使います。
- **(T$_1$)**：三重励起の振幅を反復で求める改良版です（論文では DLPNO-CCSD(T) と呼ばれていますが、区別のため (T$_1$) と書くのがふつうです）。

Guo らによると、多くの場合、相対エネルギーは (T$_0$) でも正準の (T) と同じくらい正確です。しかし、HOMO と LUMO のエネルギー差が小さい系などでは、(T$_0$) が破綻することがあります。(T$_1$) の計算時間は、(T$_0$) のおよそ 2 倍です（基底関数による）。開殻の系にも、両方の版が作られています。開殻での (T$_1$) の計算時間は、(T$_0$) の平均 4 倍でした（Guo ら 2020）。

遷移金属のスピン状態では、この違いが大きく効きます。

- Drosou らの鉄の錯体では、(T$_0$) は CASPT2/CC と比べて高スピンを 4〜17 kcal/mol 安定化しすぎました。(T$_1$) に替えると、ずれが 2〜5 kcal/mol 小さくなりました。著者らは、スピン状態のエネルギーには (T$_1$) だけを使うべきだと結論しています（Drosou ら 2022）。
- Flöser らも、スピンクロスオーバーのモデル錯体で、(T$_1$) と TightPNO を使えば正準 CCSD(T) を化学的精度（約 1 kcal/mol）で再現できることを示しました。一方で、化学反応の相対エネルギーの多くでは (T$_0$) で十分なことも述べ、迷ったら両方を比べることを勧めています（Flöser ら 2020）。

両者の差そのものも、役に立つ情報です。Semidalas と Martin は、DLPNO-CCSD(T$_1$) と DLPNO-CCSD(T$_0$) の差の大きさが、DLPNO-CCSD(T$_1$) の正準 CCSD(T) からの誤差のよい目安になることを見いだしました（Semidalas, Martin 2022）。差が大きいときは、局所近似が難しい系だと考えます。

**ここまでで分かったこと**：遷移金属のスピン状態では、半正準近似の (T$_0$) は高スピンに大きく偏るので、(T$_1$) を使います。(T$_1$) と (T$_0$) の差の大きさは、DLPNO の誤差の目安にもなります。

## 4. 3s3p の準内殻相関と凍結内殻

### 4.1 どこまで凍結するか

第 4 回で見たように、3d 金属では 3s3p 電子が 3d 電子と同じ領域にあるので、3s3p も相関させる必要があります。凍結内殻近似で凍結するのは、1s2s2p までにします。

ORCA は、この点を考えた既定の凍結内殻を使います。Altun らの説明によると、軌道エネルギーが −200 eV より低い軌道を凍結し、−80 eV より高い軌道を相関させるという基準で、各元素の凍結する軌道（「化学的な内殻」）を決めています（Altun ら 2023）。この基準は、Bistoni らが、軌道エネルギーを使って凍結してよい軌道を系統的に決めた研究に基づいています（Bistoni ら 2017）。

この既定を確かめるため、ORCA 6.1 で Fe 原子の小さな MP2 計算（UHF、def2-SVP）をしたところ、出力には次の 1 行がありました。

```text
Freezing NCore=10 chemical core electrons
```

凍結されたのは 10 個の電子、つまり 1s2s2p だけで、3s3p の 8 個の電子は相関しています。ほかのプログラムでは既定が違うことがあるので、出力で必ず確かめてください。

### 4.2 準内殻の電子対には厳しい閾値が要る

3s3p を相関させると、DLPNO には別の問題が出ます。PNO は、局所 MP2 の見積もりから作ります。3s3p を含む電子対の相関エネルギーは、価電子どうしの電子対と同じ閾値では、十分に取り込めないのです。

Bistoni らは、内殻の電子を含む電子対では、より厳しい閾値を使わないと、全電子の DLPNO-CCSD が凍結内殻の計算と同じ精度にならないことを示しました（Bistoni ら 2017）。Altun らは、さらに 3s3p の準内殻について調べ、次のことを示しました（Altun ら 2023）。

- Zn 原子で、準内殻どうし、準内殻と価電子の電子対の相関エネルギーは、$T_\mathrm{CutPNO} = 10^{-8}$ より緩い閾値では、打ち切らない極限から大きくずれました。
- 準内殻を含む電子対だけ $T_\mathrm{CutPNO}$ を 100 倍厳しくすると、このずれは消えました。
- この「準内殻の閾値を厳しくした設定」は、MOBH35 で精度を上げ（第 2 節の表の 3 列目）、しかも計算時間を短くしました。

### 4.3 基底関数も 3s3p 用のものが要る

3s3p を相関させても、基底関数がその領域を表せなければ、相関は正しく取り込めません。Semidalas と Martin は、def2 系の基底関数には、遷移金属の内殻–価電子相関用の関数が含まれていないことを指摘しています（Semidalas, Martin 2022）。Balabanov と Peterson の cc-pwCV$n$Z（第 4 回）は、そのための関数を含みます。SSE17 の正準 CCSD(T) でも、金属には cc-pwCVTZ を使っていました（第 6 回）。

Flöser らは、スピンクロスオーバーのモデル錯体で、金属の外側の内殻の相関を含めることが必要であり、さらにスカラー相対論の効果（2 次の Douglas–Kroll–Hess 法）と TightPNO の閾値が、それぞれ約 3 kcal/mol 結果に効くと報告しています（Flöser ら 2020）。

**ここまでで分かったこと**：3d 金属では 1s2s2p だけを凍結し、3s3p を相関させます。ORCA の既定はそうなっていますが、出力で確かめます。準内殻を含む電子対には厳しい PNO の閾値が、基底関数には cc-pwCV$n$Z のような内殻–価電子用の関数が必要です。

## 5. 開殻の錯体：参照の行列式

### 5.1 開殻の DLPNO

遷移金属の錯体の多くは開殻です。開殻の DLPNO-CCSD は、Saitow らが作り直した版が使われています（Saitow ら 2017）。要旨によると、この方法は高スピンの開殻（不対電子のスピンがすべてそろった状態）を対象にし、PNO の作り方を閉殻の場合と厳密に対応させています。既定の閾値で、開殻でも正準 CCSD の相関エネルギーの約 99.9% を取り込めます。

開殻の参照の行列式について、ORCA 6.1 で OH ラジカルの小さな計算（UHF、DLPNO-CCSD(T$_1$)）をしたところ、出力には「Quasi restricted orbitals will be used (this is similar to a ROHF determinant)」とありました。UHF の計算から、準制限軌道（QRO）という、ROHF に似た 1 組の軌道を作り、それを相関計算の出発点にしていることが分かります。UHF の α と β の軌道を別々に使うのではありません。

1 つの行列式で表せない状態（第 3 回の開殻一重項など）は、単一参照の CCSD(T) の対象外です。DLPNO にしても、この点は変わりません。

### 5.2 参照の軌道の選び方

参照の行列式に、HF の軌道を使うか、Kohn–Sham（KS）の軌道を使うかでも、結果が変わることがあります。

- Drosou らは、鉄の錯体の DLPNO-CCSD(T$_1$) で、BP86 の軌道と B3LYP の軌道を比べ、CASPT2/CC からの平均の符号付き誤差が −7.0 と −6.5 kcal/mol と、近い値になりました（CPS 外挿の前）。CPS 外挿により、参照の軌道への依存も小さくなりました（Drosou ら 2022）。
- Radoń らは、正準 CCSD(T) では、HF の軌道から KS の軌道に替えても精度は一貫しては上がらなかったと報告しています。そして、文献で KS の軌道のほうがよいとされた結果の多くは DLPNO によるもので、DLPNO の精度が参照の軌道に左右され、HF の軌道で落ちることがあるためだろう、と述べています（Radoń ら 2024）。

つまり、参照の軌道による違いの一部は、CCSD(T) の性質ではなく、局所近似の誤差です。CPS 外挿で局所近似の誤差を減らすと、この違いも小さくなります。

**ここまでで分かったこと**：開殻の DLPNO は、UHF から作った準制限軌道（QRO）を出発点にします。参照の軌道（HF か KS か）による違いの一部は局所近似の誤差で、CPS 外挿で小さくなります。

## 6. T1・D1 などで多参照性を確かめる

DLPNO の誤差は、静的相関が強い系で大きくなります。Altun らは、その理由を次のように説明しています（Altun ら 2023）。PNO は局所 MP2 の見積もりから作りますが、静的相関が強いと、CCSD の段階で軌道が大きく変わります（論文では「動的相関による軌道緩和」と呼んでいます）。その変化を MP2 の見積もりは捉えられないので、PNO の空間が不十分になります。

軌道の変化の大きさは、CCSD の 1 電子励起の振幅 $t_i^a$ に現れます。そこで、1 電子励起の振幅を使う診断が役に立ちます。

- **T1 診断**：T1 が大きいほど、PNO の質が落ちます（Altun ら 2023）。T1 が 0.02 を超える系では、PNO 空間の大きさに対する収束が遅いことが報告されています（[量子化学の理論 第 11 回](../qc-theory/11-dlpno-ccsd-t.md)）。
- **振幅のベクトルの長さの 2 乗** $\lVert t_1 \rVert^2$：T1 は分子の大きさによらない量ですが、DLPNO の誤差は分子が大きいほど増えます。Altun らは、誤差の大きさを見積もるには $\lVert t_1 \rVert^2$ のほうが向くとして、これに比例する補正（第 2 節の表の 4 列目）を提案しました。
- **組み合わせ**：Semidalas と Martin は、T1、D1、1 電子励起の振幅の最大値、HOMO–LUMO ギャップの逆数をあわせて見れば、問題が起こりそうな系を十分に見分けられると述べています（Semidalas, Martin 2022）。

ORCA の DLPNO の出力には、T1 診断の値が表示されます。これらの診断の定義と限界は、[第 10 回](10-multireference-tm.md)で詳しく扱います。診断の値が大きいときは、$T_\mathrm{CutPNO}$ を変えて収束を確かめるか、CPS 外挿を使います。それでも結果が安定しなければ、CCSD(T) そのものが使える範囲を超えている可能性を考えます。

**ここまでで分かったこと**：静的相関が強いと、CCSD での軌道の変化を MP2 の見積もりが捉えられず、DLPNO の誤差が大きくなります。T1、D1、振幅の最大値などで、問題が起こりそうな系を見分けます。

## 7. 基底関数の極限への外挿

スピン状態のエネルギー差は、基底関数に対する収束がとても遅い量です（第 6 回）。DLPNO を使えば、正準 CCSD(T) では届かない大きな基底関数を使えるので、外挿の精度を上げられます。

相関エネルギーの外挿には、基底関数の基数 $L$（TZ なら 3、QZ なら 4）を使った次の形の式がよく使われます。

$$
E(L) = E_\mathrm{CBS} + A\,L^{-\alpha}
$$

ここで、$E_\mathrm{CBS}$ が求めたい極限の値、$A$ と $\alpha$ は定数です。2 つの $L$ の結果から $E_\mathrm{CBS}$ と $A$ を決めます。$\alpha$ の値は基底関数の系列で違います。def2 系については、def2-SVP と def2-TZVPP の組で $\alpha = 2.40$、def2-TZVPP と def2-QZVPP の組で $\alpha = 2.97$ が使われています（Semidalas と Martin が、Neese と Valeev の表から引用した値）。HF の部分は、別の形の式で外挿します（Semidalas, Martin 2022）。

遷移金属で外挿するときの注意は、次のとおりです。

- **TZ/QZ では足りないことがある**：Flöser らのスピンクロスオーバーのモデル錯体では、QZ/5Z の外挿でおおむね収束しましたが、TZ/QZ の外挿では、すべての場合に満足な結果にはなりませんでした（Flöser ら 2020）。
- **外挿による変化は小さくない**：[Fe(H$_2$O)$_6$]$^{2+}$ では、DLPNO-CCSD(T$_1$) の一重項と五重項の差が、TZ の値から QZ/5Z の外挿値へ約 3 kcal/mol 変わりました（Mahajan ら 2024）。
- **金属だけ大きくする**：Drosou らは、Fe の原子の基底関数だけを TZ と QZ にして外挿する方法を、CPS 外挿と組み合わせました。この組み合わせで、CASPT2/CC との差は −2.4〜2.7 kcal/mol、平均の符号付き誤差は 0.1 kcal/mol になりました（Drosou ら 2022）。

**ここまでで分かったこと**：スピン状態の差は基底関数への収束が遅く、TZ/QZ の外挿では足りないことがあります。DLPNO の利点を生かして大きな基底関数を使い、金属の基底関数を重点的に大きくする方法も有効です。

## 8. 正準 CCSD(T) と比べた研究のまとめと、計算の進め方

### 8.1 研究のまとめ

DLPNO を正準 CCSD(T)（または、それに準じる基準）と比べた研究を、条件とともにまとめます。

| 研究 | 対象 | 基準 | 主な結果 |
|---|---|---|---|
| Flöser ら 2020 | スピンクロスオーバーのモデル錯体 3 個 | 正準 CCSD(T) | (T$_1$)、TightPNO、相対論効果を入れれば化学的精度で一致 |
| Drosou ら 2022 | 鉄の錯体のスピン状態 | CASPT2/CC と正準 CCSD(T) | (T$_0$) は高スピンを 4〜17 kcal/mol 有利にしすぎ。(T$_1$)、CPS 外挿、Fe の CBS 外挿で −2.4〜2.7 |
| Semidalas, Martin 2022 | MOBH35 の大部分（def2-SVP） | 正準 CCSD(T) | 静的相関が弱ければよく一致。中程度以上で DLPNO は大きく破綻しうる（反応 9）。LNO-CCSD(T) は頑健 |
| Altun ら 2023 | MOBH35（def2-SVP） | 正準 CCSD(T) | TightPNO・$10^{-7}$ で MAE 0.71。準内殻の閾値の強化と軌道緩和の補正で 0.37 |
| Mahajan ら 2024 | [Fe(H$_2$O)$_6$]$^{2+}$（Fe: cc-pwCVTZ-DK、O・H: cc-pVDZ-DK） | 正準 CCSD(T) | DLPNO-CCSD(T$_1$) と正準 CCSD(T) の差は 1.0 kcal/mol（−36.7 と −37.7） |

表から言えるのは、適切な設定を選べば DLPNO は正準 CCSD(T) を 1〜2 kcal/mol 程度で再現できる一方、既定の設定のままでは、特にスピン状態と静的相関の強い系で、数 kcal/mol 以上ずれうるということです。

### 8.2 計算の進め方

ここまでの注意を、計算の手順にまとめます。

1. **凍結内殻を確かめる**：3d 金属では 3s3p を相関させます。出力で凍結した電子の数を確かめます。
2. **(T$_1$) を使う**：スピン状態の差では必ず (T$_1$) を使います。反応エネルギーでも、(T$_1$) と (T$_0$) の差が大きければ注意します。
3. **閾値を 2 つ計算して外挿する**：TightPNO を基本にして、$T_\mathrm{CutPNO}$ だけを $10^{-6}$ と $10^{-7}$ にした 2 つの計算をし、CPS 外挿します。計算量を抑えるなら、NormalPNO を基本にした $10^{-6}$ と $3.33 \times 10^{-7}$ の組も使えます（Drosou ら 2022）。
4. **診断を見る**：T1、D1、(T$_1$) と (T$_0$) の差を確かめます。大きければ、閾値への収束をさらに確かめ、多参照性を疑います（[第 10 回](10-multireference-tm.md)）。
5. **基底関数の極限を見積もる**：金属には cc-pwCV$n$Z のような内殻–価電子用の基底関数を使い、少なくとも 2 つの大きさの基底関数から外挿します。
6. **相対論効果を入れる**：第一遷移系列でも、スカラー相対論の効果が数 kcal/mol 効くことがあります（Flöser ら 2020）。相対論効果は第 11 回（準備中）で扱います。

手順 3 の 2 つの計算は、ORCA では次のように書けます。1 つ目は TightPNO の既定（$T_\mathrm{CutPNO} = 10^{-7}$）のまま、2 つ目は `%mdci` のブロックで $T_\mathrm{CutPNO}$ だけを $10^{-6}$ に変えます。この書き方で、$T_\mathrm{CutPNO}$ 以外の閾値が TightPNO のまま残ることを、ORCA 6.1 の小さな計算の出力で確かめました。基底関数の名前は例です。

```text
! UHF DLPNO-CCSD(T1) TightPNO def2-TZVPP def2-TZVPP/C
%mdci
  TCutPNO 1e-6
end
```

2 つの計算の相関エネルギーを $E_6$、$E_7$ とすると、CPS(6/7) の値は $E_6 + 1.5\,(E_7 - E_6)$ です（[量子化学の理論 第 11 回](../qc-theory/11-dlpno-ccsd-t.md)）。外挿は相関エネルギーに対して行い、HF のエネルギーは足し戻します。

**ここまでで分かったこと**：適切な設定（3s3p の相関、(T$_1$)、CPS 外挿、内殻–価電子用の基底関数）を選べば、DLPNO は遷移金属でも正準 CCSD(T) を 1〜2 kcal/mol 程度で再現できます。既定の設定のまま使うと、数 kcal/mol 以上ずれることがあります。

## まとめ

- DLPNO-CCSD(T) の誤差は、局所近似、基底関数、CCSD(T) 自体の 3 つに分けて考えます。遷移金属では、静的相関が強いほど局所近似の誤差も大きくなります。
- スピン状態の差では、状態ごとの相関エネルギーの取り込み割合の違いが誤差になります。(T$_1$) を使い、2 つの $T_\mathrm{CutPNO}$ から CPS 外挿します。
- 3d 金属では 3s3p を相関させ、準内殻の電子対には厳しい閾値、基底関数には内殻–価電子用の関数を使います。
- 開殻では UHF から作った QRO を出発点にします。T1、D1、(T$_1$) と (T$_0$) の差で、問題が起こりそうな系を見分けます。
- スピン状態の差は基底関数への収束が遅く、TZ/QZ の外挿では足りないことがあります。

## 参考文献

- G. Bistoni, C. Riplinger, Y. Minenkov, L. Cavallo, A. A. Auer, F. Neese, "Treating Subvalence Correlation Effects in Domain Based Pair Natural Orbital Coupled Cluster Calculations: An Out-of-the-Box Approach", J. Chem. Theory Comput. 13, 3220 (2017). DOI: [10.1021/acs.jctc.7b00352](https://doi.org/10.1021/acs.jctc.7b00352)
- M. Saitow, U. Becker, C. Riplinger, E. F. Valeev, F. Neese, "A new near-linear scaling, efficient and accurate, open-shell domain-based local pair natural orbital coupled cluster singles and doubles theory", J. Chem. Phys. 146, 164105 (2017). DOI: [10.1063/1.4981521](https://doi.org/10.1063/1.4981521)
- Y. Guo, C. Riplinger, U. Becker, D. G. Liakos, Y. Minenkov, L. Cavallo, F. Neese, "Communication: An improved linear scaling perturbative triples correction for the domain based local pair-natural orbital based singles and doubles coupled cluster method [DLPNO-CCSD(T)]", J. Chem. Phys. 148, 011101 (2018). DOI: [10.1063/1.5011798](https://doi.org/10.1063/1.5011798)
- Y. Guo, C. Riplinger, D. G. Liakos, U. Becker, M. Saitow, F. Neese, "Linear scaling perturbative triples correction approximations for open-shell domain-based local pair natural orbital coupled cluster singles and doubles theory [DLPNO-CCSD(T0/T)]", J. Chem. Phys. 152, 024116 (2020). DOI: [10.1063/1.5127550](https://doi.org/10.1063/1.5127550)
- B. M. Flöser, Y. Guo, C. Riplinger, F. Tuczek, F. Neese, "Detailed Pair Natural Orbital-Based Coupled Cluster Studies of Spin Crossover Energetics", J. Chem. Theory Comput. 16, 2224 (2020). DOI: [10.1021/acs.jctc.9b01109](https://doi.org/10.1021/acs.jctc.9b01109)
- M. Drosou, C. A. Mitsopoulou, D. A. Pantazis, "Reconciling Local Coupled Cluster with Multireference Approaches for Transition Metal Spin-State Energetics", J. Chem. Theory Comput. 18, 3538 (2022). DOI: [10.1021/acs.jctc.2c00265](https://doi.org/10.1021/acs.jctc.2c00265)
- E. Semidalas, J. M. L. Martin, "The MOBH35 Metal–Organic Barrier Heights Reconsidered: Performance of Local-Orbital Coupled Cluster Approaches in Different Static Correlation Regimes", J. Chem. Theory Comput. 18, 883 (2022). DOI: [10.1021/acs.jctc.1c01126](https://doi.org/10.1021/acs.jctc.1c01126)
- A. Altun, C. Riplinger, F. Neese, G. Bistoni, "Exploring the Accuracy Limits of PNO-Based Local Coupled-Cluster Calculations for Transition-Metal Complexes", J. Chem. Theory Comput. 19, 2039 (2023). DOI: [10.1021/acs.jctc.3c00087](https://doi.org/10.1021/acs.jctc.3c00087)
- M. Radoń, G. Drabik, M. Hodorowicz, J. Szklarzewicz, "Performance of quantum chemistry methods for a benchmark set of spin-state energetics derived from experimental data of 17 transition metal complexes (SSE17)", Chem. Sci. 15, 20189 (2024). DOI: [10.1039/D4SC05471G](https://doi.org/10.1039/D4SC05471G)
- A. Mahajan, J. H. Thorpe, J. S. Kurian, D. R. Reichman, D. A. Matthews, S. Sharma, "Beyond CCSD(T) accuracy at lower scaling with auxiliary field quantum Monte Carlo", arXiv:2410.02885 (2024).
