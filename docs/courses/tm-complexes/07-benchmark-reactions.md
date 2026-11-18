# 第 7 回　ベンチマーク（2）：反応エネルギーと反応障壁

!!! abstract "この回の道筋"
    触媒反応の計算で知りたいのは、反応エネルギーと反応障壁です。この回では、遷移金属の反応の代表的なベンチマーク（MOR41、ROST61、TMC151、MOBH35）の中身と、そこから分かる汎関数の傾向を順に見ます。MOBH35 については、参照値が後の研究でどう見直されたかを追い、参照値に使う CCSD(T) そのものがどこまで信頼できるかを考えます。

    | 節 | 内容 |
    |---|---|
    | 1 | 反応のベンチマークの作り方 |
    | 2 | MOR41：閉殻の有機金属反応 |
    | 3 | ROST61：開殻の反応 |
    | 4 | TMC151：小さな系の集まり |
    | 5 | MOBH35：反応障壁と、参照値の見直し |
    | 6 | 参照値の CCSD(T) はどこまで信頼できるか |
    | 7 | データセットを通して見える傾向 |

    **前に読んでおく回**：[第 6 回 ベンチマーク（1）：スピン状態のエネルギー](06-benchmark-spin-states.md)、[量子化学の理論 第 11 回 局所相関法：DLPNO-CCSD(T)](../qc-theory/11-dlpno-ccsd-t.md)

## 1. 反応のベンチマークの作り方

第 6 回のスピン状態では、実験から参照値を作る道がありました。反応エネルギーと反応障壁では、この道は使いにくくなります。触媒反応の中間体や遷移状態のエネルギーは、ほとんど測れないからです。そのため、反応のベンチマークの参照値は、ほぼすべて計算で作られています。

現実的な触媒の錯体は、数十〜百原子を超えます。この大きさでは、正準 CCSD(T) を大きな基底関数で計算するのは困難です。そこで多くのベンチマークは、局所相関法の DLPNO-CCSD(T) を使い、2 つの基底関数の結果から基底関数の極限（CBS）へ外挿した値を参照値にしています。DLPNO-CCSD(T) の考え方は[量子化学の理論 第 11 回](../qc-theory/11-dlpno-ccsd-t.md)で扱いました。

したがって、反応のベンチマークの参照値には、次の 3 つの誤差がありえます。

1. DLPNO の近似による、正準 CCSD(T) からのずれ。
2. 基底関数の外挿の誤差。
3. CCSD(T) そのものの誤差（多参照性が強い場合に大きくなる）。

この回の後半（第 5・6 節）では、これらの誤差が実際に問題になった例を見ます。

**ここまでで分かったこと**：反応のベンチマークの参照値は、ほとんどが DLPNO-CCSD(T)/CBS などの計算値です。その参照値にも、局所近似・基底関数・CCSD(T) 自体の誤差がありえます。

## 2. MOR41：閉殻の有機金属反応

MOR41 は、Dohm らが 2018 年に発表した、41 個の閉殻の有機金属反応のデータセットです（Dohm ら 2018）。触媒反応でよく使われる変換に似た反応を集め、それまでの遷移金属のテストセットより大きな分子を含んでいます。3d・4d・5d の金属を含み、最大 120 原子の化合物があります。

- **参照値**：DLPNO-CCSD(T) を def2-TZVPP と def2-QZVPP で計算し、CBS へ外挿した値です。
- **評価した方法**：GGA 7 個、meta-GGA 3 個、混成 14 個、二重混成 17 個の計 41 個の汎関数に、2 種類の分散補正を組み合わせ、MP2 や RPA なども加えた 90 の方法です。

要旨に挙げられた平均絶対誤差（MAD）は、次のとおりです。

| 方法 | 段 | MAD（kcal/mol） |
|---|---|---|
| PWPB95-D3(BJ) | 二重混成 | 1.9 |
| ωB97X-V | 範囲分離混成 | 2.2 |
| mPW1B95-D3(BJ) | 混成 | 2.4 |
| PBE0-D3(BJ) | 混成 | 2.8 |
| SCAN-D3(BJ) | meta-GGA | 3.2 |
| TPSS-D3(BJ) | meta-GGA | 3.3 |
| BP86-D3(BJ) | GGA | 4.9 |
| B3LYP-D3(BJ) | 混成 | 4.9 |

DFT の基底関数は、要旨には書かれていません。数値は要旨から写しました。

著者らは、BP86-D3(BJ) と B3LYP-D3(BJ) は、同じ計算量でよりよい方法があるので、有機金属の熱化学には勧められないと述べています。また、汎関数の性能の傾向は、主要族元素の熱化学のベンチマーク（[量子化学の理論 第 16 回](../qc-theory/16-functionals.md)の GMTKN55）の結論と一致していました。

**ここまでで分かったこと**：MOR41 では、二重混成の PWPB95-D3(BJ) が MAD 1.9 kcal/mol で最もよく、範囲分離混成と混成がそれに続きました。B3LYP-D3(BJ) と BP86-D3(BJ) は 4.9 kcal/mol でした。

## 3. ROST61：開殻の反応

触媒反応には、不対電子をもつ開殻の錯体がよく出てきます。ROST61 は、Maurer らが MOR41 を開殻の系に広げた、61 個の反応エネルギーのデータセットです（Maurer ら 2021）。名前は「単一参照の開殻遷移金属錯体の反応」の略です。

- **中身**：13〜93 原子の錯体で、20 種類の d ブロック元素を含みます。反応エネルギーの平均は −42.8 kcal/mol です。
- **参照値**：結合クラスター法の値です。Semidalas と Martin の説明によると、DLPNO-CCSD(T$_1$)/CBS です（Semidalas, Martin 2022）。(T$_1$) は、三重励起の補正を反復的に改良した版です（[第 9 回](09-dlpno-for-tm.md)）。
- **制限**：単一参照で扱える開殻の系に限ったため、鉄や白金のように重要な元素が、ほとんど含まれていません。

31 個の汎関数を 3 種類の分散補正と組み合わせて評価した結果、要旨には次の値が挙げられています。

| 方法 | 基底関数 | MAD（kcal/mol） |
|---|---|---|
| PWPB95-D4（二重混成） | def2-QZVPP | 1.6 |
| TPSS0-D4（混成） | def2-QZVPP | 2.3 |
| r$^2$SCAN-3c（複合法） | 専用の基底関数 | 2.9 |

r$^2$SCAN-3c は、meta-GGA の r$^2$SCAN に、専用の小さめの三重ゼータ基底関数（mTZVPP）と補正を組み合わせた、計算量の小さい複合法です（Bursch ら 2022）。著者らは、汎関数の順位がおおむね梯子の段（[量子化学の理論 第 16 回](../qc-theory/16-functionals.md)）の順に並ぶこと、分散補正を加えるとほとんどすべての方法が改善することを、MOR41 と同じく見いだしました。また、厳密交換の割合が少なめの混成汎関数（TPSS0 など）がよい結果を出しています（Maurer ら 2021）。

**ここまでで分かったこと**：ROST61 は、単一参照で扱える開殻の反応のデータセットです。二重混成の PWPB95-D4 が MAD 1.6 kcal/mol で最もよく、分散補正はほぼ必須でした。鉄や白金はほとんど含まれていません。

## 4. TMC151：小さな系の集まり

TMC151 は、Chan、Gill、Kimura が、既存の参照データを集めて作った 151 個のデータの集まりです（Chan ら 2019）。二原子分子の解離エネルギーと、代表的な遷移金属の反応の反応エネルギーと障壁を含みます。

要旨によると、最も正確だった汎関数は、ωB97M-V、ωB97X-V、MN15、B97M-rV でした。これらは、主要族元素の化学でもよい結果を出すことが知られていた汎関数です。それでも、遷移金属での精度は、主要族元素よりかなり悪くなりました。また、精度の低い汎関数では、主要族元素での成績と遷移金属での成績の間に、よい相関がありませんでした。つまり、主要族元素のベンチマークだけを見て、遷移金属に使う汎関数を選ぶのは危険です。

この論文の個々の数値は、要旨に書かれていないため、ここでは挙げません。

**ここまでで分かったこと**：TMC151 では、範囲分離混成などの一部の汎関数が最もよい結果を出しましたが、遷移金属での精度は主要族元素より悪くなりました。主要族元素での成績は、遷移金属での成績を保証しません。

## 5. MOBH35：反応障壁と、参照値の見直し

### 5.1 もとのデータセット

MOBH35 は、Iron と Janes が 2019 年に発表した、35 個の有機金属反応の反応障壁のデータセットです（Iron, Janes 2019a）。各反応について、正方向と逆方向の障壁、および反応エネルギーを含みます。

- **参照値**：DLPNO-CCSD(T) を、Weizmann-1（W1）理論に似た方法で CBS へ外挿した値です。
- **結果**：def2-TZVPP 基底関数で、次の汎関数がよい結果でした。

| 方法 | 段 | MAD（kcal/mol、def2-TZVPP） |
|---|---|---|
| ωB97M-V | 範囲分離混成 meta-GGA | 1.7 |
| ωB97M-D3BJ | 同上 | 1.9 |
| ωB97X-V | 範囲分離混成 | 2.0 |
| revTPSS0-D4 | 混成 | 2.2 |
| B2K-PLYP | 二重混成 | 1.7 |
| revDOD-PBEP86-D4 | 二重混成 | 1.8 |

二重混成汎関数もよい結果でしたが、著者らは、計算量の増加と見合うかを考える必要があると述べています。

同じ年に、この論文の訂正が出ています（Iron, Janes 2019b）。Semidalas と Martin の説明によると、訂正後の参照値は、DLPNO-CCSD(T) の外挿値に、三重励起を改良した DLPNO-CCSD(T$_1$) と従来の DLPNO-CCSD(T$_0$) の差（def2-TZVP で計算）を加えたものです（Semidalas, Martin 2022）。

### 5.2 Semidalas と Martin による見直し

Semidalas と Martin は 2022 年に、MOBH35 の大部分の反応について、局所近似を使わない正準 CCSD(T) を def2-SVP 基底関数で計算し、参照値を見直しました（Semidalas, Martin 2022）。大きすぎる 4 つの二分子反応（反応 17〜20）は除いています。主な結果は次のとおりです。

- 新旧の参照値の違いは、ほとんどの反応で小さいものでした。後で述べる 2 つの反応（8 と 9）を除くと、新旧の差の二乗平均平方根（RMS）は 1.68 kcal/mol です。
- 例外は、鉄を含む連続した 2 つの反応 8 と 9 です。反応 9 では、正方向と逆方向の障壁の差がそれぞれ 5.3 と 7.0 kcal/mol に達し、反応エネルギーでは 12.3 kcal/mol も違いました。
- 反応 9 の生成物では、T1 診断が 0.05、D1 診断が 0.52 で、強い静的相関を示していました。これらの値は、反応物から遷移状態、生成物へと大きくなり、HOMO–LUMO ギャップは 10.6 eV から 9.3 eV、7.4 eV へと狭まりました。
- 局所近似の方法のうち、DLPNO-CCSD(T) と DLPNO-CCSD(T$_1$) は、この反応で大きく破綻しました。一方、LNO-CCSD(T) という別の局所相関法は、閾値を厳しくするにつれて、なめらかに正準 CCSD(T) へ近づきました。

著者らは、DFT の評価に MOBH35 を使うなら、反応 9 は静的相関が強すぎるので除くのがよい、と勧めています。反応 8 も除けば、多参照性の弱い参照値の集まりになります。また、反応 9 を除いても、多くの汎関数の MAD はわずかに下がるだけで、汎関数どうしの順位についての結論は変わりませんでした。ただし、最もよい二重混成汎関数の評価には注意が要る、としています。

この見直しの教訓は、参照値の誤差が一部の「難しい」反応に集中することです。平均の誤差が小さくても、個々の反応、特に鉄のように静的相関が強くなりやすい金属を含む反応では、参照値そのものが数 kcal/mol 以上ずれていることがあります。

**ここまでで分かったこと**：MOBH35 では、範囲分離混成と一部の二重混成汎関数が MAD 2 kcal/mol 前後でした。正準 CCSD(T) による見直しで、多くの参照値は妥当と確かめられましたが、鉄を含む反応 9 では DLPNO の参照値が大きくずれており、除くことが勧められています。

## 6. 参照値の CCSD(T) はどこまで信頼できるか

第 5 節では、DLPNO と正準 CCSD(T) の違いを見ました。では、正準 CCSD(T) そのものは、遷移金属の反応で信頼できるのでしょうか。これは、より高い精度の方法と比べなければ分かりません。ここでは、異なる角度からの研究を 3 つ紹介します。

**(1) 多参照性が強い反応**：Semidalas と Martin は、MOBH35 の反応 9 について、正準 CCSD(T) 自体が安全に使える範囲を超えている可能性を指摘しています。この大きさの系では、4 電子励起の寄与を見積もる手段がないため、確かめることもできません（Semidalas, Martin 2022）。

**(2) AFQMC との比較**：Neugebauer らは、電極触媒に関わる 28 個の 3d 金属の分子（3dTMV）の垂直イオン化エネルギーを、CCSD(T) と ph-AFQMC で計算しました（def2-SVP 基底関数、Neugebauer ら 2023）。多くの多参照性の診断を調べた結果、CCSD の波動関数や PBE0 の解でスピン対称性が破れるかどうかが、多配置の波動関数の解析とよく対応しました。この破れの大きさで「安全な領域」を決めると、その中では CCSD(T) と AFQMC の平均のずれが約 2 kcal/mol 以下になり、領域の外では CCSD(T) が失敗しうる、と結論しています。

**(3) そのまま使ったときの失敗**：Minenkov らは、11 種類の遷移金属を含む 111 個の気相の反応エンタルピーを、実験値と比べました（Minenkov ら 2016）。全体では DLPNO-CCSD(T)、B3LYP-D3、PBE0-D3 がよく、平均絶対誤差はすべて 5 kcal/mol 未満でした。しかし、個々の反応では、最もよかった汎関数でも 20〜30 kcal/mol の系統的でないずれが出ました。DLPNO-CCSD(T) も、Zr の錯体では平均絶対誤差が 14.4 kcal/mol を超えました。著者らは、CCSD(T) も「既定の設定のまま」使えば信頼できない場合があり、内殻相関、相対論効果、参照の波動関数の質を注意深く確かめれば、よい結果が得られると述べています。

これらをまとめると、次のように言えます。

- 静的相関が弱い系では、注意深く計算した CCSD(T) は、考え方の違う高精度の方法（AFQMC）と約 2 kcal/mol 以内で一致しました（3dTMV のイオン化エネルギー）。第 6 回のスピン状態でも、CCSD(T) は実験から作った参照値をよく再現しました。
- 静的相関が強い反応では、CCSD(T) 自体の誤差が大きくなりえます。その判定には、診断の値（[第 10 回](10-multireference-tm.md)）やスピン対称性の破れを確かめる必要があります。
- 局所近似（DLPNO）の誤差は、静的相関が強いほど大きくなりやすく、正準 CCSD(T) の誤差の上に重なります（[第 9 回](09-dlpno-for-tm.md)）。

**ここまでで分かったこと**：CCSD(T) は、静的相関が弱い反応では信頼できる参照値を与えますが、強い反応では自身の誤差が問題になります。そのような反応を参照値の集まりから除くか、多参照性を確かめることが必要です。

## 7. データセットを通して見える傾向

4 つのデータセットの性格をまとめます。

| データセット | 中身 | 参照値 | 最もよかった方法（要旨の値） |
|---|---|---|---|
| MOR41（2018） | 閉殻の有機金属反応 41、最大 120 原子 | DLPNO-CCSD(T)/CBS（def2-TZVPP/QZVPP） | PWPB95-D3(BJ)：MAD 1.9 |
| ROST61（2021） | 単一参照の開殻の反応 61、13〜93 原子 | DLPNO-CCSD(T$_1$)/CBS | PWPB95-D4/def2-QZVPP：MAD 1.6 |
| TMC151（2019） | 二原子の解離、代表的な反応のエネルギーと障壁 | 既存の高精度データを収集 | ωB97M-V、ωB97X-V、MN15、B97M-rV |
| MOBH35（2019、2022 見直し） | 有機金属反応 35 の正逆の障壁と反応エネルギー | DLPNO-CCSD(T)/CBS（W1 型）、のちに正準 CCSD(T) で見直し | ωB97M-V/def2-TZVPP：MAD 1.7 |

4 つを通して、次の傾向が見えます。ただし、どれも特定のデータセットの範囲での結果です。

1. **梯子の段の順序はおおむね保たれる**：平均としては、二重混成、範囲分離混成と混成、meta-GGA、GGA の順によくなります（Dohm ら 2018; Maurer ら 2021）。
2. **分散補正はほぼ必須**：分散補正を加えると、ほとんどの方法が改善します（Maurer ら 2021）。分散補正については、[量子化学の理論 第 17 回](../qc-theory/17-dispersion.md)と、この連載の第 12 回（準備中）で扱います。
3. **同じ段の中でも差が大きい**：B3LYP-D3(BJ) は、MOR41 で混成汎関数の中でも成績が悪い側でした（Dohm ら 2018）。
4. **遷移金属の誤差は主要族元素より大きい**：主要族元素で最良の汎関数も、遷移金属では誤差が大きくなります（Chan ら 2019）。
5. **難しい反応で順位が変わりうる**：平均の順位が同じでも、静的相関の強い反応では、汎関数どうしの差が大きくなります。最良の汎関数どうしを比べるときは、参照値そのものの揺れ（MOBH35 の新旧の参照値の差は、難しい 2 反応を除いても RMS で 1.7 kcal/mol 程度）より小さな差に、意味を持たせないようにします。

最近では、主要族元素と遷移金属の反応エネルギー・障壁を含む、さらに大きなデータベースも作られています。Liang と Head-Gordon の GSCDB137 は、137 個のデータセット（8377 項目）を集め、29 個の汎関数を評価しました（Liang, Head-Gordon 2025、プレプリント）。要旨によると、全体としては梯子の段の順序が確かめられ、二重混成汎関数は最良の混成汎関数より平均の誤差が約 25% 小さい一方、内殻の凍結、基底関数、多参照性の扱いに注意が要るとされています。

これらの結果を、目的に応じた汎関数の選び方にどうつなげるかは、[第 8 回](08-functional-choice.md)で扱います。

**ここまでで分かったこと**：反応のベンチマークを通して、梯子の段の順序と分散補正の大切さは共通していますが、遷移金属では誤差が大きく、難しい反応で順位が変わりえます。

## まとめ

- 反応のベンチマークの参照値は、ほとんどが DLPNO-CCSD(T)/CBS などの計算値で、局所近似・基底関数・CCSD(T) 自体の誤差がありえます。
- MOR41（閉殻）と ROST61（開殻）では、二重混成の PWPB95 が MAD 2 kcal/mol 未満で最もよく、分散補正はほぼ必須でした。
- MOBH35 では範囲分離混成と一部の二重混成がよく、正準 CCSD(T) による見直しで、鉄を含む反応 9 の参照値が大きくずれていたことが分かりました。
- CCSD(T) は、静的相関が弱い反応では信頼できますが、強い反応では自身の誤差が問題になります。
- 遷移金属の誤差は主要族元素より大きく、主要族元素のベンチマークだけで汎関数を選ぶのは危険です。

## 参考文献

- S. Dohm, A. Hansen, M. Steinmetz, S. Grimme, M. P. Checinski, "Comprehensive Thermochemical Benchmark Set of Realistic Closed-Shell Metal Organic Reactions", J. Chem. Theory Comput. 14, 2596 (2018). DOI: [10.1021/acs.jctc.7b01183](https://doi.org/10.1021/acs.jctc.7b01183)
- L. R. Maurer, M. Bursch, S. Grimme, A. Hansen, "Assessing Density Functional Theory for Chemically Relevant Open-Shell Transition Metal Reactions", J. Chem. Theory Comput. 17, 6134 (2021). DOI: [10.1021/acs.jctc.1c00659](https://doi.org/10.1021/acs.jctc.1c00659)
- B. Chan, P. M. W. Gill, M. Kimura, "Assessment of DFT Methods for Transition Metals with the TMC151 Compilation of Data Sets and Comparison with Accuracies for Main-Group Chemistry", J. Chem. Theory Comput. 15, 3610 (2019). DOI: [10.1021/acs.jctc.9b00239](https://doi.org/10.1021/acs.jctc.9b00239)
- M. A. Iron, T. Janes, "Evaluating Transition Metal Barrier Heights with the Latest Density Functional Theory Exchange–Correlation Functionals: The MOBH35 Benchmark Database", J. Phys. Chem. A 123, 3761 (2019). DOI: [10.1021/acs.jpca.9b01546](https://doi.org/10.1021/acs.jpca.9b01546)
- M. A. Iron, T. Janes, "Correction to “Evaluating Transition Metal Barrier Heights with the Latest Density Functional Theory Exchange–Correlation Functionals: The MOBH35 Benchmark Database”", J. Phys. Chem. A 123, 6379 (2019). DOI: [10.1021/acs.jpca.9b06135](https://doi.org/10.1021/acs.jpca.9b06135)
- E. Semidalas, J. M. L. Martin, "The MOBH35 Metal–Organic Barrier Heights Reconsidered: Performance of Local-Orbital Coupled Cluster Approaches in Different Static Correlation Regimes", J. Chem. Theory Comput. 18, 883 (2022). DOI: [10.1021/acs.jctc.1c01126](https://doi.org/10.1021/acs.jctc.1c01126)
- Y. Minenkov, E. Chermak, L. Cavallo, "Troubles in the Systematic Prediction of Transition Metal Thermochemistry with Contemporary Out-of-the-Box Methods", J. Chem. Theory Comput. 12, 1542 (2016). DOI: [10.1021/acs.jctc.5b01163](https://doi.org/10.1021/acs.jctc.5b01163)
- H. Neugebauer, H. T. Vuong, J. L. Weber, R. A. Friesner, J. Shee, A. Hansen, "Toward Benchmark-Quality Ab Initio Predictions for 3d Transition Metal Electrocatalysts: A Comparison of CCSD(T) and ph-AFQMC", J. Chem. Theory Comput. 19, 6208 (2023). DOI: [10.1021/acs.jctc.3c00617](https://doi.org/10.1021/acs.jctc.3c00617)
- M. Bursch, J.-M. Mewes, A. Hansen, S. Grimme, "Best-Practice DFT Protocols for Basic Molecular Computational Chemistry", Angew. Chem. Int. Ed. 61, e202205735 (2022). DOI: [10.1002/anie.202205735](https://doi.org/10.1002/anie.202205735)
- J. Liang, M. Head-Gordon, "Gold-Standard Chemical Database 137 (GSCDB137): A diverse set of accurate energy differences for assessing and developing density functionals", arXiv:2508.13468 (2025).
