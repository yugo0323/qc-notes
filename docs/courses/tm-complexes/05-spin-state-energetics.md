# 第 5 回　スピン状態のエネルギー差：なぜ計算が難しいか

!!! abstract "この回の道筋"
    高スピンと低スピンのエネルギー差は、遷移金属錯体の計算でいちばん難しい量の一つです。この回では、まず DFT の結果が汎関数の厳密交換の割合に強く依存することを、文献と自分の計算で確かめ、その理由を交換と相関の誤差の打ち消しから説明します。次に、CASPT2 や CCSD(T) といった波動関数法にも課題があることを示します。最後に、計算と実験を比べるときに必要な補正を説明します。

    | 節 | 内容 |
    |---|---|
    | 1 | 何を計算したいのか |
    | 2 | DFT の結果は厳密交換の割合で大きく変わる |
    | 3 | なぜ厳密交換に敏感なのか：交換と相関の誤差の打ち消し |
    | 4 | 波動関数法の課題：CASPT2 と CCSD(T) |
    | 5 | 実験値と比べる難しさ |

    **前に読んでおく回**：[第 3 回 スピン状態：高スピンと低スピン](03-spin-states.md)、[第 4 回 遷移金属の電子相関：静的相関と動的相関](04-electron-correlation-tm.md)

## 1. 何を計算したいのか

この回で扱う量は、高スピン（HS）と低スピン（LS）の電子エネルギーの差

$$
\Delta E_\mathrm{HL} = E_\mathrm{HS} - E_\mathrm{LS}
$$

です。$\Delta E_\mathrm{HL} < 0$ なら高スピンが基底状態、$\Delta E_\mathrm{HL} > 0$ なら低スピンが基底状態です。文献によって符号の決め方が違うので、論文を読むときは定義を確かめてください。

この差は、錯体の全エネルギー（数千ハートリー）に比べてとても小さい量です。それでも、基底状態がどちらのスピン状態かで、構造、磁性、反応性が変わります。第 3 回のスピンクロスオーバーのように、差が数 kcal/mol しかない錯体も多くあります。そのため、数 kcal/mol の精度が求められます。

**ここまでで分かったこと**：求めたいのは、高スピンと低スピンのエネルギー差です。数 kcal/mol の精度が必要です。

## 2. DFT の結果は厳密交換の割合で大きく変わる

### 2.1 Reiher らの観察と B3LYP*

混成汎関数は、交換エネルギーの一部を、Hartree–Fock の交換（厳密交換）で置き換えた汎関数です。B3LYP では、その割合が 20% です。

Reiher、Salomon、Hess は 2001 年に、第一配位圏に硫黄を多く含む鉄錯体について、スピン多重度の違う状態のエネルギー差を調べました（Reiher ら 2001）。まず、こうした錯体では、汎関数によって結果が最大で 1 eV（約 23 kcal/mol）も違うことを指摘しました。そして、Hartree–Fock 型の理論は高スピンを系統的に有利にすること、エネルギー差が厳密交換の割合に対して**直線的に**変わることを示しました。そのうえで、遷移金属の化合物で意味のあるエネルギーを得るには、B3LYP の 20% を約 15% に減らすべきだと提案しました。この 15% の汎関数は、のちに B3LYP* と呼ばれるようになりました。

### 2.2 自分の計算で確かめる

同じことを、[Fe(H$_2$O)$_6$]$^{2+}$ で確かめました。PBE 汎関数の交換の一部を厳密交換で置き換えた汎関数

$$
E_\mathrm{xc} = a\,E_\mathrm{x}^\mathrm{HF} + (1-a)\,E_\mathrm{x}^\mathrm{PBE} + E_\mathrm{c}^\mathrm{PBE}
$$

を使い、厳密交換の割合 $a$ を 0 から 0.4 まで変えました（$a = 0.25$ が PBE0 です）。高スピン（$S = 2$）は UKS、低スピン（$S = 0$）は RKS で、基底関数は def2-SVP です。比べやすくするため、両方のスピン状態に同じ構造（Fe–O 2.05 Å のモデル構造）を使っています。スクリプトは付録 A に載せます。

```text
 a     E(LS)-E(HS) / kcal mol-1   <S^2>(HS)
0.00       25.1                  6.0026   
0.10       30.4                  6.0027   
0.15       32.9                  6.0028   
0.20       35.3                  6.0028   
0.25       37.7                  6.0029   
0.30       40.0                  6.0030   
0.40       44.2                  6.0032   
```

どの $a$ でも高スピンが低いのですが、その差は $a$ によって 25 kcal/mol から 44 kcal/mol まで変わります。直線で近似すると、傾きは**厳密交換 10% あたり 4.8 kcal/mol** で、直線からのずれは最大 0.45 kcal/mol でした。厳密交換を増やすほど、高スピンが有利になります。Reiher らが見た直線的な依存が、ここでも現れています。

この傾きの大きさを、第 1 節の精度の目標と比べてください。厳密交換の割合を 5% 変えるだけで、約 2.4 kcal/mol 変わります。汎関数の選び方で、スピン状態の順序が入れ替わってもおかしくありません。

なお、この計算は同じ構造での電子エネルギーの差なので、実験と直接は比べられません（第 5 節）。ここで見たいのは、$a$ に対する依存の大きさです。

### 2.3 より広い系での傾向

Ioannidis と Kulik は、代表的な Fe(II) と Fe(III) の正八面体錯体について、厳密交換の割合を 0% から 50% まで変えて調べました（Ioannidis, Kulik 2015）。どの錯体でも、スピン状態のエネルギー差は厳密交換の割合に対して直線的に変わり、その傾きは主に金属に直接結合した配位子の原子の種類で決まりました。また、厳密交換を増やすと Fe の上の電荷はむしろ減る（電荷が配位子へ広がる）ことも見いだしています。「厳密交換は電子の非局在化を抑える」という単純な見方だけでは、スピン状態への影響を説明できないことを示す結果です。

Swart と Gruden は、初期の DFT の研究から、混成汎関数は高スピンを、初期の GGA 汎関数は低スピンを有利にする固有の傾向があることを指摘しています（Swart, Gruden 2016）。この問題に対して、鉄錯体のスピン状態をよく再現する汎関数として OPBE が見いだされ（Swart 2008）、その後も SSB-D や S12g のように、スピン状態を意識して作られた汎関数が提案されてきました（Swart, Gruden 2016）。

しかし、1 つの汎関数がすべての錯体でよい結果を出すわけではありません。Radoń は、4 つの Fe 錯体の実験から求めた参照値と比べ、すべてで釣り合いのとれた結果を出せる汎関数は少数（よかったのは B2PLYP-D3 と OPBE）しかないことを示しました（Radoń 2019）。さらに 17 個の錯体の参照値（SSE17）と比べた研究では、二重混成汎関数の PWPB95-D3(BJ) と B2PLYP-D3(BJ) が平均絶対誤差 3 kcal/mol 未満で最もよく、これまでスピン状態に推奨されてきた B3LYP*-D3(BJ) と TPSSh-D3(BJ) は、平均絶対誤差がそれぞれ 5.1 と 7.7 kcal/mol、最大誤差が 14〜15 kcal/mol でした（Radoń ら 2024）。

**ここまでで分かったこと**：スピン状態のエネルギー差は、厳密交換の割合にほぼ比例して変わり、その傾きは 10% あたり数 kcal/mol です。どの錯体にもよい汎関数は、まだありません。

## 3. なぜ厳密交換に敏感なのか：交換と相関の誤差の打ち消し

### 3.1 Fermi 相関と Coulomb 相関

第 3 回で、同じ向きのスピンの電子の対が多いほど、交換で安定になることを見ました。高スピンは、低スピンより同じスピンの対が多い状態です。

Hartree–Fock 法は、同じスピンの電子が互いを避ける効果（Pauli の原理による避け合い。**Fermi 相関**と呼びます）を、交換の項として正しく含みます。一方、逆向きのスピンの電子が避け合う効果（**Coulomb 相関**）は含みません。第 4 回の電子相関 $E_\mathrm{corr}$ は、この欠けた部分です。

低スピンの状態では、同じ空間軌道に逆スピンの電子が対になって入るので、Coulomb 相関の効果が大きくなります。HF 法はそれを欠いているので、低スピンを不利に、高スピンを有利に見積もります。HF 法が高スピンを系統的に有利にする理由は、このように説明されています（Paulsen, Trautwein 2004）。厳密交換を混ぜた汎関数が高スピンを有利にするのも、この HF 法の性質を一部受け継ぐからだと考えられます。

### 3.2 打ち消しが崩れる

この節の説明は、定量的に証明されたものではなく、結果を理解するための定性的な見方です。

近似の汎関数は、交換の部分と相関の部分を組にして使い、それぞれの誤差がある程度打ち消し合うことで、よい結果を出しています。厳密交換の割合を変えることは、この組み合わせの釣り合いを変えることです。

- 厳密交換を増やすと、Fermi 相関（同じスピンの電子の避け合い）はより正確に扱われます。しかし、それに見合うだけの Coulomb 相関は、汎関数の相関の部分から補われません。その結果、高スピンが有利になりすぎます。
- 厳密交換を減らすと、近似の交換汎関数の性質が強く出て、低スピンが有利になります（第 2 節の計算、および GGA が低スピンを有利にする傾向）。近似の交換汎関数には、電子が自分自身と反発してしまう誤差（自己相互作用誤差）が残り、これが第 4 回の非局在化誤差の原因です。この誤差も、電子の分布が違う 2 つのスピン状態で同じ大きさにはなりません。

どちらの向きの誤差も、高スピンと低スピンで電子の分布が違うために、打ち消し合わずに残ります。適切な割合は錯体ごとに違うので、1 つの割合ですべての錯体をうまく扱うことはできません。これが、第 2 節で見た「汎関数の普遍性のなさ」の理由です。

**ここまでで分かったこと**：厳密交換は同じスピンの電子の相関だけを正しく扱うので、増やすと高スピンが有利になります。交換と相関の誤差の釣り合いが錯体ごとに違うので、1 つの割合では足りません。

## 4. 波動関数法の課題：CASPT2 と CCSD(T)

それなら、汎関数に頼らない波動関数法を使えばよいように思えます。しかし、波動関数法にも課題があります。

### 4.1 CASPT2 と IPEA シフト

CASPT2 は、CASSCF の波動関数を出発点に、動的相関を 2 次の摂動論で加える多参照の方法です。静的相関を扱えるので、遷移金属錯体に広く使われてきました。

Pierloot と Vancoillie は、大きな基底関数と十分に大きな活性空間（10 電子 12 軌道）を使った CASPT2 で、[Fe(H$_2$O)$_6$]$^{2+}$、[Fe(NH$_3$)$_6$]$^{2+}$、[Fe(bpy)$_3$]$^{2+}$ の $\Delta E_\mathrm{HL}$ を計算しました（Pierloot, Vancoillie 2006）。零点エネルギーの補正を含めた値は、それぞれ −17690、−8389、+3820 cm$^{-1}$（約 −50.6、−24.0、+10.9 kcal/mol）でした。この値は、DFT の結果を較正する基準として使われてきました。

ところが、CASPT2 には、開殻の状態の相関エネルギーを大きく見積もりすぎる系統的な誤差があります。Ghigo、Roos、Malmqvist は、これを補正するため、零次のハミルトニアンに経験的なシフト（**IPEA シフト**）を入れました（Ghigo ら 2004）。その標準の値は 0.25 $E_\mathrm{h}$ です。

スピン状態の問題では、この標準の値でも足りないことが分かってきました。

- Kepenekian らは、Fe(II)N$_6$ の錯体の高スピンと低スピンの差を実験と比べ、標準の 0.25 $E_\mathrm{h}$ では不十分で、0.50〜0.70 $E_\mathrm{h}$ のほうがよいと提案しました（Kepenekian ら 2009）。
- Pierloot らは、CASPT2 の高スピンへの偏りが、価電子の相関ではなく、3s3p の準内殻相関の不正確な記述から来ていることを示しました（第 4 回）。IPEA シフトを大きくすると 3s3p の部分の誤差は減りますが、価電子の相関の記述が悪くなります（Pierloot ら 2017）。

つまり、IPEA シフトを大きくしてよい結果が得られても、それは誤差の打ち消しによる部分があります。実験と比べた検証でも、CASPT2 は高スピンを最大 5.5 kcal/mol 安定化しすぎ、同じ多参照の摂動論の NEVPT2 は最大 7 kcal/mol の誤差を出しました（Radoń 2019）。SSE17 でも、IPEA シフト 0.25 $E_\mathrm{h}$ の CASPT2 は、平均の符号付き誤差 −3.3 kcal/mol（高スピン寄り）、平均絶対誤差 4.1 kcal/mol でした（Radoń ら 2024）。3s3p の部分を CCSD(T) で置き換える CASPT2/CC 法（第 4 回）は、この偏りを一部取り除きます（Phung ら 2018; Radoń 2019）。

### 4.2 CCSD(T)

単一参照の CCSD(T) は、静的相関が強くない系では非常に正確です。実験から求めた参照値との比較では、次の結果が出ています。

- Radoń（4 つの Fe 錯体）：最大誤差 2 kcal/mol 未満、平均絶対誤差 1 kcal/mol 未満。特に、Kohn–Sham 軌道を出発点にした CCSD(T) がよい結果でした（Radoń 2019）。
- SSE17（17 個の錯体）：平均絶対誤差 1.5 kcal/mol、最大誤差 −3.5 kcal/mol で、調べた多参照の方法のどれよりもよい結果でした。ただし、Kohn–Sham 軌道を使っても、一貫した改善は見られませんでした（Radoń ら 2024）。

一方で、CCSD(T) がうまくいかないという報告もあります。Song らは、全電子の固定節拡散モンテカルロ法（DMC）で、小さな配位子をもつ 4 つの Fe(II) スピンクロスオーバー錯体の高スピンと低スピンの断熱エネルギー差の基準値を求め、半局所の DFT と CCSD(T) がどちらも大きく外れると結論しました（Song ら 2018）。このように、高い精度の方法どうしでも結論が食い違うことがあり、どの方法をどこまで信頼できるかは、まだ研究が続いている問題です。

CCSD(T) を使うときにも、第 4 回で見たように、3s3p の電子を相関させること、十分に大きな基底関数を使うことが必要です。計算量は大きく、大きな錯体には局所相関の近似（DLPNO-CCSD(T) など。第 9 回（準備中））が使われます。

**ここまでで分かったこと**：CASPT2 は高スピンを有利にしすぎる傾向があり、IPEA シフトでは根本的には直りません。CCSD(T) は実験から求めた参照値をよく再現しますが、別の基準と食い違う報告もあります。

## 5. 実験値と比べる難しさ

方法の良し悪しを決めるには、正しい参照値が必要です。しかし、実験で直接測れるのは、計算の $\Delta E_\mathrm{HL}$ そのものではありません。

### 5.1 実験で測れる量

スピン状態のエネルギー差に関わる実験の量は、主に次の 2 つです。

- **スピンクロスオーバーのエンタルピー** $\Delta H$：第 3 回で見たように、高スピンと低スピンが平衡にある錯体で、熱測定などから求めます。
- **スピン禁制の吸収帯のエネルギー**：基底状態と違うスピン多重度の状態への、弱い吸収です。構造を変えずに起こる遷移なので、垂直のエネルギー差にあたります。

### 5.2 必要な補正

計算の電子エネルギーの差と比べるには、実験の値から次の効果を取り除く（逆補正する）必要があります。

- **振動の効果**：$\Delta H$ には、零点振動のエネルギーの差と、熱による振動のエネルギーの差が含まれます。第 3 回で見たように、高スピンでは結合が弱く振動数が低いので、この差は無視できません。
- **環境の効果**：実験は結晶や溶液の中で行われます。結晶の中の分子の詰まり方や、溶媒和は、2 つのスピン状態を違う大きさで安定化します。

Radoń らは、SSE17 を作るとき、これらの補正を見積もりました（Radoń ら 2024）。スピンクロスオーバーを起こす錯体では、振動の補正は 2 kcal/mol 以内で、すべて負でした。結晶の詰まり方による環境の補正は 1.5 kcal/mol 以内でした。2 つの補正は、合わせて最大 3 kcal/mol に達し、参照値そのものの不確かさは 1〜3 kcal/mol と見積もられています。補正の大きさが、求めたい精度（数 kcal/mol）と同じ程度だということに注意してください。

### 5.3 実験の帰属そのものが問題になる場合

参照値の問題は、補正だけではありません。Reimann と Kaupp は、多参照の方法を組み合わせた精密な計算（CASPT2+δMRCI）で、3d 金属の六水和イオンのスピン反転の励起エネルギーを調べました（Reimann, Kaupp 2023）。その結果、Fe$^{3+}$ の六水和イオンの最低の吸収として 1.56 eV の値が、計算法の評価によく使われてきましたが、これは [Fe(H$_2$O)$_6$]$^{3+}$ ではなく、プロトンが 1 つ外れた [Fe(H$_2$O)$_5$OH]$^{2+}$ によるものである強い証拠を示しました。参照値に使う実験の帰属が正しいかどうかも、確かめる必要があるということです。

**ここまでで分かったこと**：実験の値を計算と比べるには、振動と環境の効果を逆補正する必要があり、その大きさは数 kcal/mol に達します。実験の帰属そのものが誤っていることもあります。

## まとめ

- スピン状態のエネルギー差は全エネルギーに比べてとても小さく、数 kcal/mol の精度が必要です。
- DFT の結果は厳密交換の割合にほぼ比例して変わり、[Fe(H$_2$O)$_6$]$^{2+}$ の計算では 10% あたり約 4.8 kcal/mol でした。厳密交換は Fermi 相関だけを正しく扱うので、増やすと高スピンが有利になります。
- どの錯体にもよい汎関数はまだなく、最近の比較では二重混成汎関数がよい結果を出しています。
- CASPT2 は高スピンを有利にしすぎる傾向があり、その一因は 3s3p 相関です。CCSD(T) は実験から求めた参照値をよく再現しますが、別の基準と食い違う報告もあります。
- 実験と比べるには振動と環境の効果の逆補正が必要で、参照値にも 1〜3 kcal/mol の不確かさがあります。

## 付録 A　[Fe(H$_2$O)$_6$]$^{2+}$ の計算のスクリプト

PySCF 2.15 で実行しました。密度フィッティング（`density_fit`）を使い、1 コアで、厳密交換の割合 1 つあたり（高スピンと低スピンの 2 つの計算で）4 分ほどかかりました。

```python
"""[Fe(H2O)6]2+ の高スピン(S=2)と低スピン(S=0)のエネルギー差を、
厳密交換の割合 a を変えた PBE0 型の汎関数 a*HF + (1-a)*PBE で計算する。
両方のスピン状態に同じ構造（Fe-O 2.05 Å のモデル構造）を使う。"""
import numpy as np
from pyscf import gto, dft
rFeO, rOH, ang = 2.05, 0.97, np.deg2rad(104.5)
atoms = [("Fe", (0, 0, 0))]
ex = np.eye(3)
plane = {0: 1, 1: 2, 2: 0}   # 軸 x の水は xy 面、y は yz 面、z は zx 面
for ax in range(3):
    for sgn in (1, -1):
        u = sgn*ex[ax]; v = ex[plane[ax]]
        O = rFeO*u
        atoms.append(("O", O))
        for t in (1, -1):
            atoms.append(("H", O + rOH*(np.cos(ang/2)*u + t*np.sin(ang/2)*v)))
def run(spin, a):
    mol = gto.M(atom=atoms, basis="def2-svp", charge=2, spin=spin, verbose=0)
    mf = (dft.RKS(mol) if spin == 0 else dft.UKS(mol)).density_fit()
    mf.xc = f"{a}*HF + {1-a}*PBE, PBE"
    mf.conv_tol = 1e-8; mf.max_cycle = 200
    e = mf.kernel()
    s2 = mf.spin_square()[0] if spin else 0.0
    return e, s2, mf.converged
print(" a     E(LS)-E(HS) / kcal mol-1   <S^2>(HS)")
for a in (0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40):
    eh, s2, c1 = run(4, a); el, _, c2 = run(0, a)
    print(f"{a:4.2f}   {(el-eh)*627.5:8.1f}                  {s2:.4f}   {'' if c1 and c2 else '未収束'}", flush=True)
```

出力は本文の 2.2 に載せたとおりです（すべての計算が収束しました）。直線の当てはめは、出力の 7 点に `numpy.polyfit` で 1 次式を当てはめて求めました（傾き 47.8 kcal/mol、すなわち 10% あたり 4.78 kcal/mol、切片 25.5 kcal/mol）。

## 参考文献

- M. Reiher, O. Salomon, B. A. Hess, "Reparameterization of hybrid functionals based on energy differences of states of different multiplicity", Theor. Chem. Acc. 107, 48 (2001). DOI: [10.1007/s00214-001-0300-3](https://doi.org/10.1007/s00214-001-0300-3)
- G. Ghigo, B. O. Roos, P.-Å. Malmqvist, "A modified definition of the zeroth-order Hamiltonian in multiconfigurational perturbation theory (CASPT2)", Chem. Phys. Lett. 396, 142 (2004). DOI: [10.1016/j.cplett.2004.08.032](https://doi.org/10.1016/j.cplett.2004.08.032)
- K. Pierloot, S. Vancoillie, "Relative energy of the high-(5T2g) and low-(1A1g) spin states of [Fe(H2O)6]2+, [Fe(NH3)6]2+, and [Fe(bpy)3]2+: CASPT2 versus density functional theory", J. Chem. Phys. 125, 124303 (2006). DOI: [10.1063/1.2353829](https://doi.org/10.1063/1.2353829)
- H. Paulsen, A. X. Trautwein, "Density Functional Theory Calculations for Spin Crossover Complexes", Top. Curr. Chem. 235, 197 (2004). DOI: [10.1007/b95428](https://doi.org/10.1007/b95428)
- M. Swart, "Accurate Spin-State Energies for Iron Complexes", J. Chem. Theory Comput. 4, 2057 (2008). DOI: [10.1021/ct800277a](https://doi.org/10.1021/ct800277a)
- M. Kepenekian, V. Robert, B. Le Guennic, "What zeroth-order Hamiltonian for CASPT2 adiabatic energetics of Fe(II)N6 architectures?", J. Chem. Phys. 131, 114702 (2009). DOI: [10.1063/1.3211020](https://doi.org/10.1063/1.3211020)
- E. I. Ioannidis, H. J. Kulik, "Towards quantifying the role of exact exchange in predictions of transition metal complex properties", J. Chem. Phys. 143, 034104 (2015). DOI: [10.1063/1.4926836](https://doi.org/10.1063/1.4926836)
- M. Swart, M. Gruden, "Spinning around in Transition-Metal Chemistry", Acc. Chem. Res. 49, 2690 (2016). DOI: [10.1021/acs.accounts.6b00271](https://doi.org/10.1021/acs.accounts.6b00271)
- K. Pierloot, Q. M. Phung, A. Domingo, "Spin State Energetics in First-Row Transition Metal Complexes: Contribution of (3s3p) Correlation and Its Description by Second-Order Perturbation Theory", J. Chem. Theory Comput. 13, 537 (2017). DOI: [10.1021/acs.jctc.6b01005](https://doi.org/10.1021/acs.jctc.6b01005)
- Q. M. Phung, M. Feldt, J. N. Harvey, K. Pierloot, "Toward Highly Accurate Spin State Energetics in First-Row Transition Metal Complexes: A Combined CASPT2/CC Approach", J. Chem. Theory Comput. 14, 2446 (2018). DOI: [10.1021/acs.jctc.8b00057](https://doi.org/10.1021/acs.jctc.8b00057)
- S. Song, M.-C. Kim, E. Sim, A. Benali, O. Heinonen, K. Burke, "Benchmarks and Reliable DFT Results for Spin Gaps of Small Ligand Fe(II) Complexes", J. Chem. Theory Comput. 14, 2304 (2018). DOI: [10.1021/acs.jctc.7b01196](https://doi.org/10.1021/acs.jctc.7b01196)
- M. Radoń, "Benchmarking quantum chemistry methods for spin-state energetics of iron complexes against quantitative experimental data", Phys. Chem. Chem. Phys. 21, 4854 (2019). DOI: [10.1039/C9CP00105K](https://doi.org/10.1039/C9CP00105K)
- M. Reimann, M. Kaupp, "Spin-State Splittings in 3d Transition-Metal Complexes Revisited: Toward a Reliable Theory Benchmark", J. Chem. Theory Comput. 19, 97 (2023). DOI: [10.1021/acs.jctc.2c00925](https://doi.org/10.1021/acs.jctc.2c00925)
- M. Radoń, G. Drabik, M. Hodorowicz, J. Szklarzewicz, "Performance of quantum chemistry methods for a benchmark set of spin-state energetics derived from experimental data of 17 transition metal complexes (SSE17)", Chem. Sci. 15, 20189 (2024). DOI: [10.1039/D4SC05471G](https://doi.org/10.1039/D4SC05471G)
