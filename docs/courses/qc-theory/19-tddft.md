# 第 19 回　線形応答と TDDFT

!!! abstract "この回の道筋"
    励起状態のエネルギーを、基底状態の計算の延長として求める方法を学びます。まず、決定子 1 個から出発する最も簡単な方法 CIS を導き、次に時間に依存する Hartree–Fock 法（TDHF）と、密度汎関数理論を時間に依存する形に広げた TDDFT を、同じ形の行列の固有値問題として見ます。ホルムアルデヒドの励起エネルギーを計算して手法を比べ、最後に TDDFT が苦手な場面（電荷移動励起、二電子励起、三重項不安定性）を整理します。

    | 節 | 内容 |
    |---|---|
    | 1 | 励起状態を基底状態から求める |
    | 2 | CIS：1 電子励起の配置だけで対角化する |
    | 3 | TDHF：時間変化への応答として求める |
    | 4 | TDDFT：Runge–Gross の定理と断熱近似 |
    | 5 | 計算例：ホルムアルデヒドの励起エネルギー |
    | 6 | TDDFT が苦手な場面 |

    **前に読んでおく回**：[量子化学の理論 第 3 回](03-hartree-fock.md)、[第 5 回](05-spin.md)、[第 8 回](08-configuration-interaction.md)、[第 15 回](15-kohn-sham.md)、[第 16 回](16-functionals.md)

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書きます。1 $E_\mathrm{h}$ = 27.211 eV です。2 電子積分は、空間軌道についての化学者の表記 $(pq|rs)$ を使います。軌道はすべて実数とします。添字は、占有軌道を $i, j$、仮想軌道（空軌道）を $a, b$ で表します。

## 1. 励起状態を基底状態から求める

光の吸収や発光、光化学反応を調べるには、励起状態のエネルギーが必要です。基底状態のエネルギーとの差を**励起エネルギー**と呼び、光の吸収のスペクトルの位置に対応します。

励起状態を求める方法は、大きく二つに分けられます。

- **状態ごとに解く方法**：励起状態ごとに、別々に波動関数を求めます。第 12 回の状態平均 CASSCF と、第 13 回の CASPT2・NEVPT2 は、この考え方です。
- **基底状態から求める方法**：基底状態の波動関数を一度求め、そこからの「励起」として、多くの励起状態をまとめて求めます。この回で扱う CIS、TDHF、TDDFT は、こちらです。

後者の方法は、基底状態が決定子 1 個でよく表せることを前提にします。そのかわり、計算が比較的安く、大きな分子の多くの励起状態を一度に求められます。

**ここまでで分かったこと**：この回の方法は、決定子 1 個で表せる基底状態から出発し、多くの励起状態をまとめて求めます。

## 2. CIS：1 電子励起の配置だけで対角化する

### 2.1 行列の形

第 8 回の CI を思い出してください。CI の波動関数は、HF の決定子 $\Phi_0$ と、そこから電子を励起した配置の線形結合でした。

HF の決定子と 1 電子励起の配置の間のハミルトニアンの行列要素は、0 になります（Brillouin の定理、第 8 回）。そこで、1 電子励起の配置だけを使って励起状態を表し、その中でハミルトニアンを対角化します。これを **CIS**（configuration interaction singles）と呼びます [Foresman 1992]。

$$
\Psi_n = \sum_{ia} X_{ia}^{n}\, \Phi_i^a
$$

$\Phi_i^a$ は、占有軌道 $i$ の電子 1 個を仮想軌道 $a$ に移した配置です。

閉殻の分子では、配置をスピンについてそろえると、一重項と三重項の問題に分かれます（第 5 回）。行列の要素は次のとおりです。導出は付録 A に回します。

$$
A_{ia,jb}^{\mathrm S} = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + 2(ia|jb) - (ij|ab)
$$

$$
A_{ia,jb}^{\mathrm T} = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) - (ij|ab)
$$

$\varepsilon_p$ は HF の軌道エネルギーです。励起エネルギー $\omega_n$ は、固有値問題

$$
\mathbf A \mathbf X^n = \omega_n \mathbf X^n
$$

の固有値として求まります。

### 2.2 各項の意味

各項の意味は、次のとおりです。

- **$\varepsilon_a - \varepsilon_i$**：軌道エネルギーの差です。これだけなら、励起エネルギーは「空軌道と占有軌道のエネルギーの差」になります。
- **$-(ij|ab)$**：励起した電子（軌道 $a$）と、電子が抜けた穴（軌道 $i$）の間のクーロン引力です。軌道エネルギーの差は、電子が抜けた後の引力を含まないので、ここで補います。
- **$2(ia|jb)$**：一重項にだけ現れる交換型の項です。三重項にはないので、同じ軌道の組の励起では、ふつう三重項のほうが一重項より低くなります。

CIS は、変分法で求めた波動関数をもち、サイズ無矛盾です [Foresman 1992]。しかし、励起状態の電子相関を取り込まないので、励起エネルギーを大きく見積もる傾向があります（5 節の計算例でも、より精度の高い方法より 0.5 eV 高くなります）。

**ここまでで分かったこと**：CIS は、1 電子励起の配置の中でハミルトニアンを対角化します。行列は、軌道エネルギーの差と、電子と穴の相互作用の項からなります。

## 3. TDHF：時間変化への応答として求める

### 3.1 応答から励起エネルギーを求める考え方

分子に、ある振動数 $\omega$ で振動する弱い電場を当てると、電子の分布が同じ振動数で揺れます。この揺れの大きさ（分極率など）は、$\omega$ が励起エネルギーに一致したところで発散します。そこで、時間に依存する摂動に対する波動関数の応答（**線形応答**）を計算し、応答が発散する $\omega$ を探せば、励起エネルギーが求まります。

HF の波動関数の時間変化をこの考え方で扱うのが、**時間に依存する HF 法**（TDHF）です。乱雑位相近似（random phase approximation, RPA）とも呼ばれます。

### 3.2 行列の形

導出は省きますが、TDHF の励起エネルギーは、次の固有値問題で求まります [Dreuw 2005]。

$$
\begin{pmatrix} \mathbf A & \mathbf B \\ \mathbf B & \mathbf A \end{pmatrix}
\begin{pmatrix} \mathbf X \\ \mathbf Y \end{pmatrix}
= \omega
\begin{pmatrix} \mathbf 1 & \mathbf 0 \\ \mathbf 0 & -\mathbf 1 \end{pmatrix}
\begin{pmatrix} \mathbf X \\ \mathbf Y \end{pmatrix}
$$

$\mathbf A$ は CIS と同じ行列です。新しく現れる $\mathbf B$ は、一重項と三重項で次のとおりです。

$$
B_{ia,jb}^{\mathrm S} = 2(ia|jb) - (ib|ja), \qquad
B_{ia,jb}^{\mathrm T} = -(ib|ja)
$$

$\mathbf X$ は占有から仮想への励起の振幅、$\mathbf Y$ は仮想から占有への「脱励起」の振幅です。$\mathbf Y$ は、基底状態が HF の決定子からずれている分（基底状態の相関の一部）を表すと考えられます。$\mathbf B = \mathbf 0$ とおけば $\mathbf Y = \mathbf 0$ となり、CIS に戻ります。

軌道が実数のとき、この問題は次の形に書き直せます。

$$
(\mathbf A - \mathbf B)^{1/2} (\mathbf A + \mathbf B) (\mathbf A - \mathbf B)^{1/2}\, \mathbf Z = \omega^2 \mathbf Z
$$

求まるのは $\omega$ ではなく $\omega^2$ です。$\mathbf A + \mathbf B$ が負の固有値をもつと、$\omega^2 < 0$、つまり励起エネルギーが虚数になります。

### 3.3 三重項不安定性とのつながり

三重項の $\mathbf A^{\mathrm T} + \mathbf B^{\mathrm T}$ は、RHF の解が、スピン対称性を破った UHF の解へ崩れるかどうかを調べる行列と同じものです。第 7 回で使った安定性解析（`stability`）は、この行列の固有値が負かどうかを調べています。RHF の解が UHF に崩れる寸前（**三重項不安定性**に近い状態）では、この行列の固有値が 0 に近づき、TDHF の三重項の励起エネルギーが 0 に向かって下がります。5 節の計算例で、実際にこれを見ます。

**ここまでで分かったこと**：TDHF は、時間に依存する摂動への応答から励起エネルギーを求めます。CIS の行列に $\mathbf B$ が加わった形で、$\mathbf B = 0$ とおくと CIS に戻ります。三重項不安定性に近いと、三重項の励起エネルギーが下がりすぎます。

## 4. TDDFT：Runge–Gross の定理と断熱近似

### 4.1 Runge–Gross の定理

第 15 回では、基底状態の電子密度が外部ポテンシャルを決める、という Hohenberg–Kohn の定理を学びました。Runge と Gross は、これを時間に依存する系に広げました [Runge 1984]。時間に依存する電子密度 $\rho(\mathbf r, t)$ は、（決まった初期状態から出発したとき）時間に依存する外部ポテンシャルを、時間だけの関数の差を除いて一つに決めます。

これにより、時間に依存する問題も、Kohn–Sham 法と同じように、相互作用のない電子の系で扱えます。この理論を、**時間に依存する密度汎関数理論**（time-dependent DFT, TDDFT）と呼びます。

### 4.2 断熱近似

厳密な理論では、ある時刻の交換相関ポテンシャルは、それより前のすべての時刻の密度に依存します（記憶をもつ）。実際の計算では、これを無視し、各時刻の密度だけで決まるとします。つまり、基底状態の汎関数（第 16 回）を、各時刻の密度にそのまま使います。これを**断熱近似**と呼びます。

### 4.3 Casida の方程式

断熱近似のもとで、TDDFT の線形応答を TDHF と同じ形の行列の問題に書いたのが、Casida です [Casida 1995]。行列は TDHF と同じ形で、要素が次のように変わります。厳密な交換を割合 $c_\mathrm x$ だけ混ぜた混成汎関数で書きます（第 16 回）。

$$
A_{ia,jb}^{\mathrm S} = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + 2(ia|jb) - c_\mathrm x (ij|ab) + (ia|f_\mathrm{xc}^{\mathrm S}|jb)
$$

$$
B_{ia,jb}^{\mathrm S} = 2(ia|jb) - c_\mathrm x (ib|ja) + (ia|f_\mathrm{xc}^{\mathrm S}|jb)
$$

ここで、$\varepsilon_p$ は Kohn–Sham の軌道エネルギーです。$f_\mathrm{xc}$ は**交換相関カーネル**と呼ばれ、交換相関ポテンシャルを密度でもう一度微分したものです。

$$
(ia|f_\mathrm{xc}|jb) = \iint \phi_i(\mathbf r)\phi_a(\mathbf r)\, \frac{\delta^2 E_\mathrm{xc}}{\delta\rho(\mathbf r)\,\delta\rho(\mathbf r')}\, \phi_j(\mathbf r')\phi_b(\mathbf r')\, d\mathbf r\, d\mathbf r'
$$

ただし、この式の $E_\mathrm{xc}$ は、汎関数のうち半局所の部分だけです。一重項の $f_\mathrm{xc}^{\mathrm S}$ は、α スピンと β スピンの密度による微分を、同じスピンどうしと違うスピンどうしで足したものです。三重項では、$2(ia|jb)$ の項が消え、カーネルは同じスピンどうしの微分から違うスピンどうしの微分を引いたものになります。

$c_\mathrm x = 1$ で $f_\mathrm{xc}$ を 0 にすると TDHF に、$c_\mathrm x = 0$ にすると純粋な半局所汎関数の TDDFT になります。$\mathbf B = 0$ とおいた近似は、TDDFT の **Tamm–Dancoff 近似**（TDA）と呼ばれます [Hirata 1999]。TDA の TDHF は、CIS そのものです。

**ここまでで分かったこと**：TDDFT は、断熱近似のもとで、TDHF と同じ形の固有値問題になります。HF の交換の項の一部が、交換相関カーネルに置き換わります。

## 5. 計算例：ホルムアルデヒドの励起エネルギー

ホルムアルデヒド（H₂CO）の低い励起状態を、CIS、TDHF、TDDFT（TDA と完全な形）で計算し、比べます。基準には、同じ基底関数の EOM-CCSD（結合クラスター法を励起状態に広げた方法、第 10 回の CCSD が出発点）の値を使います。これは、この回の手法より高い精度をもつと考えられる方法です。基底関数は aug-cc-pVDZ です。

ホルムアルデヒドの低い励起状態には、次の 2 種類があります。

- **nπ\* 状態**：酸素の非共有電子対の軌道（n）から、C=O の反結合性の π\* 軌道への励起です。
- **ππ\* 状態**：C=O の結合性の π 軌道から π\* 軌道への励起です。

手法によって状態の並び順が入れ替わることがあるので、エネルギーの低い順ではなく、励起の種類をそろえて比べます。励起の種類は、各状態で最も大きな振幅をもつ「占有軌道 → 仮想軌道」の組の対称性で判定します。分子を yz 面に置くと、n 軌道は b₂、π と π\* 軌道は b₁ の対称性をもちます。

```python
import numpy as np
from pyscf import gto, scf, dft, tdscf, cc, symm
from pyscf.cc import eom_rccsd

# ホルムアルデヒド（C=O 1.21 Å、C-H 1.11 Å、H-C-H 116°、分子は yz 面）
ang = np.radians(116 / 2)
mol = gto.M(atom=f"""C 0 0 0; O 0 0 1.21
                     H 0 {1.11*np.sin(ang):.4f} {-1.11*np.cos(ang):.4f}
                     H 0 {-1.11*np.sin(ang):.4f} {-1.11*np.cos(ang):.4f}""",
            basis="aug-cc-pvdz", symmetry=True, verbose=0)
ev = 27.2114
# この向きでは、O の非共有電子対 n は b2、π と π* は b1 の対称性をもつ
KIND = {("B2", "B1"): "nπ*", ("B1", "B1"): "ππ*"}

def kind(mf, x):
    """励起の振幅 x（占有 × 仮想）で最大の成分から、励起の種類を判定する"""
    irr = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff)
    nocc = (mf.mo_occ > 0).sum()
    i, a = np.unravel_index(np.abs(x).argmax(), x.shape)
    return KIND.get((irr[i], irr[nocc + a]), "その他")

def states(mf, cls, singlet):
    td = cls(mf); td.singlet = singlet; td.nstates = 6; td.kernel()
    out = {}
    for e, (x, _) in zip(td.e, td.xy):
        out.setdefault(kind(mf, x), e * ev)     # 種類ごとに最も低いもの
    return out

rows = []
rhf = scf.RHF(mol).run()
for name, cls in [("CIS", tdscf.TDA), ("TDHF", tdscf.TDHF)]:
    rows.append((name, states(rhf, cls, True), states(rhf, cls, False)))
for xc in ["pbe", "pbe0", "camb3lyp"]:
    ks = dft.RKS(mol, xc=xc).run()
    for label, cls in [("TDA", tdscf.TDA), ("TDDFT", tdscf.TDDFT)]:
        rows.append((f"{label}-{xc.upper()}", states(ks, cls, True), states(ks, cls, False)))
mycc = cc.CCSD(rhf).run()
eom = {}
for singlet, E in [(True, eom_rccsd.EOMEESinglet), (False, eom_rccsd.EOMEETriplet)]:
    solver = E(mycc)
    es, vs = solver.kernel(nroots=6)
    out = {}
    for e, v in zip(es, vs):
        r1 = solver.vector_to_amplitudes(v)[0]
        out.setdefault(kind(rhf, r1), e * ev)
    eom[singlet] = out
rows.append(("EOM-CCSD", eom[True], eom[False]))

print("手法             1(nπ*)  3(nπ*)  3(ππ*)   （単位 eV）")
for name, s, t in rows:
    print(f"{name:15s} {s.get('nπ*', np.nan):6.2f}  {t.get('nπ*', np.nan):6.2f}  {t.get('ππ*', np.nan):6.2f}", flush=True)
```

1 コアで 4 分ほどかかりました。出力は次のとおりです。

```text
手法             1(nπ*)  3(nπ*)  3(ππ*)   （単位 eV）
CIS               4.55    3.73    4.89
TDHF              4.37    3.40    2.02
TDA-PBE           3.79    3.09    5.98
TDDFT-PBE         3.77    3.04    5.74
TDA-PBE0          3.94    3.20    5.71
TDDFT-PBE0        3.91    3.12    5.21
TDA-CAMB3LYP      3.91    3.23    5.76
TDDFT-CAMB3LYP    3.88    3.15    5.30
EOM-CCSD          4.01    3.55    5.99
```

この表から、次のことが分かります。

1. **一重項の nπ\***：CIS は 4.55 eV で、EOM-CCSD（4.01 eV）より 0.5 eV 高くなりました。2 節で述べた、CIS が励起エネルギーを大きく見積もる傾向のとおりです。TDDFT は 3.8〜3.9 eV で、EOM-CCSD より少し低い値でした。
2. **TDA と完全な形の差**：nπ\* 状態では、TDA と完全な TDDFT の差は 0.03〜0.08 eV で、ほぼ同じです。
3. **三重項の ππ\***：TDHF は 2.02 eV で、EOM-CCSD（5.99 eV）より 4 eV も低くなりました。このため TDHF では、三重項の中で ππ\* が nπ\* より下に来て、状態の順番が入れ替わっています。同じ RHF の解で、PySCF の安定性解析（RHF から UHF への安定性）を行うと、最も小さな固有値は正でしたが 0.017 と、0 にかなり近い値でした。3.3 節で見た、三重項不安定性に近いために三重項の励起エネルギーが下がりすぎる状況です。
4. **TDDFT の三重項 ππ\***：厳密な交換を含む PBE0 と CAM-B3LYP では、完全な TDDFT の値（5.21、5.30 eV）が TDA の値（5.71、5.76 eV）より 0.5 eV ほど低く、TDA のほうが EOM-CCSD に近くなりました。

3 と 4 は、文献で報告されている傾向と一致します。Peach らは、厳密な交換を含む TDDFT の三重項の励起エネルギーが大きく下がりすぎることがあり、それが基底状態の三重項不安定性で理解できることを示しました。さらに、TDA を使うと、問題のある三重項の励起エネルギーが大きく改善することを報告しています [Peach 2011]。

ただし、この比較は 1 個の分子の、少数の状態についてのものです。汎関数の良し悪しを、この表から一般的に結論することはできません。

**ここまでで分かったこと**：CIS は励起エネルギーを大きく見積もります。TDHF の三重項は、三重項不安定性に近いと大きく下がりすぎ、TDDFT でも厳密な交換を含むと同じ傾向が出ます。TDA はこの問題を和らげます。

## 6. TDDFT が苦手な場面

### 6.1 電荷移動励起

分子の一部（電子を与える側 D）から、離れた別の部分（電子を受け取る側 A）へ電子が移る励起を、**電荷移動励起**（charge-transfer excitation）と呼びます。D と A が十分に離れていると、正しい励起エネルギーは次のようになります。

$$
\omega_\mathrm{CT} \approx \mathrm{IP}_\mathrm D - \mathrm{EA}_\mathrm A - \frac1R
$$

$\mathrm{IP}_\mathrm D$ は D のイオン化エネルギー、$\mathrm{EA}_\mathrm A$ は A の電子親和力、$R$ は D と A の距離です。最後の項は、離れた正と負の電荷の間のクーロン引力です。

4.3 節の行列で、この励起を考えます。占有軌道 $i$ が D に、仮想軌道 $a$ が A にあるとします。二つの軌道はほとんど重ならないので、積 $\phi_i\phi_a$ を含む $(ia|jb)$ とカーネルの項はほぼ 0 です。残るのは、

$$
\omega_\mathrm{CT} \approx \varepsilon_a - \varepsilon_i - c_\mathrm x\, (ii|aa), \qquad (ii|aa) \approx \frac1R
$$

です。$(ii|aa)$ は、D にある電子の分布と A にある電子の分布の間のクーロン反発なので、離れると $1/R$ に近づきます。

半局所汎関数（$c_\mathrm x = 0$）では、$-1/R$ の項がまったく現れません。さらに、半局所汎関数の軌道エネルギーの差 $\varepsilon_a - \varepsilon_i$ は、$\mathrm{IP}_\mathrm D - \mathrm{EA}_\mathrm A$ よりかなり小さいことが多いので（第 18 回の自己相互作用誤差）、電荷移動励起のエネルギーは大きく低く見積もられます。Dreuw らは、離れた電荷の間の引力が HF の非局所的な交換に由来し、局所的な汎関数の TDDFT ではこれが取り込めないことを示しました [Dreuw 2003]。

この問題への代表的な対処が、遠くで厳密な交換の割合を大きくする**範囲分離混成汎関数**です（第 16 回）。たとえば CAM-B3LYP は、厳密な交換の割合を、近距離の 19% から遠距離の 65% まで、電子間距離に応じて変えます [Yanai 2004]。

### 6.2 二電子励起

断熱近似の TDDFT で求まる励起状態は、1 電子励起の配置に対応するものだけです [Maitra 2004]。行列の大きさが 1 電子励起の数で決まるからです。二つの電子が同時に励起した性格の強い状態（**二電子励起**の状態）は、そもそも解として現れないか、正しく表せません。Maitra らは、二電子励起の状態の近くでは、厳密な交換相関カーネルが振動数に強く依存することを示しました [Maitra 2004]。断熱近似は、この振動数への依存を無視しています。

二電子励起の性格は、ポリエンの暗い励起状態や、結合が切れかけた分子の励起状態などで重要になります。

### 6.3 基底状態が多参照的な場合

この回の方法は、基底状態が決定子 1 個でよく表せることを前提にしていました。ジラジカルや、結合が切れかけた分子のように、基底状態そのものが多参照的なとき（第 7 回の静的相関）、出発点が悪いので、励起状態もうまく求まりません。

光化学で重要な、基底状態と励起状態が交わる点（円錐交差、[第 22 回](22-conical-intersections.md)）の近くでは、基底状態が強く多参照的になります。この問題に、決定子 1 個の出発点の安さを保ったまま対処しようとする方法が、次回の**スピン反転 TDDFT**です（[第 20 回](20-spin-flip-tddft.md)）。

**ここまでで分かったこと**：TDDFT は、半局所汎関数では電荷移動励起を大きく低く見積もり、断熱近似のために二電子励起を表せず、基底状態が多参照的な場合にも弱点があります。

## まとめ

- CIS は、1 電子励起の配置の中でハミルトニアンを対角化します。励起状態の相関を含まないので、励起エネルギーを大きく見積もりがちです。
- TDHF と TDDFT は、時間に依存する摂動への線形応答から励起エネルギーを求め、同じ形の行列の固有値問題になります。$\mathbf B = 0$ とおくのが Tamm–Dancoff 近似で、TDHF の TDA は CIS です。
- TDDFT は、Runge–Gross の定理を基礎にし、実際の計算では断熱近似を使います。HF の交換の項の一部が、交換相関カーネルに置き換わります。
- 三重項不安定性に近い系では、TDHF と厳密な交換を含む TDDFT の三重項が下がりすぎます。TDA はこれを和らげます。
- 半局所汎関数の TDDFT は電荷移動励起に弱く、断熱近似の TDDFT は二電子励起を表せません。

## 付録 A　CIS の行列要素の導出

閉殻の HF の決定子 $\Phi_0$ から、空間軌道 $i$ の α 電子を空間軌道 $a$ の α スピン軌道に移した決定子を $\Phi_{i\alpha}^{a\alpha}$、β 電子について同じものを $\Phi_{i\beta}^{a\beta}$ と書きます。Slater–Condon の規則（第 8 回）から、同じスピンどうしの要素は

$$
\langle \Phi_{i\alpha}^{a\alpha} | \hat H - E_0 | \Phi_{j\alpha}^{b\alpha} \rangle = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + (ia|jb) - (ij|ab)
$$

です。$E_0$ は HF のエネルギーです。違うスピンどうしでは、交換の項が消えて、

$$
\langle \Phi_{i\alpha}^{a\alpha} | \hat H | \Phi_{j\beta}^{b\beta} \rangle = (ia|jb)
$$

になります。一重項の配置は $(\Phi_{i\alpha}^{a\alpha} + \Phi_{i\beta}^{a\beta})/\sqrt2$、三重項（$M_S = 0$ の成分）の配置は $(\Phi_{i\alpha}^{a\alpha} - \Phi_{i\beta}^{a\beta})/\sqrt2$ です（第 5 回）。これらで行列要素をとると、

$$
A^{\mathrm S} = \big[\delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + (ia|jb) - (ij|ab)\big] + (ia|jb)
$$

$$
A^{\mathrm T} = \big[\delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + (ia|jb) - (ij|ab)\big] - (ia|jb)
$$

となり、本文の式が得られます。

## 参考文献

- A. Szabo, N. S. Ostlund, "Modern Quantum Chemistry" (Dover, 1996). ISBN: 9780486691862
- J. B. Foresman, M. Head-Gordon, J. A. Pople, M. J. Frisch, "Toward a systematic molecular orbital theory for excited states", J. Phys. Chem. 96, 135 (1992). DOI: [10.1021/j100180a030](https://doi.org/10.1021/j100180a030)
- E. Runge, E. K. U. Gross, "Density-Functional Theory for Time-Dependent Systems", Phys. Rev. Lett. 52, 997 (1984). DOI: [10.1103/PhysRevLett.52.997](https://doi.org/10.1103/PhysRevLett.52.997)
- M. E. Casida, "Time-Dependent Density Functional Response Theory for Molecules", in Recent Advances in Density Functional Methods, Part I (World Scientific, 1995), p. 155. DOI: [10.1142/9789812830586_0005](https://doi.org/10.1142/9789812830586_0005)
- S. Hirata, M. Head-Gordon, "Time-dependent density functional theory within the Tamm–Dancoff approximation", Chem. Phys. Lett. 314, 291 (1999). DOI: [10.1016/S0009-2614(99)01149-5](https://doi.org/10.1016/S0009-2614(99)01149-5)
- A. Dreuw, J. L. Weisman, M. Head-Gordon, "Long-range charge-transfer excited states in time-dependent density functional theory require non-local exchange", J. Chem. Phys. 119, 2943 (2003). DOI: [10.1063/1.1590951](https://doi.org/10.1063/1.1590951)
- A. Dreuw, M. Head-Gordon, "Single-Reference ab Initio Methods for the Calculation of Excited States of Large Molecules", Chem. Rev. 105, 4009 (2005). DOI: [10.1021/cr0505627](https://doi.org/10.1021/cr0505627)
- N. T. Maitra, F. Zhang, R. J. Cave, K. Burke, "Double excitations within time-dependent density functional theory linear response", J. Chem. Phys. 120, 5932 (2004). DOI: [10.1063/1.1651060](https://doi.org/10.1063/1.1651060)
- T. Yanai, D. P. Tew, N. C. Handy, "A new hybrid exchange–correlation functional using the Coulomb-attenuating method (CAM-B3LYP)", Chem. Phys. Lett. 393, 51 (2004). DOI: [10.1016/j.cplett.2004.06.011](https://doi.org/10.1016/j.cplett.2004.06.011)
- M. J. G. Peach, M. J. Williamson, D. J. Tozer, "Influence of Triplet Instabilities in TDDFT", J. Chem. Theory Comput. 7, 3578 (2011). DOI: [10.1021/ct200651r](https://doi.org/10.1021/ct200651r)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
