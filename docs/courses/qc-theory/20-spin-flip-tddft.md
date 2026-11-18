# 第 20 回　スピン反転 TDDFT

!!! abstract "この回の道筋"
    第 19 回の TDDFT は、閉殻の基底状態から出発するので、ジラジカルや結合が切れかけた分子のように基底状態そのものが多参照的な系を苦手とします。この回では、高スピンの三重項を出発点にし、電子のスピンを反転させる励起で一重項の基底状態と励起状態をまとめて求める**スピン反転法**を学びます。行列を自分で組み立てて、H₂ の解離とエチレンのねじれを計算し、利点を確かめます。最後に、スピン反転法の弱点であるスピン汚染を見ます。これが次回の MRSF-TDDFT の出発点になります。

    | 節 | 内容 |
    |---|---|
    | 1 | 閉殻から出発すると困る系 |
    | 2 | スピン反転の考え方 |
    | 3 | スピン反転 TDDFT の行列 |
    | 4 | 計算例：H₂ の解離 |
    | 5 | 計算例：エチレンのねじれ |
    | 6 | 弱点：スピン汚染 |

    **前に読んでおく回**：[量子化学の理論 第 5 回](05-spin.md)、[第 7 回](07-electron-correlation.md)、[第 12 回](12-casscf.md)、[第 19 回](19-tddft.md)

!!! tip "物理的なイメージ"
    スピン反転法が求めているのは、高スピンの三重項を**足場**にして、電子 1 個のスピンを上向き（α）から下向き（β）に返したときにできる状態のエネルギーです。開殻の 2 個の電子のスピンが同じ向きにそろった三重項は、パウリの原理で電子どうしが避け合うので、決定子 1 個でよく表せる安定な足場です。そこからスピンを 1 個返すと、2 個の電子が「両方とも下の軌道」「両方とも上の軌道」「1 個ずつ」に入った配置が、すべて同じ立場で現れます。一重項の基底状態も励起状態も、これらの混ざり方として、足場からの「励起」の一つとして求まります。閉殻の基底状態を特別な出発点にしないことが、結合が切れかけた系での強みです。

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書きます。1 $E_\mathrm{h}$ = 27.211 eV = 627.5 kcal/mol です。2 電子積分は、空間軌道についての化学者の表記 $(pq|rs)$ を使い、軌道はすべて実数とします。

## 1. 閉殻から出発すると困る系

エチレンの C=C 結合のまわりで、片方の CH₂ を 90° ねじることを考えます。平面の構造では、π 軌道に 2 個の電子が入った閉殻の配置 $\pi^2$ がよい近似です。ねじると π 結合が切れ、90° では π 軌道と π\* 軌道のエネルギーが等しくなります。このとき一重項の基底状態は、$\pi^2$ と $\pi^{*2}$ が同じ重みで混ざった状態になり、決定子 1 個では表せません（第 7 回の H₂ の解離と同じ状況です）。

このような系では、閉殻の RKS（または RHF）は、出発点として質的に間違っています。スピン対称性を破った UKS にすればエネルギーは下がりますが、一重項と三重項が混ざった状態になります（第 7 回）。第 19 回の TDDFT は閉殻の基底状態から出発するので、同じ問題を抱えます。

多配置の CASSCF（第 12 回）を使えば正しく扱えますが、活性空間を選ぶ必要があり、動的相関を足すには CASPT2 や NEVPT2（第 13 回）が必要です。決定子 1 個から出発する方法の手軽さを保ったまま、この問題に対処できないでしょうか。

**ここまでで分かったこと**：結合が切れかけた系やジラジカルでは、閉殻の決定子 1 個は出発点として質的に間違っています。

## 2. スピン反転の考え方

### 2.1 高スピンの三重項から出発する

Krylov は、次の点に注目しました [Krylov 2001a; Krylov 2001b]。ねじれたエチレンの π と π\* に電子が 1 個ずつ、同じ向きのスピン（α）で入った高スピンの三重項（$M_S = +1$）は、決定子 1 個でよく表せます。同じスピンの電子どうしはパウリの原理で避け合うので、静的相関の問題が小さいからです [Casanova 2020]。

そこで、この三重項の決定子を出発点（参照）にし、α 電子 1 個を β 電子に変える励起（**スピン反転**、spin flip）を考えます。スピン反転を 1 回すると、$M_S = 0$ の配置ができます。$M_S = 0$ は、一重項と三重項の両方がもつ成分です。

### 2.2 どんな配置ができるか

二つの開殻の軌道を $u$ と $v$（エチレンでは π と π\*）とします。参照は、閉殻の軌道に電子が 2 個ずつ入り、$u$ と $v$ に α 電子が 1 個ずつ入った決定子です。$u$ と $v$ の中だけでスピンを反転させると、次の 4 つの配置ができます。

| 反転のしかた | できる配置 |
|---|---|
| $u$ の α 電子を、$v$ の β に移す | $v^2$（閉殻） |
| $v$ の α 電子を、$u$ の β に移す | $u^2$（閉殻） |
| $u$ の α 電子を、$u$ の β に変える | $u$ に β、$v$ に α（開殻） |
| $v$ の α 電子を、$v$ の β に変える | $u$ に α、$v$ に β（開殻） |

二つの閉殻の配置 $u^2$ と $v^2$ が、両方とも入っています。これらを混ぜれば、ねじれたエチレンの一重項の基底状態（$\pi^2$ と $\pi^{*2}$ の混合）を表せます。開殻の 2 つの配置は、足すと三重項（$M_S = 0$ の成分）、引くと開殻の一重項になります（第 5 回）。

つまり、決定子 1 個から 1 電子の励起だけを考えるのに、CAS(2,2)（第 12 回）と同じ配置がすべて得られます。基底状態は、参照からの「励起」の一つとして、励起状態と同じ扱いで求まります。

**ここまでで分かったこと**：高スピンの三重項を参照にし、スピンを 1 回反転させると、二つの閉殻の配置と開殻の配置が得られます。多参照的な一重項の基底状態を、決定子 1 個からの励起として表せます。

## 3. スピン反転 TDDFT の行列

### 3.1 行列の形

第 19 回の CIS と同じように、スピン反転の配置の中でハミルトニアンを対角化します。占有軌道 $i, j$ の α 電子を、β の空いた軌道 $a, b$ に移す配置を $\Phi_{i\alpha}^{a\beta}$ と書きます。参照の決定子を $\Phi_0$ とすると、行列の要素は次のとおりです（Slater–Condon の規則から導けます）。

$$
A_{ia,jb} = \langle \Phi_{i\alpha}^{a\beta} | \hat H - E_0 | \Phi_{j\alpha}^{b\beta} \rangle = \delta_{ij} F_{ab}^\beta - \delta_{ab} F_{ji}^\alpha - (ij|ab)
$$

$F^\alpha$ と $F^\beta$ は、参照の決定子からつくった α と β の Fock 行列です。添字 $i, j$ は α の占有軌道（閉殻と開殻の軌道）、$a, b$ は β の空いた軌道（開殻と仮想の軌道）を動きます。参照が ROHF のとき、これを **SF-CIS** と呼びます [Krylov 2001b]。

第 19 回の CIS の行列と比べると、$(ia|jb)$ の項がありません。この項は、軌道 $i$ と $a$ の積の分布のクーロン相互作用ですが、$i$ は α、$a$ は β なので、スピンの積分で 0 になるからです。

### 3.2 TDDFT にすると

Shao、Head-Gordon、Krylov は、これを TDDFT の Tamm–Dancoff 近似（第 19 回）で実現しました [Shao 2003]。参照には高スピンの三重項の Kohn–Sham の決定子を使います。

ここで、交換相関カーネルが問題になります。ふつうの TDDFT で使う汎関数は、α と β の密度の関数です（**共線形**、collinear と呼びます）。スピンの向きを変える応答では、α と β の密度がどちらも変わらないので、半局所の汎関数のカーネルは寄与しません [Shao 2003; Wang 2004]。行列に残るのは、厳密な交換の分だけです。

イメージで言えば、スピンを反転させる応答は、電子のスピンの向きを、上下（α と β）の軸から横に傾けることにあたります。1 次の変化は、α と β の密度ではなく、α と β の間の「横向きのスピン」に現れます。共線形の汎関数は、上向きと下向きの電子の量しか見ていないので、スピンが横に傾いても何も感じません。そのため、傾きを元に戻す力、つまり第 19 回 3.4 節の「ばね定数」を与えません。この傾きに対するばねを与えるのは、厳密な交換だけです。

$$
A_{ia,jb} = \delta_{ij} F_{ab}^\beta - \delta_{ab} F_{ji}^\alpha - c_\mathrm x (ij|ab)
$$

$F^\alpha$、$F^\beta$ は Kohn–Sham の Fock 行列、$c_\mathrm x$ は汎関数に混ぜた厳密な交換の割合です。$c_\mathrm x$ が小さいと、配置どうしの結合がほとんどなくなり、正しい混ざり方が得られません。Shao らは、厳密な交換 50% の混成汎関数で特によい結果が得られることを示しました [Shao 2003]。この回では、厳密な交換 50% と B88 交換 50% に LYP 相関を組み合わせた **BHHLYP** を使います。

カーネルの問題に対しては、スピンの向きの変化も扱える**非共線形**（noncollinear）のカーネルを使う方法も提案されています [Wang 2004]。

### 3.3 実装

この回の計算では、上の行列を NumPy で組み立てて対角化します。参照には ROHF または ROKS を使います。プログラム `sf.py` は付録 B に載せます。各状態の $\langle \hat S^2 \rangle$ も計算します。その式は付録 A で導きます。

**ここまでで分かったこと**：スピン反転の行列は、CIS の行列から $(ia|jb)$ の項を除いた形です。共線形の TDDFT では、半局所の汎関数が反転の結合に寄与しないので、厳密な交換を多く混ぜた汎関数を使います。

## 4. 計算例：H₂ の解離

まず、2 電子の H₂ で、スピン反転の考え方を確かめます。参照は、σg と σu に α 電子が 1 個ずつ入った三重項の ROHF です。SF-CIS の一重項の基底状態を、RHF と、厳密な FCI と比べます。基底関数は cc-pVDZ です。

```python
from pyscf import gto, scf, fci
import numpy as np
from sf import sf_tda

print("R/Å     RHF        SF-CIS    (<S^2>)    FCI       （全エネルギー、E_h）")
for r in (0.74, 1.5, 2.5, 4.0, 10.0):
    atom = f"H 0 0 0; H 0 0 {r}"
    m0 = gto.M(atom=atom, basis="cc-pvdz", verbose=0)             # 一重項
    m2 = gto.M(atom=atom, basis="cc-pvdz", spin=2, verbose=0)     # 参照の三重項
    rhf = scf.RHF(m0).run()
    e_fci = fci.FCI(rhf).kernel()[0]
    e_sf, s2, e_ref = sf_tda(m2, None, nstates=4)
    k = int(np.argmax(s2 < 1.0))                                  # <S^2> で一重項を選ぶ
    print(f"{r:5.2f}  {rhf.e_tot:9.5f}  {e_sf[k]:9.5f}  ({s2[k]:.3f})  {e_fci:9.5f}")
```

```text
R/Å     RHF        SF-CIS    (<S^2>)    FCI       （全エネルギー、E_h）
 0.74   -1.12870   -1.13356  (0.016)   -1.16337
 1.50   -1.00219   -1.05191  (0.005)   -1.06153
 2.50   -0.86533   -1.00260  (0.000)   -1.00313
 4.00   -0.78220   -0.99860  (0.000)   -0.99861
10.00   -0.73384   -0.99856  (0.016)   -0.99856
```

RHF は、結合を伸ばすとエネルギーが上がり続け、10 Å では FCI より 0.26 $E_\mathrm h$（約 7 eV）も高くなります。閉殻の配置 $\sigma_g^2$ には、イオン的な配置（H⁺H⁻）が半分含まれているからです（第 7 回）。

SF-CIS は、どの距離でも RHF よりずっと低く、2.5 Å より長くなると FCI とほぼ一致しました。10 Å では、FCI と小数点以下 5 桁まで一致しています。離れた 2 個の H 原子は、それぞれ電子を 1 個しかもたないので、σg と σu の二つの配置を正しく混ぜれば厳密な答えになるからです。

一方、平衡距離の 0.74 Å では、SF-CIS は FCI より 0.03 $E_\mathrm h$ 高くなりました。SF-CIS は 1 電子の励起（反転）しか含まないので、動的相関をほとんど取り込めないからです。

また、平衡距離では、一重項の $\langle \hat S^2 \rangle$ が 0.016 と、0 からわずかにずれています。この原因は 6 節で説明します。

なお、この計算では、スピン反転で作った三重項（$M_S = 0$ の成分）のエネルギーが、どの距離でも参照の三重項（$M_S = +1$）のエネルギーと一致することも確かめました。同じ三重項の 2 つの成分なので、一致するのが正しい結果です。

**ここまでで分かったこと**：SF-CIS は、閉殻の配置と開殻の配置を混ぜることで、H₂ の解離を質的に正しく表します。ただし、動的相関はほとんど含みません。

## 5. 計算例：エチレンのねじれ

次に、エチレンの片方の CH₂ をねじったときの、一重項の基底状態のエネルギーを計算します。Shao らが SF-TDDFT の最初の論文で調べた系です [Shao 2003]。構造は C=C 1.34 Å、C–H 1.09 Å、H–C–H 117° の平面構造で、ねじるとき以外の構造は変えません。基底関数は cc-pVDZ、汎関数は BHHLYP です。

次の 4 つを比べます。

- **RKS**：閉殻の Kohn–Sham 法。
- **UKS**：スピン対称性を破った Kohn–Sham 法。安定性解析で、より低い解を探します。
- **SF-TDA**：三重項の ROKS を参照にした、スピン反転の TDA。一重項を $\langle \hat S^2 \rangle$ で選びます。90° 付近では、三重項が一重項とほぼ同じエネルギーになるからです。
- **NEVPT2**：基準として、π と π\* を活性空間にした CAS(2,2) から NEVPT2（第 13 回）で動的相関を足したもの [Angeli 2001]。

```python
import numpy as np
from pyscf import gto, scf, dft, mcscf, mrpt
from sf import sf_tda

def ethylene(theta):
    """C=C 1.34 Å、C-H 1.09 Å、H-C-H 117° の平面構造から、片方の CH2 を角度 theta だけねじる"""
    hx, hz = 1.09 * np.sin(np.radians(117 / 2)), 1.09 * np.cos(np.radians(117 / 2))
    t = np.radians(theta)
    return (f"C 0 0 0.67; C 0 0 -0.67; "
            f"H {hx} 0 {0.67 + hz}; H {-hx} 0 {0.67 + hz}; "
            f"H {hx*np.cos(t)} {hx*np.sin(t)} {-0.67 - hz}; H {-hx*np.cos(t)} {-hx*np.sin(t)} {-0.67 - hz}")

if __name__ == "__main__":
    xc = "bhandhlyp"     # 厳密な交換 50% の混成汎関数（BHHLYP）
    kcal = 627.5
    rows = []
    for theta in (0, 30, 60, 75, 85, 90):
        m0 = gto.M(atom=ethylene(theta), basis="cc-pvdz", verbose=0)
        m2 = gto.M(atom=ethylene(theta), basis="cc-pvdz", spin=2, verbose=0)
        e_rks = dft.RKS(m0, xc=xc).run().e_tot
        uks = dft.UKS(m0, xc=xc).run()                       # スピン対称性を破った解を探す
        mo, _, stable, _ = uks.stability(return_status=True)
        while not stable:
            uks.kernel(uks.make_rdm1(mo, uks.mo_occ))
            mo, _, stable, _ = uks.stability(return_status=True)
        e_sf, s2_sf, _ = sf_tda(m2, xc, nstates=4)
        k = int(np.argmax(s2_sf < 1.0))                      # <S^2> で一重項を選ぶ
        t_sf = e_sf[int(np.argmax(s2_sf > 1.0))]             # 三重項（M_S = 0 の成分）
        cas = mcscf.CASSCF(scf.RHF(m0).run(), 2, 2).run()    # 基準：CAS(2,2) + NEVPT2
        e_nev = cas.e_tot + mrpt.NEVPT(cas).kernel()
        rows.append((theta, e_rks, uks.e_tot, uks.spin_square()[0], e_sf[k], s2_sf[k], e_nev, t_sf))

    r0 = rows[0]
    print("ねじれ角  RKS    UKS (<S^2>)     SF-TDA (<S^2>)   NEVPT2   （平面からのエネルギー、kcal/mol）")
    for th, er, eu, su, es, ss, en, _ in rows:
        print(f"{th:4d}   {(er-r0[1])*kcal:6.1f}  {(eu-r0[2])*kcal:6.1f} ({su:4.2f})   {(es-r0[4])*kcal:6.1f} ({ss:4.2f})   {(en-r0[6])*kcal:6.1f}")
    print(f"90°の SF-TDA：E(三重項) - E(一重項) = {(rows[-1][7] - rows[-1][4]) * kcal:.1f} kcal/mol")
```

1 コアで 2 分ほどかかりました。出力は次のとおりです（括弧の中は $\langle \hat S^2 \rangle$）。

```text
ねじれ角  RKS    UKS (<S^2>)     SF-TDA (<S^2>)   NEVPT2   （平面からのエネルギー、kcal/mol）
   0      0.0     0.0 (0.00)      0.0 (0.01)      0.0
  30     11.8    11.8 (0.00)     11.8 (0.01)     11.3
  60     46.6    44.8 (0.37)     44.8 (0.01)     43.6
  75     72.2    58.9 (0.83)     66.4 (0.01)     64.4
  85     92.0    63.6 (0.99)     78.2 (0.01)     77.2
  90    102.6    64.2 (1.02)     80.2 (0.01)     77.0
90°の SF-TDA：E(三重項) - E(一重項) = -0.4 kcal/mol
```

この表から、次のことが分かります。

1. **RKS**：ねじれ角が大きくなるにつれて基準から離れ、90° では 102.6 kcal/mol と、NEVPT2（77.0 kcal/mol）より 26 kcal/mol 高くなりました。閉殻の配置しか使えないためです。
2. **UKS**：60° 付近からスピン対称性を破り、90° では 64.2 kcal/mol と低くなりすぎました。$\langle \hat S^2 \rangle$ は 1.02 で、一重項（0）と三重項（2）がほぼ半分ずつ混ざった状態です。ねじれたエチレンでは三重項が一重項より低いので、三重項が混ざるとエネルギーが下がりすぎます。
3. **SF-TDA**：どのねじれ角でも $\langle \hat S^2 \rangle$ は 0.01 で、ほぼ純粋な一重項です。90° の値は 80.2 kcal/mol で、NEVPT2 との差は 3 kcal/mol です。85° から 90° への増え方は 2.0 kcal/mol で、RKS（10.6 kcal/mol）のように 90° に向かって急に上がり続けることはありません。
4. 90° の SF-TDA では、三重項が一重項より 0.4 kcal/mol 低くなりました。ねじれたエチレンでは、二つの状態がほぼ縮退しています。

Shao らも、ねじれのポテンシャルで、SF-DFT が従来の Kohn–Sham 法より大きく改善し、厳密な交換 50% の汎関数で特によいことを報告しています [Shao 2003]。

**ここまでで分かったこと**：エチレンのねじれでは、RKS は高すぎ、UKS はスピン汚染のために低すぎました。SF-TDA は純粋な一重項を保ち、基準に近い滑らかな曲線を与えました。

## 6. 弱点：スピン汚染

### 6.1 励起状態の $\langle \hat S^2 \rangle$

基底状態ではうまくいきましたが、励起状態ではどうでしょうか。平面（0°）とねじれ（90°）のエチレンで、SF-TDA の低い 8 個の状態の $\langle \hat S^2 \rangle$ を見ます。

```python
from pyscf import gto
from sf import sf_tda
from ethylene import ethylene
for theta in (0, 90):
    m2 = gto.M(atom=ethylene(theta), basis="cc-pvdz", spin=2, verbose=0)
    e, s2, eref = sf_tda(m2, "bhandhlyp", nstates=8)
    print(f"ねじれ角 {theta}°：参照の三重項からのエネルギー / eV と <S^2>")
    for k in range(8):
        print(f"  {k+1}  {(e[k] - eref) * 27.2114:7.2f}   {s2[k]:5.2f}")
```

```text
ねじれ角 0°：参照の三重項からのエネルギー / eV と <S^2>
  1    -4.16    0.01
  2     0.71    2.00
  3     3.87    1.01
  4     3.95    0.01
  5     4.16    1.00
  6     4.69    1.00
  7     5.15    1.00
  8     5.27    1.00
ねじれ角 90°：参照の三重項からのエネルギー / eV と <S^2>
  1     0.62    2.00
  2     0.64    0.01
  3     2.76    0.06
  4     2.84    0.06
  5     6.89    1.01
  6     6.89    1.01
  7     7.73    1.06
  8     7.73    1.06
```

$\langle \hat S^2 \rangle$ が 0 なら一重項、2 なら三重項です。ところが、1.00 前後の状態がたくさん現れました。これらは、一重項と三重項がほぼ半分ずつ混ざった状態で、どちらとも言えません。このような状態は、エネルギーも、その状態が何であるかも信頼できません。

### 6.2 原因：スピンについて不完全な配置の空間

原因は、スピン反転でつくる配置の集まりにあります。閉殻の軌道 $k$ から仮想の軌道 $c$ への励起を考えます。スピン反転では、$k$ の α 電子を $c$ の β に移す配置はつくれます。ところが、スピンの純粋な状態をつくるには、$k$ の β 電子を $c$ の α に移し、同時に開殻の電子のスピンを変えるような配置も必要になります。この相手は、参照から 1 回のスピン反転では届かないので、配置の空間に入っていません。

そのため、閉殻から仮想への励起や、開殻から仮想への励起の性格をもつ状態は、相手を欠いたまま、一重項と三重項が混ざった状態になります。これを、配置の空間が**スピンについて不完全**（spin-incomplete）であると言います。4 節の H₂ の平衡距離で、基底状態の $\langle \hat S^2 \rangle$ がわずかに 0 からずれたのも、開殻から仮想への励起の配置が少し混ざったためです。

Sears、Sherrill、Krylov は、欠けている配置を加えてスピンについて完全にした SF-CIS を調べ、スピンの固有関数にすることで精度が大きく改善することを示しました [Sears 2003]。

### 6.3 もう一つの不整合

DFT のスピン反転には、もう一つ気になる点があります。6.1 節の出力で、$\langle \hat S^2 \rangle = 2$ の状態（スピン反転でつくった三重項の $M_S = 0$ の成分）のエネルギーは、0° で参照の三重項（$M_S = +1$）より 0.71 eV、90° で 0.62 eV 高くなりました。同じ三重項の 2 つの成分なので、本来は一致するはずです。4 節の SF-CIS（HF）では、実際に一致していました。

この違いは、3.2 節で見た共線形の扱いから生じます。参照の三重項のエネルギーには、半局所の汎関数の交換相関がすべて入っています。一方、スピンを反転させる結合には、厳密な交換の分しか入りません。そのため、二つの成分で、交換相関の扱いがそろわないのです。

**ここまでで分かったこと**：スピン反転の配置の空間はスピンについて不完全で、励起状態に一重項と三重項が混ざった状態が多く現れます。DFT では、三重項の 2 つの成分のエネルギーも一致しません。次回の MRSF-TDDFT（[第 21 回](21-mrsf-tddft.md)）は、スピン汚染の問題に対処する方法です。

## まとめ

- 結合が切れかけた系やジラジカルでは、閉殻の決定子 1 個は出発点として質的に間違っています。
- スピン反転法は、決定子 1 個でよく表せる高スピンの三重項を参照にし、スピンを反転させる 1 電子励起で、一重項の基底状態と励起状態をまとめて求めます。二つの閉殻の配置を同時に含みます。
- 共線形の SF-TDDFT では、半局所の汎関数が反転の結合に寄与しないので、BHHLYP のように厳密な交換の多い汎関数を使います。
- H₂ の解離とエチレンのねじれで、SF 法は閉殻の方法やスピン対称性を破った方法の失敗を正しました。
- 配置の空間がスピンについて不完全なため、励起状態にスピン汚染が現れます。

## 付録 A　スピン反転の状態の $\langle \hat S^2 \rangle$

参照 $\Phi_0$ を高スピンの ROHF/ROKS の決定子（$S = M_S = 1$）とし、スピン反転の状態を

$$
\Psi = \sum_{ia} c_{ia}\, \hat a_{a\beta}^\dagger \hat a_{i\alpha} \Phi_0
$$

と書きます。$\Psi$ は $M_S = 0$ なので、$\hat S^2 = \hat S_- \hat S_+ + \hat S_z^2 + \hat S_z$ から、

$$
\langle \Psi | \hat S^2 | \Psi \rangle = \langle \Psi | \hat S_- \hat S_+ | \Psi \rangle = \| \hat S_+ \Psi \|^2
$$

です。$\hat S_+ = \sum_p \hat a_{p\alpha}^\dagger \hat a_{p\beta}$ と、$\Phi_0$ が $M_S$ の最大の状態なので $\hat S_+ \Phi_0 = 0$ であることを使います。交換関係 $[\hat S_+, \hat a_{a\beta}^\dagger] = \hat a_{a\alpha}^\dagger$、$[\hat S_+, \hat a_{i\alpha}] = -\hat a_{i\beta}$ から、

$$
\hat S_+ \hat a_{a\beta}^\dagger \hat a_{i\alpha} \Phi_0 = \hat a_{a\alpha}^\dagger \hat a_{i\alpha} \Phi_0 - \hat a_{a\beta}^\dagger \hat a_{i\beta} \Phi_0
$$

となります。右辺は、スピンを変えない 1 電子励起です。閉殻の軌道を C、開殻の軌道を O、仮想の軌道を V で表すと、0 にならないのは次の 3 種類です。

1. $i = a$ がどちらも開殻の同じ軌道のとき、第 1 項は $\Phi_0$ そのものになります。
2. $a$ が仮想の軌道のとき、第 1 項は α の励起 $\Phi_{i\alpha}^{a\alpha}$ になります。
3. $i$ が閉殻の軌道のとき、第 2 項は β の励起 $\Phi_{i\beta}^{a\beta}$ になります。

これらの決定子は互いに直交するので、

$$
\langle \hat S^2 \rangle = \Big( \sum_{u \in \mathrm O} c_{uu} \Big)^2 + \sum_{i \in \mathrm{C,O}}\ \sum_{a \in \mathrm V} c_{ia}^2 + \sum_{i \in \mathrm C}\ \sum_{a \in \mathrm{O,V}} c_{ia}^2
$$

が得られます。2.2 節の 4 つの配置だけからなる状態では、第 2 項と第 3 項が 0 になります。第 1 項は、開殻の 2 つの配置の係数が等しければ 2（三重項）、符号が逆なら 0（開殻の一重項）です。閉殻の配置 $u^2$、$v^2$ の係数は現れません。第 2 項と第 3 項が、6.2 節で見た、相手を欠いた配置による汚染です。

## 付録 B　スピン反転 TDA のプログラム

本文の計算で使った `sf.py` です。PySCF 2.15 で実行しました。

```python
import numpy as np
from pyscf import scf, dft, ao2mo

def sf_tda(mol, xc=None, nstates=4):
    """高スピン三重項（M_S = +1）の ROHF/ROKS を参照にした、共線形のスピン反転 TDA。
    xc=None なら ROHF を参照にした SF-CIS。全エネルギー、<S^2>、参照のエネルギーを返す。"""
    if xc is None:
        ref, cx = scf.ROHF(mol).run(), 1.0
        dm = ref.make_rdm1()
        vhf = scf.UHF(mol).get_veff(mol, dm)
    else:
        ref = dft.ROKS(mol, xc=xc).run()
        cx = ref._numint.hybrid_coeff(xc)        # 厳密な交換の割合
        dm = ref.make_rdm1()
        u = dft.UKS(mol, xc=xc); u.grids = ref.grids
        vhf = u.get_veff(mol, dm)
    C = ref.mo_coeff
    h = mol.intor("int1e_kin") + mol.intor("int1e_nuc")
    fa = C.T @ (h + vhf[0]) @ C                   # α の Fock 行列（分子軌道の基底）
    fb = C.T @ (h + vhf[1]) @ C                   # β の Fock 行列
    na, nb = mol.nelec                            # na = nb + 2
    nmo = C.shape[1]
    o, v = slice(0, na), slice(nb, nmo)           # α の占有軌道、β の空軌道
    no, nv = na, nmo - nb
    eri = ao2mo.general(mol, (C[:, o], C[:, o], C[:, v], C[:, v]), compact=False)
    eri = eri.reshape(no, no, nv, nv)             # (ij|ab)
    A = (np.einsum("ij,ab->iajb", np.eye(no), fb[v, v])
         - np.einsum("ij,ab->iajb", fa[o, o], np.eye(nv))
         - cx * eri.transpose(0, 2, 1, 3))
    w, X = np.linalg.eigh(A.reshape(no * nv, no * nv))
    s2 = [spin_square(X[:, k].reshape(no, nv), nb, na) for k in range(nstates)]
    return ref.e_tot + w[:nstates], np.array(s2), ref.e_tot

def spin_square(c, nb, na):
    """c[i, a]：α の占有軌道 i から β の空軌道 a への反転の係数。
    閉殻の軌道は 0..nb-1、開殻の軌道は nb..na-1（a の添字は nb から始まる）。"""
    nopen = na - nb
    s_open = sum(c[nb + k, k] for k in range(nopen))    # 開殻の軌道で、同じ軌道の中での反転
    s2 = s_open**2
    s2 += np.sum(c[:, nopen:]**2)                        # 空軌道への反転（α の励起の成分）
    s2 += np.sum(c[:nb, :]**2)                           # 閉殻の軌道からの反転（β の励起の成分）
    return s2
```

## 参考文献

- A. I. Krylov, "Size-consistent wave functions for bond-breaking: the equation-of-motion spin-flip model", Chem. Phys. Lett. 338, 375 (2001). DOI: [10.1016/S0009-2614(01)00287-1](https://doi.org/10.1016/S0009-2614(01)00287-1)
- A. I. Krylov, "Spin-flip configuration interaction: an electronic structure model that is both variational and size-consistent", Chem. Phys. Lett. 350, 522 (2001). DOI: [10.1016/S0009-2614(01)01316-1](https://doi.org/10.1016/S0009-2614(01)01316-1)
- Y. Shao, M. Head-Gordon, A. I. Krylov, "The spin–flip approach within time-dependent density functional theory: Theory and applications to diradicals", J. Chem. Phys. 118, 4807 (2003). DOI: [10.1063/1.1545679](https://doi.org/10.1063/1.1545679)
- F. Wang, T. Ziegler, "Time-dependent density functional theory based on a noncollinear formulation of the exchange-correlation potential", J. Chem. Phys. 121, 12191 (2004). DOI: [10.1063/1.1821494](https://doi.org/10.1063/1.1821494)
- J. S. Sears, C. D. Sherrill, A. I. Krylov, "A spin-complete version of the spin-flip approach to bond breaking: What is the impact of obtaining spin eigenfunctions?", J. Chem. Phys. 118, 9084 (2003). DOI: [10.1063/1.1568735](https://doi.org/10.1063/1.1568735)
- D. Casanova, A. I. Krylov, "Spin-flip methods in quantum chemistry", Phys. Chem. Chem. Phys. 22, 4326 (2020). DOI: [10.1039/C9CP06507E](https://doi.org/10.1039/C9CP06507E)
- C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J.-P. Malrieu, "Introduction of n-electron valence states for multireference perturbation theory", J. Chem. Phys. 114, 10252 (2001). DOI: [10.1063/1.1361246](https://doi.org/10.1063/1.1361246)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
