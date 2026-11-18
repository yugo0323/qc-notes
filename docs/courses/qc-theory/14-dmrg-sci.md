# 第 14 回　大きな活性空間：DMRG と選択的 CI

!!! abstract "この回の道筋"
    CASSCF で扱える活性空間は、20 軌道程度が限界です。遷移金属を何個も含む錯体などでは、それでは足りません。この回では、大きな活性空間の波動関数を近似的に、しかし系統的に精度を上げられる形で求める二つの方法、密度行列繰り込み群（DMRG）と選択的 CI を学びます。どちらも、小さな計算例でしくみを確かめます。最後に、その上で動的相関を足す方法と、2023 年以降の研究から見えてきた「何ができて、何がまだ難しいか」をまとめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 活性空間が大きいと CAS が破綻する |
    | 2 | 行列積状態と結合次元 |
    | 3 | DMRG：行列積状態を最適化する |
    | 4 | 選択的 CI：大事な決定子だけを選ぶ |
    | 5 | 大きな活性空間の上での動的相関 |
    | 6 | 何ができて、何がまだ難しいか |

    **前に読んでおく回**：[量子化学の理論 第 12 回](12-casscf.md)、[第 13 回](13-caspt2-nevpt2.md)、[第 8 回](08-configuration-interaction.md)

!!! tip "物理的なイメージ"
    DMRG と選択的 CI が利用しているのは、「大きな活性空間でも、本当に大事な電子の入り方はごく一部だ」という性質です。選択的 CI は、エネルギーへの寄与が大きい決定子だけを選んで集めます。DMRG は、軌道を一列に並べ、列のどこで二つに切っても、左右の部分の間の結びつき（もつれ）のうち大きな成分だけを残します。どちらも、残す量（決定子の数、結合次元）を増やせば、系統的に FCI に近づきます。

エネルギーは原子単位（ハートリー，$E_\mathrm{h}$）で書きます。1 m$E_\mathrm{h}$ = $10^{-3}\ E_\mathrm{h}$ で、約 0.63 kcal/mol です。

## 1. 活性空間が大きいと CAS が破綻する

第 12 回の 3 節で見たように、CAS($n$,$n$) の決定子の数は、活性軌道が 4 個増えるごとに 100 倍以上に増えます。CAS(20,20) で約 $3 \times 10^{10}$ 個、CAS(24,24) で約 $7 \times 10^{12}$ 個です。超並列の計算機でも、CAS(24,24) では CI の反復 1 回がやっとでした [Vogiatzis 2017]。

ところが、化学の問題には、もっと大きな活性空間を必要とするものがあります。たとえば、Cr₂ のような遷移金属の二量体では、3d と 4s の軌道（12 軌道）だけの活性空間では不十分で、3d 軌道をもう一組加えた 28 軌道の活性空間で、ようやく解離エネルギーの過大評価の原因が説明できたと報告されています [Kurashige 2011]。複数の金属を含む錯体や、長い共役系も同じです。

そこで、活性空間の中の FCI を厳密に解くのをあきらめ、**少ない変数で FCI の波動関数を近似し、変数を増やせば FCI に近づく**方法を使います。代表的なものが、2 節と 3 節の DMRG と、4 節の選択的 CI です。どちらも、CASSCF の「活性空間の中の FCI」を置き換える部品として使えます。

**ここまでで分かったこと**：活性空間が 20 軌道を大きく超えると、CAS の FCI は解けません。そこで、系統的に精度を上げられる近似で、活性空間の FCI を置き換えます。

## 2. 行列積状態と結合次元

### 波動関数をテンソルとして見る

活性軌道が $k$ 個あるとき、決定子は「各軌道に電子がどう入っているか」で決まります。1 個の空間軌道の状態は、空、α 電子 1 個、β 電子 1 個、α と β の 2 個、の 4 通りです。これを $n_p \in \{0, \alpha, \beta, 2\}$ と書くと、FCI の波動関数は

$$
\Psi = \sum_{n_1, n_2, \dots, n_k} C_{n_1 n_2 \cdots n_k} \, |n_1 n_2 \cdots n_k\rangle
$$

と書けます。係数 $C$ は、添字を $k$ 個もつ配列（テンソル）で、成分は全部で $4^k$ 個あります。この数が指数関数的に増えることが、1 節の問題そのものです。

### 行列の積で書く

この係数を、行列の積で表すことを考えます。

$$
C_{n_1 n_2 \cdots n_k} = \mathbf A^{[1] n_1} \mathbf A^{[2] n_2} \cdots \mathbf A^{[k] n_k}
$$

ここで $\mathbf A^{[p] n_p}$ は、軌道 $p$ の占有 $n_p$ ごとに決まる行列です。最初の $\mathbf A^{[1] n_1}$ は行ベクトル、最後の $\mathbf A^{[k] n_k}$ は列ベクトルなので、積は 1 個の数になります。この形の波動関数を**行列積状態**（matrix product state, MPS）と呼びます。

行列の大きさ（の最大値）$D$ を**結合次元**（bond dimension）と呼びます。変数の数は、およそ $4 k D^2$ 個です。$D$ を固定すれば、軌道の数 $k$ に比例して増えるだけで、指数関数的には増えません。

$D$ を十分に大きくとれば、どんな $C$ も厳密に表せます。軌道を「1 番目から $p$ 番目まで」と「残り」の二つに分けたとき、必要な $D$ は、その分け目での二つのグループの**もつれ**（entanglement）の大きさで決まります。もつれが小さい、つまり二つのグループがほぼ独立に振る舞うなら、小さな $D$ で十分です。もつれは、係数 $C$ を分け目で 2 つの添字の組に分けた行列とみなし、特異値分解したときの特異値の広がりで測れます。導出は省きますが、大きな特異値から $D$ 個を残すのが、その分け目での最良の近似です。

### 計算例：H₁₀ の FCI の波動関数を行列積状態で近似する

しくみを確かめるために、H 原子 10 個を 1.5 Å 間隔で直線に並べた H₁₀（STO-3G 基底、軌道 10 個、電子 10 個）を使います。原子の間隔が広いので、静的相関の強い系です。まず FCI の係数 $C$ を求め、左の軌道から順に特異値分解して、結合次元を $D$ までに切ります。切った $C$ でエネルギーの期待値を計算し、FCI と比べます。

軌道は 2 通り試します。RHF の正準軌道（分子全体に広がった軌道）と、各 H 原子の 1s 軌道を直交化した局在化した軌道（Löwdin の直交化）です。局在化した軌道は、鎖に沿った順に並べます。FCI のエネルギーは軌道の選び方によりませんが、もつれの大きさは変わります。

```python
import numpy as np
from pyscf import gto, scf, ao2mo, fci, lo
from pyscf.fci import cistring, direct_spin1

n = 10
mol = gto.M(atom="; ".join(f"H 0 0 {1.5*k}" for k in range(n)), basis="sto-3g", verbose=0)
mf = scf.RHF(mol).run()
na = nb = n // 2
stra = cistring.make_strings(range(n), na)   # α 電子の入り方（ビット列）
strb = cistring.make_strings(range(n), nb)

def occ_and_sign(sa, sb):
    """軌道ごとの占有（0:空, 1:α, 2:β, 3:αβ）と、並べ替えの符号。"""
    occ = tuple(((sa >> p) & 1) + 2 * ((sb >> p) & 1) for p in range(n))
    # 「α 全部 → β 全部」の順から「1α 1β 2α 2β …」の順へ並べ替えるときの符号
    par = sum(((sb >> q) & 1) * bin(sa >> (q + 1)).count("1") for q in range(n))
    return occ, (-1) ** par

def to_tensor(ci):
    T = np.zeros((4,) * n)
    for ia, sa in enumerate(stra):
        for ib, sb in enumerate(strb):
            occ, sgn = occ_and_sign(sa, sb)
            T[occ] = sgn * ci[ia, ib]
    return T

def from_tensor(T):
    ci = np.zeros((len(stra), len(strb)))
    for ia, sa in enumerate(stra):
        for ib, sb in enumerate(strb):
            occ, sgn = occ_and_sign(sa, sb)
            ci[ia, ib] = sgn * T[occ]
    return ci

def truncate(T, D):
    """左から順に特異値分解し、結合次元を D 以下に切ってから、テンソルに戻す。"""
    cores, rest, left = [], T.reshape(1, -1), 1
    for p in range(n - 1):
        U, s, Vt = np.linalg.svd(rest.reshape(left * 4, -1), full_matrices=False)
        k = min(D, len(s))
        cores.append(U[:, :k].reshape(left, 4, k))   # 軌道 p の行列 A^[p]
        rest, left = s[:k, None] * Vt[:k], k
    cores.append(rest.reshape(left, 4, 1))
    out = cores[0]
    for c in cores[1:]:
        out = np.tensordot(out, c, axes=([-1], [0]))  # 行列の積をとる
    return out.reshape((4,) * n)

for label, C in [("RHF の正準軌道", mf.mo_coeff), ("局在化した軌道（Löwdin）", lo.orth_ao(mol, "lowdin"))]:
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.full(mol, C)
    e_fci, ci = fci.direct_spin1.FCI().kernel(h1, eri, n, (na, nb), ecore=mol.energy_nuc())
    h2 = direct_spin1.absorb_h1e(h1, eri, n, (na, nb), 0.5)
    T = to_tensor(ci)
    print(label, f"E_FCI = {e_fci:.6f}")
    for D in [2, 4, 8, 16, 32, 64, 1024]:
        c = from_tensor(truncate(T, D)).ravel()
        e = c @ direct_spin1.contract_2e(h2, c, n, (na, nb)).ravel() / (c @ c) + mol.energy_nuc()
        print(f"  D = {D:4d}:  E - E_FCI = {1000*(e - e_fci):9.3f} mEh")
```

```
RHF の正準軌道 E_FCI = -4.995447
  D =    2:  E - E_FCI =   392.110 mEh
  D =    4:  E - E_FCI =   353.824 mEh
  D =    8:  E - E_FCI =   294.314 mEh
  D =   16:  E - E_FCI =   143.738 mEh
  D =   32:  E - E_FCI =   101.728 mEh
  D =   64:  E - E_FCI =    45.748 mEh
  D = 1024:  E - E_FCI =     0.000 mEh
局在化した軌道（Löwdin） E_FCI = -4.995447
  D =    2:  E - E_FCI =   783.814 mEh
  D =    4:  E - E_FCI =    81.880 mEh
  D =    8:  E - E_FCI =     6.352 mEh
  D =   16:  E - E_FCI =     0.666 mEh
  D =   32:  E - E_FCI =     0.011 mEh
  D =   64:  E - E_FCI =     0.000 mEh
  D = 1024:  E - E_FCI =     0.000 mEh
```

10 個の軌道を 5 個ずつに分ける真ん中の分け目では、それぞれの側の状態は $4^5 = 1024$ 通りです。そのため $D = 1024$ では、どちらの軌道でも厳密に一致します。違いが出るのは、$D$ が小さいときです。

- **局在化した軌道**では、$D = 16$ で誤差が 0.7 m$E_\mathrm{h}$、$D = 32$ で 0.01 m$E_\mathrm{h}$ まで下がります。真ん中の分け目では、厳密に必要な 1024 の 3% の大きさで、ほぼ厳密な値が得られます。
- **正準軌道**では、$D = 64$ でも誤差が 46 m$E_\mathrm{h}$ 残ります。

鎖の上で局在化した軌道を順に並べると、離れた軌道どうしのもつれは小さくなります。一方、分子全体に広がった正準軌道では、どこで分けても強くもつれています。行列積状態の効率は、**軌道の選び方と並べ方**に強く依存するのです。このため、実際の DMRG の計算では、局在化した軌道を使い、もつれの大きい軌道どうしが近くに来るように並べ方を工夫します。

（$D = 2$ では、局在化した軌道の誤差のほうが正準軌道より大きくなっています。ここでは特異値分解で順に切っただけで、エネルギーを最適化していないため、$D$ が非常に小さいところでは、このような逆転が起こります。次節の DMRG は、与えられた $D$ の中でエネルギーを最小にします。）

**ここまでで分かったこと**：行列積状態は、FCI の係数を行列の積で表し、結合次元 $D$ で精度を調節します。必要な $D$ は軌道どうしのもつれで決まり、局在化した軌道を適切に並べると、小さな $D$ で FCI に近づけます。

## 3. DMRG：行列積状態を最適化する

2 節の計算例では、先に FCI の係数を求めてから近似しました。実際の大きな系では、FCI の係数は求められません。そこで、行列積状態の行列 $\mathbf A^{[p]}$ を変数として、エネルギーを直接最小にします。

行列全部を一度に最適化するのは難しいので、1 個（または隣り合う 2 個）の軌道の行列だけを動かし、残りを固定して最適化します。これを鎖の端から端へ順に行い、往復（スイープ）を繰り返します。1 か所の最適化は、小さな行列の固有値問題になります。隣り合う 2 個の軌道の行列をまとめて最適化したあと、特異値分解で $D$ 個まで切るので、2 節で見た「大きな特異値を残す」切り方が、ここで使われます。

この方法は、White が 1992 年に、もともとは物理の 1 次元の模型（ハイゼンベルク鎖など）のために、**密度行列繰り込み群**（density matrix renormalization group, DMRG）として提案しました [White 1992]。DMRG で得られる波動関数は、2 節の行列積状態の形をしています [Chan 2011]。White と Martin は、1999 年に DMRG を分子の量子化学計算（ガウス型基底の H₂O）に使いました [White 1999]。化学への応用の入門的な総説もあります [Chan 2011]。

DMRG の性質をまとめます。

- **変分法です。** 与えられた $D$ の中でエネルギーを最小にするので、エネルギーは FCI より下になりません。$D$ を大きくすると FCI に近づくので、いくつかの $D$ で計算して、$D \to \infty$ に外挿することもできます。
- **1 次元的なもつれに強いです。** 2 節で見たとおり、軌道を 1 列に並べたときに、離れた軌道どうしのもつれが小さい系に向いています。
- **活性空間を大きくできます。** Stein と Reiher は、DMRG で約 100 個の軌道までを活性空間に含められることを、活性軌道の自動選択に利用しました [Stein 2016]（第 12 回の 8 節）。
- **CASSCF の部品になります。** 活性空間の FCI を DMRG で置き換え、軌道も最適化する方法を DMRG-SCF と呼びます [Kurashige 2011]。

DMRG のプログラムとしては、たとえば block2 がオープンソースで公開されています。基底状態だけでなく、時間発展、応答、有限温度の計算や、5 節で扱う動的相関の理論も実装されています [Zhai 2023]。

**ここまでで分かったこと**：DMRG は、行列積状態の行列を、軌道 1〜2 個ずつ順に最適化する変分法です。結合次元を上げると FCI に近づき、約 100 軌道の活性空間を扱えます。

## 4. 選択的 CI：大事な決定子だけを選ぶ

### 考え方

もう一つの方法は、FCI の決定子のうち、係数の大きな決定子だけを選んで CI をする方法です。第 7 回の H₂ でも、FCI の波動関数は少数の決定子でほぼ表せました。大きな系でも、ほとんどの決定子の係数は非常に小さく、エネルギーに効く決定子はごく一部です。問題は、どの決定子が大事かを、FCI を解く前に見分けることです。

**選択的 CI**（selected CI）は、これを反復で行います。

1. 少数の決定子（最初は HF の決定子だけなど）の空間で CI の行列を対角化します。
2. いまの波動関数 $\Psi = \sum_i c_i \Phi_i$ と、ハミルトニアンの行列要素でつながる外側の決定子 $\Phi_k$ のうち、「大事さ」の指標が大きいものを空間に加えます。
3. 1 に戻り、エネルギーが変わらなくなるまで繰り返します。

最後に、選ばれなかった決定子の寄与を、第 9 回の 2 次の摂動論の形で足します。

### CIPSI と heat-bath CI

この考え方は、Huron、Malrieu、Rancurel が 1973 年に **CIPSI** として提案しました [Huron 1973]。CIPSI は、外側の決定子の大事さを、1 次の摂動論での係数の見積もり

$$
c_k^{(1)} = \frac{\sum_i H_{ki} c_i}{E_0 - H_{kk}}
$$

の大きさで測ります。

Holmes、Tubman、Umrigar は、指標をもっと簡単な

$$
\max_i \, |H_{ki} c_i| > \varepsilon_1
$$

に替えた **heat-bath CI**（HCI）を提案しました [Holmes 2016]。CIPSI の指標の分母は、決定子によってあまり変わらないので、分子の大きさでほぼ順位が決まります。HCI の指標なら、大きな行列要素 $H_{ki}$ をもつ決定子だけを最初から作れるので、外側のすべての決定子を調べる必要がありません。Holmes らは、C₂（cc-pVTZ 基底）と Cr₂（全電子を相関させた計算）の FCI のエネルギーを、1 個の CPU で数分のうちに 1 m$E_\mathrm{h}$ より良い精度で再現しました [Holmes 2016]。2 次の摂動の補正を確率的に計算して、メモリの壁を取り除いた版（SHCI）もあります [Sharma 2017]。同じころ、外側の決定子を探す範囲を絞って CIPSI を速くした ASCI も提案されています [Tubman 2016]。

選択的 CI では、しきい値 $\varepsilon_1$ を小さくすると、決定子が増えて FCI に近づきます。いくつかの $\varepsilon_1$ で計算して、$\varepsilon_1 \to 0$ に外挿することもよく行われます [Holmes 2016]。

### 計算例：伸ばした H₂O の選択的 CI

第 13 回で使った、O–H 結合を 2 倍に伸ばした H₂O（6-31G 基底）で、選択的 CI を試します。PySCF の `fci.selected_ci` は、HCI に近い基準（係数の大きさと 2 電子積分の大きさの積がしきい値を超えるもの）で、α 電子と β 電子の入り方を別々に選ぶ簡単な実装です。PySCF のソースコードの説明にも、効率のよい実装ではないと書かれています。ここでは、2 次の摂動の補正を足さない、変分のエネルギーだけを見ます。

```python
import numpy as np
from pyscf import gto, scf, ao2mo, fci
from pyscf.fci import selected_ci

r, th = 0.96 * 2.0, np.deg2rad(104.5) / 2      # O–H を 2 倍に伸ばした H2O
mol = gto.M(atom=f"O 0 0 0; H 0 {r*np.sin(th):.6f} {r*np.cos(th):.6f}; H 0 {-r*np.sin(th):.6f} {r*np.cos(th):.6f}",
            basis="6-31g", verbose=0)
mf = scf.RHF(mol).run()
C = mf.mo_coeff
norb, nelec = C.shape[1], (5, 5)
h1 = C.T @ mf.get_hcore() @ C
eri = ao2mo.full(mol, C)

e_fci = fci.FCI(mf).kernel()[0]
n_fci = fci.cistring.num_strings(norb, 5) ** 2
print(f"FCI: 決定子 {n_fci:,d} 個, E = {e_fci:.6f}")
print("しきい値   決定子の数   FCI に対する割合   E - E_FCI / mEh")
for cut in [3e-2, 1e-2, 3e-3, 1e-3, 3e-4]:
    sci = selected_ci.SCI(mol)
    sci.select_cutoff = sci.ci_coeff_cutoff = cut
    e, civec = sci.kernel(h1, eri, norb, nelec, ecore=mol.energy_nuc())
    ndet = len(civec._strs[0]) * len(civec._strs[1])
    print(f"{cut:7.0e}  {ndet:10,d}       {100*ndet/n_fci:6.2f} %       {1000*(e - e_fci):8.3f}")
```

```
FCI: 決定子 1,656,369 個, E = -75.879753
しきい値   決定子の数   FCI に対する割合   E - E_FCI / mEh
  3e-02         625         0.04 %         40.290
  1e-02       5,776         0.35 %         12.744
  3e-03      32,041         1.93 %          2.558
  1e-03      96,100         5.80 %          0.400
  3e-04     336,400        20.31 %          0.022
```

FCI の計算に 2 分ほどかかり、選択的 CI はすべて合わせて 30 秒ほどです。しきい値 $10^{-3}$ では、FCI の決定子の 6% 弱で、誤差が 0.4 m$E_\mathrm{h}$ まで下がります。

この実装は α 電子と β 電子の入り方を別々に選び、その組み合わせをすべて使うので、本当に必要な数より多くの決定子を含みます。決定子を 1 個ずつ選ぶ HCI などの実装では、必要な数はもっと少なくなります。また、この H₂O は FCI が解ける小さな例です。Holmes らの Cr₂ の計算では、FCI の決定子の数は約 $2 \times 10^{22}$ 個と見積もられており、選ばれた決定子は、そのごくわずかな割合にすぎません [Holmes 2016]。

**ここまでで分かったこと**：選択的 CI は、ハミルトニアンの行列要素と係数から大事な決定子を選び、残りを摂動で足す方法です。しきい値を下げると、FCI の決定子のごく一部で FCI のエネルギーに近づきます。

## 5. 大きな活性空間の上での動的相関

DMRG や選択的 CI で活性空間を大きくしても、第 12 回の 9 節で見た動的相関の不足は残ります。そこで、その上に第 13 回の多参照摂動論を足します。

- **DMRG-CASPT2**：Kurashige と Yanai は、DMRG-SCF の波動関数を出発点とする CASPT2 を開発し、Cr₂ のポテンシャルエネルギー曲線を、28 軌道の活性空間で計算しました [Kurashige 2011]。12 軌道の活性空間の CASPT2 が解離エネルギーを大きく見積もりすぎ、0 次のハミルトニアンの選び方にも敏感すぎる理由を、活性空間の大きさから説明しています。
- **DMRG-NEVPT2**：Guo らは、DMRG の波動関数から高次の密度行列を効率よく計算する方法を開発し、SC-NEVPT2 と組み合わせました [Guo 2016]。
- **選択的 CI と多参照摂動論**：2025 年には、選択的 CI に近い方法で得た近似的な CASSCF の波動関数の上で、CASPT2 や NEVPT2 を計算した報告が続いています。活性空間は、最大で CAS(34,34) です [Min 2025; Guo 2025]。

これらに共通する壁は、第 13 回の 2 節と 6 節で見た、**3 電子・4 電子の密度行列**です。活性軌道の数を $k$ とすると、4 電子の密度行列の成分の数は $k^8$ に比例します。$k = 30$ なら $6.6 \times 10^{11}$ 個で、そのまま保存するのは現実的ではありません。2025 年の総説は、高次の密度行列を避けながら大きな活性空間の外側の動的相関を扱う方法を、7 つの種類に分けて整理しています [Song 2025]。上の CAS(34,34) の報告も、高次の密度行列をそのまま保存せずに、必要な中間量を直接つくる工夫をしています [Min 2025; Guo 2025]。

**ここまでで分かったこと**：DMRG や選択的 CI の上に CASPT2 や NEVPT2 を足すことで、大きな活性空間と動的相関を両立できます。その主な壁は、高次の密度行列の大きさです。

## 6. 何ができて、何がまだ難しいか

最後に、確かめられた文献の範囲で、現状を整理します。分かっていることと、まだ分かっていないことを分けて書きます。

### 分かっていること・できるようになったこと

- **数十軌道の活性空間の近似的な FCI**：選択的 CI を使った近似的な CASSCF で、50 電子・50 軌道を超える活性空間が、相対エネルギーで 1 kcal/mol 以内の精度で扱えると報告されています [Levine 2020]。DMRG では約 100 軌道までが扱えます [Stein 2016]。
- **難しい遷移金属二量体の決着**：Cr₂ は、長いあいだ理論と実験が合わない問題でした。2022 年の大規模な計算は、ポテンシャルエネルギー曲線と振動スペクトルを大きく改善し、実験の振動の帰属を一つずらすことで、すべての結合距離で理論と実験が定量的に一致することを示しました [Larsson 2022]。
- **ソフトウェアの整備**：DMRG の多くの手法を実装したオープンソースのプログラムが公開され、使える人が増えています [Zhai 2023]。

### まだ難しいこと

- **「ほぼ厳密」どうしの食い違い**：ベンゼン（cc-pVDZ 基底、内殻を凍結）の FCI のエネルギーを、多くのグループがさまざまな「ほぼ厳密」な方法で、互いの結果を知らずに求める試みが行われました [Eriksen 2020]。相関エネルギーはおおむね $-863$ m$E_\mathrm{h}$ 付近で一致したものの、方法の間のばらつき（二乗平均平方根）は 1.3 m$E_\mathrm{h}$ ありました。小さな分子で実績のある方法でも、分子が大きくなると、1 m$E_\mathrm{h}$ の精度を保証するのは難しいことを示しています。
- **2 次元・3 次元的なもつれ**：2 節で見たとおり、行列積状態は軌道を 1 列に並べたときのもつれに強く、そうでない系では必要な $D$ が大きくなります。2025 年にも、ふつうの DMRG にはサイズ無矛盾性がなく、2 次元・3 次元的な系では精度を保ちにくいことを課題として、埋め込み法と組み合わせる提案が出ています [Weatherly 2025]。
- **選択的 CI の軌道最適化**：選択的 CI を使った近似的な CASSCF では、活性軌道どうしの回転に対してエネルギーが不変でなくなります（第 12 回の 2 節で見た、CAS の便利な性質が失われます）。そのため、軌道の最適化で多くの極小点に落ちやすいことが指摘されています [Levine 2020]。
- **動的相関**：5 節で見たとおり、大きな活性空間の上で動的相関を足す方法は、高次の密度行列の扱いを中心に、まだ発展の途中です [Song 2025]。
- **量子コンピュータとの関係**：量子コンピュータが基底状態のエネルギーの計算で古典計算に対して指数関数的に有利だという主張について、2023 年の研究は、古典計算機での近似的な解法の性能も含めて証拠を調べ、化学の問題全般でそのような有利さを示す証拠はまだ見つかっていない、と結論しました [Lee 2023]。

**ここまでで分かったこと**：DMRG と選択的 CI により、数十〜100 軌道程度の活性空間が扱えるようになり、Cr₂ のような難問にも決着がつきました。一方、大きな分子での m$E_\mathrm{h}$ の精度、2・3 次元的なもつれ、軌道の最適化、動的相関の追加は、まだ研究が続いています。

## まとめ

- 活性空間が 20 軌道を大きく超えると、CAS の FCI は解けないので、系統的に改善できる近似で置き換えます。
- 行列積状態は、FCI の係数を行列の積で表し、結合次元 $D$ で精度を調節します。必要な $D$ は軌道どうしのもつれで決まり、軌道の選び方と並べ方が効きます。
- DMRG は、行列積状態を軌道 1〜2 個ずつ順に最適化する変分法で、約 100 軌道の活性空間を扱えます。
- 選択的 CI（CIPSI、HCI など）は、大事な決定子だけを選んで対角化し、残りを摂動で足します。
- これらの上に CASPT2 や NEVPT2 を足せますが、高次の密度行列が壁です。大きな分子での高精度、2・3 次元的な系、軌道の最適化などは、まだ研究の途中です。

## 参考文献

- K. D. Vogiatzis, D. Ma, J. Olsen, L. Gagliardi, W. A. de Jong, "Pushing configuration-interaction to the limit: Towards massively parallel MCSCF calculations", J. Chem. Phys. 147, 184111 (2017). DOI: [10.1063/1.4989858](https://doi.org/10.1063/1.4989858)
- Y. Kurashige, T. Yanai, "Second-order perturbation theory with a density matrix renormalization group self-consistent field reference function: Theory and application to the study of chromium dimer", J. Chem. Phys. 135, 094104 (2011). DOI: [10.1063/1.3629454](https://doi.org/10.1063/1.3629454)
- S. R. White, "Density matrix formulation for quantum renormalization groups", Phys. Rev. Lett. 69, 2863 (1992). DOI: [10.1103/PhysRevLett.69.2863](https://doi.org/10.1103/PhysRevLett.69.2863)
- S. R. White, R. L. Martin, "Ab initio quantum chemistry using the density matrix renormalization group", J. Chem. Phys. 110, 4127 (1999). DOI: [10.1063/1.478295](https://doi.org/10.1063/1.478295)
- G. K.-L. Chan, S. Sharma, "The Density Matrix Renormalization Group in Quantum Chemistry", Annu. Rev. Phys. Chem. 62, 465 (2011). DOI: [10.1146/annurev-physchem-032210-103338](https://doi.org/10.1146/annurev-physchem-032210-103338)
- C. J. Stein, M. Reiher, "Automated Selection of Active Orbital Spaces", J. Chem. Theory Comput. 12, 1760 (2016). DOI: [10.1021/acs.jctc.6b00156](https://doi.org/10.1021/acs.jctc.6b00156)
- H. Zhai et al., "Block2: A comprehensive open source framework to develop and apply state-of-the-art DMRG algorithms in electronic structure and beyond", J. Chem. Phys. 159, 234801 (2023). DOI: [10.1063/5.0180424](https://doi.org/10.1063/5.0180424)
- B. Huron, J. P. Malrieu, P. Rancurel, "Iterative perturbation calculations of ground and excited state energies from multiconfigurational zeroth-order wavefunctions", J. Chem. Phys. 58, 5745 (1973). DOI: [10.1063/1.1679199](https://doi.org/10.1063/1.1679199)
- A. A. Holmes, N. M. Tubman, C. J. Umrigar, "Heat-Bath Configuration Interaction: An Efficient Selected Configuration Interaction Algorithm Inspired by Heat-Bath Sampling", J. Chem. Theory Comput. 12, 3674 (2016). DOI: [10.1021/acs.jctc.6b00407](https://doi.org/10.1021/acs.jctc.6b00407)
- S. Sharma, A. A. Holmes, G. Jeanmairet, A. Alavi, C. J. Umrigar, "Semistochastic Heat-Bath Configuration Interaction Method: Selected Configuration Interaction with Semistochastic Perturbation Theory", J. Chem. Theory Comput. 13, 1595 (2017). DOI: [10.1021/acs.jctc.6b01028](https://doi.org/10.1021/acs.jctc.6b01028)
- N. M. Tubman, J. Lee, T. Y. Takeshita, M. Head-Gordon, K. B. Whaley, "A deterministic alternative to the full configuration interaction quantum Monte Carlo method", J. Chem. Phys. 145, 044112 (2016). DOI: [10.1063/1.4955109](https://doi.org/10.1063/1.4955109)
- S. Guo, M. A. Watson, W. Hu, Q. Sun, G. K.-L. Chan, "N-Electron Valence State Perturbation Theory Based on a Density Matrix Renormalization Group Reference Function, with Applications to the Chromium Dimer and a Trimer Model of Poly(p-Phenylenevinylene)", J. Chem. Theory Comput. 12, 1583 (2016). DOI: [10.1021/acs.jctc.5b01225](https://doi.org/10.1021/acs.jctc.5b01225)
- K. S. Min, J. W. Park, "Second-Order Complete Active Space Perturbation Theory (CASPT2) and N-Electron Valence State Perturbation Theory (NEVPT2) Based on Adaptive Sampling Configuration Interaction Self-Consistent Field (ASCI-SCF)", J. Chem. Theory Comput. 21, 5425 (2025). DOI: [10.1021/acs.jctc.5c00190](https://doi.org/10.1021/acs.jctc.5c00190)
- Y. Guo, K. Sivalingam, V. G. Chilkuri, F. Neese, "Approximations of density matrices in N-electron valence state second-order perturbation theory (NEVPT2). III. Large active space calculations with selected configuration interaction reference", J. Chem. Phys. 162, 144110 (2025). DOI: [10.1063/5.0262473](https://doi.org/10.1063/5.0262473)
- Y. Song, Y. Cheng, H. Ma, "Describing dynamic electron correlation beyond a large active space", Phys. Chem. Chem. Phys. 27, 12622 (2025). DOI: [10.1039/d5cp00998g](https://doi.org/10.1039/d5cp00998g)
- D. S. Levine, D. Hait, N. M. Tubman, S. Lehtola, K. B. Whaley, M. Head-Gordon, "CASSCF with Extremely Large Active Spaces Using the Adaptive Sampling Configuration Interaction Method", J. Chem. Theory Comput. 16, 2340 (2020). DOI: [10.1021/acs.jctc.9b01255](https://doi.org/10.1021/acs.jctc.9b01255)
- H. R. Larsson, H. Zhai, C. J. Umrigar, G. K.-L. Chan, "The Chromium Dimer: Closing a Chapter of Quantum Chemistry", J. Am. Chem. Soc. 144, 15932 (2022). DOI: [10.1021/jacs.2c06357](https://doi.org/10.1021/jacs.2c06357)
- J. J. Eriksen et al., "The Ground State Electronic Energy of Benzene", J. Phys. Chem. Lett. 11, 8922 (2020). DOI: [10.1021/acs.jpclett.0c02621](https://doi.org/10.1021/acs.jpclett.0c02621)
- S. Weatherly, T. Van Voorhis, "The Embedded Density Matrix Renormalization Group: Size-Extensive and Quasi-Exact for Nonlinear Quantum Chemistry", J. Chem. Theory Comput. 21, 7818 (2025). DOI: [10.1021/acs.jctc.5c00725](https://doi.org/10.1021/acs.jctc.5c00725)
- S. Lee et al., "Evaluating the evidence for exponential quantum advantage in ground-state quantum chemistry", Nat. Commun. 14, 1952 (2023). DOI: [10.1038/s41467-023-37587-6](https://doi.org/10.1038/s41467-023-37587-6)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## 関連リンク

- block2（DMRG のプログラム）：https://github.com/block-hczhai/block2-preview
