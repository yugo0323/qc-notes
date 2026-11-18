# 第 4 回　第二量子化：生成・消滅演算子

!!! abstract "この回の道筋"
    決定子をたくさん使う方法（電子相関の方法）では、行列式を毎回書き下すのは大変です。第二量子化では、決定子を「どのスピン軌道に電子がいるか」の並びで表し、電子を出し入れする演算子で計算します。この回では、生成・消滅演算子とその反交換関係を導入し、ハミルトニアンを書き直します。最後に numpy で H₂ のハミルトニアン行列を組み立て、Slater–Condon 則と一致することを確かめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 占有数表示 |
    | 2 | 生成・消滅演算子 |
    | 3 | 反交換関係 |
    | 4 | ハミルトニアンを演算子で書く |
    | 5 | Slater–Condon 則との対応 |
    | 6 | numpy で確かめる：H₂ の 2 軌道 2 電子 |

    **前に読んでおく回**：量子化学の理論 [第 2 回](02-determinants.md)、[第 3 回](03-hartree-fock.md)

!!! tip "物理的なイメージ"
    第二量子化は、電子の状態を「どの席（スピン軌道）が埋まっているか」の名簿として見る見方です。生成演算子は空いた席に電子を 1 個座らせ、消滅演算子は座っている電子を 1 個立たせます。電子は区別できないので、「どの電子がどの席にいるか」ではなく、「どの席が埋まっているか」だけが意味をもちます。演算子の反交換関係が、電子を入れ替えると符号が変わること（反対称性）と、同じ席に 2 個座れないこと（パウリの原理）を自動的に守ります。ハミルトニアンは、電子を 1 個または 2 個、席から席へ動かす操作の和として書けます。

## 1. 占有数表示

規格直交なスピン軌道を $M$ 個用意し、番号 $\chi_1, \chi_2, \dots, \chi_M$ を付けます。
たとえば第 3 回の Hartree–Fock 法で得た軌道（占有軌道と仮想軌道）がこれに当たります。
この中から $N$ 個を選んで電子を入れると、決定子が 1 つ決まります。

第 2 回で見たように、決定子は「どのスピン軌道を占有したか」で決まり、軌道の並べ方は符号を変えるだけです。
そこで、軌道の番号の小さい順に並べると決めておけば、決定子は各スピン軌道の占有数 $n_p$（0 か 1）の並びで表せます。

$$
|\chi_1\chi_3\rangle \;\longleftrightarrow\; |n_1 n_2 n_3 n_4\rangle = |1\,0\,1\,0\rangle
\qquad(M = 4,\ N = 2\ \text{の例})
$$

これを占有数表示（occupation-number representation）と呼びます。
電子が 1 個もない状態 $|0\,0\cdots 0\rangle$ を真空（vacuum）と呼び、$|\mathrm{vac}\rangle$ と書きます。

占有数の並びが違う決定子どうしは直交し、同じなら内積は 1 です。

$$
\langle \mathbf n|\mathbf m\rangle = \prod_{p=1}^{M}\delta_{n_p m_p}
$$

**ここまでで分かったこと**：軌道の順番を決めておけば、決定子は占有数の並び $|n_1 n_2\cdots n_M\rangle$ で表せます。

## 2. 生成・消滅演算子

### 2.1 定義

スピン軌道 $\chi_p$ に電子を 1 個加える演算子 $\hat a_p^\dagger$ を、生成演算子（creation operator）と呼びます。
次のように定義します。

$$
\hat a_p^\dagger\,|\cdots n_p\cdots\rangle =
\begin{cases}
\Gamma_p\,|\cdots 1_p\cdots\rangle & (n_p = 0)\\
0 & (n_p = 1)
\end{cases}
\qquad
\Gamma_p = (-1)^{\sum_{q<p} n_q}
$$

すでに電子がいる軌道に加えようとすると 0 になります。これがパウリの排他原理です。
符号 $\Gamma_p$ は、$p$ より番号の小さい占有軌道の数が奇数なら $-1$ になります。
この符号の意味は次の 2.2 で説明します。

同じように、$\chi_p$ から電子を 1 個取り除く演算子 $\hat a_p$ を、消滅演算子（annihilation operator）と呼びます。

$$
\hat a_p\,|\cdots n_p\cdots\rangle =
\begin{cases}
\Gamma_p\,|\cdots 0_p\cdots\rangle & (n_p = 1)\\
0 & (n_p = 0)
\end{cases}
$$

$\hat a_p$ は $\hat a_p^\dagger$ のエルミート共役です。

### 2.2 符号の意味

符号 $\Gamma_p$ は、決定子の反対称性を占有数表示に持ち込むためのものです。
占有数の並び $|\mathbf n\rangle$ を、生成演算子を番号の小さい順に左から並べて真空にかけたもの

$$
|\mathbf n\rangle = \big(\hat a_1^\dagger\big)^{n_1}\big(\hat a_2^\dagger\big)^{n_2}\cdots\big(\hat a_M^\dagger\big)^{n_M}|\mathrm{vac}\rangle
$$

と約束します。
たとえば $|1\,0\,1\,0\rangle = \hat a_1^\dagger\hat a_3^\dagger|\mathrm{vac}\rangle$ です。
ここに $\hat a_2^\dagger$ を左からかけると、$\hat a_2^\dagger\hat a_1^\dagger\hat a_3^\dagger|\mathrm{vac}\rangle$ になります。
約束の順番に直すには $\hat a_2^\dagger$ を $\hat a_1^\dagger$ の右へ 1 回入れ替える必要があり、その分だけ符号が変わります（次の節の反交換関係）。
したがって $\hat a_2^\dagger|1\,0\,1\,0\rangle = -|1\,1\,1\,0\rangle$ で、これが $\Gamma_2 = (-1)^{n_1} = -1$ です。

**ここまでで分かったこと**：生成演算子 $\hat a_p^\dagger$ は電子を加え、消滅演算子 $\hat a_p$ は電子を取り除きます。決定子は、生成演算子を順に真空にかけたものです。

## 3. 反交換関係

生成・消滅演算子の性質は、次の 3 つの反交換関係にまとめられます。

$$
\{\hat a_p, \hat a_q^\dagger\} = \delta_{pq},\qquad
\{\hat a_p, \hat a_q\} = 0,\qquad
\{\hat a_p^\dagger, \hat a_q^\dagger\} = 0
$$

ここで $\{\hat A, \hat B\} = \hat A\hat B + \hat B\hat A$ は反交換子（anticommutator）です。
定義から直接確かめられますが、この連載では第 6 節で数値的に確かめます。

この関係から、次のことが分かります。

- $\hat a_p^\dagger\hat a_q^\dagger = -\hat a_q^\dagger\hat a_p^\dagger$：電子を加える順番を入れ替えると符号が変わります。決定子の反対称性です。
- $p = q$ とおくと $\hat a_p^\dagger\hat a_p^\dagger = 0$：同じ軌道に 2 個の電子は加えられません。パウリの排他原理です。
- $\hat n_p = \hat a_p^\dagger\hat a_p$ は、$|\mathbf n\rangle$ にかけると $n_p$ を返します。これを数演算子（number operator）と呼びます。$\hat N = \sum_p \hat n_p$ は電子の数を返します。

フェルミ粒子の演算子に反交換関係を使うことは Jordan と Wigner が 1928 年に導入し [Jordan 1928]、Fock が 1932 年に多粒子の波動関数との関係を整理しました [Fock 1932]。
教科書では、Szabo–Ostlund の 2.4 節と、Helgaker らの第 1 章で詳しく扱われています [Szabo 1996; Helgaker 2000]。

第一量子化（行列式を書き下すやり方）では、反対称性は波動関数の側にありました。
第二量子化では、反対称性は演算子の反交換関係の側に移ります。
このため、波動関数は占有数の並びという単純なものになります。

**ここまでで分かったこと**：決定子の反対称性とパウリの排他原理は、生成・消滅演算子の反交換関係にすべて含まれています。

## 4. ハミルトニアンを演算子で書く

### 4.1 1 電子演算子

第 3 回で使った $\hat h$ のように、各電子に別々に働く演算子の和 $\sum_{i=1}^{N}\hat h(i)$ を、1 電子演算子と呼びます。
1 電子演算子は、第二量子化では次のように書けます。

$$
\sum_{i=1}^{N}\hat h(i) \;\longrightarrow\; \sum_{p,q=1}^{M} h_{pq}\,\hat a_p^\dagger\hat a_q,
\qquad h_{pq} = \langle\chi_p|\hat h|\chi_q\rangle
$$

$\hat a_p^\dagger\hat a_q$ は、電子を $\chi_q$ から $\chi_p$ へ移す演算子です。
その「重み」が積分 $h_{pq}$ です。
決定子の間の行列要素が、元の形と同じになるように決めたものです。

### 4.2 2 電子演算子

電子反発 $\sum_{i<j}1/r_{ij}$ のように、電子の組に働く演算子を 2 電子演算子と呼びます。
2 電子演算子は、2 個の電子を同時に移す演算子で書けます。
スピン軌道についての化学者の表記の積分

$$
[pq|rs] = \iint \chi_p^*(\mathbf x_1)\chi_q(\mathbf x_1)\,\frac{1}{r_{12}}\,\chi_r^*(\mathbf x_2)\chi_s(\mathbf x_2)\,d\mathbf x_1\,d\mathbf x_2
$$

を使うと、

$$
\sum_{i<j}\frac{1}{r_{ij}} \;\longrightarrow\; \frac12\sum_{p,q,r,s=1}^{M}[pq|rs]\;\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q
$$

です。
添字の並びに注意してください。
電子 1 が $\chi_q$ から $\chi_p$ へ、電子 2 が $\chi_s$ から $\chi_r$ へ移ります。
消滅演算子の順番 $\hat a_s\hat a_q$ は、生成演算子の順番 $\hat a_p^\dagger\hat a_r^\dagger$ と逆になります。
係数 $\frac12$ は、電子の組を 2 回数えないためのものです。

### 4.3 電子のハミルトニアン

まとめると、第 1 回の電子のハミルトニアンは次のように書けます。

$$
\hat H_\mathrm{el} = \sum_{pq} h_{pq}\,\hat a_p^\dagger\hat a_q + \frac12\sum_{pqrs}[pq|rs]\,\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q
$$

この形には、2 つの特徴があります。

1. 電子の数 $N$ が式に出てきません。分子の情報はすべて積分 $h_{pq}$ と $[pq|rs]$ に入っています。
2. 演算子が軌道の基底の上だけで定義されています。$M$ 個のスピン軌道で作れる決定子の範囲での、ハミルトニアンの「射影」です。$M$ を無限に大きくすれば、元のハミルトニアンに一致します。

スピン軌道を空間軌道 $\phi_P$ とスピン $\sigma$ に分けると、積分はスピンが合うときだけ 0 でなくなります。
そのため、次の演算子でまとめて書くことがよくあります。

$$
\hat E_{PQ} = \sum_{\sigma=\alpha,\beta}\hat a_{P\sigma}^\dagger\hat a_{Q\sigma}
$$

$\hat E_{PQ}$ は、スピンを変えずに電子を空間軌道 $\phi_Q$ から $\phi_P$ へ移す演算子です。
これを使うと $\hat H_\mathrm{el} = \sum_{PQ}h_{PQ}\hat E_{PQ} + \frac12\sum_{PQRS}(PQ|RS)(\hat E_{PQ}\hat E_{RS} - \delta_{QR}\hat E_{PS})$ と書けます [Helgaker 2000]。
$(PQ|RS)$ は空間軌道についての化学者の表記の積分です。
この形は、第 5 回でスピンを扱うときに役立ちます。

**ここまでで分かったこと**：1 電子演算子は $\sum h_{pq}\hat a_p^\dagger\hat a_q$、2 電子演算子は $\frac12\sum[pq|rs]\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q$ と書けます。分子の情報は積分の値に入っています。

## 5. Slater–Condon 則との対応

決定子 $|\Phi\rangle$、$|\Phi'\rangle$ の間のハミルトニアンの行列要素 $\langle\Phi'|\hat H_\mathrm{el}|\Phi\rangle$ を考えます。
これを行列式から直接求める規則を、Slater–Condon 則と呼びます [Slater 1929; Condon 1930]。
第二量子化では、この規則は反交換関係から機械的に出てきます。

$\hat H_\mathrm{el}$ には、消滅演算子が多くても 2 個しかありません。
したがって、$|\Phi\rangle$ に $\hat H_\mathrm{el}$ をかけても、電子は多くても 2 個しか移りません。
このことから、2 つの決定子の違いの大きさで場合分けができます。
以下、$i, j$ は $|\Phi\rangle$ で占有されたスピン軌道、$a, b$ は占有されていないスピン軌道です。
$|\Phi_i^a\rangle$ は $\chi_i$ の電子を $\chi_a$ へ移した決定子、$|\Phi_{ij}^{ab}\rangle$ は 2 個の電子を移した決定子です。

| 2 つの決定子の違い | 行列要素 |
|---|---|
| 同じ | $\langle\Phi\vert\hat H_\mathrm{el}\vert\Phi\rangle = \sum_i h_{ii} + \frac12\sum_{ij}\big([ii\vert jj] - [ij\vert ji]\big)$ |
| スピン軌道 1 個 | $\langle\Phi_i^a\vert\hat H_\mathrm{el}\vert\Phi\rangle = h_{ai} + \sum_j\big([ai\vert jj] - [aj\vert ji]\big)$ |
| スピン軌道 2 個 | $\langle\Phi_{ij}^{ab}\vert\hat H_\mathrm{el}\vert\Phi\rangle = [ai\vert bj] - [aj\vert bi]$ |
| スピン軌道 3 個以上 | 0 |

和の $i, j$ は、$|\Phi\rangle$ の占有スピン軌道すべてについてとります。
1 行目は、第 3 回の付録 B で使った決定子のエネルギーそのものです。
2 行目は、第 3 回の Fock 演算子の行列要素 $\langle\chi_a|\hat f|\chi_i\rangle$ に等しくなります。
Hartree–Fock の軌道では、この値が 0 になります。これを Brillouin の定理と呼びます。

1 行目の導出は付録 A に回しますが、手順は次のとおりです。
$\hat a_p^\dagger\hat a_q$ などを $|\Phi\rangle$ にかけ、反交換関係で消滅演算子を右へ動かします。
真空に消滅演算子がかかると 0 になるので、残る項は数えるだけで求まります。
2 行目以降の符号は、移した電子の軌道の位置関係（$\Gamma_p$）で決まります。
上の表は、$|\Phi_i^a\rangle$ を「$|\Phi\rangle$ の中の $\chi_i$ を、同じ位置で $\chi_a$ に置き換えた決定子」と定義したときの符号です。
占有数表示（番号の小さい順）に並べ直すと、符号が変わることがあります。
次の節で、この符号の扱いを数値で確かめます。

**ここまでで分かったこと**：Slater–Condon 則は、反交換関係から出てきます。ハミルトニアンは、電子が 3 個以上違う決定子どうしを直接つなぎません。

## 6. numpy で確かめる：H₂ の 2 軌道 2 電子

ここまでの話を、いちばん小さな分子で確かめます。
H₂（核間距離 0.74 Å）を STO-3G 基底 [Hehre 1969] で解くと、空間軌道は結合性の $\phi_1$ と反結合性の $\phi_2$ の 2 個です。
スピン軌道は 4 個で、番号を $0 = \phi_1\alpha$、$1 = \phi_1\beta$、$2 = \phi_2\alpha$、$3 = \phi_2\beta$ とします。

占有数の並びは $2^4 = 16$ 通りなので、演算子は $16\times 16$ の行列で表せます。
消滅演算子 $\hat a_p$ は、各スピン軌道の占有数を表す 2 次元の空間のテンソル積として作ります。
$p$ より前の軌道には $(-1)^{n}$ を返す行列を置き、これで符号 $\Gamma_p$ を表します。
（この作り方は Jordan–Wigner 変換と呼ばれます [Jordan 1928]。）
積分は PySCF [Sun 2020] で計算します。

```python
import numpy as np
from functools import reduce
from pyscf import gto, scf, ao2mo, fci

# --- 1. H2/STO-3G の分子軌道の積分（空間軌道 2 個）---
mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)
mf = scf.RHF(mol).run()
C = mf.mo_coeff
h_mo = C.T @ mf.get_hcore() @ C                       # h_PQ
eri_mo = ao2mo.restore(1, ao2mo.kernel(mol, C), 2)   # (PQ|RS) 化学者の表記

# --- 2. スピン軌道の積分 0=1α, 1=1β, 2=2α, 3=2β ---
M = 4
sp = lambda p: p // 2          # 空間軌道の番号
sg = lambda p: p % 2           # スピン（0=α, 1=β）
h = np.zeros((M, M)); g = np.zeros((M, M, M, M))
for p in range(M):
    for q in range(M):
        if sg(p) == sg(q):
            h[p, q] = h_mo[sp(p), sp(q)]
        for r in range(M):
            for s in range(M):
                if sg(p) == sg(q) and sg(r) == sg(s):
                    g[p, q, r, s] = eri_mo[sp(p), sp(q), sp(r), sp(s)]

# --- 3. 消滅演算子 a_p を 16×16 行列で作る ---
I2 = np.eye(2)
Z = np.diag([1.0, -1.0])                  # (-1)^n
low = np.array([[0.0, 1.0], [0.0, 0.0]])  # |0><1|：電子を 1 個消す
def annihilate(p):
    ops = [Z] * p + [low] + [I2] * (M - p - 1)
    return reduce(np.kron, ops)
a = [annihilate(p) for p in range(M)]
ad = [x.T for x in a]                     # 生成演算子 a_p^†

# 反交換関係を確かめる
ok = all(np.allclose(a[p] @ ad[q] + ad[q] @ a[p], np.eye(16) * (p == q))
         and np.allclose(a[p] @ a[q] + a[q] @ a[p], 0)
         for p in range(M) for q in range(M))
print("反交換関係が成り立つ:", ok)

# --- 4. 第二量子化のハミルトニアン ---
H = sum(h[p, q] * ad[p] @ a[q] for p in range(M) for q in range(M))
H = H + 0.5 * sum(g[p, q, r, s] * ad[p] @ ad[r] @ a[s] @ a[q]
                  for p in range(M) for q in range(M)
                  for r in range(M) for s in range(M))

# --- 5. 2 電子・Sz = 0 の決定子だけを取り出す ---
def index(occ):                         # 占有数 (n0,n1,n2,n3) → 行列の番号
    return int("".join(map(str, occ)), 2)
dets = {"|1α 1β>": (1, 1, 0, 0), "|1α 2β>": (1, 0, 0, 1),
        "|1β 2α>": (0, 1, 1, 0), "|2α 2β>": (0, 0, 1, 1)}
idx = [index(o) for o in dets.values()]
Hsub = H[np.ix_(idx, idx)]
np.set_printoptions(precision=5, suppress=True)
print("決定子:", list(dets.keys()))
print("ハミルトニアン行列（電子部分、Eh）:\n", Hsub)

# --- 6. 対角化して PySCF の FCI と比べる ---
w = np.linalg.eigvalsh(Hsub)
print("最低固有値 + V_NN :", w[0] + mol.energy_nuc())
print("PySCF の FCI       :", fci.FCI(mf).kernel()[0])
print("RHF                :", mf.e_tot)
print("K12 = (12|21) =", eri_mo[0, 1, 1, 0])
```

出力は次のとおりです。

```text
反交換関係が成り立つ: True
決定子: ['|1α 1β>', '|1α 2β>', '|1β 2α>', '|2α 2β>']
ハミルトニアン行列（電子部分、Eh）:
 [[-1.83186 -0.       0.       0.18121]
 [ 0.      -1.06467 -0.18121 -0.     ]
 [-0.      -0.18121 -1.06467  0.     ]
 [ 0.18121  0.      -0.      -0.25249]]
最低固有値 + V_NN : -1.137283834488502
PySCF の FCI       : -1.1372838344885023
RHF                : -1.1167593073964255
K12 = (12|21) = 0.18121046201519703
```

$\hat S_z = 0$（$\alpha$ と $\beta$ の電子が 1 個ずつ）の決定子は 4 個です。
決定子はすべて、軌道の番号の小さい順に並べたものです（$|1\beta\,2\alpha\rangle$ は $\hat a_1^\dagger\hat a_2^\dagger|\mathrm{vac}\rangle$）。
行列の各要素を、第 5 節の Slater–Condon 則と照らし合わせます。

- **対角要素**：$|1\alpha\,1\beta\rangle$ のエネルギー $-1.83186$ は、第 3 回の閉殻の式で $2h_{11} + J_{11}$ です。これに核の反発を足したものが RHF のエネルギーです。
- **$|1\alpha\,1\beta\rangle$ と $|2\alpha\,2\beta\rangle$ の間（2 個違い）**：表の 3 行目で $i = 1\alpha$、$j = 1\beta$、$a = 2\alpha$、$b = 2\beta$ とおくと、$[2\alpha\,1\alpha|2\beta\,1\beta] - [2\alpha\,1\beta|2\beta\,1\alpha]$ です。第 2 項はスピンが合わないので 0、第 1 項は空間の積分 $(21|21)$ です。実数の軌道では $(21|21) = (12|21) = K_{12}$ なので、値は最後の行の $K_{12} = 0.18121$ に一致します。
- **$|1\alpha\,2\beta\rangle$ と $|1\beta\,2\alpha\rangle$ の間（2 個違い）**：値は $-K_{12}$ です。$|1\beta\,2\alpha\rangle = -|2\alpha\,1\beta\rangle$ と並べ替えてから表の 3 行目を使うと $+K_{12}$ になり、並べ替えの符号で $-K_{12}$ になります。演算子の行列は、この符号を自動的に扱っています。
- **$|1\alpha\,1\beta\rangle$ と 1 個違いの決定子の間**：0 です。Brillouin の定理のとおりです（この分子では、結合性と反結合性の軌道の対称性が違うことからも 0 になります）。

最後に、この $4\times 4$ 行列を対角化した最低の固有値は、PySCF の完全 CI のエネルギーと一致しました。
決定子 1 個の RHF より $0.0205\,E_\mathrm h$ 低くなっています。
この差が、この基底での電子相関のエネルギーです。
完全 CI は[第 8 回](08-configuration-interaction.md)で詳しく扱います。

中央の $2\times 2$ のブロックは、第 2 回で扱った「$\phi_1$ と $\phi_2$ に 1 個ずつ、反対のスピン」の決定子です。
このブロックを対角化すると、エネルギーは $h_{11} + h_{22} + J_{12} \mp K_{12}$ に分かれます。
これが一重項と三重項に対応することは、第 5 回で確かめます。

**ここまでで分かったこと**：生成・消滅演算子を行列で作り、積分をかけて足すだけで、ハミルトニアン行列ができます。その要素は Slater–Condon 則に一致し、対角化すると完全 CI のエネルギーが得られます。

## まとめ

- 決定子は、スピン軌道の占有数の並び $|n_1 n_2\cdots n_M\rangle$ で表せます。生成演算子を順に真空にかけたものです。
- 生成・消滅演算子の反交換関係 $\{\hat a_p, \hat a_q^\dagger\} = \delta_{pq}$ などに、決定子の反対称性とパウリの排他原理がすべて含まれます。
- ハミルトニアンは $\sum h_{pq}\hat a_p^\dagger\hat a_q + \frac12\sum[pq|rs]\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q$ と書け、分子の情報は積分の値に入ります。
- 決定子の間の行列要素（Slater–Condon 則）は、反交換関係から機械的に出てきます。numpy で作った H₂ のハミルトニアン行列は、この規則と完全 CI のエネルギーに一致しました。

## 付録 A　決定子のエネルギーの導出

$|\Phi\rangle$ を占有数 $n_p$ の決定子とします。
まず 1 電子の部分です。
$\hat a_q|\Phi\rangle$ は $\chi_q$ の電子を取り除いた決定子（$n_q = 0$ なら 0）で、$\hat a_p^\dagger$ でそれを別の決定子に移します。
$p \neq q$ なら結果は $|\Phi\rangle$ と直交するので、

$$
\langle\Phi|\hat a_p^\dagger\hat a_q|\Phi\rangle = \delta_{pq}\,n_p
$$

です。

次に 2 電子の部分です。
$\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q|\Phi\rangle$ が $|\Phi\rangle$ に戻るには、取り除いた 2 つの軌道 $\{q, s\}$ に、同じ 2 つの軌道 $\{p, r\}$ で電子を戻す必要があります。
その方法は $p = q, r = s$ と $p = s, r = q$ の 2 通りです。
前者では $\hat a_p^\dagger\hat a_r^\dagger\hat a_r\hat a_p$ となり、反交換関係で並べ替えると $\hat n_p\hat n_r$（$p\neq r$）になります。
後者では $\hat a_p^\dagger\hat a_r^\dagger\hat a_p\hat a_r = -\hat a_p^\dagger\hat a_r^\dagger\hat a_r\hat a_p$ なので、符号が逆になります。
したがって

$$
\langle\Phi|\hat a_p^\dagger\hat a_r^\dagger\hat a_s\hat a_q|\Phi\rangle = n_p n_r\big(\delta_{pq}\delta_{rs} - \delta_{ps}\delta_{rq}\big)
$$

です（$p = r$ のときは $\hat a_r\hat a_p = 0$ ですが、右辺の括弧も 0 になるので、この式はそのまま成り立ちます）。
これらを $\hat H_\mathrm{el}$ に入れると、

$$
\langle\Phi|\hat H_\mathrm{el}|\Phi\rangle = \sum_p n_p h_{pp} + \frac12\sum_{pr}n_p n_r\big([pp|rr] - [pr|rp]\big)
$$

となり、第 5 節の表の 1 行目が得られます。

## 参考文献

- A. Szabo, N. S. Ostlund, "Modern Quantum Chemistry" (Dover, 1996). ISBN: 9780486691862
- T. Helgaker, P. Jørgensen, J. Olsen, "Molecular Electronic-Structure Theory" (Wiley, 2000). ISBN: 9780471967552
- J. C. Slater, "The Theory of Complex Spectra", Phys. Rev. 34, 1293 (1929). DOI: [10.1103/PhysRev.34.1293](https://doi.org/10.1103/PhysRev.34.1293)
- E. U. Condon, "The Theory of Complex Spectra", Phys. Rev. 36, 1121 (1930). DOI: [10.1103/PhysRev.36.1121](https://doi.org/10.1103/PhysRev.36.1121)
- V. Fock, "Konfigurationsraum und zweite Quantelung", Z. Phys. 75, 622 (1932). DOI: [10.1007/BF01344458](https://doi.org/10.1007/BF01344458)
- P. Jordan, E. Wigner, "Über das Paulische Äquivalenzverbot", Z. Phys. 47, 631 (1928). DOI: [10.1007/BF01331938](https://doi.org/10.1007/BF01331938)
- W. J. Hehre, R. F. Stewart, J. A. Pople, "Self-Consistent Molecular-Orbital Methods. I. Use of Gaussian Expansions of Slater-Type Atomic Orbitals", J. Chem. Phys. 51, 2657 (1969). DOI: [10.1063/1.1672392](https://doi.org/10.1063/1.1672392)
- Q. Sun et al., "Recent developments in the PySCF program package", J. Chem. Phys. 153, 024109 (2020). DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)
