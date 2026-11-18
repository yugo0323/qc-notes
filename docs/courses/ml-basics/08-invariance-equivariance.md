# 第 8 回　対称性を組み込む：不変と等変

!!! abstract "この回の道筋"
    分子を回しても、ずらしても、同じ元素の原子の番号を付け替えても、エネルギーは変わりません。一方、力や双極子モーメントは、分子と一緒に回ります。
    この回では、こうした性質を「群」と「表現」の言葉で最小限に整理し、不変と等変を定義します。
    そのうえで、対称性をモデルに組み込む方法と、データを回して増やす方法（データ拡張）の違いを説明し、スカラー・ベクトル・テンソルの特徴量の作り方を、numpy で回転行列を使って確かめます。

    | 節 | 内容 |
    |---|---|
    | 1 | 分子の対称性：何をしても答えが変わらないか |
    | 2 | 群：変換の集まり |
    | 3 | 表現：量がどう変わるかの規則 |
    | 4 | 不変と等変 |
    | 5 | スカラー・ベクトル・テンソルの特徴量 |
    | 6 | 計算例：回転行列で等変性を確かめる |
    | 7 | 対称性を組み込むか、データで教えるか |

    **前に読んでおく回**：[機械学習の基礎 第 5 回](05-neural-networks.md)、[第 7 回](07-graph-neural-networks.md)

!!! tip "物理的なイメージ"
    不変と等変が表しているのは、分子を回したりずらしたりしたときに、予測がどう振る舞うべきかという約束です。エネルギーは、分子をどう回しても同じ値（不変）でなければなりません。力の矢印は、分子と一緒に回らなければなりません（等変）。この約束をモデルの形に組み込んでおけば、データで覚えさせる必要がなくなり、少ないデータで学習できます。

## 1. 分子の対称性：何をしても答えが変わらないか

分子の構造を、$N$ 個の原子の位置 $\boldsymbol{r}_1, \dots, \boldsymbol{r}_N$（それぞれ 3 次元のベクトル）で表します。
孤立した分子（外から電場などがかかっていない分子）を考えると、次の 4 種類の操作をしても、分子のエネルギーは変わりません。

- **並進**：すべての原子を同じだけずらす。$\boldsymbol{r}_i \to \boldsymbol{r}_i + \boldsymbol{t}$。
- **回転**：分子全体を回す。$\boldsymbol{r}_i \to \mathbf{Q}\boldsymbol{r}_i$（$\mathbf{Q}$ は回転行列）。
- **反転**：すべての位置の符号を変える。$\boldsymbol{r}_i \to -\boldsymbol{r}_i$。鏡に映す操作（鏡映）は、反転と回転を組み合わせたものです。
- **置換**：同じ元素の原子の番号を付け替える。

反転でエネルギーが変わらないのは、鏡像異性体（エナンチオマー）のエネルギーが等しいことに対応します。
置換は、第 7 回で分子グラフについて考えた操作と同じです。

一方、エネルギー以外の量には、これらの操作で変わるものがあります。

- **力** $\boldsymbol{F}_i$：分子を回すと、各原子にかかる力も同じように回ります。原子の番号を付け替えると、どの原子の力かも付け替わります。
- **双極子モーメント** $\boldsymbol{\mu}$：分子を回すと同じように回り、反転すると符号が変わります。
- **分極率** $\boldsymbol{\alpha}$：3×3 の行列で、回転すると $\mathbf{Q}\boldsymbol{\alpha}\mathbf{Q}^\mathsf{T}$ に変わります。

機械学習のモデルで、たとえば力を予測するなら、モデルの出力もこれと同じ規則で変わってほしいわけです。
この「同じ規則で変わる」ことを正確に言うために、2〜4 節で群と表現の言葉を準備します。

**ここまでで分かったこと**：孤立した分子のエネルギーは、並進・回転・反転・同種原子の置換で変わりません。力、双極子モーメント、分極率は、決まった規則で一緒に変わります。

## 2. 群：変換の集まり

### 2.1 定義

1 節の操作には、共通の性質があります。
たとえば、2 つの回転を続けて行うと、また 1 つの回転になります。
何もしない操作（角度 0 の回転）があり、どの回転にも元に戻す逆の回転があります。
このような性質を持つ操作の集まりを、**群**（group）と呼びます。
正確には、集合 $G$ と、その 2 つの要素から 1 つの要素を作る演算（ここでは「続けて行う」こと、$gh$ と書く）が、次の 4 つを満たすとき、$G$ を群と呼びます [Bronstein 2021]。

1. **閉じている**：$g, h \in G$ なら $gh \in G$。
2. **結合法則**：$(gh)k = g(hk)$。
3. **単位元がある**：何もしない要素 $e$ があり、$eg = ge = g$。
4. **逆元がある**：どの $g$ にも $g^{-1}$ があり、$gg^{-1} = g^{-1}g = e$。

化学では、点群（水分子の $C_{2v}$ など）で群に出会っていると思います。
ただし、この回で考える群は、1 つの分子が持つ対称性ではなく、どんな構造の入力にも作用する変換の群です。

### 2.2 分子の機械学習に出てくる群

| 群 | 要素 | 記号 |
|---|---|---|
| 並進群 | 並進 $\boldsymbol{r} \to \boldsymbol{r} + \boldsymbol{t}$ | $T(3)$ |
| 回転群 | 回転行列 $\mathbf{Q}$（$\mathbf{Q}^\mathsf{T}\mathbf{Q} = \mathbf{I}$、$\det\mathbf{Q} = +1$） | $SO(3)$ |
| 直交群 | 回転と反転（$\mathbf{Q}^\mathsf{T}\mathbf{Q} = \mathbf{I}$、$\det\mathbf{Q} = \pm 1$） | $O(3)$ |
| ユークリッド群 | 並進・回転・反転のすべての組み合わせ | $E(3)$ |
| 対称群 | $n$ 個のものの並べ替え | $S_n$ |

$\det\mathbf{Q} = -1$ の直交行列は、反転 $-\mathbf{I}$ と回転の積で書けます。
等変なモデルの論文の題によく出てくる「$E(3)$ 等変」「$SE(3)$ 等変」は、この表の群に対する等変性のことです（$SE(3)$ は並進と回転だけで、反転を含みません）。

**ここまでで分かったこと**：群は、続けて行えて、元に戻せる変換の集まりです。分子では、並進・回転・反転（$E(3)$）と、同種原子の置換（$S_n$）が大事です。

## 3. 表現：量がどう変わるかの規則

### 3.1 定義

1 節で見たように、同じ回転 $\mathbf{Q}$ でも、エネルギーは変わらず、双極子モーメントは $\mathbf{Q}\boldsymbol{\mu}$ に、分極率は $\mathbf{Q}\boldsymbol{\alpha}\mathbf{Q}^\mathsf{T}$ になります。
群の 1 つの要素が、量の種類ごとに違う形で作用するわけです。
これを表すのが**表現**（representation）です。

群 $G$ の $n$ 次元の表現とは、各要素 $g$ に $n \times n$ の正則な行列 $\rho(g)$ を対応させる規則で、次を満たすものです [Bronstein 2021]。

$$
\rho(gh) = \rho(g)\,\rho(h) \qquad (\text{すべての } g, h \in G)
$$

「$g$ と $h$ を続けて行う」ことが、「行列 $\rho(g)$ と $\rho(h)$ を掛ける」ことに対応する、という条件です。
量 $\boldsymbol{a}$（$n$ 個の数を並べたもの）は、$g$ によって $\boldsymbol{a} \to \rho(g)\boldsymbol{a}$ と変わります。

### 3.2 回転と反転の表現の例

$O(3)$ の要素 $\mathbf{Q}$ について、よく出てくる表現を挙げます。

| 量の種類 | $\rho(\mathbf{Q})$ | 化学での例 |
|---|---|---|
| スカラー | $1$（変わらない） | エネルギー、原子間距離 |
| ベクトル | $\mathbf{Q}$ | 位置の差 $\boldsymbol{r}_i - \boldsymbol{r}_j$、力、双極子モーメント |
| 擬ベクトル | $(\det\mathbf{Q})\,\mathbf{Q}$ | 角運動量、磁気モーメント、外積（ベクトル積）$\boldsymbol{u} \times \boldsymbol{v}$ |
| 2 階のテンソル | $\boldsymbol{T} \to \mathbf{Q}\boldsymbol{T}\mathbf{Q}^\mathsf{T}$ | 分極率、四極子モーメント |

ベクトルと擬ベクトルは、回転（$\det\mathbf{Q} = 1$）では同じように変わり、反転（$\mathbf{Q} = -\mathbf{I}$）で違いが出ます。
ベクトル（極性ベクトルとも呼ぶ）は反転で $-1$ 倍になり、擬ベクトルは変わりません。
外積が擬ベクトルになるのは、$(-\boldsymbol{u}) \times (-\boldsymbol{v}) = \boldsymbol{u} \times \boldsymbol{v}$ だからです。
6 節の計算例で、このことを確かめます。

2 階のテンソルの規則も表現です。
$\mathbf{Q}_1\mathbf{Q}_2$ で変換した $\mathbf{Q}_1\mathbf{Q}_2\boldsymbol{T}\mathbf{Q}_2^\mathsf{T}\mathbf{Q}_1^\mathsf{T}$ は、先に $\mathbf{Q}_2$、次に $\mathbf{Q}_1$ で変換したものと同じだからです（9 個の成分を 1 列に並べれば、$\rho(\mathbf{Q})$ は $9 \times 9$ の行列になります）。

### 3.3 置換の表現

$n$ 個の原子の番号の付け替え $\pi \in S_n$ は、第 7 回の置換行列 $\mathbf{P}$ で表せます。
原子ごとの量を行に並べた行列 $\mathbf{H}$（$n$ 行）は、$\mathbf{H} \to \mathbf{P}\mathbf{H}$ と変わります。
エネルギーのような分子全体の量は、置換で変わりません（$\rho = 1$）。

**ここまでで分かったこと**：表現は、群の要素を行列に対応させ、量がどう変わるかを決める規則です。回転と反転では、スカラー、ベクトル、擬ベクトル、2 階のテンソルが代表的な表現です。

## 4. 不変と等変

### 4.1 定義

入力 $\boldsymbol{x}$ から出力 $f(\boldsymbol{x})$ を作るモデルを考えます。
群 $G$ が、入力には表現 $\rho_\text{in}$ で、出力には表現 $\rho_\text{out}$ で作用するとします。
次が成り立つとき、$f$ は $G$ について**等変**（equivariant）であると言います [Bronstein 2021]。

$$
f\bigl(\rho_\text{in}(g)\,\boldsymbol{x}\bigr) = \rho_\text{out}(g)\, f(\boldsymbol{x}) \qquad (\text{すべての } g \in G)
$$

言葉で言えば、「入力を変換してからモデルに通しても、モデルに通してから出力を変換しても、同じ結果になる」ということです。
$\rho_\text{out}(g) = 1$ の特別な場合、つまり $f(\rho_\text{in}(g)\boldsymbol{x}) = f(\boldsymbol{x})$ を、**不変**（invariant）と呼びます。

1 節の例に当てはめると、次のようになります。

- エネルギーを予測するモデル：$E(3)$ と $S_n$ について不変。
- 力を予測するモデル：回転と反転について等変（$\rho_\text{out}$ は各原子に $\mathbf{Q}$）、並進について不変、置換について等変（$\rho_\text{out} = \mathbf{P}$）。
- 双極子モーメントを予測するモデル：回転と反転について等変（$\rho_\text{out} = \mathbf{Q}$）。

### 4.2 不変なエネルギーから、等変な力が出てくる

力は、エネルギーの座標についての微分の符号を変えたもの $\boldsymbol{F}_i = -\nabla_{\boldsymbol{r}_i} E$ です。
エネルギーが並進・回転・反転で不変なら、そこから計算した力は自動的に等変になります。
導出は付録 A に回しますが、結果は次のとおりです。

$$
E(\mathbf{Q}\boldsymbol{r}_1 + \boldsymbol{t}, \dots) = E(\boldsymbol{r}_1, \dots)
\quad\Longrightarrow\quad
\boldsymbol{F}_i(\mathbf{Q}\boldsymbol{r}_1 + \boldsymbol{t}, \dots) = \mathbf{Q}\,\boldsymbol{F}_i(\boldsymbol{r}_1, \dots)
$$

多くの機械学習ポテンシャルが、エネルギーだけを不変なモデルで予測し、力をその微分として求めるのは、このためです。

### 4.3 等変な関数をつなげても等変

モデルを作るうえで最も大事な性質は、次のことです。
$f$ が $\rho_1$ から $\rho_2$ への等変な関数、$h$ が $\rho_2$ から $\rho_3$ への等変な関数なら、それらをつなげた $h \circ f$ は $\rho_1$ から $\rho_3$ への等変な関数です。
実際、次のようになります。

$$
h\bigl(f(\rho_1(g)\boldsymbol{x})\bigr) = h\bigl(\rho_2(g) f(\boldsymbol{x})\bigr) = \rho_3(g)\, h\bigl(f(\boldsymbol{x})\bigr)
$$

したがって、等変な層を何層重ねても全体は等変で、最後の層だけを不変にすれば、全体は不変なモデルになります。
第 7 回のメッセージパッシングの「層の中は置換について等変、読み出しで不変」という作りは、この一例です。

**ここまでで分かったこと**：等変は「変換してから通しても、通してから変換しても同じ」、不変はその特別な場合です。等変な関数をつなげても等変なので、等変な層を重ねて最後に不変な量を取り出せば、不変なモデルが作れます。

## 5. スカラー・ベクトル・テンソルの特徴量

### 5.1 不変な特徴量だけを使う方法

対称性を組み込む最も簡単な方法は、入力の段階で不変な量だけを使うことです。
原子間の距離 $r_{ij} = |\boldsymbol{r}_i - \boldsymbol{r}_j|$ は、並進・回転・反転のどれでも変わりません。
3 つの原子の作る角度も同じです。
距離や角度から作った記述子（[量子化学のための機械学習 第 2 回](../ml-chem/02-descriptors.md)）を入力にすれば、どんなモデルを使っても、出力は不変になります。

ただし、この方法では、モデルの途中の量もすべてスカラーです。
向きを持った量（ベクトル）をモデルの中で扱えず、ベクトルを直接出力することもできません。

### 5.2 等変な特徴量を作る規則

モデルの途中でベクトルやテンソルを扱うには、それらを等変性を壊さずに組み合わせる規則が必要です。
回転と反転（$O(3)$）について、確かめやすい規則を挙げます。
$s$ はスカラー、$\boldsymbol{u}, \boldsymbol{v}$ はベクトルです。

1. **相対位置はベクトル**：$\boldsymbol{r}_i - \boldsymbol{r}_j$ は回転・反転で $\mathbf{Q}$ 倍になり、並進では変わりません。
2. **内積と大きさはスカラー**：$\boldsymbol{u}^\mathsf{T}\boldsymbol{v} \to \boldsymbol{u}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{Q}\boldsymbol{v} = \boldsymbol{u}^\mathsf{T}\boldsymbol{v}$。
3. **スカラー × ベクトルはベクトル**：$s\boldsymbol{u} \to s\mathbf{Q}\boldsymbol{u} = \mathbf{Q}(s\boldsymbol{u})$。ベクトルどうしの和もベクトルです。
4. **直積 $\boldsymbol{u}\boldsymbol{v}^\mathsf{T}$（3×3 の行列）は 2 階のテンソル**：$\mathbf{Q}\boldsymbol{u}\boldsymbol{v}^\mathsf{T}\mathbf{Q}^\mathsf{T}$ と変わります。
5. **外積（ベクトル積）$\boldsymbol{u} \times \boldsymbol{v}$ は擬ベクトル**：回転だけなら等変ですが、反転を含めると $(\det\mathbf{Q})\mathbf{Q}$ で変わります。

規則 2 と 3 を組み合わせると、たとえば次の量は、原子 $i$ にかかるベクトルとして等変です。

$$
\boldsymbol{v}_i = \sum_{j \ne i} \phi(r_{ij})\,(\boldsymbol{r}_i - \boldsymbol{r}_j)
$$

$\phi$ は距離の任意の関数で、ニューラルネットワークで学習させてもかまいません。
Satorras らの EGNN は、原子の位置そのものを、この形（相対位置にスカラーを掛けて足したもの）で更新するモデルです [Satorras 2021]。

### 5.3 活性化関数に注意する

第 5 回の活性化関数は、ベクトルの成分ごとに作用させました。
しかし、ベクトルの量に成分ごとの $\tanh$ を作用させると、等変性が壊れます。
$\tanh(\mathbf{Q}\boldsymbol{u}) \ne \mathbf{Q}\tanh(\boldsymbol{u})$ だからです（6 節で確かめます）。
そのため、等変なモデルでは、非線形な関数はスカラー（ベクトルの大きさなど）にだけ作用させ、その結果をベクトルに掛ける（規則 3）、といった工夫をします。

### 5.4 2 階のテンソルを分ける

2 階のテンソル $\boldsymbol{T}$（9 個の成分）は、回転で互いに混ざらない 3 つの部分に分けられます。

$$
\boldsymbol{T} = \underbrace{\tfrac{1}{3}(\operatorname{tr}\boldsymbol{T})\,\mathbf{I}}_{1\ \text{成分}}
+ \underbrace{\tfrac{1}{2}(\boldsymbol{T} - \boldsymbol{T}^\mathsf{T})}_{3\ \text{成分}}
+ \underbrace{\Bigl[\tfrac{1}{2}(\boldsymbol{T} + \boldsymbol{T}^\mathsf{T}) - \tfrac{1}{3}(\operatorname{tr}\boldsymbol{T})\,\mathbf{I}\Bigr]}_{5\ \text{成分}}
$$

1 つ目はトレースに比例する部分（スカラー）、2 つ目は反対称な部分、3 つ目はトレースが 0 の対称な部分です。
どの部分も、変換すると同じ部分の中だけで混ざり、ほかの部分には移りません。
$\mathbf{Q}\mathbf{I}\mathbf{Q}^\mathsf{T} = \mathbf{I}$ であり、$\mathbf{Q}\boldsymbol{T}\mathbf{Q}^\mathsf{T}$ の転置は $\mathbf{Q}\boldsymbol{T}^\mathsf{T}\mathbf{Q}^\mathsf{T}$ で、トレースも変わらないからです。
反対称な部分の独立な成分は 3 個で、たとえば $\boldsymbol{T} = \boldsymbol{u}\boldsymbol{v}^\mathsf{T}$ なら、外積 $\boldsymbol{u} \times \boldsymbol{v}$ の成分の 1/2 になります。

回転だけを考えると、3 つの部分は、角運動量の量子数 $l = 0, 1, 2$ の量と同じように変わります（成分の数 1、3、5 は、s、p、d 軌道の数と同じです）。
ただし、反転では違いが出ます。
$\boldsymbol{T}$ は反転で $(-\mathbf{I})\boldsymbol{T}(-\mathbf{I})^\mathsf{T} = \boldsymbol{T}$ となり、3 つの部分はどれも変わりません。
したがって、$l = 1$ の反対称な部分は、反転で符号が変わる普通のベクトル（p 軌道と同じ）ではなく、反転で変わらない擬ベクトルです（外積と同じです）。
このように $l$ が同じでも反転での変わり方が違う量があるので、反転まで含めて等変なモデル（NequIP など）では、$l$ に加えて反転での符号（パリティ）でも量を区別します [Batzner 2021]。

量を $l$ ごとに分け、球面調和関数とその組み合わせ（テンソル積）で等変な特徴量を作る方法は、Thomas らの tensor field network [Thomas 2018] などで導入されました。
NequIP などで使われているこの方法は、[量子化学のための機械学習 第 5 回](../ml-chem/05-equivariant-mlips.md)で扱っています。

**ここまでで分かったこと**：相対位置、内積、スカラー倍、直積の規則を使うと、スカラー・ベクトル・テンソルの等変な特徴量が作れます。非線形な関数はスカラーにだけ作用させます。2 階のテンソルは、$l = 0, 1, 2$ の部分に分けられます。

## 6. 計算例：回転行列で等変性を確かめる

### 6.1 回転で確かめる

水分子の座標を使い、5 節の規則で作った量が、本当に決まった規則で変わるかを数値で確かめます。
でたらめな回転行列は、でたらめな行列を QR 分解して作ります。
エネルギーは、距離だけから作った適当な関数です（単位は任意）。
力は、そのエネルギーを中心差分で微分して求めます。

```python
import numpy as np

rng = np.random.default_rng(0)

def random_rotation(rng):
    """一様にでたらめな回転行列（det = +1）"""
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q * np.sign(np.diag(R))                 # QR 分解の符号の任意性を取り除く
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q

# 水分子の座標（Å）。行が原子 O, H, H
R = np.array([[0.0,  0.0,    0.1173],
              [0.0,  0.7572, -0.4692],
              [0.0, -0.7572, -0.4692]])

def pair_vectors(R):
    d = R[:, None, :] - R[None, :, :]           # d[i, j] = r_i - r_j
    r = np.linalg.norm(d, axis=-1)
    return d, r

def energy(R):
    """スカラー（不変）：距離だけの関数"""
    d, r = pair_vectors(R)
    iu = np.triu_indices(len(R), 1)
    return np.sum((1 - np.exp(-(r[iu] - 1.0)))**2)

def vec_feature(R):
    """ベクトル（等変）：原子ごとに Σ_j φ(r_ij) (r_i - r_j)"""
    d, r = pair_vectors(R)
    phi = np.exp(-r**2); np.fill_diagonal(phi, 0.0)
    return np.einsum("ij,ijk->ik", phi, d)

def tensor_feature(R):
    """2 階のテンソル（等変）：原子ごとに Σ_j φ(r_ij) (r_i - r_j)(r_i - r_j)^T"""
    d, r = pair_vectors(R)
    phi = np.exp(-r**2); np.fill_diagonal(phi, 0.0)
    return np.einsum("ij,ijk,ijl->ikl", phi, d, d)

def forces(R, h=1e-5):
    """力 = エネルギーの座標についての微分の符号を変えたもの（中心差分）"""
    F = np.zeros_like(R)
    for i in range(R.shape[0]):
        for k in range(3):
            Rp, Rm = R.copy(), R.copy()
            Rp[i, k] += h; Rm[i, k] -= h
            F[i, k] = -(energy(Rp) - energy(Rm)) / (2 * h)
    return F

Q = random_rotation(rng)
t = rng.normal(size=3)                          # 並進
RQ = R @ Q.T + t                                # 各原子の位置 r_i を Q r_i + t に移す
print("det Q =", round(np.linalg.det(Q), 6))
print("エネルギー（不変）   :", np.isclose(energy(RQ), energy(R)))
print("ベクトル（等変）     :", np.allclose(vec_feature(RQ), vec_feature(R) @ Q.T))
print("テンソル（等変）     :", np.allclose(tensor_feature(RQ), Q @ tensor_feature(R) @ Q.T))
print("力（等変）           :", np.allclose(forces(RQ), forces(R) @ Q.T, atol=1e-8))
```

```text
det Q = 1.0
エネルギー（不変）   : True
ベクトル（等変）     : True
テンソル（等変）     : True
力（等変）           : True
```

座標は原子ごとに行に並べているので、各行に $\mathbf{Q}$ を掛ける操作は `R @ Q.T` と書けます。
回転と並進を同時に行っても、エネルギーは変わらず、ベクトル・テンソル・力は 3.2 節の規則どおりに変わりました（力については 4.2 節の結果の確認です）。

### 6.2 反転で確かめる

次に、反転 $-\mathbf{I}$ と、回転と反転を組み合わせた $\mathbf{M} = \mathbf{Q}(-\mathbf{I})$（$\det\mathbf{M} = -1$）で調べます。
擬ベクトルの例として、O から見た 2 つの O–H ベクトルの外積 $\boldsymbol{n}$（分子面に垂直なベクトル）を使います。

```python
# 上の続き：反転（すべての位置 r を -r に移す）
Inv = -np.eye(3)
RI = R @ Inv.T

def normal_vector(R):
    """O から見た 2 つの O-H ベクトルの外積（分子面に垂直なベクトル）"""
    return np.cross(R[1] - R[0], R[2] - R[0])

print("ベクトル特徴は -1 倍になるか:", np.allclose(vec_feature(RI), -vec_feature(R)))
print("外積:", np.round(normal_vector(R), 4) + 0.0, "→ 反転後", np.round(normal_vector(RI), 4) + 0.0)

# 回転と反転を組み合わせた変換（det = -1）でも確かめる
M = Q @ Inv
print("det M =", round(np.linalg.det(M), 6))
print("外積は M n になるか      :", np.allclose(normal_vector(R @ M.T), M @ normal_vector(R)))
print("外積は det(M) M n になるか:", np.allclose(normal_vector(R @ M.T), np.linalg.det(M) * M @ normal_vector(R)))
```

```text
ベクトル特徴は -1 倍になるか: True
外積: [-0.8882  0.      0.    ] → 反転後 [-0.8882  0.      0.    ]
det M = -1.0
外積は M n になるか      : False
外積は det(M) M n になるか: True
```

相対位置から作ったベクトル特徴は、反転で $-1$ 倍になりましたが、外積は変わりませんでした。
$\det\mathbf{M} = -1$ の変換では、外積は $\mathbf{M}\boldsymbol{n}$ ではなく $(\det\mathbf{M})\mathbf{M}\boldsymbol{n}$ に変わりました。
3.2 節の表のとおり、外積は擬ベクトルです。

### 6.3 等変にならない例

最後に、等変にならない例を 2 つ見ます。
1 つは、座標の 9 個の数をそのまま多層パーセプトロンに入れる「エネルギー」です。
もう 1 つは、5.3 節の、ベクトルに成分ごとの $\tanh$ を作用させる操作です。

```python
# 上の続き：座標をそのまま入れる多層パーセプトロン（重みはでたらめ）は不変にならない
W1 = rng.normal(size=(16, 9)); W2 = rng.normal(size=16)
def mlp_energy(R):
    return W2 @ np.tanh(W1 @ R.reshape(-1))

E_rot = [mlp_energy(R @ random_rotation(rng).T) for _ in range(1000)]
print(f"もとの向き: {mlp_energy(R):.3f}   1000 通りに回したとき: 最小 {min(E_rot):.3f}, 最大 {max(E_rot):.3f}")

# ベクトルの成分ごとに tanh を作用させると、等変でなくなる
v = vec_feature(R)[1]                           # 原子 1 のベクトル特徴
print("tanh(Q v) = Q tanh(v) か:", np.allclose(np.tanh(Q @ v), Q @ np.tanh(v)))
print("|Q v| = |v| か          :", np.isclose(np.linalg.norm(Q @ v), np.linalg.norm(v)))
```

```text
もとの向き: 2.829   1000 通りに回したとき: 最小 -7.250, 最大 7.137
tanh(Q v) = Q tanh(v) か: False
|Q v| = |v| か          : True
```

座標をそのまま入れたモデルは、同じ分子を回しただけで、出力が $-7.3$ から $7.1$ まで変わりました。
成分ごとの $\tanh$ は回転と順番を入れ替えられません。一方、ベクトルの大きさ（スカラー）は回転で変わらないので、非線形な関数を作用させてもかまいません（5.3 節）。

**ここまでで分かったこと**：距離から作ったエネルギーは不変、相対位置から作ったベクトルとテンソル、エネルギーを微分した力は等変であることを、数値で確かめました。外積は反転で擬ベクトルとして変わります。座標をそのまま入れたモデルや、成分ごとの活性化関数は、等変性を壊します。

## 7. 対称性を組み込むか、データで教えるか

### 7.1 データ拡張

対称性を扱うもう 1 つの方法は、モデルは普通のものを使い、訓練データを変換して増やすことです。
分子なら、訓練データの各構造をでたらめに回したものを、同じエネルギーのデータとして使います。
この方法を**データ拡張**（data augmentation）と呼びます [Goodfellow 2016, 7.4 節]。

データ拡張は、モデルの形を変えずに使えるので手軽で、画像の分野などで広く使われています [Bronstein 2021]。
しかし、モデルが学ぶのは「訓練で見た向きの範囲で、おおよそ不変」なことだけで、対称性が厳密に成り立つ保証はありません。
また、Bronstein らの総説は、データ拡張は必要なデータの量の点で最適でないことが示されている、と述べています [Bronstein 2021]。

### 7.2 計算で比べる

6 節の続きで、O–H 距離 0.9〜1.1 Å、HOH 角 95〜115° の水分子 500 個で、6.1 節のエネルギーをカーネルリッジ回帰（第 4 回）で学習させます。
次の 3 つを比べます。

- 座標の 9 個の数を入力にし、訓練データは決まった向きのまま（拡張なし）
- 座標の 9 個の数を入力にし、訓練データの各分子をでたらめな向きに回したものに置き換える（データ拡張）
- 2 つの O–H 距離（小さい順に並べる）と H–H 距離を入力にする（不変な記述子）

ガウスカーネルの幅 $\gamma$ は、第 3 回と同じく、訓練にもテストにも使わない検証データ（でたらめな向きの 100 個の分子）での誤差が最も小さくなる値を、$10^{-4}$〜$100$ の中から選びます。
選んだ $\gamma$ で、別に用意したテストデータ（でたらめな向きの 200 個の分子）の誤差を 1 回だけ計算します。
さらに、1 つの分子を 200 通りの向きに回して予測し、その標準偏差（予測が向きでどれだけ揺れるか）も調べます。

```python
# 上の続き：データ拡張と、不変な記述子を比べる
from sklearn.kernel_ridge import KernelRidge

def water(rng):
    """O-H 距離 0.9〜1.1 Å、HOH 角 95〜115° の水分子（決まった向き）"""
    r1, r2 = rng.uniform(0.9, 1.1, size=2)
    th = np.deg2rad(rng.uniform(95, 115))
    return np.array([[0, 0, 0],
                     [0, r1 * np.sin(th / 2), r1 * np.cos(th / 2)],
                     [0, -r2 * np.sin(th / 2), r2 * np.cos(th / 2)]])

def distances(R):
    """不変な記述子：2 つの O-H 距離（小さい順）と H-H 距離"""
    d, r = pair_vectors(R)
    return np.array(sorted([r[0, 1], r[0, 2]]) + [r[1, 2]])

train = [water(rng) for _ in range(500)]
train_rot = [Ri @ random_rotation(rng).T for Ri in train]   # データ拡張：でたらめに回す
val = [water(rng) @ random_rotation(rng).T for _ in range(100)]    # 検証データ（γ を選ぶ）
test = [water(rng) @ random_rotation(rng).T for _ in range(200)]   # テストデータ（最後に 1 回だけ使う）
y_tr, y_val, y_te = (np.array([energy(Ri) for Ri in S]) for S in (train, val, test))

def fit(feat, train_set, gamma):
    krr = KernelRidge(kernel="rbf", alpha=1e-8, gamma=gamma)
    return krr.fit(np.array([feat(Ri) for Ri in train_set]), y_tr)

def rmse(krr, feat, data, y):
    return np.sqrt(np.mean((krr.predict(np.array([feat(Ri) for Ri in data])) - y)**2))

flat = lambda Ri: Ri.reshape(-1)
gammas = [1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]   # 試す γ
R0 = water(rng)                                 # 1 つの分子を 200 通りに回して、予測の揺れを見る
R0_rots = [R0 @ random_rotation(rng).T for _ in range(200)]
print(f"テストのエネルギーの標準偏差 {y_te.std():.1e}")
for name, feat, train_set in [("座標そのまま、拡張なし  ", flat, train),
                              ("座標そのまま、回転で拡張", flat, train_rot),
                              ("距離の記述子（不変）    ", distances, train)]:
    gamma = min(gammas, key=lambda g: rmse(fit(feat, train_set, g), feat, val, y_val))  # 検証データで選ぶ
    krr = fit(feat, train_set, gamma)
    err = rmse(krr, feat, test, y_te)
    spread = np.std(krr.predict(np.array([feat(Ri) for Ri in R0_rots])))
    print(f"{name}: γ = {gamma:<5g} テスト RMSE {err:.1e}, 回したときの予測の標準偏差 {spread:.1e}")
```

```text
テストのエネルギーの標準偏差 4.4e-02
座標そのまま、拡張なし  : γ = 1     テスト RMSE 1.3e-01, 回したときの予測の標準偏差 1.0e-01
座標そのまま、回転で拡張: γ = 0.001 テスト RMSE 4.1e-03, 回したときの予測の標準偏差 1.4e-03
距離の記述子（不変）    : γ = 10    テスト RMSE 5.3e-06, 回したときの予測の標準偏差 9.1e-12
```

### 7.3 結果を読む

- **拡張なし**では、テスト RMSE（$1.3 \times 10^{-1}$）が、テストのエネルギーそのもののばらつき（$4.4 \times 10^{-2}$）より大きくなりました。訓練で見たことのない向きの分子は、まったく予測できていません。
- **データ拡張**で、誤差は $4.1 \times 10^{-3}$ まで下がりました。しかし、同じ分子を回しただけで、予測が $1.4 \times 10^{-3}$ ほど揺れます。不変性は近似的にしか成り立っていません。
- **不変な記述子**では、誤差は $5.3 \times 10^{-6}$ で、データ拡張よりおよそ 3 桁小さくなりました。回したときの揺れ（$9 \times 10^{-12}$）は、数値計算の誤差の程度です。不変な記述子では、データを向きの違いを覚えるのに使わずに済むためと考えられます。

この例は、不変な記述子が簡単に作れる 3 原子の場合です。
大きな分子や、ベクトル・テンソルを出力するモデルでは、5 節の等変な特徴量をモデルの中で使って、同じ利点を得ようとします。
等変なモデルが少ないデータで高い精度を示したという報告は、[量子化学のための機械学習 第 5 回](../ml-chem/05-equivariant-mlips.md)で紹介しています。

**ここまでで分かったこと**：データ拡張は手軽ですが、対称性は近似的にしか成り立たず、データの一部を向きの違いを覚えるのに使います。不変・等変な形をモデルに組み込むと、対称性は厳密に成り立ち、同じデータ量でずっと小さな誤差になりました。

## まとめ

- 孤立した分子のエネルギーは並進・回転・反転・同種原子の置換で不変で、力や双極子モーメントは等変です。
- 群は「続けて行えて、元に戻せる変換の集まり」、表現は「群の要素が量にどう作用するか」を行列で表す規則です。
- 等変は $f(\rho_\text{in}(g)\boldsymbol{x}) = \rho_\text{out}(g)f(\boldsymbol{x})$、不変は $\rho_\text{out} = 1$ の場合です。等変な層を重ね、最後に不変な量を取り出せば、不変なモデルになります。
- 相対位置・内積・スカラー倍・直積で、スカラー・ベクトル・テンソルの等変な特徴量が作れます。外積は擬ベクトルで、反転では普通のベクトルと違う変わり方をします。
- データ拡張は対称性を近似的にしか満たしません。対称性をモデルに組み込むと厳密に満たされ、この回の例では同じデータ量でずっと小さな誤差になりました。

## 付録 A　不変なエネルギーから等変な力が出てくること

エネルギー $E(\boldsymbol{r}_1, \dots, \boldsymbol{r}_N)$ が、すべての直交行列 $\mathbf{Q}$ と並進 $\boldsymbol{t}$ について次を満たすとします。

$$
E(\mathbf{Q}\boldsymbol{r}_1 + \boldsymbol{t}, \dots, \mathbf{Q}\boldsymbol{r}_N + \boldsymbol{t}) = E(\boldsymbol{r}_1, \dots, \boldsymbol{r}_N)
$$

$\boldsymbol{r}_i' = \mathbf{Q}\boldsymbol{r}_i + \boldsymbol{t}$ と書きます。
両辺を $\boldsymbol{r}_i$ の成分 $r_{i,a}$（$a = x, y, z$）で微分します。
左辺は、連鎖律と $\partial r'_{i,b}/\partial r_{i,a} = Q_{ba}$ を使って、次のようになります。

$$
\sum_b \frac{\partial E}{\partial r'_{i,b}}(\boldsymbol{r}') \, Q_{ba} = \frac{\partial E}{\partial r_{i,a}}(\boldsymbol{r})
$$

ベクトルで書くと $\mathbf{Q}^\mathsf{T}\nabla_{\boldsymbol{r}_i'}E(\boldsymbol{r}') = \nabla_{\boldsymbol{r}_i}E(\boldsymbol{r})$ です。
両辺に左から $\mathbf{Q}$ を掛け、$\mathbf{Q}\mathbf{Q}^\mathsf{T} = \mathbf{I}$ を使い、符号を変えると、次が得られます。

$$
\boldsymbol{F}_i(\boldsymbol{r}') = \mathbf{Q}\,\boldsymbol{F}_i(\boldsymbol{r})
$$

並進 $\boldsymbol{t}$ は微分で消えるので、力は並進で変わりません。
回転（$\det\mathbf{Q} = 1$）でも反転（$\det\mathbf{Q} = -1$）でも同じ式になるので、力は擬ベクトルではなく、普通のベクトルとして変わります。

## 参考文献

- M. M. Bronstein, J. Bruna, T. Cohen, P. Veličković, "Geometric Deep Learning: Grids, Groups, Graphs, Geodesics, and Gauges", arXiv:2104.13478 (2021).
- I. Goodfellow, Y. Bengio, A. Courville, "Deep Learning" (MIT Press, 2016). ISBN: 9780262035613
- V. G. Satorras, E. Hoogeboom, M. Welling, "E(n) Equivariant Graph Neural Networks", arXiv:2102.09844 (2021).
- S. Batzner, A. Musaelian, L. Sun, M. Geiger, J. P. Mailoa, M. Kornbluth, N. Molinari, T. E. Smidt, B. Kozinsky, "E(3)-Equivariant Graph Neural Networks for Data-Efficient and Accurate Interatomic Potentials", arXiv:2101.03164 (2021).
- N. Thomas, T. Smidt, S. Kearnes, L. Yang, L. Li, K. Kohlhoff, P. Riley, "Tensor field networks: Rotation- and translation-equivariant neural networks for 3D point clouds", arXiv:1802.08219 (2018).
