# 第 4 回　初期の機械学習ポテンシャル：Behler–Parrinello、ANI、SchNet

!!! abstract "この回の道筋"
    第 2 回と第 3 回で準備した考え方（原子ごとの記述子、原子エネルギーの和）が、実際のモデルの中でどう使われてきたかを、3 つの代表的なモデルで見ます。Behler–Parrinello のモデルが枠組みを作り、ANI がそれを有機分子に広げ、SchNet が記述子そのものをデータから学ぶ方向に進めました。それぞれ何が新しかったのかを、原論文に沿って説明します。

    | 節 | 内容 |
    |---|---|
    | 1 | Behler–Parrinello（2007 年）：高次元ニューラルネットワークポテンシャル |
    | 2 | ANI-1（2017 年）：有機分子への拡張と大量のデータ |
    | 3 | SchNet（2017〜2018 年）：記述子を学習する |
    | 4 | 3 つのモデルの比較と、残った課題 |

    **前に読んでおく回**：[第 2 回](02-descriptors.md)、[第 3 回](03-mlip-basics.md)、[機械学習の基礎 第 5 回](../ml-basics/05-neural-networks.md)、[機械学習の基礎 第 7 回](../ml-basics/07-graph-neural-networks.md)

!!! tip "物理的なイメージ"
    初期の機械学習ポテンシャルが行っていたのは、第 2 回の指紋（記述子）を原子ごとに計算し、それをニューラルネットワークに入れて原子のエネルギーを出し、足し合わせることです。Behler–Parrinello のモデルがこの枠組みをつくり、ANI はそれを多くの有機分子に広げました。SchNet は、指紋を人が設計する代わりに、原子どうしの距離から指紋そのものをデータで学ぶようにしました。

## 1. Behler–Parrinello（2007 年）：高次元ニューラルネットワークポテンシャル

### 1.1 それ以前の問題

ニューラルネットワークで PES を表す試みは、2007 年より前からありました。
初期のモデルは、原子間の距離や角度などの内部座標を、決まった順番で並べて入力にしていました [Unke 2021]。
この方法には、第 2 回で見た問題がそのまま現れます。
内部座標の並べ方が原子の番号の付け方に依存するので、対称的に同じ構造に違うエネルギーを返すことがありえます。
また、入力の長さが原子の数で決まるので、学習した系と同じ大きさの系にしか使えません [Unke 2021]。

### 1.2 新しかったこと

Behler と Parrinello は、2007 年の Physical Review Letters の論文で、この 2 つの問題を同時に解決しました [Behler 2007]。
その中身は、第 2 回と第 3 回で説明した 2 つの考え方の組み合わせです。

1. 各原子のまわりの環境を、原子中心対称関数 $\mathbf G_i$ で表す（対称関数の性質は、Behler の 2011 年の論文で詳しく調べられています [Behler 2011]）。
2. 全エネルギーを原子エネルギーの和 $E = \sum_i E_i$ で表し、$E_i$ を同じ元素で共通のニューラルネットワークで計算する。

この形のモデルは、高次元ニューラルネットワークポテンシャル（high-dimensional neural network potential、HDNNP）と呼ばれます。
Unke らの総説は、Behler と Parrinello がこの HDNNP を初めて提案したとしています [Unke 2021]。
原論文は、密度汎関数理論（DFT）で計算したケイ素の固体を例に、エネルギーと力を DFT より何桁も速く計算できることを示しました [Behler 2007]。

### 1.3 原子ごとのニューラルネットワーク

$E_i$ を計算するニューラルネットワークは、全結合の層を重ねたものです。
たとえば隠れ層が 2 つなら、次のように書けます。

$$
E_i = \mathbf w^{(3)} \cdot f\!\left( \mathsf W^{(2)} f\!\left( \mathsf W^{(1)} \mathbf G_i + \mathbf b^{(1)} \right) + \mathbf b^{(2)} \right) + b^{(3)} .
$$

$\mathsf W^{(l)}$ と $\mathbf b^{(l)}$ は学習で決めるパラメーター、$f$ は $\tanh$ のような非線形の関数（活性化関数）で、ベクトルの成分ごとに作用します。
ニューラルネットワークの一般論は、[機械学習の基礎 第 5 回](../ml-basics/05-neural-networks.md)で扱っています。

ここで大事なのは、ネットワークの入力が記述子 $\mathbf G_i$ で、その長さが原子の数によらないことです。
同じ元素の原子は同じパラメーターを使うので、原子が何個あっても、同じネットワークを原子の数だけ使えば済みます。
複数の元素を含む系では、記述子を元素の組み合わせごとに用意し、元素ごとに別のネットワークを使うのが普通です [Unke 2021]。

**ここまでで分かったこと**：Behler–Parrinello のモデルは、原子中心対称関数と、原子エネルギーの和を組み合わせて、どんな大きさの系にも使えるニューラルネットワークポテンシャルを初めて作りました。

## 2. ANI-1（2017 年）：有機分子への拡張と大量のデータ

### 2.1 目標

Behler–Parrinello の原論文が扱ったのは、ケイ素という 1 つの元素の系でした。
Smith、Isayev、Roitberg は、2017 年の Chemical Science の論文で、H、C、N、O からなる有機分子一般に使えるポテンシャル ANI-1 を発表しました [Smith 2017]。
「一般に使える」とは、学習に使っていない分子、それも学習に使った分子より大きな分子にも使える、という意味です。

### 2.2 記述子の工夫

ANI-1 は、Behler–Parrinello の対称関数に手を加えた記述子を使っています [Smith 2017]。
論文では、これを原子環境ベクトル（atomic environment vector、AEV）と呼んでいます。

動径の部分は、第 2 回の動径の対称関数と同じ形で、幅 $\eta$ を 1 つに固定し、中心 $R_\mathrm s$ を多数並べます。
角度の部分には、2 つの変更が加えられています。

$$
G_i^{\text{角度}} = 2^{1-\zeta} \sum_{j, k \neq i}
\left(1 + \cos(\theta_{ijk} - \theta_\mathrm s)\right)^{\zeta}
\exp\!\left[-\eta \left(\frac{R_{ij} + R_{ik}}{2} - R_\mathrm s\right)^2\right]
f_\mathrm c(R_{ij})\, f_\mathrm c(R_{ik}).
$$

1 つ目の変更は、角度のずらし $\theta_\mathrm s$ です。
第 2 回の角度の関数は、$\lambda = \pm 1$ で 0° 側か 180° 側かを見るだけでした。
$\theta_\mathrm s$ をいくつも用意すると、特定の角度の付近を選んで調べられます。
2 つ目の変更は、2 つの距離の平均を $R_\mathrm s$ のまわりで見るガウス関数です。
これで、角度の情報を「中心からどれくらいの距離にある原子の組か」で分けて調べられます。

ANI-1 では、動径の部分のカットオフ半径を 4.6 Å、角度の部分を 3.1 Å とし、4 種類の元素の組み合わせを合わせて、1 原子あたり 768 個の数からなる記述子にしています [Smith 2017]。

### 2.3 データの規模

ANI-1 で新しかったもう 1 つの点は、学習データの量と作り方です。
GDB-11 というデータベースから、重原子（H 以外の原子）が 8 個以下の分子、約 5 万 8 千個を選びました [Smith 2017]。
各分子について、エネルギー最小の構造で基準振動を求め、基準振動の座標に沿ってランダムにずらした構造を多数作りました。
ずらす大きさは、調和振動子のエネルギーがある上限を超えないように決めています。
論文ではこの方法を基準振動サンプリング（normal mode sampling）と呼んでいます。
こうして作った約 1720 万個の構造について、ωB97X/6-31G(d) でエネルギーを計算し、学習データにしています [Smith 2017]。

### 2.4 何が示されたか

論文は、学習に使っていない重原子 10 個の分子 134 個（8245 構造）で精度を調べています [Smith 2017]。
二乗平均平方根誤差は、全体で 1.8 kcal/mol、各分子のエネルギー最小の構造から 30 kcal/mol 以内の構造に限ると 0.6 kcal/mol でした。
要旨では、学習データの分子よりずっと大きな分子（54 原子まで）でも、参照の DFT の計算と比べて化学的精度が得られたと述べています。
第 1 回で注意したように、これらは学習に使った ωB97X/6-31G(d) に対する誤差です。
ANI-1 の題にある「DFT の精度」も、この意味です。

その後、同じグループは、DFT で学習したモデルを CCSD(T)/CBS のデータで学習し直す ANI-1ccx を発表しています（第 1 回で紹介しました）[Smith 2019]。

**ここまでで分かったこと**：ANI-1 は、対称関数に角度と距離のずらしを加えた記述子と、約 1720 万個の DFT のデータで、有機分子に広く使えるポテンシャルを作りました。

## 3. SchNet（2017〜2018 年）：記述子を学習する

### 3.1 記述子を人が決めることの限界

Behler–Parrinello と ANI では、記述子の形（$\eta$、$R_\mathrm s$、$\theta_\mathrm s$ などの組）を人が前もって決めます。
選び方がまずく、まったく違う 2 つの構造がほとんど同じ記述子になってしまうと、その後のネットワークをいくら大きくしても区別できません [Unke 2021]。
また、元素の数が増えると、元素の組み合わせごとに記述子を用意するので、記述子の長さがどんどん増えます [Unke 2021]。

そこで、記述子そのものもデータから学習する方向が生まれました。
Unke らの総説は、前者を「記述子にもとづく」モデル、後者を「端から端まで学習する」（end-to-end）モデルと分けています [Unke 2021]。
その代表が SchNet です。

### 3.2 SchNet の流れ

SchNet は、Schütt らが 2017 年に NIPS（現在の NeurIPS）という会議で発表し [Schütt 2017]、2018 年に Journal of Chemical Physics で詳しく報告したモデルです [Schütt 2018]。
同じグループの DTNN（deep tensor neural network）の考え方を発展させたものです [Schütt 2017b]。
計算は次の順に進みます [Schütt 2018]。

1. **原子の初期ベクトル**：原子 $i$ に、元素 $Z_i$ ごとに決まるベクトル $\mathbf x_i^{(0)}$ を割り当てます。このベクトルは、最初はランダムで、学習で決まります（埋め込み、embedding）。
2. **相互作用ブロック**：各原子のベクトルを、まわりの原子のベクトルと距離を使って更新します。これを $T$ 回くり返します。
3. **原子エネルギー**：最後のベクトル $\mathbf x_i^{(T)}$ から、小さなネットワークで原子エネルギー $E_i$ を計算し、和を取ります。

### 3.3 連続フィルター畳み込み

相互作用ブロックの中心は、次の更新です。

$$
\mathbf x_i^{(t+1)} = \mathbf x_i^{(t)} + \sum_{j \neq i} \mathbf x_j^{(t)} \circ \mathbf W^{(t)}(r_{ij}) .
$$

（実際のモデルでは、この前後に原子ごとの全結合の層が入ります。）
$\circ$ は成分ごとの積です。
$\mathbf W^{(t)}(r)$ は、原子間の距離 $r$ からベクトルを作る小さなニューラルネットワークで、フィルター生成ネットワークと呼ばれます。
画像の畳み込みでは、フィルターは格子点ごとの値の表です。
原子は格子の上にないので、SchNet ではフィルターを距離の連続関数にしました。
これが「連続フィルター畳み込み」（continuous-filter convolution）という名前の由来です [Schütt 2017]。

フィルターの入力には、原子間のベクトルではなく距離を使います。
そのため、更新されたベクトルは回転で変わりません [Schütt 2018]。
距離はそのまま使わず、0 からカットオフ半径まで 0.1 Å 間隔で中心を並べたガウス関数

$$
e_k(r) = \exp\!\left[-\gamma (r - \mu_k)^2\right]
$$

の値の並びに直してから、ネットワークに入れます [Schütt 2018]。
これは、第 2 回の動径の対称関数とよく似た形です。

力は、エネルギーを座標で微分して求めます。
こうすると、回転で変わらないエネルギーから、回転と一緒に回る力が得られ、力は保存力になります [Schütt 2018]。

### 3.4 メッセージパッシングとして見る

SchNet の更新は、「隣の原子からメッセージを受け取り、自分のベクトルを書き換える」操作と見ることができます。
Gilmer らは 2017 年に、分子のためのいくつかのニューラルネットワークを、この形の共通の枠組みにまとめ、メッセージパッシングニューラルネットワーク（message passing neural network、MPNN）と名付けました [Gilmer 2017]。
グラフニューラルネットワークの一般論は、[機械学習の基礎 第 7 回](../ml-basics/07-graph-neural-networks.md)で扱います。

メッセージパッシングには、記述子にもとづくモデルにない性質があります。
1 回の更新で、原子 $i$ はカットオフ半径の中の原子から情報を受け取ります。
2 回目の更新では、その原子たちが 1 回目に受け取った情報も届きます。
したがって、更新を $T$ 回くり返すと、原子 $i$ には、最大でカットオフ半径の $T$ 倍の距離にある原子の影響が届きます。
この範囲を受容野（receptive field）と呼びます [Batatia 2022]。

### 3.5 計算例：2 回の更新で情報が届く距離

SchNet と同じ形の更新を、ランダムな重みで 2 回行う小さなモデルを numpy で作ります。
学習はしていません。
C 原子 3 個を 3 Å 間隔で一直線に並べ、カットオフ半径を 4 Å とします。
原子 0 と原子 2 は 6 Å 離れているので、原子 2 は原子 0 のカットオフ半径の外にあります。
原子 2 だけを少し動かして、原子 0 のベクトルが変わるかを調べます。

```python
import numpy as np

rng = np.random.default_rng(0)
RC = 4.0                                   # 各層のカットオフ半径（Å）
MU = np.arange(0.0, RC, 0.1)               # ガウス関数の中心（0.1 Å 間隔）
GAMMA = 10.0                               # ガウス関数の幅（Å^-2）
F = 4                                      # 原子の特徴ベクトルの長さ
A = [rng.normal(size=(F, len(MU))) * 0.3 for _ in range(2)]   # 層ごとの重み
embed = {1: rng.normal(size=F), 6: rng.normal(size=F)}        # 元素ごとの初期ベクトル

def fc(r):
    return 0.5 * (np.cos(np.pi * r / RC) + 1.0) * (r < RC)

def filt(r, layer):
    """距離 r から、長さ F のフィルターを作る（距離をガウス関数で展開してから変換）"""
    e = np.exp(-GAMMA * (r - MU) ** 2)
    return np.tanh(A[layer] @ e) * fc(r)

def run(Z, R, n_layers):
    x = np.array([embed[z] for z in Z])
    for t in range(n_layers):
        new = x.copy()
        for i in range(len(Z)):
            for j in range(len(Z)):
                if i != j:
                    r = np.linalg.norm(R[i] - R[j])
                    new[i] += x[j] * filt(r, t)      # 隣の原子からの「メッセージ」
        x = new
    return x

# C 原子 3 個を 3 Å 間隔で一直線に並べる。原子 0 と原子 2 は 6 Å 離れている
Z = [6, 6, 6]
R = np.array([[0.0, 0, 0], [3.0, 0, 0], [6.0, 0, 0]])
R2 = R.copy(); R2[2, 0] += 0.2                     # 原子 2 だけを少し動かす
np.set_printoptions(precision=5, suppress=True)
for n in (1, 2):
    d = run(Z, R2, n)[0] - run(Z, R, n)[0]
    print(f"{n} 層：原子 2 を動かしたときの、原子 0 の特徴ベクトルの変化 {d}")
```

```text
1 層：原子 2 を動かしたときの、原子 0 の特徴ベクトルの変化 [0. 0. 0. 0.]
2 層：原子 2 を動かしたときの、原子 0 の特徴ベクトルの変化 [-0.00414  0.00004  0.00067 -0.     ]
```

1 回の更新では、原子 0 のベクトルはまったく変わりません。
2 回の更新では変わります。
原子 2 の動きが、原子 1 を経由して原子 0 に届いたのです。

これは、第 3 回の「$E_i$ はカットオフ半径の中の原子だけで決まる」という局所性が、メッセージパッシングでは「カットオフ半径の $T$ 倍の中の原子で決まる」に広がることを意味します。
遠くの影響を取り込めるという利点がある一方、計算を複数の GPU に分けるのが難しくなるという指摘もあります [Batatia 2022]。

**ここまでで分かったこと**：SchNet は、元素ごとの初期ベクトルを、距離の連続関数のフィルターを使った更新で書き換えることで、記述子をデータから学習します。更新をくり返すと、情報はカットオフ半径の外まで届きます。

## 4. 3 つのモデルの比較と、残った課題

### 4.1 比較

3 つのモデルを並べます。

| | Behler–Parrinello | ANI-1 | SchNet |
|---|---|---|---|
| 発表 | 2007 年 | 2017 年 | 2017 年（会議）、2018 年（論文誌） |
| 原子環境の表し方 | 原子中心対称関数（人が決める） | 対称関数を改良した AEV（人が決める） | 距離のガウス展開から更新で学習 |
| 元素の扱い | 元素ごとに記述子とネットワーク | 元素ごとに記述子とネットワーク | 元素ごとの初期ベクトル（学習） |
| エネルギー | 原子エネルギーの和 | 原子エネルギーの和 | 原子エネルギーの和 |
| 力 | エネルギーの微分 | エネルギーの微分 | エネルギーの微分 |

3 つとも、全エネルギーを原子エネルギーの和で表し、力をエネルギーの微分で求める点は同じです。
違うのは、原子環境の表し方です。

### 4.2 不変な量だけを使うことの限界

3 つのモデルに共通するもう 1 つの点は、モデルの中で扱う量が、すべて回転で変わらない量（距離、角度、それらから作ったベクトル）だということです。
NequIP の論文は、当時の対称性を考慮したモデルの多くが、回転で変わらない量（スカラー）だけを扱っていたと指摘しています [Batzner 2022]。

不変な量だけで原子環境を表すと、2 つの困りごとがあります。

- SchNet の更新では、メッセージは相手の原子のベクトルと、2 つの原子の距離だけで決まります。角度の情報は、更新をくり返すことで間接的にしか入りません。MACE の論文は、このような 2 体の不変なメッセージによる原子環境の表し方が完全ではない（違う環境を区別できないことがある）と指摘しています [Batatia 2022]。
- 第 2 回で見たように、距離と角度から作った不変な記述子でも、違う原子環境を区別できない例があります [Pozdnyakov 2020]。

この限界を越えるために、モデルの中で「回転と一緒に回る量」（ベクトルやテンソル）を扱う等変なモデルが生まれました。
それが第 5 回の NequIP と MACE です。

**ここまでで分かったこと**：3 つのモデルは、原子エネルギーの和という枠組みを共有し、原子環境の表し方で違います。モデルの中で不変な量だけを使うことには限界があります。

## まとめ

- Behler–Parrinello（2007 年）は、原子中心対称関数と原子エネルギーの和を組み合わせ、どんな大きさの系にも使えるニューラルネットワークポテンシャルの枠組みを作りました。
- ANI-1（2017 年）は、記述子を改良し、約 1720 万個の DFT の構造で学習して、H、C、N、O からなる有機分子に広く使えるポテンシャルを作りました。
- SchNet（2017〜2018 年）は、距離の連続関数のフィルターを使った更新（メッセージパッシング）で、原子環境の表し方そのものを学習しました。
- 3 つとも、モデルの中で回転で変わらない量だけを使っています。この限界が、第 5 回の等変なモデルにつながります。

## 参考文献

- J. Behler, M. Parrinello, "Generalized Neural-Network Representation of High-Dimensional Potential-Energy Surfaces", Phys. Rev. Lett. 98, 146401 (2007). DOI: [10.1103/PhysRevLett.98.146401](https://doi.org/10.1103/PhysRevLett.98.146401)
- J. Behler, "Atom-centered symmetry functions for constructing high-dimensional neural network potentials", J. Chem. Phys. 134, 074106 (2011). DOI: [10.1063/1.3553717](https://doi.org/10.1063/1.3553717)
- J. S. Smith, O. Isayev, A. E. Roitberg, "ANI-1: an extensible neural network potential with DFT accuracy at force field computational cost", Chem. Sci. 8, 3192 (2017). DOI: [10.1039/C6SC05720A](https://doi.org/10.1039/C6SC05720A)
- J. S. Smith, B. T. Nebgen, R. Zubatyuk, N. Lubbers, C. Devereux, K. Barros, S. Tretiak, O. Isayev, A. E. Roitberg, "Approaching coupled cluster accuracy with a general-purpose neural network potential through transfer learning", Nat. Commun. 10, 2903 (2019). DOI: [10.1038/s41467-019-10827-4](https://doi.org/10.1038/s41467-019-10827-4)
- K. T. Schütt, F. Arbabzadah, S. Chmiela, K.-R. Müller, A. Tkatchenko, "Quantum-chemical insights from deep tensor neural networks", Nat. Commun. 8, 13890 (2017). DOI: [10.1038/ncomms13890](https://doi.org/10.1038/ncomms13890)
- K. T. Schütt, P.-J. Kindermans, H. E. Sauceda, S. Chmiela, A. Tkatchenko, K.-R. Müller, "SchNet: A continuous-filter convolutional neural network for modeling quantum interactions", arXiv:1706.08566 (2017).
- K. T. Schütt, H. E. Sauceda, P.-J. Kindermans, A. Tkatchenko, K.-R. Müller, "SchNet – A deep learning architecture for molecules and materials", J. Chem. Phys. 148, 241722 (2018). DOI: [10.1063/1.5019779](https://doi.org/10.1063/1.5019779)
- J. Gilmer, S. S. Schoenholz, P. F. Riley, O. Vinyals, G. E. Dahl, "Neural Message Passing for Quantum Chemistry", arXiv:1704.01212 (2017).
- O. T. Unke, S. Chmiela, H. E. Sauceda, M. Gastegger, I. Poltavsky, K. T. Schütt, A. Tkatchenko, K.-R. Müller, "Machine Learning Force Fields", Chem. Rev. 121, 10142 (2021). DOI: [10.1021/acs.chemrev.0c01111](https://doi.org/10.1021/acs.chemrev.0c01111)
- S. Batzner, A. Musaelian, L. Sun, M. Geiger, J. P. Mailoa, M. Kornbluth, N. Molinari, T. E. Smidt, B. Kozinsky, "E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials", Nat. Commun. 13, 2453 (2022). DOI: [10.1038/s41467-022-29939-5](https://doi.org/10.1038/s41467-022-29939-5)
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, G. Csányi, "MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields", arXiv:2206.07697 (2022).
- S. N. Pozdnyakov, M. J. Willatt, A. P. Bartók, C. Ortner, G. Csányi, M. Ceriotti, "Incompleteness of Atomic Structure Representations", Phys. Rev. Lett. 125, 166001 (2020). DOI: [10.1103/PhysRevLett.125.166001](https://doi.org/10.1103/PhysRevLett.125.166001)
