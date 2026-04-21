# Naive Gauss-Jordan Solver

線形代数演算のための、軽量な純粋Python実装によるガウス・ジョルダン消去法ライブラリです。このライブラリは、教育目的、および**厳密な有理数演算**が必要なシナリオ向けに設計されています。

## ✨ 特徴

- **厳密な演算サポート**: Pythonの `fractions.Fraction` オブジェクトを使用することで、浮動小数点の丸め誤差を回避し、数学的に正確な結果を得ることができます。
- **連立一次方程式の解法**: 連立一次方程式 ($Ax = b$) を解きます。
- **逆行列の計算**: 正方行列の逆行列を計算します。
- **行列演算**: 
  - 行列の積 (Dot product)
  - 行列の転置
  - 行列式の計算
  - 行列のランク (Rank) 計算
  - ムーア・ペンローズ擬似逆行列
- **依存関係なし**: Pythonの標準ライブラリのみを使用して構築されています。`numpy` は不要です。
- **整形表示**: 行列を読みやすい形式に整形して表示するための組み込みユーティリティを備えています。

## 🚀 インストール

このライブラリは単一のファイルで構成されているため、`gaussjordan.py` をダウンロードしてプロジェクトにインポートするだけで使用できます。

```bash
# ファイルをダウンロード
curl -O https://raw.githubusercontent.com/kitanokitsune/naive_gauss_jordan/main/gaussjordan.py
```

## 🛠 使い方

### 1. 連立一次方程式を解く
以下を解きます：
$1x + 2y + 5z = 12$
$1x + 3y + 1z = 10$
$2x + 3y + 1z = 14$

```python
from gaussjordan import solve

M = [[1, 2, 5], 
     [1, 3, 1], 
     [2, 3, 1]]
V = [12, 10, 14]

ans = solve(M, V)
print(f"解: {ans}") 
# 出力: [4.0, 1.6923076923076916, 0.9230769230769231]
```

### 2. 分数を用いた厳密な演算
浮動小数点の誤差を避けるには、処理の前に行列を `Fraction` オブジェクトに変換してください。

```python
from gaussjordan import solve, toFraction

M = [[1, 2, 5], [1, 3, 1], [2, 3, 1]]
V = [12, 10, 14]

# 厳密な有理数に変換
Mf = toFraction(M)
Vf = toFraction(V)

# 厳密に解く
ans_exact = solve(Mf, Vf)

print(f"厳密な解: {ans_exact}")
# 出力: [Fraction(4, 1), Fraction(22, 13), Fraction(12, 13)]
```

### 3. 逆行列と整形表示
```python
from gaussjordan import invert, toFraction, pprint_mat

M = toFraction([ 
     [1, 2], 
     [3, 4]])

inv_M = invert(M)
pprint_mat("逆行列=", inv_M)
# 出力: 
# 逆行列={  -2,    1 |
#        | 3/2, -1/2 }
```

### 4. ムーア・ペンローズ擬似逆行列
```python
from gaussjordan import rank, moore_penrose, toFraction, pprint_mat

M = toFraction([
     [1, 2, 3, 1],
     [2, 4, 6, 4],
     [0, 0, 0, 2]])

r = rank(M)
print(f"ランク: {r}")
# 出力: 2

pinv_M = moore_penrose(M)
print("ムーア・ペンローズ擬似逆行列:")
pprint_mat("", pinv_M)
# 出力: 
# { 1/28, 1/56, -3/56 |
# | 1/14, 1/28, -3/28 |
# | 3/28, 3/56, -9/56 |
# | -1/6, 1/12,  5/12 }
```

## 📖 APIリファレンス

| 関数 | 説明 |
| :--- | :--- |
| `solve(A, b)` | 連立一次方程式 $Ax = b$ を解きます。 |
| `invert(mat)` | 正方行列の逆行列を返します。 |
| `dot(A, B)` | 行列の積（ドット積）を行います。 |
| `det(A)` | 行列式を計算します。 |
| `rank(mat)` | 行列のランクを返します。 |
| `moore_penrose(mat)` | ムーア・ペンローズ擬似逆行列を計算します。 |
| `transpose(A)` | 行列の転置を返します。 |
| `toFraction(M)` | 行列の要素を再帰的に `fractions.Fraction` に変換します。 |
| `toReal(M)` | 行列の要素を再帰的に `float` に変換します。 |
| `pprint_mat(prefix="", mat=[])` | 行列を整形して読みやすい形式で表示します。 |
| `spprint_mat(prefix="", mat=[])` | `pprint_mat()` の出力を文字列として返します。 |

## ⚠️ 注意事項

この実装は「**ナイーブ（素朴な）**」ソルバーです。学習や厳密な有理数演算には非常に適していますが、`NumPy` のような極めて大きな行列の処理や、ハイパフォーマンス・コンピューティング向けに最適化されているわけではありません。

## 📜 ライセンス

このプロジェクトは [MIT License](LICENSE) の下で公開されています。

