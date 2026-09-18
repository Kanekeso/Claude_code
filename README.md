# CCE から任意圧力の Bg を計算する

CCE（定組成膨張試験, Constant Composition Expansion）の測定結果から、
**表に載っていない任意の圧力**におけるガス容積係数 Bg を求めるための
Python モジュールです。標準ライブラリのみで動作します
（numpy / scipy / pandas 不要）。

```
Bg(p) = (p_sc / T_sc) · z(p) · T / p
```

## 考え方：Bg ではなく Z を補間する

Bg は `1/p` に比例して双曲線的に変化するため、CCE 表の点と点の間を
**Bg のまま線形補間すると系統誤差が出ます**（低圧側ほど大きい）。
一方 Z は圧力に対してなだらかに変化するので、

1. CCE 表から `Z vs p` を作る（相対体積や Bg しかない場合は Z に換算する）
2. `Z(p)` を**単調保存 3 次補間**（PCHIP 相当。点の間で振動しない）で求める
3. 解析式で Bg に戻す

という順序を採っています。同じ補間の微分から `cg = 1/p − (1/z)(dz/dp)` も
得られるので、等温圧縮率も一緒に返します。

## 入力形式

CCE レポートの体裁に合わせて 3 通りの入力に対応しています。

| 手元にあるもの | 使うコンストラクタ |
| --- | --- |
| 圧力 と Z | `GasCCE(...)` |
| 圧力 と 相対体積 V/Vsat | `GasCCE.from_relative_volume(..., z_ref=...)` |
| 圧力 と Bg | `GasCCE.from_bg(...)` |
| 上記いずれかの CSV | `GasCCE.from_csv(...)` |

相対体積からは `z(p) = z_ref · (p·V) / (p_ref·V_ref)` で Z を復元します
（等温セル内で `V = zNRT/p` が成り立つため）。基準点 `z_ref` は
レポートの露点における Z を使うのが普通です。`p_ref` を省略すると
`V/Vsat = 1.0` に最も近い点（＝飽和圧）が自動的に基準になります。

## 使い方（Python）

```python
from cce_bg import GasCCE

cce = GasCCE(
    pressures=[5000, 4600, 4200, 3800, 3400, 3000, 2600, 2200, 1800, 1400, 1000],
    z_factors=[0.988, 0.956, 0.933, 0.914, 0.899, 0.889, 0.884, 0.885, 0.892, 0.906, 0.926],
    reservoir_temperature=250.0,   # degF
    units="field",                 # psia / degF -> Bg [ft3/scf]
    dew_point=4200.0,              # 任意。露点未満の結果にフラグが付く
)

cce.z_at(3250)      # 任意圧力の Z
cce.bg_at(3250)     # 任意圧力の Bg  [ft3/scf]
cce.eg_at(3250)     # 1/Bg
cce.cg_at(3250)     # 等温圧縮率 [1/psia]

print(cce.format_table(range(5000, 999, -500)))
cce.write_csv("bg_table.csv", range(5000, 999, -500))
```

CSV から読む場合（ヘッダー名は自動判定。`Pressure (psia)` のように
単位が付いていても認識します）：

```python
cce = GasCCE.from_csv("examples/cce_gas_field.csv", reservoir_temperature=250.0,
                      units="field", dew_point=4200.0)
```

## 使い方（コマンドライン）

```bash
# 指定した圧力だけ
python cce_bg.py examples/cce_gas_field.csv -T 250 -u field -p 3250 -p 2750

# 等間隔の圧力列 + CSV 出力
python cce_bg.py examples/cce_gas_field.csv -T 250 --range 5000 1000 -500 \
    --dew-point 4200 -o bg_table.csv

# 相対体積しかない CSV（基準の Z が必要）
python cce_bg.py examples/cce_relative_volume_only.csv -T 250 \
    --z-ref 0.9330 --p-ref 4200 --range 5000 1000 -500
```

出力例：

```
        p [psia]     Z [-]        Bg [ft3/scf]     Eg [1/Bg]       cg [1/psia]  note
----------------------------------------------------------------------------------------
            5000    0.9880          0.00396565       252.165       0.000107642
            4500    0.9496          0.00423511       236.122       0.000158012
            4000    0.9230          0.00463092        215.94       0.000198249  露点未満(2相Z)
```

## 単位系

`-u / units` で選びます。標準状態は `--p-std` / `--t-std`、
Python からは `UnitSystem.with_standard()` で上書きできます。

| 名前 | 圧力 | 温度 | Bg | 標準状態 |
| --- | --- | --- | --- | --- |
| `field` | psia | degF | ft³/scf | 14.696 psia / 60 °F |
| `metric` | kPa(a) | degC | m³/sm³ | 101.325 kPa / 15 °C |
| `bar` | bara | degC | m³/sm³ | 1.01325 bara / 15 °C |
| `mpa` | MPa(a) | degC | m³/sm³ | 0.101325 MPa / 15 °C |

field と metric で Bg の値が 0.2 % ほどずれるのは、標準状態の温度が
60 °F と 15 °C で異なるためです（バグではありません）。

## 注意点

- **データ範囲外**は既定でエラーになります。外挿が必要なときは
  `extrapolation="linear"`（端の勾配で直線外挿）または `"clamp"`
  （端の値で一定）を指定してください。結果には `extrapolated=True` の
  フラグが付きます。
- **露点未満**の CCE の Z は通常 *two-phase Z* です。そこから計算した値は
  ガス単相の Bg ではなく 2 相の容積係数になるため、`dew_point` を渡すと
  該当する行に警告フラグが付きます。露点未満の乾きガス Bg が必要な場合は
  CVD（定容枯渇試験）の single-phase Z を使ってください。
- CSV の空欄行は読み飛ばします（露点未満で single-phase Z が空欄の
  レポートをそのまま読めます）。

## 環境構築（VS Code / ローカル）

`cce_bg.py` 本体は標準ライブラリだけで動くので、**何もインストールしなくても実行できます**。
numpy などを併用したい場合や、テスト・リンターを動かす場合だけ以下を行ってください。

### 1. リポジトリを取得する

```bash
git clone https://github.com/Kanekeso/Claude_code.git
cd Claude_code
git checkout claude/cce-bg-pressure-calculation-ef9xp3
```

### 2. 仮想環境を作って有効化する

プロジェクトごとに `.venv` を作るのが VS Code の標準的な進め方です。

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate

# Windows (PowerShell)
py -m venv .venv
.venv\Scripts\Activate.ps1
```

VS Code では `Ctrl+Shift+P` →「Python: Select Interpreter」で `.venv` を選びます
（`.vscode/settings.json` で既定の解釈系を `.venv` に向けてあるので、通常は自動で選ばれます）。

### 3. パッケージをインストールする

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt   # numpy / pandas / matplotlib / pytest / ruff
```

numpy だけで十分なら `python -m pip install numpy` でも構いません。
新しくパッケージを追加したら `requirements.txt` にも書き足しておくと、
他の環境やこのリポジトリの Web セッションでも同じ構成を再現できます。

> **注意**: `pip install` は必ず有効化した `.venv` の中で実行してください。
> `python -m pip ...` の形で呼ぶと、いま選んでいる Python に確実に入ります。

### 4. 動作確認

```bash
python -m pytest              # テスト（28 件）
python -m ruff check .        # リンター
python -m ruff format .       # フォーマッタ
```

### Claude Code on the web で使う場合

`.claude/hooks/session-start.sh` を用意してあるので、Web セッションの開始時に
`requirements-dev.txt` が自動でインストールされます（ローカル実行時は何もしません）。
このフックが効くのは、**変更をリポジトリの既定ブランチにマージした後**の
セッションからです。

## 変更をコミットする

```bash
git add -A
git commit -m "変更内容の要約"
git push -u origin claude/cce-bg-pressure-calculation-ef9xp3
```

`.gitignore` で `__pycache__/` と `.venv/` は除外済みなので、
`git add -A` しても仮想環境そのものはコミットされません。

## テスト

```bash
python -m pytest                        # pytest を入れている場合
python -m unittest discover -s tests    # 追加インストールなしで実行する場合
```

Bg の式の手計算照合、補間がデータ点を通ること、オーバーシュートしないこと、
相対体積 / Bg からの往復変換、単位系の整合、CSV 読み込み、CLI、docstring の
使用例（doctest）までを 28 件で検証しています。
