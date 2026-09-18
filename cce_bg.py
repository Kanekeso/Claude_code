"""CCE（定組成膨張試験）の結果から、任意の圧力における Bg を計算する。

CCE (Constant Composition Expansion) レポートには通常、圧力ごとに

    * ガス偏差係数 Z（single-phase Z / two-phase Z）
    * 相対体積 V/Vsat（relative volume）
    * ガス容積係数 Bg そのもの

のいずれかが表として与えられる。本モジュールはそれらを読み込み、
**Z を圧力の関数として単調保存の 3 次補間（PCHIP 相当）**してから

    Bg(p) = (p_sc / T_sc) * z(p) * T / p

で任意圧力の Bg を求める。

なぜ Bg ではなく Z を補間するのか
---------------------------------
Bg は 1/p に比例して双曲線的に変化するため、表の点の間を Bg で直接
線形補間すると系統的な誤差が出る（特に低圧側で大きい）。一方 Z は
圧力に対してなだらかに変化するので、Z を補間してから解析式で Bg に
戻すほうが精度が高い。本モジュールはこの方針を採る。

単位系
------
`UNIT_SYSTEMS` に主要な単位系を用意している。標準状態（p_sc, T_sc）は
`UnitSystem.with_standard()` で差し替えられる。

使用例
------
>>> cce = GasCCE(
...     pressures=[5000, 4000, 3000, 2000, 1000],
...     z_factors=[0.985, 0.920, 0.885, 0.890, 0.930],
...     reservoir_temperature=250.0,   # degF
...     units="field",                 # psia / degF / ft3/scf
... )
>>> round(cce.bg_at(3500), 6)
0.005142
"""

from __future__ import annotations

import argparse
import bisect
import csv
import dataclasses
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

__all__ = [
    "UnitSystem",
    "UNIT_SYSTEMS",
    "MonotoneCubic",
    "BgRow",
    "GasCCE",
    "main",
]


# --------------------------------------------------------------------------
# 単位系
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class UnitSystem:
    """圧力・温度の単位と標準状態の定義。

    Attributes
    ----------
    p_std, t_std:
        標準状態の圧力・温度（この単位系の単位で表したもの）。
    absolute_zero:
        この単位系の温度目盛りにおける絶対零度（degF なら -459.67）。
    """

    name: str
    pressure_unit: str
    temperature_unit: str
    bg_unit: str
    p_std: float
    t_std: float
    absolute_zero: float

    def absolute_temperature(self, temperature: float) -> float:
        """温度を絶対温度（degR / K）に変換する。"""
        t_abs = temperature - self.absolute_zero
        if t_abs <= 0.0:
            raise ValueError(
                f"絶対温度が正になりません: {temperature} {self.temperature_unit}"
            )
        return t_abs

    @property
    def t_std_absolute(self) -> float:
        return self.absolute_temperature(self.t_std)

    @property
    def bg_coefficient(self) -> float:
        """Bg = coefficient * z * T_abs / p の係数 (p_sc / T_sc)。"""
        return self.p_std / self.t_std_absolute

    def with_standard(
        self, p_std: float | None = None, t_std: float | None = None
    ) -> "UnitSystem":
        """標準状態だけを差し替えた新しい単位系を返す。"""
        return dataclasses.replace(
            self,
            name=f"{self.name}*",
            p_std=self.p_std if p_std is None else p_std,
            t_std=self.t_std if t_std is None else t_std,
        )


#: よく使う単位系。キーは大文字小文字を区別しない（`GasCCE` 側で lower される）。
UNIT_SYSTEMS: dict[str, UnitSystem] = {
    # 米国石油業界の慣用単位。標準状態 14.696 psia / 60 degF
    "field": UnitSystem(
        name="field",
        pressure_unit="psia",
        temperature_unit="degF",
        bg_unit="ft3/scf",
        p_std=14.696,
        t_std=60.0,
        absolute_zero=-459.67,
    ),
    # SI。標準状態 101.325 kPa / 15 degC
    "metric": UnitSystem(
        name="metric",
        pressure_unit="kPa(a)",
        temperature_unit="degC",
        bg_unit="m3/sm3",
        p_std=101.325,
        t_std=15.0,
        absolute_zero=-273.15,
    ),
    "bar": UnitSystem(
        name="bar",
        pressure_unit="bara",
        temperature_unit="degC",
        bg_unit="m3/sm3",
        p_std=1.01325,
        t_std=15.0,
        absolute_zero=-273.15,
    ),
    "mpa": UnitSystem(
        name="mpa",
        pressure_unit="MPa(a)",
        temperature_unit="degC",
        bg_unit="m3/sm3",
        p_std=0.101325,
        t_std=15.0,
        absolute_zero=-273.15,
    ),
}


def resolve_units(units: "str | UnitSystem") -> UnitSystem:
    """単位系の名前または `UnitSystem` を `UnitSystem` に解決する。"""
    if isinstance(units, UnitSystem):
        return units
    try:
        return UNIT_SYSTEMS[units.strip().lower()]
    except KeyError:
        raise ValueError(
            f"未知の単位系 {units!r}。使用可能: {', '.join(sorted(UNIT_SYSTEMS))}"
        ) from None


# --------------------------------------------------------------------------
# 単調保存 3 次補間（Fritsch-Carlson / PCHIP）
# --------------------------------------------------------------------------
def _end_slope(h0: float, h1: float, d0: float, d1: float) -> float:
    """PCHIP の端点勾配（片側 3 点公式 + 単調性のクリップ）。"""
    slope = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if slope * d0 <= 0.0:
        return 0.0
    if d0 * d1 < 0.0 and abs(slope) > abs(3.0 * d0):
        return 3.0 * d0
    return slope


class MonotoneCubic:
    """単調性を保存する 3 次 Hermite 補間（scipy.PchipInterpolator 相当）。

    表の点の間でオーバーシュート（物理的にありえない Z の振動）を起こさない
    ため、PVT 表の補間に適している。scipy に依存しない純 Python 実装。

    Parameters
    ----------
    x, y:
        データ点。x は重複不可（順序は自動で昇順ソートされる）。
    extrapolation:
        ``"error"``（既定）／``"linear"``（端の勾配で直線外挿）／
        ``"clamp"``（端の値で一定）。
    """

    def __init__(
        self,
        x: Sequence[float],
        y: Sequence[float],
        extrapolation: str = "error",
    ) -> None:
        if len(x) != len(y):
            raise ValueError("x と y の長さが一致しません")
        if len(x) < 2:
            raise ValueError("補間には最低 2 点必要です")
        if extrapolation not in ("error", "linear", "clamp"):
            raise ValueError(
                f"extrapolation は 'error' / 'linear' / 'clamp' のいずれか: {extrapolation!r}"
            )

        pairs = sorted(zip(x, y))
        self.x = [float(a) for a, _ in pairs]
        self.y = [float(b) for _, b in pairs]
        self.extrapolation = extrapolation

        for i in range(len(self.x) - 1):
            if self.x[i + 1] == self.x[i]:
                raise ValueError(f"x に重複があります: {self.x[i]}")

        self._h = [self.x[i + 1] - self.x[i] for i in range(len(self.x) - 1)]
        self._delta = [
            (self.y[i + 1] - self.y[i]) / self._h[i] for i in range(len(self.x) - 1)
        ]
        self._d = self._slopes()

    def _slopes(self) -> list[float]:
        n = len(self.x)
        h, delta = self._h, self._delta
        if n == 2:
            return [delta[0], delta[0]]

        d = [0.0] * n
        for i in range(1, n - 1):
            if delta[i - 1] * delta[i] <= 0.0:
                d[i] = 0.0  # 極値点では勾配 0（オーバーシュート防止）
            else:
                w1 = 2.0 * h[i] + h[i - 1]
                w2 = h[i] + 2.0 * h[i - 1]
                d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
        d[0] = _end_slope(h[0], h[1], delta[0], delta[1])
        d[-1] = _end_slope(h[-1], h[-2], delta[-1], delta[-2])
        return d

    # ---------------------------------------------------------------- 評価
    def _locate(self, xv: float) -> int:
        i = bisect.bisect_right(self.x, xv) - 1
        return min(max(i, 0), len(self.x) - 2)

    def _outside(self, xv: float) -> int:
        """範囲外なら -1（左）／+1（右）、範囲内なら 0。"""
        if xv < self.x[0]:
            return -1
        if xv > self.x[-1]:
            return 1
        return 0

    def __call__(self, xv: float) -> float:
        side = self._outside(xv)
        if side:
            return self._extrapolate(xv, side, derivative=False)
        i = self._locate(xv)
        h = self._h[i]
        t = (xv - self.x[i]) / h
        t2, t3 = t * t, t * t * t
        return (
            (2.0 * t3 - 3.0 * t2 + 1.0) * self.y[i]
            + (t3 - 2.0 * t2 + t) * h * self._d[i]
            + (-2.0 * t3 + 3.0 * t2) * self.y[i + 1]
            + (t3 - t2) * h * self._d[i + 1]
        )

    def derivative(self, xv: float) -> float:
        """dy/dx を返す。"""
        side = self._outside(xv)
        if side:
            return self._extrapolate(xv, side, derivative=True)
        i = self._locate(xv)
        h = self._h[i]
        t = (xv - self.x[i]) / h
        t2 = t * t
        return (
            (6.0 * t2 - 6.0 * t) / h * self.y[i]
            + (3.0 * t2 - 4.0 * t + 1.0) * self._d[i]
            + (-6.0 * t2 + 6.0 * t) / h * self.y[i + 1]
            + (3.0 * t2 - 2.0 * t) * self._d[i + 1]
        )

    def _extrapolate(self, xv: float, side: int, derivative: bool) -> float:
        if self.extrapolation == "error":
            raise ValueError(
                f"{xv:g} はデータ範囲 [{self.x[0]:g}, {self.x[-1]:g}] の外です。"
                " 外挿する場合は extrapolation='linear' または 'clamp' を指定してください。"
            )
        idx = 0 if side < 0 else -1
        slope = 0.0 if self.extrapolation == "clamp" else self._d[idx]
        if derivative:
            return slope
        return self.y[idx] + slope * (xv - self.x[idx])


# --------------------------------------------------------------------------
# 計算結果の 1 行
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class BgRow:
    """ある圧力における計算結果。"""

    pressure: float
    z: float
    bg: float
    eg: float          # 1/Bg（ガス膨張係数 scf/ft3 など）
    cg: float          # 等温圧縮率 1/p - (1/z)(dz/dp)
    extrapolated: bool = False
    two_phase: bool = False

    def as_dict(self) -> dict[str, float | bool]:
        return dataclasses.asdict(self)


# --------------------------------------------------------------------------
# 本体
# --------------------------------------------------------------------------
class GasCCE:
    """CCE 表を保持し、任意圧力の Z / Bg を返すオブジェクト。

    Parameters
    ----------
    pressures:
        CCE の圧力点（昇順・降順どちらでもよい）。
    z_factors:
        各圧力における偏差係数 Z。
    reservoir_temperature:
        CCE の試験温度（= 貯留層温度）。単位系の温度単位で与える。
    units:
        ``"field"`` / ``"metric"`` / ``"bar"`` / ``"mpa"`` または `UnitSystem`。
    dew_point:
        露点圧力。与えると、それ未満の圧力の結果に ``two_phase=True`` の
        フラグが立つ（CCE の Z が two-phase Z の領域であることの注意喚起）。
    extrapolation:
        データ範囲外の扱い。``"error"`` / ``"linear"`` / ``"clamp"``。
    """

    def __init__(
        self,
        pressures: Sequence[float],
        z_factors: Sequence[float],
        reservoir_temperature: float,
        units: "str | UnitSystem" = "field",
        dew_point: float | None = None,
        extrapolation: str = "error",
        name: str = "",
    ) -> None:
        self.units = resolve_units(units)
        self.temperature = float(reservoir_temperature)
        self.t_absolute = self.units.absolute_temperature(self.temperature)
        self.dew_point = None if dew_point is None else float(dew_point)
        self.name = name

        if len(pressures) != len(z_factors):
            raise ValueError("pressures と z_factors の長さが一致しません")
        if len(pressures) < 2:
            raise ValueError("CCE 表には最低 2 点必要です")
        for p, z in zip(pressures, z_factors):
            if p <= 0.0:
                raise ValueError(f"圧力は正である必要があります: {p}")
            if z <= 0.0:
                raise ValueError(f"Z は正である必要があります: {z}")

        self._interp = MonotoneCubic(pressures, z_factors, extrapolation=extrapolation)

    # ------------------------------------------------------------- 別入力形式
    @classmethod
    def from_relative_volume(
        cls,
        pressures: Sequence[float],
        relative_volumes: Sequence[float],
        reservoir_temperature: float,
        z_ref: float,
        p_ref: float | None = None,
        **kwargs,
    ) -> "GasCCE":
        """相対体積 V/Vsat から Z を復元して構築する。

        V = z n R T / p より、等温の CCE セル内では

            z(p) = z_ref * (p * V) / (p_ref * V_ref)

        が成り立つ（V は相対体積でよい。比なので基準は相殺する）。

        Parameters
        ----------
        z_ref, p_ref:
            基準点の Z と圧力。`p_ref` 省略時は相対体積が 1.0 に最も近い点
            （通常は露点＝飽和圧）を基準にする。
        """
        if len(pressures) != len(relative_volumes):
            raise ValueError("pressures と relative_volumes の長さが一致しません")
        if z_ref <= 0.0:
            raise ValueError("z_ref は正である必要があります")

        if p_ref is None:
            idx = min(
                range(len(relative_volumes)),
                key=lambda i: abs(relative_volumes[i] - 1.0),
            )
        else:
            try:
                idx = [float(p) for p in pressures].index(float(p_ref))
            except ValueError:
                raise ValueError(
                    f"p_ref={p_ref} が CCE 表の圧力点に見つかりません"
                ) from None

        p_ref_v = float(pressures[idx])
        v_ref = float(relative_volumes[idx])
        if v_ref <= 0.0:
            raise ValueError("基準点の相対体積は正である必要があります")

        z = [
            z_ref * (float(p) * float(v)) / (p_ref_v * v_ref)
            for p, v in zip(pressures, relative_volumes)
        ]
        return cls(pressures, z, reservoir_temperature, **kwargs)

    @classmethod
    def from_bg(
        cls,
        pressures: Sequence[float],
        bg_values: Sequence[float],
        reservoir_temperature: float,
        units: "str | UnitSystem" = "field",
        **kwargs,
    ) -> "GasCCE":
        """すでに Bg が与えられている表から構築する（内部では Z に戻す）。"""
        u = resolve_units(units)
        t_abs = u.absolute_temperature(float(reservoir_temperature))
        z = [
            float(bg) * float(p) / (u.bg_coefficient * t_abs)
            for p, bg in zip(pressures, bg_values)
        ]
        return cls(pressures, z, reservoir_temperature, units=u, **kwargs)

    @classmethod
    def from_csv(
        cls,
        path: "str | Path",
        reservoir_temperature: float,
        units: "str | UnitSystem" = "field",
        pressure_column: str | None = None,
        value_column: str | None = None,
        value_kind: str | None = None,
        z_ref: float | None = None,
        p_ref: float | None = None,
        **kwargs,
    ) -> "GasCCE":
        """CSV から構築する。

        ヘッダー名から列を自動判定する（大文字小文字・空白・記号は無視）。
        認識する列:

        ==========  ==================================================
        圧力        p, pressure, press, 圧力
        Z           z, z-factor, zfactor, gas z, deviation factor
        相対体積    vrel, v/vsat, relative volume, relvol, 相対体積
        Bg          bg, gas fvf, gas formation volume factor
        ==========  ==================================================

        `value_column` と `value_kind`（``"z"`` / ``"vrel"`` / ``"bg"``）を
        明示すれば自動判定を上書きできる。相対体積の場合は `z_ref` が必須。
        """
        pressures, values, kind = _read_cce_csv(
            path, pressure_column, value_column, value_kind
        )
        if kind == "z":
            return cls(pressures, values, reservoir_temperature, units=units, **kwargs)
        if kind == "bg":
            return cls.from_bg(
                pressures, values, reservoir_temperature, units=units, **kwargs
            )
        if z_ref is None:
            raise ValueError(
                "相対体積の列しかありません。基準点の Z を z_ref で指定してください。"
            )
        return cls.from_relative_volume(
            pressures,
            values,
            reservoir_temperature,
            z_ref=z_ref,
            p_ref=p_ref,
            units=units,
            **kwargs,
        )

    # ------------------------------------------------------------------ 計算
    @property
    def pressures(self) -> list[float]:
        """内部に保持している圧力点（昇順）。"""
        return list(self._interp.x)

    @property
    def z_factors(self) -> list[float]:
        return list(self._interp.y)

    @property
    def pressure_range(self) -> tuple[float, float]:
        return self._interp.x[0], self._interp.x[-1]

    def z_at(self, pressure: float) -> float:
        """任意圧力の Z（単調保存 3 次補間）。"""
        return self._interp(float(pressure))

    def dz_dp(self, pressure: float) -> float:
        """dZ/dp。"""
        return self._interp.derivative(float(pressure))

    def bg_at(self, pressure: float) -> float:
        """任意圧力のガス容積係数 Bg [貯留層体積 / 標準体積]。

        Bg = (p_sc / T_sc) * z T / p
        """
        p = float(pressure)
        if p <= 0.0:
            raise ValueError(f"圧力は正である必要があります: {pressure}")
        return self.units.bg_coefficient * self.z_at(p) * self.t_absolute / p

    def eg_at(self, pressure: float) -> float:
        """ガス膨張係数 Eg = 1/Bg。"""
        return 1.0 / self.bg_at(pressure)

    def cg_at(self, pressure: float) -> float:
        """等温ガス圧縮率 cg = 1/p - (1/z)(dz/dp)。"""
        p = float(pressure)
        return 1.0 / p - self.dz_dp(p) / self.z_at(p)

    def is_extrapolated(self, pressure: float) -> bool:
        lo, hi = self.pressure_range
        return not (lo <= float(pressure) <= hi)

    def is_two_phase(self, pressure: float) -> bool:
        """露点未満か（露点が与えられていない場合は常に False）。"""
        return self.dew_point is not None and float(pressure) < self.dew_point

    def row_at(self, pressure: float) -> BgRow:
        p = float(pressure)
        return BgRow(
            pressure=p,
            z=self.z_at(p),
            bg=self.bg_at(p),
            eg=self.eg_at(p),
            cg=self.cg_at(p),
            extrapolated=self.is_extrapolated(p),
            two_phase=self.is_two_phase(p),
        )

    def table(self, pressures: Iterable[float]) -> list[BgRow]:
        """複数圧力をまとめて計算する。"""
        return [self.row_at(p) for p in pressures]

    # ------------------------------------------------------------------ 出力
    def format_table(self, pressures: Iterable[float]) -> str:
        """人間が読める表に整形する。"""
        u = self.units
        rows = self.table(pressures)
        header = (
            f"{'p [' + u.pressure_unit + ']':>16}"
            f"{'Z [-]':>10}"
            f"{'Bg [' + u.bg_unit + ']':>20}"
            f"{'Eg [1/Bg]':>14}"
            f"{'cg [1/' + u.pressure_unit + ']':>18}"
            "  note"
        )
        lines = [header, "-" * (len(header) + 4)]
        for r in rows:
            notes = []
            if r.extrapolated:
                notes.append("外挿")
            if r.two_phase:
                notes.append("露点未満(2相Z)")
            lines.append(
                f"{r.pressure:>16.6g}"
                f"{r.z:>10.4f}"
                f"{r.bg:>20.6g}"
                f"{r.eg:>14.6g}"
                f"{r.cg:>18.6g}"
                f"  {' '.join(notes)}"
            )
        return "\n".join(lines)

    def write_csv(self, path: "str | Path", pressures: Iterable[float]) -> None:
        """計算結果を CSV に書き出す。"""
        u = self.units
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    f"pressure [{u.pressure_unit}]",
                    "z [-]",
                    f"Bg [{u.bg_unit}]",
                    f"Eg [{u.bg_unit}^-1]",
                    f"cg [1/{u.pressure_unit}]",
                    "extrapolated",
                    "two_phase",
                ]
            )
            for r in self.table(pressures):
                writer.writerow(
                    [
                        f"{r.pressure:.6g}",
                        f"{r.z:.5f}",
                        f"{r.bg:.6g}",
                        f"{r.eg:.6g}",
                        f"{r.cg:.6g}",
                        int(r.extrapolated),
                        int(r.two_phase),
                    ]
                )

    def __repr__(self) -> str:  # pragma: no cover - デバッグ用
        lo, hi = self.pressure_range
        return (
            f"GasCCE(name={self.name!r}, n={len(self._interp.x)}, "
            f"p=[{lo:g}, {hi:g}] {self.units.pressure_unit}, "
            f"T={self.temperature:g} {self.units.temperature_unit})"
        )


# --------------------------------------------------------------------------
# CSV 読み込み
# --------------------------------------------------------------------------
_PRESSURE_KEYS = {"p", "pressure", "press", "pres", "圧力"}
_Z_KEYS = {
    "z", "zfactor", "zfactors", "gasz", "zgas", "twophasez", "singlephasez",
    "deviationfactor", "gasdeviationfactor", "compressibilityfactor",
}
_VREL_KEYS = {
    "vrel", "vvsat", "vvd", "relativevolume", "relvol", "relativevol",
    "relvolume", "相対体積",
}
_BG_KEYS = {"bg", "gasfvf", "fvf", "gasformationvolumefactor"}


def _normalize(key: str) -> str:
    """ヘッダー名を照合用に正規化する。

    単位の注記（``Pressure (psia)``、``p [psia]``、``Z, -``）は切り落とし、
    残りから英数字と非 ASCII 文字だけを取り出す。
    """
    head = key
    for sep in ("(", "[", "{", ",", ":"):
        head = head.split(sep, 1)[0]
    return "".join(ch for ch in head.lower() if ch.isalnum() or ord(ch) > 127)


def _match_column(normalized: dict[str, str], keys: set[str]) -> str | None:
    """完全一致を優先し、なければ前方一致でヘッダー名を探す。"""
    for norm, original in normalized.items():
        if norm in keys:
            return original
    for norm, original in normalized.items():
        if any(norm.startswith(k) for k in keys if len(k) > 1):
            return original
    return None


def _read_cce_csv(
    path: "str | Path",
    pressure_column: str | None,
    value_column: str | None,
    value_kind: str | None,
) -> tuple[list[float], list[float], str]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh) if any(v for v in r.values())]
    if not rows:
        raise ValueError(f"{path} にデータ行がありません")

    fields = [f for f in rows[0].keys() if f is not None]
    normalized = {_normalize(f): f for f in fields}

    if pressure_column is None:
        p_col = _match_column(normalized, _PRESSURE_KEYS)
        if p_col is None:
            raise ValueError(f"圧力の列が見つかりません。列: {fields}")
    else:
        p_col = pressure_column

    if value_column is not None:
        if value_kind is None:
            raise ValueError("value_column を指定する場合は value_kind も必要です")
        v_col, kind = value_column, value_kind.lower()
    else:
        for keys, kind_name in (
            (_Z_KEYS, "z"),
            (_BG_KEYS, "bg"),
            (_VREL_KEYS, "vrel"),
        ):
            col = _match_column(normalized, keys)
            if col is not None:
                v_col, kind = col, kind_name
                break
        else:
            raise ValueError(
                f"Z / Bg / 相対体積 のいずれの列も見つかりません。列: {fields}"
            )

    if kind not in ("z", "bg", "vrel"):
        raise ValueError(f"value_kind は 'z' / 'bg' / 'vrel': {kind!r}")

    pressures: list[float] = []
    values: list[float] = []
    for i, row in enumerate(rows, start=2):
        raw_p, raw_v = row.get(p_col), row.get(v_col)
        if raw_p in (None, "") or raw_v in (None, ""):
            continue  # 空欄（露点未満の single-phase Z など）は読み飛ばす
        try:
            pressures.append(float(str(raw_p).replace(",", "")))
            values.append(float(str(raw_v).replace(",", "")))
        except ValueError:
            raise ValueError(f"{path} の {i} 行目を数値に変換できません: {raw_p}, {raw_v}") from None

    if len(pressures) < 2:
        raise ValueError(f"{path} から有効なデータ点を 2 点以上読み取れませんでした")
    return pressures, values, kind


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def _parse_pressure_args(
    values: list[float] | None, rng: list[float] | None
) -> list[float]:
    pressures: list[float] = list(values or [])
    if rng:
        start, stop, step = rng
        if step == 0:
            raise ValueError("--range の step に 0 は指定できません")
        n = int(math.floor((stop - start) / step + 1e-9)) + 1
        if n <= 0:
            raise ValueError("--range の start / stop / step の向きが矛盾しています")
        pressures.extend(start + step * i for i in range(n))
    if not pressures:
        raise ValueError("--pressure か --range のどちらかを指定してください")
    return pressures


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="cce_bg",
        description="CCE 試験結果から任意圧力の Bg を計算する",
    )
    parser.add_argument("csv", help="CCE 表の CSV（p と Z / 相対体積 / Bg の列）")
    parser.add_argument(
        "-T", "--temperature", type=float, required=True, help="貯留層温度（CCE 試験温度）"
    )
    parser.add_argument(
        "-u", "--units", default="field",
        choices=sorted(UNIT_SYSTEMS), help="単位系（既定: field）",
    )
    parser.add_argument(
        "-p", "--pressure", type=float, action="append", help="計算する圧力（複数指定可）"
    )
    parser.add_argument(
        "--range", type=float, nargs=3, metavar=("START", "STOP", "STEP"),
        help="等間隔の圧力列（例: --range 5000 1000 -500）",
    )
    parser.add_argument("--dew-point", type=float, help="露点圧力（2 相領域のフラグ用）")
    parser.add_argument(
        "--extrapolation", default="error", choices=["error", "linear", "clamp"],
        help="データ範囲外の扱い（既定: error）",
    )
    parser.add_argument("--z-ref", type=float, help="相対体積入力のときの基準 Z")
    parser.add_argument("--p-ref", type=float, help="相対体積入力のときの基準圧力")
    parser.add_argument("--p-std", type=float, help="標準圧力を上書き")
    parser.add_argument("--t-std", type=float, help="標準温度を上書き")
    parser.add_argument("-o", "--output", help="結果を書き出す CSV のパス")
    args = parser.parse_args(argv)

    try:
        units = resolve_units(args.units)
        if args.p_std is not None or args.t_std is not None:
            units = units.with_standard(p_std=args.p_std, t_std=args.t_std)

        cce = GasCCE.from_csv(
            args.csv,
            reservoir_temperature=args.temperature,
            units=units,
            z_ref=args.z_ref,
            p_ref=args.p_ref,
            dew_point=args.dew_point,
            extrapolation=args.extrapolation,
        )
        pressures = _parse_pressure_args(args.pressure, args.range)
    except (ValueError, OSError) as exc:
        print(f"エラー: {exc}", file=sys.stderr)
        return 2

    try:
        print(cce.format_table(pressures))
        if args.output:
            cce.write_csv(args.output, pressures)
            print(f"\n{args.output} に書き出しました。")
    except ValueError as exc:
        print(f"エラー: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
