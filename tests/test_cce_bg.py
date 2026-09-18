"""cce_bg のテスト（標準ライブラリの unittest のみ使用）。

実行:  python -m unittest discover -s tests -v
"""

import csv
import doctest
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cce_bg import (  # noqa: E402
    UNIT_SYSTEMS,
    GasCCE,
    MonotoneCubic,
    main,
    resolve_units,
)

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

P = [5000, 4600, 4200, 3800, 3400, 3000, 2600, 2200, 1800, 1400, 1000]
Z = [0.9880, 0.9560, 0.9330, 0.9140, 0.8990, 0.8890, 0.8840, 0.8850, 0.8920, 0.9060, 0.9260]
T_F = 250.0


def field_cce(**kwargs) -> GasCCE:
    return GasCCE(P, Z, T_F, units="field", **kwargs)


class TestBgFormula(unittest.TestCase):
    def test_matches_hand_calculation(self):
        """Bg = 0.0282793 * z T / p（field 単位）と一致すること。"""
        cce = field_cce()
        p = 3000.0
        expected = 14.696 / (60.0 + 459.67) * 0.8890 * (250.0 + 459.67) / p
        self.assertAlmostEqual(cce.bg_at(p), expected, places=10)
        # 教科書の係数 0.02827 [ft3/scf] とも一致
        self.assertAlmostEqual(UNIT_SYSTEMS["field"].bg_coefficient, 0.02827, delta=1e-5)

    def test_eg_is_reciprocal_of_bg(self):
        cce = field_cce()
        self.assertAlmostEqual(cce.eg_at(2750.0), 1.0 / cce.bg_at(2750.0), places=12)

    def test_ideal_gas_compressibility(self):
        """Z が一定なら cg = 1/p になること。"""
        cce = GasCCE([1000, 3000, 5000], [1.0, 1.0, 1.0], T_F, units="field")
        for p in (1500.0, 2500.0, 4000.0):
            self.assertAlmostEqual(cce.cg_at(p), 1.0 / p, places=12)

    def test_bg_decreases_with_pressure(self):
        cce = field_cce()
        bg = [cce.bg_at(p) for p in range(1000, 5001, 250)]
        self.assertTrue(all(a > b for a, b in zip(bg, bg[1:], strict=False)))

    def test_unit_systems_agree(self):
        """field と metric で同じ物理状態なら Bg が 0.3% 以内で一致すること。

        標準状態が 60 degF/14.696 psia と 15 degC/101.325 kPa で異なるため
        完全一致はしない。
        """
        p_psia, t_f, z = 3000.0, 250.0, 0.889
        p_kpa = p_psia * 6.894757
        t_c = (t_f - 32.0) / 1.8

        field = GasCCE([p_psia * 0.5, p_psia * 1.5], [z, z], t_f, units="field")
        metric = GasCCE([p_kpa * 0.5, p_kpa * 1.5], [z, z], t_c, units="metric")
        rel = abs(field.bg_at(p_psia) - metric.bg_at(p_kpa)) / field.bg_at(p_psia)
        self.assertLess(rel, 0.003)

    def test_custom_standard_conditions(self):
        """標準状態を差し替えると Bg が比例して変わること。"""
        units = UNIT_SYSTEMS["field"].with_standard(p_std=14.73)
        base = field_cce()
        custom = GasCCE(P, Z, T_F, units=units)
        self.assertAlmostEqual(custom.bg_at(3000.0) / base.bg_at(3000.0), 14.73 / 14.696, places=10)


class TestInterpolation(unittest.TestCase):
    def test_passes_through_data_points(self):
        cce = field_cce()
        for p, z in zip(P, Z, strict=True):
            self.assertAlmostEqual(cce.z_at(p), z, places=12)

    def test_no_overshoot_between_points(self):
        """単調保存補間なので、隣り合う 2 点の値の外に出ないこと。"""
        cce = field_cce()
        pairs = sorted(zip(P, Z, strict=True))
        for (p0, z0), (p1, z1) in zip(pairs, pairs[1:], strict=False):
            lo, hi = min(z0, z1), max(z0, z1)
            for k in range(1, 20):
                z = cce.z_at(p0 + (p1 - p0) * k / 20.0)
                self.assertGreaterEqual(z, lo - 1e-12)
                self.assertLessEqual(z, hi + 1e-12)

    def test_derivative_matches_finite_difference(self):
        cce = field_cce()
        p, h = 3100.0, 1e-3
        fd = (cce.z_at(p + h) - cce.z_at(p - h)) / (2.0 * h)
        self.assertAlmostEqual(cce.dz_dp(p), fd, places=8)

    def test_linear_interpolation_with_two_points(self):
        interp = MonotoneCubic([1000.0, 2000.0], [0.9, 0.8])
        self.assertAlmostEqual(interp(1500.0), 0.85, places=12)

    def test_extrapolation_error_by_default(self):
        cce = field_cce()
        with self.assertRaises(ValueError):
            cce.bg_at(6000.0)
        with self.assertRaises(ValueError):
            cce.bg_at(500.0)

    def test_extrapolation_modes(self):
        clamped = field_cce(extrapolation="clamp")
        self.assertAlmostEqual(clamped.z_at(6000.0), Z[0], places=12)
        linear = field_cce(extrapolation="linear")
        self.assertNotAlmostEqual(linear.z_at(6000.0), Z[0], places=6)
        self.assertTrue(linear.is_extrapolated(6000.0))
        self.assertFalse(linear.is_extrapolated(3000.0))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            GasCCE([1000], [0.9], T_F)
        with self.assertRaises(ValueError):
            GasCCE([1000, 2000], [0.9], T_F)
        with self.assertRaises(ValueError):
            GasCCE([1000, -2000], [0.9, 0.8], T_F)
        with self.assertRaises(ValueError):
            GasCCE([1000, 2000], [0.9, 0.0], T_F)
        with self.assertRaises(ValueError):
            resolve_units("imperial")


class TestAlternativeInputs(unittest.TestCase):
    def test_relative_volume_round_trip(self):
        """相対体積から Z を復元できること。"""
        p_ref, z_ref = 4200.0, 0.9330
        vrel = [(z / p) / (z_ref / p_ref) for p, z in zip(P, Z, strict=True)]
        cce = GasCCE.from_relative_volume(P, vrel, T_F, z_ref=z_ref, p_ref=p_ref, units="field")
        for p, z in zip(P, Z, strict=True):
            self.assertAlmostEqual(cce.z_at(p), z, places=10)

    def test_relative_volume_default_reference(self):
        """p_ref 省略時は V/Vsat = 1 の点が基準になること。"""
        p_ref, z_ref = 4200.0, 0.9330
        vrel = [(z / p) / (z_ref / p_ref) for p, z in zip(P, Z, strict=True)]
        cce = GasCCE.from_relative_volume(P, vrel, T_F, z_ref=z_ref, units="field")
        self.assertAlmostEqual(cce.z_at(3000.0), field_cce().z_at(3000.0), places=10)

    def test_bg_round_trip(self):
        base = field_cce()
        bg = [base.bg_at(p) for p in P]
        cce = GasCCE.from_bg(P, bg, T_F, units="field")
        for p, z in zip(P, Z, strict=True):
            self.assertAlmostEqual(cce.z_at(p), z, places=10)
            self.assertAlmostEqual(cce.bg_at(p), base.bg_at(p), places=12)

    def test_relative_volume_requires_valid_reference(self):
        vrel = [1.0] * len(P)
        with self.assertRaises(ValueError):
            GasCCE.from_relative_volume(P, vrel, T_F, z_ref=0.9, p_ref=1234.0)


class TestCsv(unittest.TestCase):
    def test_reads_z_column(self):
        cce = GasCCE.from_csv(EXAMPLES / "cce_gas_field.csv", T_F, units="field")
        self.assertEqual(len(cce.pressures), len(P))
        self.assertAlmostEqual(cce.z_at(3000.0), 0.8890, places=6)

    def test_reads_relative_volume_column(self):
        cce = GasCCE.from_csv(
            EXAMPLES / "cce_relative_volume_only.csv",
            T_F,
            units="field",
            z_ref=0.9330,
            p_ref=4200.0,
        )
        self.assertAlmostEqual(cce.z_at(3000.0), 0.8890, places=3)

    def test_relative_volume_without_z_ref_raises(self):
        with self.assertRaises(ValueError):
            GasCCE.from_csv(EXAMPLES / "cce_relative_volume_only.csv", T_F)

    def test_skips_blank_rows_and_thousands_separators(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cce.csv"
            path.write_text(
                'Pressure,Z\n"5,000",0.988\n\n3000,\n"1,000",0.926\n',
                encoding="utf-8",
            )
            cce = GasCCE.from_csv(path, T_F, units="field")
            self.assertEqual(cce.pressures, [1000.0, 5000.0])

    def test_write_csv(self):
        cce = field_cce(dew_point=4200.0)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "bg.csv"
            cce.write_csv(out, [4400.0, 3000.0])
            with out.open(encoding="utf-8") as fh:
                rows = list(csv.reader(fh))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[1][-1], "0")  # 4400 psia は露点以上
            self.assertEqual(rows[2][-1], "1")  # 3000 psia は露点未満

    def test_missing_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("foo,bar\n1,2\n3,4\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                GasCCE.from_csv(path, T_F)


class TestRowsAndCli(unittest.TestCase):
    def test_two_phase_flag(self):
        cce = field_cce(dew_point=4200.0)
        self.assertFalse(cce.row_at(4500.0).two_phase)
        self.assertTrue(cce.row_at(3000.0).two_phase)

    def test_table(self):
        rows = field_cce().table([4000.0, 3000.0, 2000.0])
        self.assertEqual([r.pressure for r in rows], [4000.0, 3000.0, 2000.0])
        self.assertTrue(all(r.bg > 0 for r in rows))

    def test_cli_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.csv"
            code = main(
                [
                    str(EXAMPLES / "cce_gas_field.csv"),
                    "-T",
                    "250",
                    "-u",
                    "field",
                    "--range",
                    "5000",
                    "1000",
                    "-500",
                    "--dew-point",
                    "4200",
                    "-o",
                    str(out),
                ]
            )
            self.assertEqual(code, 0)
            with out.open(encoding="utf-8") as fh:
                rows = list(csv.reader(fh))
            self.assertEqual(len(rows), 10)  # ヘッダー + 9 点

    def test_cli_rejects_out_of_range_pressure(self):
        code = main(
            [
                str(EXAMPLES / "cce_gas_field.csv"),
                "-T",
                "250",
                "-p",
                "9000",
            ]
        )
        self.assertEqual(code, 2)


def load_tests(loader, tests, ignore):
    """モジュール docstring の使用例（doctest）も一緒に実行する。"""
    import cce_bg

    tests.addTests(doctest.DocTestSuite(cce_bg))
    return tests


if __name__ == "__main__":
    unittest.main()
