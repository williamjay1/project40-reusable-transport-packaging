"""Small deterministic semantic tests; no data download or repository writes."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
import run_lifecycle_v38 as producer
import lifecycle_v38_path_verify as independent


def jobs(starts, owners=None, delay=1, lost=()):
    return [dict(id=i, owner=0 if owners is None else owners[i], start=t,
                 ready=t + delay, lost=i in lost) for i, t in enumerate(starts)]


class LifecycleSemantics(unittest.TestCase):
    def minimum(self, js, K, L, mode="pooled", initial=(), endpoint="ready"):
        return independent.assignment_optimum(independent.job_key(js), K, L,
                                               mode, initial, endpoint)[0]

    def test_ready_expiry_is_inclusive(self):
        js = jobs([0], delay=3)
        self.assertEqual(self.minimum(js, 1, 3), 1)
        self.assertIsNone(self.minimum(js, 1, 2))

    def test_delivery_endpoint_not_ready_endpoint(self):
        js = jobs([0], delay=4)
        self.assertIsNone(self.minimum(js, 1, 3))
        self.assertEqual(self.minimum(js, 1, 3, endpoint="start"), 1)

    def test_loss_is_absorbing(self):
        self.assertEqual(self.minimum(jobs([0, 1], lost=(0,)), 4, 10), 2)
        self.assertEqual(self.minimum(jobs([0, 1]), 4, 10), 1)

    def test_ownership_sharing_has_boundary(self):
        js = jobs([0, 1], [0, 1])
        self.assertEqual(self.minimum(js, 4, 10, "dedicated"), 2)
        self.assertEqual(self.minimum(js, 4, 10), 1)
        self.assertEqual(self.minimum(js, 1, 10, "dedicated"), 2)
        self.assertEqual(self.minimum(js, 1, 10), 2)

    def test_old_age_counts_waiting(self):
        old = ((0, 0, Fraction(3), 0),)
        self.assertEqual(self.minimum(jobs([4], delay=2), 3, 8, initial=old), 1)
        self.assertEqual(self.minimum(jobs([3], delay=2), 3, 8, initial=old), 0)

    def test_previous_services_consume_quota(self):
        old = ((0, 0, Fraction(0), 2),)
        self.assertEqual(self.minimum(jobs([0, 1]), 3, 8, initial=old), 1)

    def test_explicit_old_origin_sentinel(self):
        checked = independent.semantic_probes(producer.path_masks)
        self.assertTrue(checked["loss_is_absorbing"])
        self.assertFalse(checked["owner_none_helper_probe"]["affects_any_of_578"])

    def test_solver_dp_and_independent_assignment_agree(self):
        js = jobs([0, 1, 4, 5], [0, 0, 1, 1])
        for mode in ("dedicated", "pooled"):
            result = producer.solve_case("public_smoke", [producer.Job(**j) for j in js],
                                         4, 101, mode, "none")
            self.assertEqual(result["introductions"], self.minimum(js, 4, 101, mode))
            self.assertTrue(result["certificate"]["all_path_constraints_exactly_verified"])
            self.assertEqual(result["service_coverage"], [1] * 4)
            independent.simulate_witness(result)


if __name__ == "__main__":
    unittest.main()
