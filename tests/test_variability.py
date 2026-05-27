from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")

from vgsot_sim.analysis.variability import (
    PDKBudgetInputs,
    cv_delta_budget,
    sample_macrospin_process_constants,
)
from vgsot_sim.configs import PhysicalConstantsConfig


def test_cv_delta_budget_components_partition_total_variance():
    budget = cv_delta_budget(PDKBudgetInputs())

    assert sum(budget.components.values()) == pytest.approx(1.0)
    assert budget.components["V_mag_area"] < budget.CV_V ** 2 / budget.CV_Delta ** 2


def test_macrospin_process_sample_is_reproducible_and_maps_scales():
    base = PhysicalConstantsConfig()
    sample_a = sample_macrospin_process_constants(
        base, PDKBudgetInputs(), rng=np.random.default_rng(123)
    )
    sample_b = sample_macrospin_process_constants(
        base, PDKBudgetInputs(), rng=np.random.default_rng(123)
    )

    assert sample_a.scales == pytest.approx(sample_b.scales)
    assert sample_a.constants.Ki == pytest.approx(base.Ki * sample_a.scales["Ki_proxy"])
    assert sample_a.constants.rho == pytest.approx(base.rho * sample_a.scales["Rsot"])
    assert sample_a.constants.D_elec == pytest.approx(base.D_elec * sample_a.scales["D"])
    assert sample_a.constants.R_W > 0.0
