#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geometric Phase-Coincidence Mechanism of Radioactive Decay on M_4 = S^3 x S^1:
Quadrupole S^2 Gate Pinning in Superheavy alpha-Emitters and 229Thm Nuclear Clocks.
Author: Evgeniy A. Tarasenko (2026)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

hbar_c = 197.3269804   # MeV * fm
e2_4pi = 1.439964548   # MeV * fm
m_u    = 931.4941024   # MeV / c^2
c_fm_s = 2.99792458e23 # fm / s
r_0    = 1.215         # fm


def gamow_half_life(Z_parent, A_parent, Q_alpha_MeV, beta2, c_def=0.0):
    Z_d = Z_parent - 2
    A_d = A_parent - 4
    mu_MeV = (4.001506 * A_d / (4.001506 + A_d)) * m_u

    R_0 = r_0 * (A_d**(1.0 / 3.0) + 4.0**(1.0 / 3.0))
    R_b = (2.0 * Z_d * e2_4pi) / Q_alpha_MeV
    x = min(0.999, R_0 / R_b)

    G = (
        np.sqrt(2.0 * mu_MeV / (hbar_c**2))
        * (2.0 * Z_d * e2_4pi / np.sqrt(Q_alpha_MeV))
        * (np.arccos(np.sqrt(x)) - np.sqrt(x * (1.0 - x)))
    )

    v_alpha = np.sqrt(2.0 * Q_alpha_MeV / mu_MeV) * c_fm_s
    nu_0 = v_alpha / (2.0 * R_0)

    mu_0_gate = 0.025
    F_geom = np.exp(c_def * beta2 * np.sqrt(Z_d))

    lam = nu_0 * mu_0_gate * F_geom * np.exp(-2.0 * G)
    return float(np.log(2.0) / lam)


def main():
    train_nuclei = (
        ("212Po (sph.)",  84, 212, 8.954,  0.000, 2.99e-7),
        ("226Ra (def.)",  88, 226, 4.871,  0.151, 5.049e10),
        ("238U (act.)",   92, 238, 4.270,  0.215, 1.41e17),
        ("240Pu (act.)",  94, 240, 5.256,  0.223, 2.07e11),
    )

    blind_superheavy = (
        ("270Ds (Z=110)", 110, 270, 11.117, 0.221, 1.00e-4),
        ("285Cn (Z=112)", 112, 285,  9.320, 0.140, 2.90e1),
        ("288Fl (Z=114)", 114, 288, 10.076, 0.110, 6.60e-1),
        ("292Lv (Z=116)", 116, 292, 10.791, 0.080, 1.80e-2),
        ("294Og (Z=118)", 118, 294, 11.810, 0.050, 6.90e-4),
    )

    def loss_fn(c_d):
        err = 0.0
        for _, Z, A, Q, b2, t_exp in train_nuclei:
            t_mod = gamow_half_life(Z, A, Q, b2, c_def=c_d)
            err += (np.log10(t_mod) - np.log10(t_exp))**2
        return err

    res = minimize_scalar(loss_fn, bounds=(-2.0, 2.0), method="bounded")
    c_def_opt = float(res.x)

    print(f"Calibrated universal quadrupole pinning constant: c_def = {c_def_opt:.4f}
")
    print("=== 1. CALIBRATION SET (NUBASE2020) ===")
    for name, Z, A, Q, b2, t_exp in train_nuclei:
        t_mod = gamow_half_life(Z, A, Q, b2, c_def=c_def_opt)
        print(f"{name:16s} | Q_a={Q:6.3f} MeV | beta2={b2:5.3f} | T1/2_mod={t_mod:9.2e} s | T1/2_exp={t_exp:9.2e} s")

    print("
=== 2. BLIND PREDICTIONS FOR SUPERHEAVY NUCLEI (JINR DUBNA, Z=110-118) ===")
    sq_err_wkb = []
    sq_err_s2 = []
    for name, Z, A, Q, b2, t_exp in blind_superheavy:
        t_wkb = gamow_half_life(Z, A, Q, b2, c_def=0.0)
        t_s2 = gamow_half_life(Z, A, Q, b2, c_def=c_def_opt)
        sq_err_wkb.append((np.log10(t_wkb) - np.log10(t_exp))**2)
        sq_err_s2.append((np.log10(t_s2) - np.log10(t_exp))**2)
        print(f"{name:16s} | WKB={t_wkb:9.2e} s | S2_Gate={t_s2:9.2e} s | Exp={t_exp:9.2e} s")

    print(f"
Blind RMS log10 error (Z=110-118): WKB = {np.sqrt(np.mean(sq_err_wkb)):.3f} dec -> S2 Gate = {np.sqrt(np.mean(sq_err_s2)):.3f} dec")

    # 229Thm nuclear clock EFG pinning shift
    E_isom_eV = 8.3557
    Q_quad_m2 = 3.11e-28
    V_zz_caf2 = 1.5e22
    V_zz_trap = 1.0e14

    dE_caf2 = 0.25 * Q_quad_m2 * V_zz_caf2
    dnz2_caf2 = (4.0 / 45.0) * (dE_caf2 / E_isom_eV)
    dlam_caf2 = 1.5 * dnz2_caf2

    dE_trap = 0.25 * Q_quad_m2 * V_zz_trap
    dnz2_trap = (4.0 / 45.0) * (dE_trap / E_isom_eV)
    dlam_trap = 1.5 * dnz2_trap

    print("
=== 3. 229Thm NUCLEAR CLOCK RADIATIVE M1 DECAY SHIFT ===")
    print(f"CaF2 axial vs cubic site: dE_quad = {dE_caf2*1e6:.2f} ueV | d<nz^2> = {dnz2_caf2:.2e} | d(lambda)/lambda = {dlam_caf2:.2e}")
    print(f"Isolated 229Th3+ ion trap: d<nz^2> = {dnz2_trap:.2e} | d(lambda)/lambda = {dlam_trap:.2e}")

    # Plot Figure 1
    all_nuclei = train_nuclei + blind_superheavy
    labels = [n for (n, _, _, _, _, _) in all_nuclei]
    log_exp = [np.log10(t) for (_, _, _, _, _, t) in all_nuclei]
    log_mod = [np.log10(gamow_half_life(Z, A, Q, b2, c_def=c_def_opt)) for (_, Z, A, Q, b2, _) in all_nuclei]
    log_raw = [np.log10(gamow_half_life(Z, A, Q, b2, c_def=0.0)) for (_, Z, A, Q, b2, _) in all_nuclei]

    fig, ax = plt.subplots(figsize=(7.6, 4.2), dpi=300)
    x_idx = np.arange(len(labels))
    ax.plot(x_idx, log_exp, "ko-", lw=2.0, markersize=6, label="Experimental Half-Life (NUBASE / JINR)")
    ax.plot(x_idx, log_mod, "rs--", lw=1.8, markersize=5.5, label=rf"Model with $S^2$ Gate Pinning ($c_{{m def}} = {c_def_opt:.4f}$)")
    ax.plot(x_idx, log_raw, "b^:", lw=1.2, markersize=5, alpha=0.75, label=r"Isotropic WKB without deformation pinning ($c_{m def} = 0$)")

    ax.axvline(3.5, color="gray", ls="-.", lw=1.2)
    ax.text(1.5, 14, "Calibration
(4 nuclei)", ha="center", fontsize=8.5, color="#333333")
    ax.text(6.0, 14, "Blind Prediction
(Superheavy Z = 110-118)", ha="center", fontsize=8.5, color="#900c3f", fontweight="bold")

    ax.set_xticks(x_idx)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel(r"$\log_{10}(T_{1/2} / \mathrm{s})$", fontsize=10)
    ax.set_title(r"$lpha$-Decay Half-Lives Across 24 Orders of Magnitude via $S^2$ Gate Pinning", fontsize=10.5)
    ax.legend(fontsize=8, loc="lower left", framealpha=0.95)
    ax.grid(True, ls="--", alpha=0.35)

    plt.tight_layout()
    plt.savefig("fig_alpha_decay_s2_gate.png", dpi=300)
    plt.close()
    print("Saved figure: fig_alpha_decay_s2_gate.png")


if __name__ == "__main__":
    main()
