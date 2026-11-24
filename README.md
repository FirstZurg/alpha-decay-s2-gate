# Geometric Phase-Coincidence Nuclear Decay Model (`alpha_decay_s2_gate_sim.py`)

Reproducible Python implementation accompanying the manuscript:  
**"Geometric phase-coincidence mechanism of radioactive decay on $M_4 = S^3 \times S^1$: Quadrupole $S^2$ gate pinning in superheavy $\alpha$-emitters and $^{229}\mathrm{Th}^m$ nuclear clocks"**  
*Author:* Evgeniy A. Tarasenko (Independent Researcher, Chita, Russia)

---

## Overview

This self-contained script reproduces all numerical calibrations, blind predictions, and figures presented in the study:

1. **Universal Constant Calibration:** Calibrates the single dimensionless quadrupole gate-pinning parameter $c_{\mathrm{def}} = 1.5097$ against four reference nuclei ($^{212}\mathrm{Po}$, $^{226}\mathrm{Ra}$, $^{238}\mathrm{U}$, and $^{240}\mathrm{Pu}$) spanning 24 orders of magnitude in half-life.
2. **Blind Superheavy Half-Life Predictions ($Z = 110\text{--}118$):** Evaluates $\alpha$-decay half-lives for five superheavy isotopes synthesized at JINR Dubna ($^{270}\mathrm{Ds}$, $^{285}\mathrm{Cn}$, $^{288}\mathrm{Fl}$, $^{292}\mathrm{Lv}$, and $^{294}\mathrm{Og}$) with $c_{\mathrm{def}}$ frozen, reducing the RMS logarithmic error from $0.525$ to $0.487$ decades.
3. **$^{229}\mathrm{Th}^m$ Nuclear Clock EFG Shift:** Computes the quadrupole normal-vector pinning shift $\delta\lambda/\lambda = 1.86 \times 10^{-8}$ between axial $\mathrm{F}^-$-compensated sites ($V_{zz} \approx 1.5 \times 10^{22}~\text{V/m}^2$) and cubic zero-EFG sites in solid-state $\mathrm{CaF}_2$ hosts, as well as the shift in isolated vacuum ion traps ($1.24 \times 10^{-16}$).

---

## Requirements & Installation

- Python 3.8+
- `numpy`
- `scipy`
- `matplotlib`

Install dependencies via `pip`:
```bash
pip install numpy scipy matplotlib
```

---

## Usage

Run the script directly from the terminal:
```bash
python alpha_decay_s2_gate_sim.py
```

### Expected Output
- Prints the calibrated coupling constant $c_{\mathrm{def}} = 1.5097$, calibration table, blind superheavy $\alpha$-decay predictions, and $^{229}\mathrm{Th}^m$ radiative rate shifts to standard output.
- Generates the Geiger--Nuttall comparison plot `fig_alpha_decay_s2_gate.png` (300 dpi) in the working directory.

---

## License

This code is released under the **MIT License**.
