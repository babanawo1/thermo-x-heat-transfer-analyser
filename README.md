# 🔥 THERMO-X | Heat Transfer Analyser

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/)
[![React 18+](https://img.shields.io/badge/React-18.3.1-61DAFB.svg?style=flat&logo=react&logoColor=black)](https://reactjs.org/)
[![TypeScript](https://img.shields.io/badge/TypeScript-5.5+-3178C6.svg?style=flat&logo=typescript&logoColor=white)](https://www.typescriptlang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**THERMO-X** is an engineering-grade thermal analysis and simulation tool designed for mechanical engineers, researchers, educators, and students. It performs deterministic, validated calculations across the five fundamental domains of heat transfer, paired with interactive visualisations, step-by-step mathematical breakdowns, and exportable engineering reports in PDF and CSV formats.

All calculations strictly use standard **SI units** and adhere to rigorous thermodynamic boundary constraints.

---

## 🌟 Key Features

### 1. 🧱 Conduction Analysis (Fourier's Law)
- **1D Steady-State Conduction**: Plane wall heat transfer $Q = \frac{k A (T_{hot} - T_{cold})}{L}$.
- **Material Preset Library**: Copper, Aluminum, Structural Steel, Glass, Concrete, Brick, Wood, Fiberglass, and Aerogel.
- **Dynamic Gradient Profile**: Interactive $T(x)$ temperature profile through the solid layer.

### 2. 💨 Convection Analysis (Newton's Law of Cooling)
- **Fluid Boundary Exchange**: $Q = h A (T_{surface} - T_{\infty})$.
- **Convective Regimes**: Free convection in air, forced convection in air, forced convection in water, and boiling/condensation heat transfer.
- **Cooling vs. Heating Detection**: Automatically detects heat flow direction and computes convective boundary resistance $R_{conv} = \frac{1}{hA}$.

### 3. ☀️ Radiation Analysis (Stefan-Boltzmann Law)
- **Net Radiative Exchange**: $Q_{rad} = \varepsilon \sigma A (T_{surface}^4 - T_{surroundings}^4)$ using NIST standard $\sigma = 5.670374419 \times 10^{-8} \text{ W/m}^2\cdot\text{K}^4$.
- **Surface Presets**: Polished Silver, Polished Aluminum, Oxidized Copper, White Paint, Black Matte Paint, and Ideal Blackbody ($\varepsilon = 1.0$).
- **Non-Linear Characteristic Curve**: $T^4$ temperature sweep chart illustrating radiative power acceleration.

### 4. 🔄 Heat Exchanger Analysis (LMTD & $\varepsilon$-NTU)
- **Flow Configurations**: Counter-Current and Parallel Flow.
- **Singularity-Protected LMTD**: Numerical singularity handling when $\Delta T_1 \approx \Delta T_2$ via analytical limit evaluation:
  $$\lim_{\Delta T_1 \to \Delta T_2} \text{LMTD} = \Delta T_1$$
- **First Law Energy Balances**: Hot stream duty $Q_h = \dot{m}_h C_{p,h} (T_{h,in} - T_{h,out})$ vs. Cold stream duty $Q_c = \dot{m}_c C_{p,c} (T_{c,out} - T_{c,in})$ with percentage discrepancy tracking.
- **Effectiveness-NTU Rating**: Computes $C_{min}$, $C_{max}$, capacity ratio $C_r$, maximum possible heat transfer $Q_{max}$, thermal effectiveness $\varepsilon$, and NTU.
- **Second Law Violation Warnings**: Traps temperature crossovers before computing invalid thermodynamic states.

### 5. 🧱 Thermal Resistance Networks
- **Series Composite Walls**: Multi-layer combined solid conduction and fluid boundary convection.
- **Interface Temperature Computation**: Determines intermediate plane temperatures $T_i$ across all layers.
- **Resistance Distribution**: Identifies dominant insulating resistances and overall heat-transfer coefficient $U$.

### 6. 📄 Professional Reports & Exporting
- **PDF Reports**: Comprehensive engineering calculation sheets with governing formulas, input parameters, result tables, governing assumptions, and engineering interpretations.
- **CSV Data Logs**: Complete machine-readable data logs for spreadsheet integration and audit trails.

---

## 📐 Governing Equations Summary

| Domain | Governing Equation | Primary SI Unit |
| :--- | :--- | :--- |
| **Conduction** | $Q = \frac{k A \Delta T}{L} = \frac{\Delta T}{R_{cond}}$ | $\text{W}$ |
| **Convection** | $Q = h A (T_s - T_\infty) = \frac{\Delta T}{R_{conv}}$ | $\text{W}$ |
| **Radiation** | $Q_{rad} = \varepsilon \sigma A (T_s^4 - T_{surr}^4)$ | $\text{W}$ |
| **Heat Exchanger** | $Q = U A \text{LMTD}, \quad \text{LMTD} = \frac{\Delta T_1 - \Delta T_2}{\ln(\Delta T_1 / \Delta T_2)}$ | $\text{W}$ |
| **Resistance Network** | $R_{total} = \sum R_i, \quad Q = \frac{\Delta T_{total}}{R_{total}}$ | $\text{K/W}, \text{W}$ |

---

## 🚀 Quickstart & Local Installation

### Prerequisites
- Python 3.9+ (for the Streamlit application)
- Node.js 18+ & npm / bun (for the React / TypeScript web client)

---

### Option A: Running the Python / Streamlit App

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/thermo-x-heat-transfer-analyser.git
   cd thermo-x-heat-transfer-analyser
## Live Application

[Launch THERMO-X](https://thermo-x-heat-transfer-analyser-3zttxzncjadtkusswxtxhs.streamlit.app/)

## GitHub Repository

[View the source code](https://github.com/babanawo1/thermo-x-heat-transfer-analyser)
