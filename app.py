"""
THERMO-X | Heat Transfer Analyser (SI PRO)
Complete standalone Streamlit Application

AI DEVELOPMENT DOCUMENTATION:
AI tools used:
- Google AI Studio Build (Gemini 3.7 Flash & Antigravity Agent Engine)

Key prompts used:
1. Build a Streamlit engineering application for heat-transfer analysis with interactive inputs, charts, Pandas results tables, validation, GitHub documentation and Streamlit deployment.
2. Implement and validate conduction, convection, radiation, heat exchanger and thermal resistance calculations using SI units.
3. Improve the application's UI, error handling, engineering validation, downloadable CSV/PDF reports, GitHub structure and deployment readiness.

Manual verification/fix:
The engineering equations, units, numerical results, heat-exchanger energy balances, physical input constraints, singularity edge cases (e.g. LMTD ΔT1 ≈ ΔT2), and Stefan-Boltzmann T⁴ radiation limits were manually reviewed, cross-tested against analytical benchmark test cases, and verified for production-grade engineering deployment.
"""

import math
import io
from datetime import datetime
from typing import Dict, Any, List

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

# ReportLab imports for PDF reporting
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle

# -------------------------------------------------------------
# CONSTANTS & THERMODYNAMIC PARAMETERS
# -------------------------------------------------------------
STEFAN_BOLTZMANN_CONSTANT = 5.670374419e-8  # W/m²·K⁴


# -------------------------------------------------------------
# 1. CORE PHYSICS CALCULATION ENGINES
# -------------------------------------------------------------
def calculate_conduction(k: float, area: float, length: float, t_hot: float, t_cold: float) -> Dict[str, Any]:
    """1D Steady-State Heat Conduction (Fourier's Law)"""
    warnings = []
    if area <= 0:
        warnings.append("Area must be strictly greater than 0 m².")
    if length <= 0:
        warnings.append("Wall thickness (L) must be strictly greater than 0 m.")
    if k <= 0:
        warnings.append("Thermal conductivity (k) must be strictly greater than 0 W/m·K.")
    if t_hot <= 0 or t_cold <= 0:
        warnings.append("Temperatures must be strictly positive on absolute Kelvin scale (T > 0 K).")

    delta_t = t_hot - t_cold
    if delta_t == 0:
        warnings.append("Isothermal boundary condition: Heat transfer rate is zero.")

    is_valid = (len(warnings) == 0 or (len(warnings) == 1 and delta_t == 0)) and (area > 0 and length > 0 and k > 0)
    if not is_valid:
        return {"delta_t": delta_t, "r_cond": 0.0, "q": 0.0, "q_flux": 0.0, "is_valid": False, "warnings": warnings}

    r_cond = length / (k * area)
    q = delta_t / r_cond
    q_flux = q / area

    return {"delta_t": delta_t, "r_cond": r_cond, "q": q, "q_flux": q_flux, "is_valid": True, "warnings": warnings}


def calculate_convection(h: float, area: float, t_surface: float, t_fluid: float) -> Dict[str, Any]:
    """Newton's Law of Cooling"""
    warnings = []
    if h <= 0:
        warnings.append("Convection coefficient (h) must be strictly greater than 0 W/m²·K.")
    if area <= 0:
        warnings.append("Surface area must be strictly greater than 0 m².")
    if t_surface <= 0 or t_fluid <= 0:
        warnings.append("Temperatures must be strictly positive on absolute Kelvin scale (T > 0 K).")

    delta_t = t_surface - t_fluid
    mode = "equilibrium"
    if delta_t > 0:
        mode = "cooling"
    elif delta_t < 0:
        mode = "heating"
    else:
        warnings.append("Thermal equilibrium (T_surface = T_fluid): Heat transfer rate is zero.")

    is_valid = (len(warnings) == 0 or (len(warnings) == 1 and delta_t == 0)) and (h > 0 and area > 0)
    if not is_valid:
        return {"delta_t": delta_t, "r_conv": 0.0, "q": 0.0, "q_flux": 0.0, "mode": mode, "is_valid": False, "warnings": warnings}

    r_conv = 1.0 / (h * area)
    q = delta_t / r_conv
    q_flux = h * delta_t

    return {"delta_t": delta_t, "r_conv": r_conv, "q": q, "q_flux": q_flux, "mode": mode, "is_valid": True, "warnings": warnings}


def calculate_radiation(emissivity: float, area: float, t_surface: float, t_surroundings: float) -> Dict[str, Any]:
    """Stefan-Boltzmann Law for Net Diffuse-Gray Thermal Radiation"""
    warnings = []
    if emissivity < 0.0 or emissivity > 1.0:
        warnings.append("Surface emissivity (ε) must be within physical bounds 0.0 ≤ ε ≤ 1.0.")
    if area <= 0:
        warnings.append("Surface area must be strictly greater than 0 m².")
    if t_surface <= 0 or t_surroundings <= 0:
        warnings.append("Absolute temperatures must be strictly greater than 0 K.")

    delta_t = t_surface - t_surroundings
    mode = "equilibrium"
    if delta_t > 0:
        mode = "emission"
    elif delta_t < 0:
        mode = "absorption"
    else:
        warnings.append("Radiative equilibrium: Net radiation is zero.")

    is_valid = (
        (len(warnings) == 0 or (len(warnings) == 1 and delta_t == 0))
        and (0 <= emissivity <= 1.0) and (area > 0) and (t_surface > 0) and (t_surroundings > 0)
    )
    if not is_valid:
        return {"q": 0.0, "q_flux": 0.0, "hr": 0.0, "r_rad": 0.0, "ts4": 0.0, "tsurr4": 0.0, "mode": mode, "is_valid": False, "warnings": warnings}

    ts4 = t_surface**4
    tsurr4 = t_surroundings**4
    q = emissivity * STEFAN_BOLTZMANN_CONSTANT * area * (ts4 - tsurr4)
    q_flux = emissivity * STEFAN_BOLTZMANN_CONSTANT * (ts4 - tsurr4)
    hr = emissivity * STEFAN_BOLTZMANN_CONSTANT * (t_surface + t_surroundings) * (t_surface**2 + t_surroundings**2)
    r_rad = 1.0 / (hr * area) if hr > 0 else 0.0

    return {"q": q, "q_flux": q_flux, "hr": hr, "r_rad": r_rad, "ts4": ts4, "tsurr4": tsurr4, "mode": mode, "is_valid": True, "warnings": warnings}


def calculate_lmtd(delta_t1: float, delta_t2: float) -> float:
    """Log Mean Temperature Difference with analytical singularity handling"""
    if delta_t1 <= 0 or delta_t2 <= 0:
        return 0.0
    if abs(delta_t1 - delta_t2) < 1e-6:
        return delta_t1
    return (delta_t1 - delta_t2) / math.log(delta_t1 / delta_t2)


def calculate_heat_exchanger(
    flow_arrangement: str,
    m_dot_hot: float,
    cp_hot: float,
    t_hot_in: float,
    t_hot_out: float,
    m_dot_cold: float,
    cp_cold: float,
    t_cold_in: float,
    t_cold_out: float,
    u: float,
    area: float,
) -> Dict[str, Any]:
    """Heat Exchanger LMTD Rating, First Law Balance & ε-NTU Method"""
    warnings = []
    if m_dot_hot <= 0 or m_dot_cold <= 0:
        warnings.append("Mass flow rates must be strictly positive (ṁ > 0 kg/s).")
    if cp_hot <= 0 or cp_cold <= 0:
        warnings.append("Specific heat capacities must be strictly positive (Cp > 0 J/kg·K).")
    if t_hot_in <= 0 or t_hot_out <= 0 or t_cold_in <= 0 or t_cold_out <= 0:
        warnings.append("Inlet and outlet temperatures must be strictly greater than 0 K.")
    if t_hot_in <= t_cold_in:
        warnings.append("Hot stream inlet (T_h,in) must be hotter than cold stream inlet (T_c,in).")
    if t_hot_in <= t_hot_out:
        warnings.append("Hot stream must reject heat: T_h,in must exceed T_h,out.")
    if t_cold_out <= t_cold_in:
        warnings.append("Cold stream must absorb heat: T_c,out must exceed T_c,in.")
    if u <= 0 or area <= 0:
        warnings.append("Overall heat transfer coefficient U and area A must be strictly positive.")

    c_hot = m_dot_hot * cp_hot
    c_cold = m_dot_cold * cp_cold
    q_hot = c_hot * (t_hot_in - t_hot_out)
    q_cold = c_cold * (t_cold_out - t_cold_in)
    q_avg = (q_hot + q_cold) / 2.0 if (q_hot + q_cold) > 0 else 0.0

    energy_discrepancy = abs(q_hot - q_cold)
    energy_discrepancy_pct = (energy_discrepancy / q_avg * 100.0) if q_avg > 0 else 0.0
    if energy_discrepancy_pct > 5.0 and q_avg > 0:
        warnings.append(f"Energy balance discrepancy: Q_hot ({q_hot/1e3:.2f} kW) and Q_cold ({q_cold/1e3:.2f} kW) differ by {energy_discrepancy_pct:.1f}%.")

    if flow_arrangement == "counter_current":
        delta_t1 = t_hot_in - t_cold_out
        delta_t2 = t_hot_out - t_cold_in
    else:
        delta_t1 = t_hot_in - t_cold_in
        delta_t2 = t_hot_out - t_cold_out

    has_crossover = False
    if delta_t1 <= 0 or delta_t2 <= 0:
        has_crossover = True
        warnings.append(f"Temperature Crossover Detected: ΔT1 = {delta_t1:.2f} K, ΔT2 = {delta_t2:.2f} K violates Second Law of Thermodynamics.")

    if flow_arrangement == "parallel_flow" and t_cold_out >= t_hot_out:
        has_crossover = True
        warnings.append("In parallel flow, T_c,out cannot exceed T_h,out.")

    lmtd = 0.0
    if not has_crossover and delta_t1 > 0 and delta_t2 > 0:
        lmtd = calculate_lmtd(delta_t1, delta_t2)

    q_sizing = u * area * lmtd
    c_min = min(c_hot, c_cold)
    c_max = max(c_hot, c_cold)
    cr = c_min / c_max if c_max > 0 else 0.0
    q_max = c_min * (t_hot_in - t_cold_in)
    effectiveness = q_avg / q_max if q_max > 0 else 0.0
    ntu = (u * area) / c_min if c_min > 0 else 0.0

    is_valid = not has_crossover and delta_t1 > 0 and delta_t2 > 0 and m_dot_hot > 0 and m_dot_cold > 0 and u > 0 and area > 0

    return {
        "c_hot": c_hot, "c_cold": c_cold, "q_hot": q_hot, "q_cold": q_cold, "q_avg": q_avg,
        "energy_discrepancy": energy_discrepancy, "energy_discrepancy_pct": energy_discrepancy_pct,
        "delta_t1": delta_t1, "delta_t2": delta_t2, "lmtd": lmtd, "q_sizing": q_sizing,
        "c_min": c_min, "c_max": c_max, "cr": cr, "q_max": q_max, "effectiveness": effectiveness,
        "ntu": ntu, "has_crossover": has_crossover, "is_valid": is_valid, "warnings": warnings,
    }


def calculate_thermal_resistance(t_inlet_fluid: float, t_outlet_fluid: float, layers: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Series Composite Network Analysis"""
    warnings = []
    if t_inlet_fluid <= 0 or t_outlet_fluid <= 0:
        warnings.append("Fluid temperatures must be strictly greater than 0 K.")
    if len(layers) == 0:
        warnings.append("At least one layer must be specified in the network.")

    delta_t_total = t_inlet_fluid - t_outlet_fluid
    r_total = 0.0
    calculated_layers = []

    for i, layer in enumerate(layers):
        name = layer.get("name", f"Layer {i+1}")
        l_type = layer.get("type", "conduction")
        area = layer.get("area", 1.0)
        r_val = 0.0

        if area <= 0:
            warnings.append(f"{name}: Area must be strictly positive.")
        elif l_type == "convection":
            h = layer.get("h", 0.0)
            if h <= 0:
                warnings.append(f"{name}: Convective coefficient h must be positive.")
            else:
                r_val = 1.0 / (h * area)
        elif l_type == "conduction":
            k = layer.get("k", 0.0)
            thickness = layer.get("thickness", 0.0)
            if k <= 0:
                warnings.append(f"{name}: Thermal conductivity k must be positive.")
            elif thickness <= 0:
                warnings.append(f"{name}: Thickness L must be positive.")
            else:
                r_val = thickness / (k * area)

        r_total += r_val
        calculated_layers.append({"name": name, "type": l_type, "r_value": r_val, "area": area})

    is_valid = len(warnings) == 0 and r_total > 0
    if not is_valid:
        return {"layers": [], "r_total": 0.0, "u_overall": 0.0, "q": 0.0, "q_flux": 0.0, "delta_t_total": delta_t_total, "is_valid": False, "warnings": warnings}

    q = delta_t_total / r_total
    primary_area = layers[0].get("area", 1.0)
    q_flux = q / primary_area if primary_area > 0 else 0.0
    u_overall = 1.0 / (r_total * primary_area) if primary_area > 0 else 0.0

    current_t = t_inlet_fluid
    final_layers = []
    for item in calculated_layers:
        delta_t_layer = q * item["r_value"]
        t_in = current_t
        t_out = current_t - delta_t_layer
        current_t = t_out
        pct = (item["r_value"] / r_total * 100.0) if r_total > 0 else 0.0
        final_layers.append({**item, "t_in": t_in, "t_out": t_out, "delta_t": delta_t_layer, "percent_of_total": pct})

    return {"layers": final_layers, "r_total": r_total, "u_overall": u_overall, "q": q, "q_flux": q_flux, "delta_t_total": delta_t_total, "is_valid": True, "warnings": warnings}


# -------------------------------------------------------------
# 2. PDF REPORT GENERATOR (ReportLab)
# -------------------------------------------------------------
def generate_pdf_report_bytes(
    analysis_title: str,
    analysis_type: str,
    governing_equation: str,
    inputs: List[Dict[str, str]],
    results: List[Dict[str, str]],
    assumptions: List[str],
    interpretation: str,
) -> bytes:
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle("HeaderTitle", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=colors.HexColor("#0f172a"))
    subtitle_style = ParagraphStyle("HeaderSub", parent=styles["Normal"], fontName="Helvetica", fontSize=9, leading=12, textColor=colors.HexColor("#64748b"))
    section_style = ParagraphStyle("SectionHeading", parent=styles["Normal"], fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=colors.HexColor("#ea580c"), spaceBefore=10, spaceAfter=4)
    body_style = ParagraphStyle("Body", parent=styles["Normal"], fontName="Helvetica", fontSize=8.5, leading=11, textColor=colors.HexColor("#334155"))
    eq_style = ParagraphStyle("Equation", parent=styles["Normal"], fontName="Courier-Bold", fontSize=9, leading=11, textColor=colors.HexColor("#0f172a"))

    story = []
    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    story.append(Paragraph("THERMO-X | Heat Transfer Analyser", title_style))
    story.append(Paragraph(f"Engineering-Grade Thermal Analysis Report &bull; {analysis_title} &bull; Generated: {timestamp}", subtitle_style))
    story.append(Spacer(1, 10))

    eq_data = [[Paragraph("<b>Governing Equation:</b>", body_style), Paragraph(governing_equation, eq_style)]]
    eq_table = Table(eq_data, colWidths=[120, 420])
    eq_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f8fafc")),
        ("BOX", (0, 0), (-1, -1), 1, colors.HexColor("#cbd5e1")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
    ]))
    story.append(eq_table)
    story.append(Spacer(1, 12))

    story.append(Paragraph("1. Operating & Boundary Parameters", section_style))
    input_data = [["Parameter Description", "Symbol", "Entered Value", "SI Unit"]]
    for inp in inputs:
        input_data.append([inp.get("name", ""), inp.get("symbol", ""), str(inp.get("value", "")), inp.get("unit", "")])
    inp_table = Table(input_data, colWidths=[200, 70, 150, 120])
    inp_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#334155")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#ffffff"), colors.HexColor("#f1f5f9")]),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd5e1")),
        ("ALIGN", (2, 1), (2, -1), "RIGHT"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(inp_table)
    story.append(Spacer(1, 12))

    story.append(Paragraph("2. Computed Engineering Results", section_style))
    res_data = [["Parameter Description", "Symbol", "Calculated Value", "SI Unit", "Meaning"]]
    for res in results:
        res_data.append([res.get("parameter", ""), res.get("symbol", ""), str(res.get("value", "")), res.get("unit", ""), res.get("description", "")])
    res_table = Table(res_data, colWidths=[150, 60, 110, 60, 160])
    res_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#ea580c")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#ffffff"), colors.HexColor("#fff7ed")]),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#fed7aa")),
        ("ALIGN", (2, 1), (2, -1), "RIGHT"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(res_table)
    story.append(Spacer(1, 12))

    story.append(Paragraph("3. Governing Assumptions", section_style))
    for a in assumptions:
        story.append(Paragraph(f"&bull; {a}", body_style))
    story.append(Spacer(1, 10))

    story.append(Paragraph("4. Engineering Interpretation & Insights", section_style))
    story.append(Paragraph(interpretation, body_style))
    story.append(Spacer(1, 14))

    disclaimer_text = "<b>DISCLAIMER:</b> THERMO-X is an educational and engineering analysis tool. Results should be independently verified before use in safety-critical applications."
    story.append(Paragraph(disclaimer_text, subtitle_style))

    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()


# -------------------------------------------------------------
# 3. STREAMLIT APP CONFIGURATION & UI
# -------------------------------------------------------------
st.set_page_config(
    page_title="THERMO-X | Heat Transfer Analyser",
    page_icon="🔥",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .main { background-color: #f8fafc; }
    .stMetric { background-color: #ffffff; padding: 14px; border-radius: 10px; border: 1px solid #e2e8f0; }
    .stMetric label { font-size: 0.8rem; font-weight: 700; color: #475569; }
    .stMetric div[data-testid="stMetricValue"] { font-size: 1.4rem; font-weight: 800; color: #0f172a; }
    .block-container { padding-top: 2rem; padding-bottom: 2rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

with st.sidebar:
    st.markdown("### 🔥 **THERMO-X**")
    st.caption("**Heat Transfer Analyser &bull; SI PRO**")
    st.markdown("---")

    analysis_mode = st.selectbox(
        "Select Analysis Type",
        [
            "Home / Dashboard",
            "Conduction",
            "Convection",
            "Radiation",
            "Heat Exchanger",
            "Thermal Resistance",
        ],
        index=0,
    )

    st.markdown("---")
    st.markdown("#### **SI Engineering Standards**")
    st.markdown(
        """
        - Temperature: **K**
        - Heat Rate: **W**
        - Heat Flux: **W/m²**
        - Conductivity: **W/m·K**
        - Convection Coeff: **W/m²·K**
        - Resistance: **K/W**
        """
    )
    st.markdown("---")
    st.caption("Streamlit Community Cloud Ready &bull; v2.4")


# -------------------------------------------------------------
# 4. VIEW ROUTER
# -------------------------------------------------------------
if analysis_mode == "Home / Dashboard":
    st.title("THERMO-X")
    st.subheader("Heat Transfer Analyser")
    st.markdown("#### **Engineering-grade thermal analysis for students, engineers and researchers.**")
    st.info("THERMO-X is an interactive engineering tool for analysing heat-transfer systems using established thermal engineering equations. Select an analysis method from the sidebar to begin.")
    st.markdown("---")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("### 🧱 Conduction\nAnalyse 1D steady-state conductive transport using Fourier's Law ($Q = kA\\Delta T/L$).")
    with c2:
        st.markdown("### 💨 Convection\nEvaluate boundary convective heat transfer using Newton's Law of Cooling ($Q = hA\\Delta T$).")
    with c3:
        st.markdown("### ☀️ Radiation\nCalculate net thermal radiation with Stefan-Boltzmann Law ($Q = \\varepsilon \\sigma A (T^4 - T_{surr}^4)$).")
    c4, c5, _ = st.columns(3)
    with c4:
        st.markdown("### 🔄 Heat Exchanger\nCounter-Current & Parallel Flow heat exchangers with LMTD rating and $\\varepsilon$-NTU sizing.")
    with c5:
        st.markdown("### 🧱 Thermal Resistance\nMulti-layer series composite walls with internal/external convection boundaries.")

elif analysis_mode == "Conduction":
    st.title("Conduction Analysis")
    st.caption("One-dimensional steady-state heat conduction through a plane wall (Fourier's Law)")
    col_in, col_res = st.columns([1, 1])
    with col_in:
        st.markdown("### Operating Parameters")
        k = st.number_input("Thermal Conductivity k [W/m·K]", min_value=0.001, value=0.50, step=0.05, format="%.3f")
        area = st.number_input("Cross-Sectional Area A [m²]", min_value=0.001, value=10.0, step=0.5, format="%.2f")
        length = st.number_input("Wall Thickness L [m]", min_value=0.0001, value=0.10, step=0.01, format="%.4f")
        t_hot = st.number_input("Hot Surface Temp T_hot [K]", min_value=1.0, value=400.0, step=5.0)
        t_cold = st.number_input("Cold Surface Temp T_cold [K]", min_value=1.0, value=300.0, step=5.0)

    res = calculate_conduction(k=k, area=area, length=length, t_hot=t_hot, t_cold=t_cold)
    for w in res["warnings"]:
        st.warning(w) if res["is_valid"] else st.error(w)

    with col_res:
        st.markdown("### Principal Results")
        m1, m2 = st.columns(2)
        m1.metric("Heat Rate Q", f"{res['q']:,.2f} W", f"{res['q']/1000:,.3f} kW")
        m2.metric("Heat Flux q''", f"{res['q_flux']:,.2f} W/m²")
        m3, m4 = st.columns(2)
        m3.metric("Thermal Resistance R_cond", f"{res['r_cond']:.4f} K/W")
        m4.metric("Temperature Potential ΔT", f"{res['delta_t']:.2f} K")

    st.markdown("---")
    st.subheader("📈 Temperature Profile T(x) Through Solid Wall")
    if res["is_valid"]:
        x_vals = np.linspace(0, length, 50)
        t_vals = t_hot - (res["q_flux"] / k) * x_vals
        fig = go.Figure(go.Scatter(x=x_vals, y=t_vals, mode="lines", name="T(x) Profile", line=dict(color="#ea580c", width=3)))
        fig.update_layout(xaxis_title="Position x [m]", yaxis_title="Temperature T [K]", height=320, margin=dict(l=40, r=20, t=20, b=40))
        st.plotly_chart(fig, use_container_width=True)

    df_results = pd.DataFrame([
        {"Parameter": "Temperature Difference", "Symbol": "ΔT", "Value": f"{res['delta_t']:.2f}", "Unit": "K", "Description": "Thermal potential across plane wall"},
        {"Parameter": "Thermal Resistance", "Symbol": "R_cond", "Value": f"{res['r_cond']:.4f}", "Unit": "K/W", "Description": "Conduction resistance L/(kA)"},
        {"Parameter": "Heat-Transfer Rate", "Symbol": "Q", "Value": f"{res['q']:,.2f}", "Unit": "W", "Description": "Total steady heat dissipation rate"},
        {"Parameter": "Heat Flux", "Symbol": "q''", "Value": f"{res['q_flux']:,.2f}", "Unit": "W/m²", "Description": "Heat rate per unit area (Q/A)"},
    ])
    st.dataframe(df_results, use_container_width=True, hide_index=True)

    assumptions = ["Steady-state 1D heat flow.", "Constant conductivity k.", "No internal generation.", "Constant area A."]
    interpretation = f"Steady-state heat rate is {res['q']:,.2f} W with conductive resistance {res['r_cond']:.4f} K/W."
    
    st.markdown("---")
    c_csv, c_pdf = st.columns(2)
    with c_csv:
        st.download_button("📥 Download Results as CSV", df_results.to_csv(index=False).encode("utf-8"), "THERMO_X_Conduction.csv", "text/csv")
    with c_pdf:
        pdf_bytes = generate_pdf_report_bytes("1D Plane Wall Conduction Analysis", "Conduction", "Q = kA(T_hot - T_cold)/L = ΔT/R_cond",
            [{"name": "Thermal Conductivity", "symbol": "k", "value": str(k), "unit": "W/m·K"},
             {"name": "Cross-Sectional Area", "symbol": "A", "value": str(area), "unit": "m²"},
             {"name": "Wall Thickness", "symbol": "L", "value": str(length), "unit": "m"},
             {"name": "Hot Surface Temp", "symbol": "T_hot", "value": str(t_hot), "unit": "K"},
             {"name": "Cold Surface Temp", "symbol": "T_cold", "value": str(t_cold), "unit": "K"}],
            [{"parameter": row["Parameter"], "symbol": row["Symbol"], "value": row["Value"], "unit": row["Unit"], "description": row["Description"]} for _, row in df_results.iterrows()],
            assumptions, interpretation)
        st.download_button("📄 Generate Engineering Report (PDF)", pdf_bytes, "THERMO_X_Conduction_Report.pdf", "application/pdf")

elif analysis_mode == "Convection":
    st.title("Convection Analysis")
    st.caption("Newton's Law of Cooling for fluid boundary convective exchange")
    col_in, col_res = st.columns([1, 1])
    with col_in:
        st.markdown("### Operating Parameters")
        h = st.number_input("Heat Transfer Coefficient h [W/m²·K]", min_value=0.1, value=20.0, step=2.0)
        area = st.number_input("Surface Area A [m²]", min_value=0.001, value=5.0, step=0.5)
        t_surface = st.number_input("Surface Temp T_surface [K]", min_value=1.0, value=350.0, step=5.0)
        t_fluid = st.number_input("Fluid Temp T_fluid [K]", min_value=1.0, value=300.0, step=5.0)

    res = calculate_convection(h=h, area=area, t_surface=t_surface, t_fluid=t_fluid)
    for w in res["warnings"]:
        st.warning(w) if res["is_valid"] else st.error(w)

    with col_res:
        st.markdown("### Principal Results")
        m1, m2 = st.columns(2)
        m1.metric("Heat Rate Q", f"{res['q']:,.2f} W", f"{res['q']/1000:,.3f} kW")
        m2.metric("Heat Flux q''", f"{res['q_flux']:,.2f} W/m²")
        m3, m4 = st.columns(2)
        m3.metric("Convective Resistance R_conv", f"{res['r_conv']:.5f} K/W")
        m4.metric("Temperature Potential ΔT", f"{res['delta_t']:.2f} K")

    df_conv = pd.DataFrame([
        {"Parameter": "Temperature Difference", "Symbol": "ΔT", "Value": f"{res['delta_t']:.2f}", "Unit": "K", "Description": "Surface to fluid temperature potential"},
        {"Parameter": "Convective Resistance", "Symbol": "R_conv", "Value": f"{res['r_conv']:.5f}", "Unit": "K/W", "Description": "Boundary layer resistance 1/(hA)"},
        {"Parameter": "Heat-Transfer Rate", "Symbol": "Q", "Value": f"{res['q']:,.2f}", "Unit": "W", "Description": f"Convective heat transfer ({res['mode']} mode)"},
        {"Parameter": "Heat Flux", "Symbol": "q''", "Value": f"{res['q_flux']:,.2f}", "Unit": "W/m²", "Description": "Heat flux h*(Ts - T_inf)"},
    ])
    st.dataframe(df_conv, use_container_width=True, hide_index=True)

    st.markdown("---")
    c_csv, c_pdf = st.columns(2)
    with c_csv:
        st.download_button("📥 Download Results as CSV", df_conv.to_csv(index=False).encode("utf-8"), "THERMO_X_Convection.csv", "text/csv")
    with c_pdf:
        pdf_bytes = generate_pdf_report_bytes("Convective Heat Transfer Analysis", "Convection", "Q = hA(T_surface - T_fluid) = ΔT/R_conv",
            [{"name": "Heat Transfer Coefficient", "symbol": "h", "value": str(h), "unit": "W/m²·K"},
             {"name": "Surface Area", "symbol": "A", "value": str(area), "unit": "m²"},
             {"name": "Surface Temp", "symbol": "T_surface", "value": str(t_surface), "unit": "K"},
             {"name": "Fluid Temp", "symbol": "T_fluid", "value": str(t_fluid), "unit": "K"}],
            [{"parameter": row["Parameter"], "symbol": row["Symbol"], "value": row["Value"], "unit": row["Unit"], "description": row["Description"]} for _, row in df_conv.iterrows()],
            ["Steady-state conditions.", "Uniform convective coefficient h.", "Constant bulk fluid temperature."],
            f"Newton's law indicates a net convective heat-transfer rate of {res['q']:,.2f} W in {res['mode']} mode.")
        st.download_button("📄 Generate Engineering Report (PDF)", pdf_bytes, "THERMO_X_Convection_Report.pdf", "application/pdf")

elif analysis_mode == "Radiation":
    st.title("Radiation Analysis")
    st.caption("Stefan-Boltzmann T⁴ Law for diffuse-gray net thermal radiation exchange")
    col_in, col_res = st.columns([1, 1])
    with col_in:
        st.markdown("### Operating Parameters")
        emissivity = st.slider("Surface Emissivity ε [-]", min_value=0.0, max_value=1.0, value=0.80, step=0.01)
        area = st.number_input("Surface Area A [m²]", min_value=0.001, value=2.0, step=0.2)
        t_surface = st.number_input("Surface Temp T_surface [K]", min_value=1.0, value=500.0, step=10.0)
        t_surroundings = st.number_input("Surroundings Temp T_surr [K]", min_value=1.0, value=300.0, step=10.0)

    res = calculate_radiation(emissivity=emissivity, area=area, t_surface=t_surface, t_surroundings=t_surroundings)
    for w in res["warnings"]:
        st.warning(w) if res["is_valid"] else st.error(w)

    with col_res:
        st.markdown("### Principal Results")
        m1, m2 = st.columns(2)
        m1.metric("Net Radiation Q_rad", f"{res['q']:,.2f} W", f"{res['q']/1000:,.3f} kW")
        m2.metric("Radiative Flux q''_rad", f"{res['q_flux']:,.2f} W/m²")
        m3, m4 = st.columns(2)
        m3.metric("Linearized Coeff h_r", f"{res['hr']:.3f} W/m²·K")
        m4.metric("Radiative Resistance R_rad", f"{res['r_rad']:.5f} K/W")

    df_rad = pd.DataFrame([
        {"Parameter": "Net Radiative Heat Rate", "Symbol": "Q_rad", "Value": f"{res['q']:,.2f}", "Unit": "W", "Description": "Net radiative energy exchange"},
        {"Parameter": "Radiative Heat Flux", "Symbol": "q''_rad", "Value": f"{res['q_flux']:,.2f}", "Unit": "W/m²", "Description": "Radiative flux per unit area"},
        {"Parameter": "Linearized Radiation Coeff", "Symbol": "h_r", "Value": f"{res['hr']:.3f}", "Unit": "W/m²·K", "Description": "Effective radiative coefficient"},
        {"Parameter": "Radiative Resistance", "Symbol": "R_rad", "Value": f"{res['r_rad']:.5f}", "Unit": "K/W", "Description": "Linearized equivalent resistance 1/(h_r*A)"},
    ])
    st.dataframe(df_rad, use_container_width=True, hide_index=True)

    st.markdown("---")
    c_csv, c_pdf = st.columns(2)
    with c_csv:
        st.download_button("📥 Download Results as CSV", df_rad.to_csv(index=False).encode("utf-8"), "THERMO_X_Radiation.csv", "text/csv")
    with c_pdf:
        pdf_bytes = generate_pdf_report_bytes("Thermal Radiation Analysis", "Radiation", "Q_rad = ε·σ·A·(T_s⁴ - T_surr⁴)",
            [{"name": "Emissivity", "symbol": "ε", "value": str(emissivity), "unit": "-"},
             {"name": "Surface Area", "symbol": "A", "value": str(area), "unit": "m²"},
             {"name": "Surface Temp", "symbol": "T_surface", "value": str(t_surface), "unit": "K"},
             {"name": "Surroundings Temp", "symbol": "T_surr", "value": str(t_surroundings), "unit": "K"}],
            [{"parameter": row["Parameter"], "symbol": row["Symbol"], "value": row["Value"], "unit": row["Unit"], "description": row["Description"]} for _, row in df_rad.iterrows()],
            ["Diffuse and gray surface (α = ε).", "Large blackbody surroundings.", "Transparent non-participating medium."],
            f"The net radiative heat-transfer rate is {res['q']:,.2f} W based on the Stefan-Boltzmann fourth-power law.")
        st.download_button("📄 Generate Engineering Report (PDF)", pdf_bytes, "THERMO_X_Radiation_Report.pdf", "application/pdf")

elif analysis_mode == "Heat Exchanger":
    st.title("Heat Exchanger Analysis")
    st.caption("Counter-Current & Parallel Flow thermal rating, LMTD singularity handling, and ε-NTU evaluation")
    flow = st.radio("Flow Arrangement", ["Counter-Current", "Parallel Flow"], horizontal=True)
    flow_key = "counter_current" if flow == "Counter-Current" else "parallel_flow"

    col_h, col_c, col_s = st.columns(3)
    with col_h:
        st.markdown("#### Hot Fluid Stream")
        m_dot_hot = st.number_input("Mass Flow Rate ṁ_h [kg/s]", min_value=0.001, value=1.50, step=0.1)
        cp_hot = st.number_input("Specific Heat Cp_h [J/kg·K]", min_value=1.0, value=4184.0, step=50.0)
        t_hot_in = st.number_input("Inlet Temp T_h,in [K]", min_value=1.0, value=363.15, step=5.0)
        t_hot_out = st.number_input("Outlet Temp T_h,out [K]", min_value=1.0, value=323.15, step=5.0)
    with col_c:
        st.markdown("#### Cold Fluid Stream")
        m_dot_cold = st.number_input("Mass Flow Rate ṁ_c [kg/s]", min_value=0.001, value=2.00, step=0.1)
        cp_cold = st.number_input("Specific Heat Cp_c [J/kg·K]", min_value=1.0, value=4184.0, step=50.0)
        t_cold_in = st.number_input("Inlet Temp T_c,in [K]", min_value=1.0, value=293.15, step=5.0)
        t_cold_out = st.number_input("Outlet Temp T_c,out [K]", min_value=1.0, value=323.15, step=5.0)
    with col_s:
        st.markdown("#### Exchanger Sizing")
        u = st.number_input("Overall Coeff U [W/m²·K]", min_value=1.0, value=850.0, step=50.0)
        area = st.number_input("Exchanger Area A [m²]", min_value=0.01, value=12.0, step=1.0)

    res = calculate_heat_exchanger(flow_key, m_dot_hot, cp_hot, t_hot_in, t_hot_out, m_dot_cold, cp_cold, t_cold_in, t_cold_out, u, area)
    for w in res["warnings"]:
        st.warning(w) if res["is_valid"] else st.error(w)

    st.markdown("---")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("LMTD", f"{res['lmtd']:.2f} K")
    m2.metric("Sizing Heat Duty Q", f"{res['q_sizing']:,.1f} W", f"{res['q_sizing']/1e3:,.2f} kW")
    m3.metric("Effectiveness ε", f"{res['effectiveness']*100:.1f} %")
    m4.metric("NTU", f"{res['ntu']:.3f}")

    df_he = pd.DataFrame([
        {"Parameter": "Hot Stream Capacity Rate", "Symbol": "C_hot", "Value": f"{res['c_hot']:,.1f}", "Unit": "W/K", "Description": "ṁ_h * Cp_h"},
        {"Parameter": "Cold Stream Capacity Rate", "Symbol": "C_cold", "Value": f"{res['c_cold']:,.1f}", "Unit": "W/K", "Description": "ṁ_c * Cp_c"},
        {"Parameter": "Log Mean Temp Difference", "Symbol": "LMTD", "Value": f"{res['lmtd']:.2f}", "Unit": "K", "Description": "Mean driving temperature difference"},
        {"Parameter": "Sizing Heat Rate", "Symbol": "Q_sizing", "Value": f"{res['q_sizing']:,.1f}", "Unit": "W", "Description": "U * A * LMTD"},
        {"Parameter": "Thermal Effectiveness", "Symbol": "ε", "Value": f"{res['effectiveness']*100:.2f}", "Unit": "%", "Description": "Q / Q_max"},
        {"Parameter": "Number of Transfer Units", "Symbol": "NTU", "Value": f"{res['ntu']:.3f}", "Unit": "-", "Description": "UA / C_min"},
    ])
    st.dataframe(df_he, use_container_width=True, hide_index=True)

    st.markdown("---")
    c_csv, c_pdf = st.columns(2)
    with c_csv:
        st.download_button("📥 Download Results as CSV", df_he.to_csv(index=False).encode("utf-8"), "THERMO_X_HeatExchanger.csv", "text/csv")
    with c_pdf:
        pdf_bytes = generate_pdf_report_bytes(f"Heat Exchanger Analysis ({flow})", "HeatExchanger", "Q = U·A·LMTD,  ε = Q/Q_max",
            [{"name": "Flow Arrangement", "symbol": "Flow", "value": flow, "unit": "-"},
             {"name": "Hot Mass Flow", "symbol": "ṁ_h", "value": str(m_dot_hot), "unit": "kg/s"},
             {"name": "Hot In/Out Temp", "symbol": "Th,in/out", "value": f"{t_hot_in}/{t_hot_out}", "unit": "K"},
             {"name": "Cold Mass Flow", "symbol": "ṁ_c", "value": str(m_dot_cold), "unit": "kg/s"},
             {"name": "Cold In/Out Temp", "symbol": "Tc,in/out", "value": f"{t_cold_in}/{t_cold_out}", "unit": "K"},
             {"name": "Overall U", "symbol": "U", "value": str(u), "unit": "W/m²·K"},
             {"name": "Area", "symbol": "A", "value": str(area), "unit": "m²"}],
            [{"parameter": row["Parameter"], "symbol": row["Symbol"], "value": row["Value"], "unit": row["Unit"], "description": row["Description"]} for _, row in df_he.iterrows()],
            ["Steady-state adiabatic operation.", "Constant fluid properties Cp and overall U.", "No phase change occurring."],
            f"Exchanger operates at LMTD = {res['lmtd']:.2f} K with thermal effectiveness ε = {res['effectiveness']*100:.1f}%.")
        st.download_button("📄 Generate Engineering Report (PDF)", pdf_bytes, "THERMO_X_HeatExchanger_Report.pdf", "application/pdf")

elif analysis_mode == "Thermal Resistance":
    st.title("Thermal Resistance Network")
    st.caption("Series composite wall and multi-layer thermal resistance network")
    col_in, col_res = st.columns([1, 1])
    with col_in:
        st.markdown("### Boundary Temperatures & Layer Properties")
        t_in = st.number_input("Inside Fluid Temp T_inlet [K]", min_value=1.0, value=370.0, step=5.0)
        t_out = st.number_input("Outside Fluid Temp T_outlet [K]", min_value=1.0, value=290.0, step=5.0)
        area = st.number_input("Surface Area A [m²]", min_value=0.01, value=4.0, step=0.5)

        st.markdown("#### Series Layers")
        h_in = st.number_input("Inside Convection h_in [W/m²·K]", min_value=0.1, value=25.0, step=2.0)
        k_wall = st.number_input("Wall Conductivity k_wall [W/m·K]", min_value=0.001, value=0.72, step=0.05)
        l_wall = st.number_input("Wall Thickness L_wall [m]", min_value=0.001, value=0.12, step=0.01)
        k_ins = st.number_input("Insulation Conductivity k_ins [W/m·K]", min_value=0.001, value=0.038, step=0.005)
        l_ins = st.number_input("Insulation Thickness L_ins [m]", min_value=0.001, value=0.05, step=0.01)
        h_out = st.number_input("Outside Convection h_out [W/m²·K]", min_value=0.1, value=50.0, step=5.0)

    layers_config = [
        {"name": "Inside Convection", "type": "convection", "h": h_in, "area": area},
        {"name": "Brick Wall", "type": "conduction", "k": k_wall, "thickness": l_wall, "area": area},
        {"name": "Fiberglass Insulation", "type": "conduction", "k": k_ins, "thickness": l_ins, "area": area},
        {"name": "Outside Convection", "type": "convection", "h": h_out, "area": area},
    ]

    res = calculate_thermal_resistance(t_inlet_fluid=t_in, t_outlet_fluid=t_out, layers=layers_config)
    for w in res["warnings"]:
        st.warning(w) if res["is_valid"] else st.error(w)

    with col_res:
        st.markdown("### Principal Results")
        m1, m2 = st.columns(2)
        m1.metric("Total Heat Rate Q", f"{res['q']:,.2f} W", f"{res['q']/1000:,.3f} kW")
        m2.metric("Total Resistance R_tot", f"{res['r_total']:.4f} K/W")
        m3, m4 = st.columns(2)
        m3.metric("Overall U-Value", f"{res['u_overall']:.3f} W/m²·K")
        m4.metric("Heat Flux q''", f"{res['q_flux']:,.2f} W/m²")

    table_data = [
        {"Parameter": "Total Heat Rate", "Symbol": "Q", "Value": f"{res['q']:,.2f}", "Unit": "W", "Description": "ΔT_total / R_total"},
        {"Parameter": "Total Resistance", "Symbol": "R_total", "Value": f"{res['r_total']:.4f}", "Unit": "K/W", "Description": "Σ R_i"},
        {"Parameter": "Overall U-Value", "Symbol": "U", "Value": f"{res['u_overall']:.3f}", "Unit": "W/m²·K", "Description": "1 / (R_total * A)"},
    ]
    for l in res["layers"]:
        table_data.append({
            "Parameter": l["name"],
            "Symbol": "R_i",
            "Value": f"{l['r_value']:.4f} ({l['percent_of_total']:.1f}%)",
            "Unit": "K/W",
            "Description": f"T_in = {l['t_in']:.1f}K → T_out = {l['t_out']:.1f}K (ΔT = {l['delta_t']:.1f}K)",
        })
    df_tr = pd.DataFrame(table_data)
    st.dataframe(df_tr, use_container_width=True, hide_index=True)

    st.markdown("---")
    c_csv, c_pdf = st.columns(2)
    with c_csv:
        st.download_button("📥 Download Results as CSV", df_tr.to_csv(index=False).encode("utf-8"), "THERMO_X_ThermalResistance.csv", "text/csv")
    with c_pdf:
        pdf_bytes = generate_pdf_report_bytes("Thermal Resistance Network Analysis", "ThermalResistance", "R_total = Σ R_i,  Q = ΔT / R_total",
            [{"name": "Inside Fluid Temp", "symbol": "T_in", "value": str(t_in), "unit": "K"},
             {"name": "Outside Fluid Temp", "symbol": "T_out", "value": str(t_out), "unit": "K"},
             {"name": "Surface Area", "symbol": "A", "value": str(area), "unit": "m²"}],
            [{"parameter": row["Parameter"], "symbol": row["Symbol"], "value": row["Value"], "unit": row["Unit"], "description": row["Description"]} for _, row in df_tr.iterrows()],
            ["1D steady-state heat conduction and convection.", "Negligible contact resistance.", "Uniform heat transfer coefficients."],
            f"Total thermal resistance is {res['r_total']:.4f} K/W across ΔT = {res['delta_t_total']:.1f} K.")
        st.download_button("📄 Generate Engineering Report (PDF)", pdf_bytes, "THERMO_X_ThermalResistance_Report.pdf", "application/pdf")
