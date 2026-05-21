import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

# --- 1. CORE KINETIC MATH ---

def calculate_threshold_cq(cycles, rfu, threshold):
    """Interpolates the fractional cycle where signal crosses threshold."""
    baseline = np.mean(rfu[:5])
    f_norm = rfu - baseline
    for i in range(len(f_norm) - 1):
        if f_norm[i] < threshold <= f_norm[i+1]:
            m = (f_norm[i+1] - f_norm[i]) / (cycles[i+1] - cycles[i])
            if m == 0: return None
            b = f_norm[i] - m * cycles[i]
            return (threshold - b) / m
    return None

def calculate_linreg_efficiency(cycles, rfu, cq, window=3):
    """Calculates instantaneous efficiency at the Cq point."""
    if cq is None or np.isnan(cq): return None
    baseline = np.mean(rfu[:5])
    f_norm = rfu - baseline
    center_idx = np.abs(cycles - cq).argmin()
    start_idx = max(0, center_idx - (window // 2))
    end_idx = min(len(cycles), start_idx + window)
    x_win, y_win = cycles[start_idx:end_idx], f_norm[start_idx:end_idx]
    if np.any(y_win <= 0): return None
    slope, _ = np.polyfit(x_win, np.log10(y_win), 1)
    return 10**slope

# --- 2. UI CONFIGURATION ---

st.set_page_config(page_title="qPCR Diagnostic Suite", layout="wide")
st.title("🛡️ qPCR Diagnostic Suite: Global Kinetics & Local Precision")

with st.sidebar:
    st.header("Input & Settings")
    raw_file = st.file_uploader("Upload Raw RFU Values in Excel (Cycle in Column A, Well RFU A1, A2, A3.... in Subsequent Columns)", type=["xlsx"])
    global_threshold = st.number_input("Global Threshold (RFU)", value=50.0, step=5.0)
    st.divider()
    st.info("Approach: This step calculates Global Master Efficiency, which is a key kinetic parameter of gene expression assays. Select all wells and adjust Global Threshold in the exponential phase of amplification to minimize replicate standard deviation while keeping Global Master Efficiency near 2.0. An ideal assay efficiency lands between 1.9 and 2.1.")
    st.divider()
    st.info("Application: Assay Global Master Efficiency for target and reference genes, and utilize these efficiency parameters to calculate relative gene expression using the Pfaffl equation (See Pfaffl Application).")
    st.divider()
    st.info("Pfaffl Method: A mathematical model for calculating relative gene expression, which factors in specific PCR amplification efficiencies for both target and reference genes thereby increasing accuracy of gene expression measurements.")

tab_map, tab_qc, tab_health, tab_help = st.tabs([
    "📍 Mapping", "📈 Efficiency Computation and Quality Control", "🏥 Plate Health Dashboard", "ℹ️ Help & Documentation"
])

if raw_file:
    df_raw = pd.read_excel(raw_file)
    
    # Case-insensitive look for any cycle identifier label 
    cycle_candidates = [col for col in df_raw.columns if 'cycle' in str(col).lower().strip()]
    chosen_cycle_col = cycle_candidates[0] if cycle_candidates else df_raw.columns[0]
    
    well_names = [col for col in df_raw.columns if col != chosen_cycle_col]
    cycles_arr = df_raw[chosen_cycle_col].to_numpy()

    if 'map_df' not in st.session_state:
        st.session_state.map_df = pd.DataFrame({
            "Well_ID": well_names,
            "Sample_ID": ["Control"] * len(well_names),
            "Gene_ID": ["IL33"] * len(well_names),
            "Include": [True] * len(well_names)
        })

    # --- TAB 1: MAPPING ---
    with tab_map:
        st.subheader("Experimental Mapping")
        st.session_state.map_df = st.data_editor(
            st.session_state.map_df,
            num_rows="fixed", width="stretch", hide_index=True,
            column_config={
                "Well_ID": st.column_config.Column(disabled=True),
                "Include": st.column_config.CheckboxColumn("Include in Analysis")
            },
            key="map_editor"
        )

    # --- TAB 2: EFFICIENCY COMPUTATION AND QUALITY CONTROL ---
    with tab_qc:
        st.subheader("Selection Diagnostics")
        
        mapping = st.session_state.map_df
        well_options = [f"{row.Well_ID} ({row.Sample_ID} | {row.Gene_ID})" for row in mapping.itertuples()]
        
        col_sel_a, col_sel_b = st.columns([1, 4])
        select_all = col_sel_a.checkbox("Select All Wells")
        
        selected_options = well_options if select_all else st.multiselect(
            "Select Wells to Assess", 
            options=well_options, 
            default=well_options[:3] if len(well_options) >= 3 else well_options
        )

        if selected_options:
            subset_data = []
            plot_data = []

            for opt in selected_options:
                w_id = opt.split(" ")[0]
                y_raw = df_raw[w_id].to_numpy()
                cq = calculate_threshold_cq(cycles_arr, y_raw, global_threshold)
                eff = calculate_linreg_efficiency(cycles_arr, y_raw, cq)
                
                meta = mapping[mapping['Well_ID'] == w_id].iloc[0]
                
                if cq:
                    subset_data.append({
                        "Well_ID": w_id, "Cq": cq, "E": eff, 
                        "Sample_ID": meta.Sample_ID, "Gene_ID": meta.Gene_ID
                    })
                
                plot_data.append(pd.DataFrame({
                    'Cycle': cycles_arr, 
                    'RFU': y_raw - np.mean(y_raw[:5]), 
                    'Well': opt
                }))

            if subset_data:
                res_df = pd.DataFrame(subset_data)
                
                # Metrics Header
                m1, m2, m3 = st.columns(3)
                
                master_e = res_df['E'].mean()
                e_sd = res_df['E'].std() if len(res_df) > 1 else 0.0
                e_cv = (e_sd / master_e) * 100 if master_e != 0 else 0.0
                m1.metric("Global Master Efficiency", f"{master_e:.3f}", f"Standard Deviation: {e_sd:.4f}")
                
                m2.metric("Efficiency Coefficient of Variance (%)", f"{e_cv:.2f}%")
                
                mean_cq = res_df['Cq'].mean()
                cq_sd = res_df['Cq'].std() if len(res_df) > 1 else 0.0
                m3.metric("Selection Mean Cq", f"{mean_cq:.2f}", f"Standard Deviation: {cq_sd:.3f}")

                # Overlay Plot with RED Threshold Line
                combined_plot_df = pd.concat(plot_data)
                fig = px.line(combined_plot_df, x='Cycle', y='RFU', color='Well', title="Kinetic Overlay")
                fig.add_hline(y=global_threshold, line_dash="dash", line_color="red", 
                             annotation_text=f"Threshold: {global_threshold}", 
                             annotation_position="top left")
                st.plotly_chart(fig, use_container_width=True)

                # Sample-Level Summary Table Below Graph
                st.markdown("#### Selection Summary (Local Variance)")
                selection_summary = res_df.groupby(['Sample_ID', 'Gene_ID']).agg(
                    Mean_Cq=('Cq', 'mean'),
                    SD_Cq=('Cq', 'std'),
                    Mean_E=('E', 'mean'),
                    SD_E=('E', 'std')
                ).reset_index()
                
                st.dataframe(selection_summary.style.format({
                    'Mean_Cq': '{:.2f}', 'SD_Cq': '{:.3f}', 
                    'Mean_E': '{:.3f}', 'SD_E': '{:.4f}'
                }), width="stretch")
            else:
                st.warning("Selected wells did not cross threshold.")
        else:
            st.info("Select wells to begin assessment.")

    # --- TAB 3: PLATE HEALTH DASHBOARD ---
    with tab_health:
        if st.button("🚀 Run Plate Health Audit"):
            results = []
            for w in well_names:
                rfu = df_raw[w].to_numpy()
                c = calculate_threshold_cq(cycles_arr, rfu, global_threshold)
                e = calculate_linreg_efficiency(cycles_arr, rfu, c)
                results.append({"Well_ID": w, "Cq": c, "Indiv_E": e})
            
            full_df = pd.merge(st.session_state.map_df, pd.DataFrame(results), on="Well_ID")
            included_df = full_df[full_df['Include'] == True].copy()
            
            global_metrics = included_df.groupby('Gene_ID').agg(
                Master_E=('Indiv_E', 'mean'), E_SD=('Indiv_E', 'std')
            ).reset_index()
            
            summary = pd.merge(included_df, global_metrics, on="Gene_ID")
            final_summary = summary.groupby(['Sample_ID', 'Gene_ID']).agg(
                Mean_Cq=('Cq', 'mean'), SD_Cq=('Cq', 'std'), Group_E=('Indiv_E', 'mean')
            ).reset_index()

            st.subheader("🏥 Diagnostic Health Checks")
            h1, h2, h3 = st.columns(3)
            avg_local_sd = final_summary['SD_Cq'].mean()
            h1.metric("Mean Replicate Standard Deviation", f"{avg_local_sd:.3f}", delta="PASS" if avg_local_sd < 0.3 else "HIGH VAR")
            plate_e = global_metrics['Master_E'].mean()
            h2.metric("Mean Global Efficiency", f"{plate_e:.3f}", delta="OPTIMAL" if 1.9 <= plate_e <= 2.1 else "SUBOPTIMAL")
            h3.metric("Wells Scanned", len(included_df))

            st.divider()
            c_left, c_right = st.columns(2)
            with c_left:
                st.write("**Target Kinetic Consistency**")
                st.plotly_chart(px.box(included_df, x="Gene_ID", y="Indiv_E", color="Gene_ID", points="all"), use_container_width=True)
            with c_right:
                st.write("**Technical Replicate Jitter**")
                st.plotly_chart(px.strip(included_df, x="Sample_ID", y="Cq", color="Gene_ID"), use_container_width=True)

            st.subheader("🚩 Variance Flag Table")
            def flag_sd(val):
                return 'background-color: #ffcccc' if val > 0.3 else 'background-color: #ccffcc'
            st.dataframe(final_summary.style.map(flag_sd, subset=['SD_Cq']).format({
                'Mean_Cq': '{:.2f}', 'SD_Cq': '{:.3f}', 'Group_E': '{:.3f}'
            }), width="stretch")

    # --- TAB 4: HELP & DOCUMENTATION ---
    with tab_help:
        st.header("📘 Help & Documentation Suite")
        st.markdown("---")
        
        doc_col_1, doc_col_2 = st.columns([2, 1])
        
        with doc_col_1:
            st.subheader("💡 Application Overview")
            st.markdown("""
            This suite is an interactive tool designed to optimize and audit Real-Time qPCR data using kinetic principles. By calculating a **Global Master Efficiency** per gene target rather than assuming an arbitrary 100% efficiency ($E = 2.0$), this application eliminates systemic mathematical bias before your dataset proceeds to relative expression analysis (e.g., via the Pfaffl equation).
            
            ### 🛠️ Core Optimization Workflow
            1. **Map Your Plate:** Organize your layout using your custom `Sample_ID` and `Gene_ID` labels in the **Mapping** tab.
            2. **Adjust the Threshold:** Use the sidebar control to position the red threshold line securely within the *log-linear exponential phase* of your amplification curves.
            3. **Audit Variance:** Watch the green **Standard Deviation** and **Coefficient of Variance (CV%)** metrics. Your optimal threshold resides where technical replicate variance is minimized while keeping the Global Efficiency parameter stable.
            4. **Examine Plate Health:** Switch to the **Plate Health Dashboard** to spot systemic outlier wells or potential tissue-derived amplification inhibition.
            """)
            
            st.markdown("""
            ### 🧮 Mathematical Model for Efficiency Computation
            
            The application computes individual well amplification efficiency ($E$) by applying an **Ordinary Least Squares (OLS) Linear Regression** model to the log-transformed fluorescence values within a localized cycle window centered around the calculated $Cq$ point.
            
            #### 1. Baseline Correction and Normalization
            Before kinetic modeling, raw fluorescence data ($RFU_{\text{raw}}$) undergoes background baseline subtraction to yield the true normalized target fluorescence ($F_n$):
            $$F_n = RFU_{\text{raw}} - \mu(RFU_{1\dots5})$$
            where $\mu(RFU_{1\dots5})$ represents the mean background noise calculated from cycles 1 through 5.
            
            #### 2. Localized Log-Linear Regression
            A sliding analytical window of 3 cycles is centered dynamically on the integer cycle nearest to the calculated fractional $Cq$ value. Within this window—where the reaction kinetics adhere strictly to exponential doubling—the data is log-transformed to establish a linear relationship:
            $$\log_{10}(F_n) = m \cdot \text{Cycle} + b$$
            
            Using the calculated slope ($m$) from the linear regression, the instantaneous **Amplification Efficiency ($E$)** parameter is derived via the inverse log transformation:
            $$E = 10^m$$
            
            #### 3. Mathematical Interpretation and Application
            The fundamental kinetic equation for exponential PCR amplification is modeled as:
            $$N_n = N_0 \cdot E^n$$

            Where:
            * $N_n$ = Total template accumulation amount (normalized fluorescence) at cycle $n$
            * $N_0$ = **Starting copy number** (initial target template quantity at cycle 0)
            * $E$ = **Amplification Efficiency parameter** (where $E = 2.000$ indicates perfect 100% template duplication per cycle)
            * $n$ = **Cycle number**
            
            The **Global Master Efficiency** displayed in the app is the arithmetic mean ($\bar{E}$) of these individual slopes across all selected wells:
            $$\bar{E} = \frac{1}{k}\sum_{i=1}^{k} E_i$$
            
            This validated master parameter directly replaces the theoretical value of $2.0$ within downstream quantification equations (e.g., the Pfaffl model: *Pfaffl, M. W. (2001). A new mathematical model for relative quantification in real-time RT–PCR. Nucleic Acids Research, 29(9), e45–e45. https://doi.org/10.1093/nar/29.9.e45*), correcting for assay-specific kinetic biases.
            """)

            st.markdown("---")
            st.subheader("🖋️ How to Cite This Tool")
            st.markdown("""
            If you utilize this diagnostic software suite to analyze, optimize, or audit qPCR kinetics for a peer-reviewed publication, please use the following citation format:
            
            > **Frank, M. G.** *(2026)*. *qPCR Diagnostic Suite: Global Kinetics & Local Precision Analytics for Relative Expression Optimization.* Open-source software framework for kinetic parameter evaluation.
            """)

        with doc_col_2:
            st.subheader("⚖️ Open-Source License")
            st.markdown("""
            This software is released under the standard **MIT License**. You are free to modify, distribute, and utilize it for both private academic research and commercial pipelines.
            """)
            
            st.code("""
MIT License

Copyright (c) 2026 M.G. Frank

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
""", language="text")
else:
    st.info("Awaiting Excel file upload in the sidebar.")