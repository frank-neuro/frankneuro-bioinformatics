import streamlit as st
import primer3
import re
import os
import signal

st.set_page_config(page_title="Multi-Junction Primer Designer", layout="wide")

st.markdown("""
    <style>
    div.stButton > button:first-child { background-color: #28a745; color: white; border: none; }
    div.stButton > button:first-child:hover { background-color: #218838; color: white; }
    .blast-button { 
        display: inline-block; width: 100%; text-align: center; 
        background-color: #007bff; color: white !important; 
        padding: 10px; text-decoration: none; border-radius: 5px; 
        font-weight: bold; margin-bottom: 10px; 
    }
    .summary-card {
        background-color: #2d3436;
        border: 2px solid #2ecc71;
        border-radius: 12px;
        padding: 25px;
        margin-bottom: 25px;
        color: #dfe6e9;
    }
    .legend-box {
        display: inline-block;
        padding: 2px 10px;
        border-radius: 3px;
        color: black;
        font-weight: bold;
        margin-right: 10px;
    }
    </style>
""", unsafe_allow_html=True)

if 'junction_results' not in st.session_state:
    st.session_state.junction_results = []
if 'full_sequence' not in st.session_state:
    st.session_state.full_sequence = None

def get_dimer_info(seq1, seq2=None):
    try:
        if seq2 is None: res = primer3.bindings.calc_homodimer(seq1)
        else: res = primer3.bindings.calc_heterodimer(seq1, seq2)
        return res.dg / 1000 
    except: return 0.0

def calculate_gc(seq):
    if not seq: return 0.0
    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100

def get_rev_comp(seq):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))

st.sidebar.header("🔍 Evaluation Strategy")
num_per_junction = st.sidebar.slider("Pairs to find per Junction", 1, 5, 2)

st.sidebar.header("Global Primer Settings")
opt_size = st.sidebar.slider("Optimal Size", 15, 30, 20)
min_size = st.sidebar.number_input("Min Size", value=18)
max_size = st.sidebar.number_input("Max Size", value=25)

st.sidebar.markdown("---")
opt_tm = st.sidebar.slider("Optimal Tm (°C)", 50.0, 70.0, 60.0)
min_tm = st.sidebar.number_input("Min Tm", value=55.0)
max_tm = st.sidebar.number_input("Max Tm", value=65.0)

st.sidebar.markdown("---")
prod_min = st.sidebar.number_input("Min Product Size", value=100)
prod_max = st.sidebar.number_input("Max Product Size", value=500)

tab1, tab2 = st.tabs(["🧬 Multi-Junction RT-PCR Primer Designer", "📘 Help & Documentation"])

with tab1:
    st.title("🧬 Multi-Junction RT-PCR Primer Designer")
    st.info("**Strategy:** To determine optimal primer pairs for RT-PCR that generate amplicons spanning exon-exon boundaries across an entire genomic sequence resulting in the exclusion of genomic DNA amplification.")
    
    raw_input = st.text_area("Paste sequence here with **EXONS** in UPPERCASE and **introns** in lowercase::", height=200)
    
    if st.button("Analyze All Junctions", type="primary"):
        if raw_input:
            lines = [line.strip() for line in raw_input.splitlines() if not line.startswith(">")]
            full_seq = "".join(lines)
            st.session_state.full_sequence = full_seq
            
            s_template = ""
            junctions = []
            for i, char in enumerate(full_seq):
                if char.isupper():
                    s_template += char
                    if i > 0 and full_seq[i-1].islower():
                        junctions.append(len(s_template) - 1)
            
            all_found = []
            global_args = {
                'PRIMER_NUM_RETURN': num_per_junction,
                'PRIMER_OPT_SIZE': opt_size, 'PRIMER_MIN_SIZE': min_size, 'PRIMER_MAX_SIZE': max_size,
                'PRIMER_OPT_TM': opt_tm, 'PRIMER_MIN_TM': min_tm, 'PRIMER_MAX_TM': max_tm,
                'PRIMER_PRODUCT_SIZE_RANGE': [[prod_min, prod_max]],
            }

            for idx, j_pos in enumerate(junctions):
                seq_args = {
                    'SEQUENCE_ID': f'Junction_{idx+1}',
                    'SEQUENCE_TEMPLATE': s_template.upper(),
                    'SEQUENCE_TARGET': [j_pos, 1]
                }
                res = primer3.bindings.design_primers(seq_args, global_args)
                if res.get('PRIMER_PAIR_NUM_RETURNED', 0) > 0:
                    all_found.append({'id': idx+1, 'pos': j_pos, 'results': res})
            
            st.session_state.junction_results = all_found

    if st.session_state.junction_results:
        results_list = st.session_state.junction_results
        
        st.divider()
        st.header("📋 Summary Report: Recommended Optimal Pair")
        
        master_data = []
        for j_entry in results_list:
            j_id = j_entry['id']
            res = j_entry['results']
            for i in range(res.get('PRIMER_PAIR_NUM_RETURNED')):
                f_seq = res.get(f'PRIMER_LEFT_{i}_SEQUENCE')
                r_seq = res.get(f'PRIMER_RIGHT_{i}_SEQUENCE')
                
                f_tm = round(res.get(f'PRIMER_LEFT_{i}_TM'), 2)
                r_tm = round(res.get(f'PRIMER_RIGHT_{i}_TM'), 2)
                f_gc = round(calculate_gc(f_seq), 1)
                r_gc = round(calculate_gc(r_seq), 1)
                f_sdg = round(get_dimer_info(f_seq), 2)
                r_sdg = round(get_dimer_info(r_seq), 2)
                c_dg = round(get_dimer_info(f_seq, r_seq), 2)
                
                master_data.append({
                    "Junction ID": j_id,
                    "Penalty Score": round(res.get(f'PRIMER_PAIR_{i}_PENALTY'), 2),
                    "Cross-ΔG": c_dg,
                    "Forward": f_seq,
                    "Reverse": r_seq,
                    "Fwd Tm": f_tm,
                    "Rev Tm": r_tm,
                    "Fwd GC": f_gc,
                    "Rev GC": r_gc,
                    "Fwd Self-ΔG": f_sdg,
                    "Rev Self-ΔG": r_sdg,
                    "Product Size": res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE')
                })
        
        master_data.sort(key=lambda x: (x['Penalty Score'], -x['Cross-ΔG']))
        
        for idx, item in enumerate(master_data):
            item["Rank"] = idx + 1
        
        best = master_data[0]
        
        f_status = "Low Risk" if best["Fwd Self-ΔG"] > -9 else "High Risk"
        r_status = "Low Risk" if best["Rev Self-ΔG"] > -9 else "High Risk"
        c_status = "Low Risk" if best["Cross-ΔG"] > -9 else "High Risk"

        st.markdown(f"""
        <div class="summary-card">
            <h3 style="color: #2ecc71; margin-top:0;">🌟 Top Recommendation: Junction {best['Junction ID']}</h3>
            <table style="width:100%; color: #dfe6e9; border-collapse: collapse;">
                <tr style="border-bottom: 1px solid #444;">
                    <td><b>Forward:</b> <code>{best["Forward"]}</code></td>
                    <td><b>Tm:</b> {best["Fwd Tm"]}°C | <b>GC:</b> {best["Fwd GC"]}%</td>
                    <td><b>Self-ΔG:</b> {best["Fwd Self-ΔG"]}</td>
                </tr>
                <tr>
                    <td><b>Reverse:</b> <code>{best["Reverse"]}</code></td>
                    <td><b>Tm:</b> {best["Rev Tm"]}°C | <b>GC:</b> {best["Rev GC"]}%</td>
                    <td><b>Self-ΔG:</b> {best["Rev Self-ΔG"]}</td>
                </tr>
            </table>
        </div>
        """, unsafe_allow_html=True)

        report_text = f"PRIMER DESIGN SUMMARY REPORT\n"
        report_text += f"============================\n"
        report_text += f"Target: Junction {best['Junction ID']}\n\n"
        report_text += f"FORWARD: {best['Forward']}\n"
        report_text += f"  - Tm: {best['Fwd Tm']}C | GC: {best['Fwd GC']}%\n"
        report_text += f"  - Self-Dimer ΔG: {best['Fwd Self-ΔG']} (Status: {f_status})\n\n"
        report_text += f"REVERSE: {best['Reverse']}\n"
        report_text += f"  - Tm: {best['Rev Tm']}C | GC: {best['Rev GC']}%\n"
        report_text += f"  - Self-Dimer ΔG: {best['Rev Self-ΔG']} (Status: {r_status})\n\n"
        report_text += f"PAIR METRICS:\n"
        report_text += f"  - Product Size: {best['Product Size']} bp\n"
        report_text += f"  - Tm Difference: {round(abs(best['Fwd Tm'] - best['Rev Tm']), 2)}C\n"
        report_text += f"  - Cross-Dimer ΔG: {best['Cross-ΔG']}\n"
        report_text += f"  - Cross-Dimer Status: {c_status}\n"
        
        st.download_button(
            label="📥 Download Detailed Summary (.txt)", 
            data=report_text, 
            file_name=f"junction_{best['Junction ID']}_report.txt",
            mime="text/plain"
        )

        st.subheader("📊 Global Junction Evaluation")
        st.caption(
            "⚠️ **Note on Specificity:** High-efficiency primers may occasionally hit multiple genes in "
            "multigenic genomic sequences. If non-specific hits occur, evaluate alternative junctions where "
            "selectivity is prioritized over peak thermodynamic scores."
        )
        display_cols = ["Rank", "Junction ID", "Penalty Score", "Cross-ΔG", "Forward", "Reverse", "Fwd Tm", "Rev Tm", "Product Size"]
        st.dataframe(master_data, column_order=display_cols, use_container_width=True)

        with st.expander("📖 Terminology & Scoring Guide"):
            st.markdown("""
            * **Penalty Score:** The 'cost' assigned to a primer pair. Points are added when a primer deviates from your 'Optimal' settings. **0.00 is a perfect match.**
            * **Cross- and Self-ΔG:** The change in **Gibbs Free Energy** (kcal/mol) required for the Forward and Reverse primers to form a dimer. 
                * **Values closer to 0:** Indicate an **unstable** dimer (Desired).
                * **More negative values:** Indicate a **more stable** dimer, increasing the risk of primer-dimer artifacts.
                * **Threshold:** ΔG < -9.0 is generally considered high risk for primer dimer formation in RT-PCR.
            * **Tm (Melting Temperature):** The temperature at which 50% of the DNA duplex dissociates. Optimal RT-PCR usually requires pairs within 2°C of each other.
            """)

        st.divider()
        st.subheader("📍 Junction-Specific Results")
        j_choice = st.selectbox("Select Junction:", [f"Junction {r['id']}" for r in results_list])
        target_id = int(j_choice.split()[-1])
	# 1. Get the ID number from the string (e.g., "Junction 2" -> 2)
        target_id = int(j_choice.split()[-1])

        # 2. Find the entry where the 'id' matches
        j_entry = next((item for item in results_list if item.get('id') == target_id), None)

        if j_entry:
            j_res = j_entry['results']
            st.success(f"Displaying results for Junction {target_id}")
        else:
            st.error(f"Could not find data for Junction {target_id}.")
            j_res = []
        p_idx = st.radio(f"Select Pair for {j_choice}:", range(j_res.get('PRIMER_PAIR_NUM_RETURNED')), horizontal=True)
        
        f_s = j_res.get(f'PRIMER_LEFT_{p_idx}_SEQUENCE')
        r_s = j_res.get(f'PRIMER_RIGHT_{p_idx}_SEQUENCE')
        f_tm_val = j_res.get(f'PRIMER_LEFT_{p_idx}_TM')
        r_tm_val = j_res.get(f'PRIMER_RIGHT_{p_idx}_TM')
        
        f_gc = calculate_gc(f_s); r_gc = calculate_gc(r_s)
        f_dg = get_dimer_info(f_s); r_dg = get_dimer_info(r_s); c_dg = get_dimer_info(f_s, r_s)
        p_size = j_res.get(f'PRIMER_PAIR_{p_idx}_PRODUCT_SIZE')

        c1, c2 = st.columns(2)
        with c1:
            st.write(f"**Forward:** `{f_s}`"); st.write(f"**Tm:** `{f_tm_val:.2f}°C` | **GC:** `{f_gc:.1f}%` ")
            st.write(f"**Self-Dimer ΔG:** `{f_dg:.2f}`"); st.caption("✅ Status: Low Risk" if f_dg > -9 else "⚠️ Status: High Risk")
        with c2:
            st.write(f"**Reverse:** `{r_s}`"); st.write(f"**Tm:** `{r_tm_val:.2f}°C` | **GC:** `{r_gc:.1f}%` ")
            st.write(f"**Self-Dimer ΔG:** `{r_dg:.2f}`"); st.caption("✅ Status: Low Risk" if r_dg > -9 else "⚠️ Status: High Risk")

        st.divider()
        st.subheader("📊 Primer Pair Analysis")
        m1, m2 = st.columns(2)
        with m1: st.write(f"**Expected Product Size:** `{p_size} bp`")
        with m2: 
            st.write(f"**Cross-Dimer ΔG:** `{c_dg:.2f}`")
            if c_dg > -9: st.caption("✅ Status: Low Risk")
            else: st.warning("⚠️ Status: High Risk (Cross-Dimer)")

        st.divider()
        st.subheader("📍 Transcript-Level Primer Binding Map")
        st.caption("""
            This map visualizes the primer binding sites within the spliced mRNA (cDNA) sequence. 
            **Note:** The reported 'Reverse' sequence in tables and reports is the **5'→3' reverse complement**, 
            which is complementary to the cyan highlighted binding site shown below.
        """)
        st.markdown('<span class="legend-box" style="background-color: yellow;">Forward Primer</span> <span class="legend-box" style="background-color: cyan;">Reverse Primer</span>', unsafe_allow_html=True)

        full_gen = st.session_state.full_sequence
        r_rc = get_rev_comp(r_s)
        
        genomic_to_mrna = []
        for i, char in enumerate(full_gen):
            if char.isupper():
                genomic_to_mrna.append(i)
        
        mrna_temp = "".join([c for c in full_gen if c.isupper()]).upper()
        
        def get_exact_exon_indices(query, mrna_str, mapping):
            hit_indices = set()
            if not query: return hit_indices
            start = 0
            while True:
                idx = mrna_str.find(query.upper(), start)
                if idx == -1: break
                for char_offset in range(len(query)):
                    hit_indices.add(mapping[idx + char_offset])
                start = idx + 1
            return hit_indices

        f_exon_set = get_exact_exon_indices(f_s, mrna_temp, genomic_to_mrna)
        r_exon_set = get_exact_exon_indices(r_rc, mrna_temp, genomic_to_mrna)

        html_out = ""
        for i, char in enumerate(full_gen):
            if i in f_exon_set:
                html_out += f"<mark style='background-color: yellow; color: black;'><b>{char}</b></mark>"
            elif i in r_exon_set:
                html_out += f"<mark style='background-color: cyan; color: black;'><b>{char}</b></mark>"
            else:
                html_out += char

        st.markdown(f"""<div style="font-family:monospace; background-color:#0e1117; color:#d1d1d1; padding:20px; border-radius:10px; word-wrap:break-word; line-height:2.2;">{html_out}</div>""", unsafe_allow_html=True)
        


with tab2:
    st.header("📘 Help & Documentation")
    st.subheader("1. How the Designer Works")
    st.write("This tool designs primers to generate amplicons that span exon-exon junctions, thereby minimizing amplification of genomic DNA in RT-PCR applications.")
    
    st.subheader("2. Step-by-Step Instructions")
    st.markdown("""
    1. **Enter**: Paste your converted FASTA sequence in the gray box. Ensure Exons are UPPERCASE and introns are lowercase. Use the "FASTA Converter" program to convert Genbank sequences into this FASTA format. 
    2. **Settings**: Select the number of primer pairs per junction under Evaluation Strategy and Adjust Global Primer Settings to define the thermodynamic and structural constraints of a primer pair. 
    4. **Process**: Click the **'Analyze All Junctions'** button to generate primer pairs for each exon-exon junction.
    5. **Review**: Check primer pairs to assess the differential in Tm (< 2 °C is optimal) and risk of primer dimer formation (A ΔG > -5 kcal/mol is low risk, ΔG between -5 to -9 kcal/mol is moderate risk, and ΔG =/< -9 kcal/mol is high risk).
    6. **Validate**: Click the 'Primer BLAST Validation' button in the sidebar, which will take you to Primer-BLAST at the National Center for Biotechnology Information (NCBI). Enter the primer sequences under 'Primer Parameters'. Select organism under 'Primer Pair Specificity Parameters' and click the 'Get Primers' button. Primer-BLAST results will validate the gene specificity of your primer pair and size of the amplicon.
    """)
    
    st.subheader("3. How to Cite")
    st.write("If you use this tool in your research, please use the following citation:")
    st.code("Frank, M.G. (2026). Primer Designer (v1.0). [Python Software].")
    
    st.subheader("4. MIT License")
    with st.expander("View License Terms"):
        st.markdown("""
        **Copyright (c) 2026 M.G. Frank**

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
        """)
    
    st.divider()
    st.caption("Developed by M.G. Frank | February 2026")

st.sidebar.markdown("---")
st.sidebar.markdown('<a href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/" target="_blank" class="blast-button">🚀 Primer BLAST Validation</a>', unsafe_allow_html=True)
if st.sidebar.button("Shut Down App", use_container_width=True):
    os.kill(os.getpid(), signal.SIGINT)