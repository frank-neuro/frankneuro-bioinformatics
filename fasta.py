# fasta.py - GenBank to FASTA Converter
# Copyright (c) 2026 Matthew G. Frank
# Licensed under the MIT License

import streamlit as st
from Bio import SeqIO
import io


# --- PAGE CONFIG ---
st.set_page_config(page_title="GenBank Converter", layout="centered")

tab1, tab2 = st.tabs(["🧬 Converter", "❓ Help & Documentation"])

with tab1:
    st.title("GenBank to FASTA Converter")
    uploaded_file = st.file_uploader("Upload a GenBank file (.gb, .gbk)", type=["gb", "gbk",    	"genbank"])
    if uploaded_file is not None:
        if st.button("Convert and Show"):
            try:
                # --- YOUR CODE LIVES HERE NOW ---
                results = []
                stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
                
                for record in SeqIO.parse(stringio, "genbank"):
                    sequence = list(str(record.seq).lower())
                    for feature in record.features:
                        if feature.type == "mRNA":
                            for part in feature.location.parts:
                                start, end = int(part.start), int(part.end)
                                sequence[start:end] = list("".join(sequence[start:end]).upper())
                    
                    fasta_entry = f">{record.id} {record.description}\n" + "".join(sequence)
                    results.append(fasta_entry)
                
                full_fasta = "\n\n".join(results)
                st.success("Successfully processed!")
                st.text_area("FASTA Result:", value=full_fasta, height=400)
                
                st.download_button(
                    label="Download FASTA File",
                    data=full_fasta,
                    file_name=f"{uploaded_file.name.split('.')[0]}.fasta",
                    mime="text/plain"
                )
            except Exception as e:
                st.error(f"An error occurred: {e}")
    else:
        st.info("Please upload a GenBank file to begin.")

with tab2:
    st.header("📘 Help & Documentation")
    
    st.subheader("1. How the Conversion Works")
    st.write("""
    This tool converts GenBank (.gb) records into FASTA format while preserving 
    structural information through **Letter Casing**:
    """)
    
    # Visual Legend for the user
    st.info("""
    - **UPPERCASE (ATGC):** Exons / Coding Regions (derived from `mRNA` features).
    - **lowercase (atgc):** Introns / Non-coding regions.
    """)

    

    st.subheader("2. Step-by-Step Instructions")
    st.markdown("""
    1. **Upload**: Drag your `.gb` or `.gbk` file into the uploader in the **Converter** tab.
    2. **Process**: Click the **'Convert and Show'** button to parse and highlight exons and introns.
    3. **Review**: Check the 'FASTA Result' text area to verify the sequence looks correct.
    4. **Download**: Click the download button. The file will automatically be named after your original file but with a `.fasta` extension.
    """)

    st.subheader("3. Troubleshooting")
    with st.expander("Why is my entire sequence in lowercase?"):
        st.write("""
        If the output is entirely lowercase, it usually means the GenBank file 
        does not contain any features labeled as **'mRNA'**. The script relies 
        specifically on the `mRNA` tag to determine which parts to uppercase.
        """)
    st.subheader("4. Technical Specifications")
    st.markdown("""
    1. **Port**: This application runs on port `8520`.
    2. **Environment**: Uses the `primer_env` virtual environment.
    3. **Input Format**: Requires standard GenBank records with `mRNA` feature tags for proper case-delineation.
    4. **Core Dependencies**: Built with `Biopython` and `Streamlit`.
    5. **Memory**: Large genomic records (>50MB) may experience slower processing times.
    6. **Troubleshooting Port Conflicts**: If the application fails to start because the port is "already in use," you can force-close the background process by typing this into your Mac Terminal: lsof -ti:8520 | xargs kill -9
    """)

    st.subheader("5. How to Cite")
    st.write("If you use this tool in your research or lab reports, please use the following citation:")
    st.code("Frank, M.G. (2026). GenBank to FASTA Converter (v1.0). [Python Software].")

    st.subheader("6. MIT License")
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

