**Referee Report --- v3 (commit main@HEAD)**

Your latest revision tightens structure, adds references, and inserts a
**Methodological** section. CI is green and the manuscript now passes
the build hook. Below is a focused third‑pass review.

**1 · Progress since v2**

| **Area**                  | **Improvement**                                                                                   | **Evidence**                                                                                                                        |
|---------------------------|---------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Section heads**         | Added numbered subsections; *Methodological Considerations* and *Solar‑System Tests* now present. |  Lines 24--30 ([GitHub](https://raw.githubusercontent.com/charlesrotter/UDT/main/Universal_Distance_Dilation_Theory_Manuscript.md)) |
| **References list**       | 8 peer‑reviewed citations inserted.                                                               |  Lines 38‑41 ([GitHub](https://raw.githubusercontent.com/charlesrotter/UDT/main/Universal_Distance_Dilation_Theory_Manuscript.md))  |
| **Limitations paragraph** | Explicitly notes data‑contamination issues (bullets 1--4).                                        |  Lines 25‑27 ([GitHub](https://raw.githubusercontent.com/charlesrotter/UDT/main/Universal_Distance_Dilation_Theory_Manuscript.md))  |

*Nice momentum---layout and scholarly apparatus now resemble a real
pre‑print.*

**2 · Show‑stopper scientific problems (must fix)**

| **\#** | **Issue**                                                                                                                                               | **Location**              | **Required action**                                                                                               |
|--------|---------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|-------------------------------------------------------------------------------------------------------------------|
| **A**  | **Quantum‑scale radius still mis‑computed**: Eq. 5 gives R0,q=1.319 × 107 mR\_{0,q}=1.319 × 10\^{7}\\m; correct value is 8.4 × 10−5 m8.4 × 10\^{-5}\\m. |  L16                      | Replace numeric, propagate through force table & any derived masses.                                              |
| **B**  | **F (τ) singularity** unchanged: τ ≈ 0.01 → F≈ 2 × 10¹ (!). Your "galactic 1.02--1.09" claim conflicts.                                                 |  L13 & force table L16‑17 | Provide regulator term *or* prove observed τ ≥ 0.1 for all SPARC radii; update text & table.                      |
| **C**  | **"Perfect 5 / 5 validation" & χ²/DOF ≈ 3.13 called "excellent"**---χ² \> 2 is normally *poor*.                                                         |  Abstract L2 & L20        | Drop superlatives; quote median χ², RMS, and p‑value. Include a table or figure.                                  |
| **D**  | **Instantaneous c_fundamental** breaks local Lorentz invariance.                                                                                        |  L8‑L9                    | Add a paragraph outlining the foliation / preferred‑frame mechanism or cite a formalism (e.g., Horava--Lifshitz). |
| **E**  | **Muon g‑2 "418 % agreement"** is not agreement.                                                                                                        |  L22                      | Supply numeric derivation; express deviation in σ, not %. If outside 2σ, classify as *open problem*.              |

Any journal referee will halt the process on A--C alone.

**3 · Repository integration gaps**

1.  **Figures & Tables**  
    *Repository contains no fig/ outputs.*  
    ‑ Modify scripts/analyze_sparc_galaxies.py to write:  
    fig/rotation_examples.png, fig/chi2_cdf.png,
    results/sparc_full.csv  
    Add Markdown links: \![Figure 1\](fig/rotation_examples.png).

2.  **Manuscript build in CI**  
    Add job build-paper to .github/workflows/audit.yml:

3.  \- name: Build manuscript

4.  run: \|

5.  pip install pandoc

6.  pandoc Universal_Distance_Dilation_Theory_Manuscript.md \\

7.  \--citeproc -o udt_preprint.pdf

8.  \- uses: actions/upload-artifact@v4

9.  with:

10. name: udt_preprint

11. path: udt_preprint.pdf

12. **Data manifest citation**  
    In *Data Availability* paragraph add:  
    "Integrity verified by sha256sum -c data/manifest_sha256.txt."

**4 · Minor editorial notes**

- Remove marketing terms ("first successful", "paradigm shift") unless
  backed by citations or mainstream recognition.

- Use SI spacing: 4.74 km s⁻¹, 10 ms, 8.4 × 10⁻⁵ m.

- Cross‑link equations: "...as shown in Eq. (3)".

**5 · Action checklist for next revision**

1.  **Fix Eq. 5 numeric and cascade changes.**

2.  **Resolve F(τ) inconsistency** (regulator or empirical bound).

3.  **Replace validation rhetoric with quantitative stats + figure.**

4.  **Add figure PNGs & table CSV to repo; link in manuscript.**

5.  **Implement build-paper CI job; artifact PDF appears.**

6.  **Clarify Lorentz‑symmetry mechanism or cite literature.**

When these items land and CI remains green, I'll focus on stylistic
polish and final compliance with journal templates.
