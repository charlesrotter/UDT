Below is a **self‑contained blueprint** you can hand to whatever agent
(human or AI) is maintaining the codebase.  
It condenses everything we've discussed into a standards‑driven
*Definition of Done* for a professional, collaborative, and
fully‑reproducible scientific repository.

**0.  Big‑Picture Goal --- "Reproducible by Strangers"**

**Anyone** with Python ≥ 3.9, Git, and Internet access can:  
1. clone the repo,  
2. install one deterministic environment command,  
3. run one target (make audit-all) and  
4. receive the *same* tables, plots, and pass/fail badges you do.

**1.  Repository Layout & Naming**

UDT/ \# repo root

│

├── src/udt/ \# importable package; no loose .py at top level

│ ├── \_\_init\_\_.py

│ ├── models/ \# theory equations, likelihoods

│ ├── data_loaders/

│ ├── diagnostics/

│ ├── cli.py \# entry‑point wrapper

│ └── \...

│

├── data/ \# \*small\* files only (≤ 100 KB) + manifest

│ ├── sample/ \# 1‑2 example files browsable on GitHub

│ └── manifest_sha256.txt

│

├── docs/ \# rendered by mkdocs or Sphinx

│ ├── index.md

│ ├── theory.md

│ └── api/

│

├── tests/ \# pytest suite; mirror src/udt/ tree

│

├── examples/ \# Jupyter / MyST notebooks

│

├── scripts/ \# thin CLI wrappers; no core logic

│

└── .github/workflows/ \# CI pipelines

*Follow the
[cookiecutter‑scientific‑python](https://github.com/cookiecutter/cookiecutter)
layout so external reviewers know where to look.*

**2.  Environment & Dependency Management**

| **Requirement**                | **Implementation**                                                          |
|--------------------------------|-----------------------------------------------------------------------------|
| **Single‑command setup**       | pip install \--no-build-isolation -e .\[dev\]                               |
| **Deterministic versions**     | Use pyproject.toml + requirements.lock *(pinned hashes via pip‑tools)*      |
| **Heavy wheels pre‑installed** | CI step: pip install \--prefer-binary numpy scipy *before* editable install |
| **Cross‑platform**             | Test matrix: ubuntu‑latest, macos‑latest, optional windows‑latest           |

**3.  Data Management & Integrity**

| **Item**                  | **Standard**                                                                                                  |
|---------------------------|---------------------------------------------------------------------------------------------------------------|
| **Large files**           | Keep in a release asset or externally hosted tarball; never in Git LFS unless reviewers explicitly need it.   |
| **Immutable fingerprint** | data/manifest_sha256.txt (one SHA‑256 per file, sorted). CI job sha256sum -c gates every push.                |
| **Sample data**           | ≤ 10 KB per file, committed in‑tree for quick notebook rendering.                                             |
| **Data licence**          | Include data/README.md citing original surveys (SPARC, CSP DR3, Pantheon +, Planck) and redistribution terms. |

**4.  Testing & Validation**

| **Layer**               | **What to test**                                             | **Target coverage**                             |
|-------------------------|--------------------------------------------------------------|-------------------------------------------------|
| **Unit**                | Pure functions: field eqs, likelihoods, loaders              | ≥ 70 %                                          |
| **Property‑based**      | SymPy vs. numeric equality across random inputs (hypothesis) | Catch algebraic slip‑ups                        |
| **Integration / Smoke** | 3‑galaxy SPARC fit; asserts χ²/d.o.f \< 2                    | Runs in \< 60 s                                 |
| **Full bench**          | 175‑galaxy sweep, Pantheon +, Planck hook                    | Manual or nightly CI (can be behind make bench) |

**CI fails on *any* test error or checksum mismatch.**

**5.  Continuous Integration (GitHub Actions)**

name: UDT Validation Tests

on:

push: {branches: \[main\]}

pull_request:

jobs:

test:

strategy:

matrix: {pyver: \[\"3.9\", \"3.10\", \"3.11\"\]}

runs-on: ubuntu-latest

steps:

\- uses: actions/checkout@v4

\- name: Install deps

run: \|

python -m pip install \--upgrade pip wheel

pip install \--prefer-binary numpy scipy matplotlib

pip install \--no-build-isolation -e .\[dev\]

\- name: Verify data integrity

run: sha256sum -c data/manifest_sha256.txt

\- name: Run test suite

run: pytest -v \--durations=10

\- name: Upload artefacts (smoke results)

if: success()

uses: actions/upload-artifact@v4

with:

name: sparc_smoke

path: results/sparc_smoke\_\*.\*

Add a **nightly cron** job bench that fetches the big tarball and runs
the full suite.

**6.  Documentation**

| **Section**     | **Tooling / Content**                                            |
|-----------------|------------------------------------------------------------------|
| **README.md**   | 1‑minute install, 1‑minute run, project status badge             |
| **Theory docs** | docs/theory.md (Derivation → Assumptions → Limits)               |
| **API docs**    | auto‑generated by Sphinx (sphinx‑autosummary)                    |
| **Repro guide** | docs/reproduce.md -- exact commands to replicate full figures    |
| **CHANGELOG**   | Keep with [Keep a Changelog](https://keepachangelog.com/) format |
| **Citation**    | CITATION.cff + Zenodo DOI on releases                            |

**7.  Collaboration & Governance**

- **CONTRIBUTING.md** -- coding style (black + ruff), PR checklist, how
  to run tests locally.

- **CODE_OF_CONDUCT.md** -- adapt Contributor Covenant v2.

- **Issue & PR templates** -- bug report, feature request, data update.

- **Semantic versioning** -- bump minor on new results, major on
  breaking math changes.

**8.  Licensing**

| **Asset** | **Licence**                                                                          |
|-----------|--------------------------------------------------------------------------------------|
| **Code**  | MIT or BSD‑3 (OSI‑approved, permissive)                                              |
| **Docs**  | CC‑BY‑4.0                                                                            |
| **Data**  | Mirror upstream licences; if not compatible, store only checksums + download script. |

**9.  External Standards & References**

| **Purpose**               | **Reference**                                                   |
|---------------------------|-----------------------------------------------------------------|
| **FAIR data principles**  | Wilkinson *et al.*, *Sci. Data* (2016)                          |
| **Reproducible research** | Sandve *et al.*, *PLoS Comp Bio* 8 Rules (2013)                 |
| **Software citation**     | Smith *et al.*, JOSS guidelines                                 |
| **Python packaging**      | PyPA specifications (pyproject.toml, PEP 517/518)               |
| **CI best practice**      | Wilson *et al.* "Good Enough Practices in Scientific Computing" |

**10.  Makefile / CLI targets (developer ergonomics)**

\# Create and sync dev environment

env:

pip install \--prefer-binary numpy scipy matplotlib

pip install \--no-build-isolation -e .\[dev\]

\# Run fast tests & smoke fits

test:

pytest -q

python scripts/analyze_sparc_galaxies.py \--max-galaxies 3

\# Full reproducibility package (heavy)

bench: env

python scripts/analyze_sparc_galaxies.py \--max-galaxies 175

python scripts/run_supernova_likelihood.py

python scripts/run_cmb_hook.py

\# Build & serve docs

docs:

mkdocs serve

**Deliver this guide to the implementing AI**

*The bullet‑proof outcome:* a repo that any reviewer can clone, test,
and audit without hand‑holding, meeting modern open‑science and
engineering norms.
