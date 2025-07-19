**I received suggestions to make it easier to get external validation.
Is it possible to do these things?**

**Suggested next commits**

| **Priority** | **Action**              | **One‑liner implementation**                                                                                      |
|--------------|-------------------------|-------------------------------------------------------------------------------------------------------------------|
| 🔑           | **Add a tiny manifest** | python -m hashlib sha256 data/sparc_database/\*.dat \> data/manifest_sha256.txt and commit that text file.        |
| 🔑           | **Document LFS pull**   | In README's *Installation* section, add: bash git lfs install && git lfs pull                                     |
| ⚙️           | **Pin environment**     | pip‑compile \--generate‑hashes requirements.in → commit requirements.txt with hashes, or export conda env export. |
| 🏗️           | **Bootstrap CI**        | Drop .github/workflows/audit.yml that runs pytest (or your audit scripts) on every push.                          |
| 🧪           | **Sample notebook**     | Provide examples/reproduce_one_galaxy.ipynb that loads a single SPARC file and plots the legacy vs rebuilt fit.   |
