Got it. I'll search for publicly available experimental
datasets---preferably in their rawest or minimally processed form---from
any type of quantum or high-energy physics experiment. I'll focus on
finding data that has minimal reliance on Standard Model-based
calibrations or theoretical assumptions, including legacy and
educational/test beam data where possible. I'll let you know once I've
gathered and assessed the most promising sources.

**Searching for Truly Raw Physics Data**

**Standard Model Bias in Experiment Data**

Modern physics experiments **heavily rely on the Standard Model (SM)**
at virtually every stage -- from calibrating instruments to analyzing
results. This means most publicly available data are *"contaminated"* by
SM-based processing and assumptions. For example, in the Muon *g--2*
experiment, the magnetic field was calibrated using known constants and
proton magnetic properties derived from the SM. The measured field (via
proton NMR probes) had to be converted to an equivalent **free-proton
Larmor frequency** using CODATA values. Equation (22) in the E821 final
report shows how the **proton precession in a water sample** was related
to a free-proton precession using reference data (Mohr's constants). In
short, the *absolute B-field scale* is tied to well-established SM
quantities like the proton's magnetic moment. Any anomaly measurement
aμa\_\mu then uses this field plus known muon/proton magnetic moment
ratios (from muonium spectroscopy) to compute
aμ:contentReference\[oaicite:2\]index=2:contentReference\[oaicite:3\]index=3a\_\mu:contentReference\[oaicite:2\]{index=2}:contentReference\[oaicite:3\]{index=3}.
All these inputs -- \$e\$, \$\hbar\$, \$m\_\mu/m_p\$, etc. -- come from
prior experiments interpreted via standard physics.

Likewise, **detector systems** are tuned using SM processes. High-energy
collider detectors (ATLAS/CMS at LHC, for instance) perform *in-situ
calibrations* with well-known SM signals. For example, the ATLAS
electromagnetic calorimeter calibration uses test beams and Monte Carlo
simulations *and then is refined using known physics decays (invariant
masses of particles)*. They explicitly mention using "**well-known
physics processes, such as the invariant mass of particle decays**" to
adjust calibrations -- e.g. using the \$Z\$ boson mass peak to set the
energy scale for electrons/photons. Similarly, tracking detectors and
triggers are optimized for SM-like events: triggers are designed to
capture signatures of known particles (high-\$p_T\$ leptons, jets from
QCD, etc.), which means *data that don't fit SM patterns might be
down-weighted or even discarded at the trigger level*. Even *background
subtraction* relies on SM predictions -- e.g. pile-up corrections assume
known interaction rates and cross-sections, beam backgrounds are
estimated from known decay chains, etc. All these steps bake in the
assumption that **Standard Model physics is the baseline**.

**Data Analysis and Theory Bias**

The **analysis methods and theoretical definitions** also incorporate SM
expectations. For instance, likelihood fits often assume the probability
distributions predicted by SM-based simulations. If one is searching for
new effects but using an SM-shaped likelihood, subtle deviations might
be glossed over. Bayesian analyses might use **priors based on previous
SM measurements** or theoretical constraints. Even fundamental constants
used in data reduction are part of the established framework -- \$e\$
(elementary charge), \$\hbar\$, \$c\$ are all defined or measured in
contexts assuming standard physics. (In fact, many "fundamental
constants" are measured by experiments that validated QED, general
relativity, etc., so their values inherently assume those theories are
correct to high precision.)

A telling example is how the **muon \$g-2\$ anomaly is defined**. The
quantity \$a\_\mu \equiv (g\_\mu - 2)/2\$ is measured by comparing the
muon's spin precession frequency to its cyclotron frequency in a
magnetic field. *However, the reference point "2" in \$g-2\$ comes from
Dirac's theory.* We say "anomalous" magnetic moment because Dirac's
relativistic QM predicts \$g=2\$ exactly for a pointlike spin-½
particle. The deviation arises from quantum loops (QED, etc.). So by
definition **\$a\_\mu\$ is measured relative to the Standard Model's
lowest-order expectation**. Furthermore, extracting \$a\_\mu\$ from raw
observables requires inserting SM constants: one formula rewrites
\$a\_\mu\$ in terms of measurable frequencies and known ratios as :

aμ  =  ge/2    μp/μe    mμ/me    ωa⟨ωp⟩ ,a\_\mu \\=\\ \frac{g_e/2 \\\\
\mu_p/\mu_e \\\\ m\_\mu/m_e \\\\ \omega_a}{\langle\omega_p\rangle}\\,

where \$\omega_a\$ is the muon's anomalous precession frequency, and
\$\langle \omega_p \rangle\$ is the average proton Larmor frequency in
the ring. Here \$g_e/2\$, \$\mu_p/\mu_e\$, \$m\_\mu/m_e\$ are
*independently measured* constants, but **each of those measurements
assumed "ordinary" physics** (QED, atomic theory, etc.) in their
methodology. For example, \$\mu_p/\mu_e\$ (the proton-electron magnetic
moment ratio) is known to 3 ppb from spectroscopy of muonium and
hydrogen, *relying on QED calculations to extract the magnetic moments*.
In essence, no number in that equation is truly theory-agnostic. Even
the **Bohr magneton** (\$\mu_B = e\hbar/2m_e\$) used as a unit is a
product of SM quantities. The collaboration then compares \$a\_\mu\$ to
the SM *predicted* value, which includes electroweak and hadronic loop
corrections. The "discrepancy" gets headlines, but it's *by construction
a difference between experiment and the Standard Model calculation*.

**Is "Truly Raw" Data Available?**

Given the above, finding *completely raw, SM-unbiased data* is extremely
challenging. By the time data is made public or published, it's usually
been calibrated, normalized, or processed with standard physics
assumptions. **Large collider experiments**, for instance, do not
release **Level-0 raw detector signals** to outsiders -- and even
internally, individual physicists rarely access totally raw data without
reconstruction software. The CERN Open Data policy explicitly states
that making full raw datasets public is *"not practically possible"* due
to the complexity and required context. Instead, they release **Level-3
"reconstructed" data**, where calibrations and basic event
reconstruction have already been applied. (Even within CMS/ATLAS, only
central production teams handle raw detector hits; end-users get
calibrated data.) As the policy notes, using raw data would demand
detailed knowledge of detector geometry, electronics response, trigger
conditions, etc., which is why *"general direct access to raw data is
not even available to individuals within the collaboration"*. In short,
*by the time data is human-friendly, it's been through the SM
machinery*.

That said, there **are some datasets and measurements closer to "raw"**
that you might leverage for a new theory test:

- **Muon \$g-2\$ Time-Series (wiggle) Data:** While the final \$a\_\mu\$
  involves many SM inputs, the *direct observables* in the Muon \$g-2\$
  experiment are more model-independent. Specifically, they measure the
  **decay positron counts vs time** in the storage ring. This famous
  **"wiggle plot"** (see below) is essentially a histogram of
  high-energy positron counts as a function of time modulo the muon
  cyclotron period. The raw time spectrum will show an exponential decay
  (muon lifetime) modulated by a cosine oscillation (spin precession).
  *That* dataset -- time stamps or binned counts of decay events -- is
  largely free of theoretical bias. It's just what the detectors saw
  (after basic noise subtraction and event selection). The frequency of
  the oscillation \$\omega_a\$ can be extracted by fitting this time
  distribution. If one had access to the **un-blinded time-series
  data**, one could apply a new physics fit to see if the oscillation
  frequency or decay curve deviates from the SM expectation. In
  practice, the Muon \$g-2\$ team hasn't publicly posted the full raw
  time-series for each fill (they have many TB of raw records). However,
  their papers often show summary plots. For example, **David Sweigart's
  PhD thesis** (FNAL, 2020) contains a blinded wiggle plot with the fit
  overlay. An excerpt is shown in **Figure 1** below. This represents a
  *nearly raw* dataset: just **time vs count**. If you could get the
  underlying numbers (perhaps by contacting the collaboration or via a
  data release), you could reanalyze with your own model (e.g. testing
  for any deviation in modulation form, alternative decay law, etc.).
  Keep in mind, even here the **time axis is calibrated** to an atomic
  clock and the **selection of "high-energy" positrons** depends on
  energy calibration of the calorimeters (which was done using standard
  electromagnetic shower physics). So it's *not 100% theory-free*, but
  it's close to the direct experimental observation.

*Figure 1: Example "wiggle plot" from the Muon \$g!-!2\$ experiment. It
shows the number of decay positrons detected (above an energy threshold)
versus time, modulo a certain period. The exponential fall-off is due to
muon decays, and the oscillation is the muon spin precessing (frequency
\$\omega_a\$). This particular plot is from a blinded analysis (for
illustrative purposes), but it reflects the raw measurement concept.*

- **Frequency Ratio Measurements:** Instead of relying on SM constants
  directly, one can look at **dimensionless ratios** that experiments
  measure. In \$g-2\$, a key intermediate result is the ratio of two
  frequencies: R=ωaωpR = \frac{\omega_a}{\omega_p} (muon spin precession
  frequency to proton Larmor frequency in the same magnetic field). This
  ratio *itself* is a directly measured quantity with minimal theory
  input. The E821 Brookhaven experiment reported
  Rμ=0.0037072064(20)R\_{\mu} = 0.0037072064(20) (for \$\mu\^+\$ and
  \$\mu\^-\$ combined). To get \$a\_\mu\$, they then used the
  independently measured \$\mu\_\mu/\mu_p\$ (denoted \$\lambda\$) and
  the formula \$a\_\mu = \frac{R}{\lambda - R}\$. If your new framework
  alters the relationship between magnetic moment and mass or other
  constants, you might predict a different \$\lambda\$ or interpret
  \$R\$ differently. The **quantity \$R\$ is an almost "raw"
  experimental invariant** -- it's the ratio of two precession
  frequencies, which cancels a lot of systematics. (Still, \$\omega_p\$
  is measured via NMR probes that were absolutely calibrated using the
  SM value of the proton magnetic moment in water. So even \$R\$ isn't
  fully free of SM input; the *shielding correction* and water probe
  calibration come from conventional physics.)

- **Electron \$g-2\$ (Penning Trap) Data:** The most precise test of QED
  is the electron magnetic moment. In that experiment, one electron is
  stored in a Penning trap and both its cyclotron frequency and spin
  flip (Larmor) frequency are measured. The **ratio of those two
  frequencies** yields \$g/2\$ directly (to first order). Gabrielse et
  al. have measured \$g_e\$ to better than 1 part in \$10\^{13}\$. If
  you could access their *raw frequency measurements*, you'd have an
  extremely clean comparison for a new QED theory. The nice thing is
  that a single trapped electron's frequencies are measured in terms of
  an oscillator clock -- not requiring an external particle's
  properties. In practice, though, these experiments still rely on
  calibrations (the trap's magnetic field is measured via a current or a
  proton NMR probe, etc.). The data itself (a series of resonance curves
  and frequency counts) might not be publicly archived, but their
  published result **\$a_e\$** and uncertainty are available. One could
  in principle rederive the fine-structure constant \$\alpha\$ from
  those data using a new theory's formula and compare to other
  \$\alpha\$ measurements.

- **Open Data Repositories:** Many physics collaborations have started
  releasing **open data** after an embargo period. For example, **CERN's
  Open Data Portal** provides **LHC data** (CMS and ATLAS primarily) at
  various processing levels. In 2014--2016, CMS released around 0.5 fb⁻¹
  of 7 TeV collision data (2010--2011 runs), and more in subsequent
  releases -- on the order of **hundreds of terabytes** of events. These
  are in ROOT format, with reconstructed physics objects (muons,
  electrons, jets, etc.) and the necessary calibration constants. While
  not *truly raw*, this gives you a chance to apply alternative analysis
  techniques. For instance, one could try to **reconstruct certain
  observables without certain SM assumptions** -- e.g. look at track
  curvature without applying the nominal magnetic field value (testing
  if momentum calibration has a systematic bias), or summing energies in
  calorimeters without the default calibration to search for anomalies.
  Be aware that **even this data has gone through the full standard
  reconstruction pipeline** (alignment, calibration, particle
  identification trained on SM signatures). So any *new physics effect*
  would likely manifest as an *unexpected distribution of these
  reconstructed objects*, rather than new objects outright. Still, you
  can do **reinterpretations**: many published searches provide
  higher-level data (histograms, likelihoods on HEPData) that you could
  compare against new-model predictions. If your question is
  specifically about *spacetime geometry*, you might, for example, look
  at whether kinematic distributions might fit better under a different
  invariant metric -- but that's a very involved reanalysis.

- **Neutrino and Cosmic-Ray Data:** Some "raw-ish" data can be found in
  areas like neutrino physics or cosmic-ray observations. For instance,
  the **Super-Kamiokande** experiment has released event datasets of
  atmospheric neutrino interactions, and **IceCube** shares data for
  certain high-energy events. These are typically counts of events
  versus reconstructed energies/angles. The reconstruction (e.g.
  inferring a neutrino energy from photomultiplier signals in water/ice)
  uses standard particle interaction models, though, so again not
  theory-free. **Cosmic ray air shower arrays** (Auger, Telescope Array)
  sometimes publish detailed data about extensive air showers -- number
  of secondary particles, timing, etc. Most of their analysis assumes
  well-known physics of air interactions (QED/QCD cascades). If your new
  framework predicts, say, different cross-sections at ultra-high
  energy, one could try to compare with these observations. But you
  would need to propagate your model through a simulation of the
  detector anyway.

- **Astrophysical Observations:** Perhaps data least "touched" by the
  Standard Model of *particle physics* are observations in astronomy.
  These still require calibration (telescope sensitivities, etc.), but
  they're not **calibrated to SM particle processes** in the same way.
  For example, the **COSMIC MICROWAVE BACKGROUND (CMB)** data from
  Planck or WMAP satellites -- the raw maps of temperature anisotropies
  -- are calibrated using thermodynamic temperature references and
  instrument models, not directly by SM particle interactions. If your
  new spacetime geometry or quantum framework predicts subtle
  differences in these anisotropy patterns or spectral distortions, you
  could use the publicly available CMB maps (Planck data is online) to
  test that. However, extracting meaningful results still requires the
  standard cosmological model or an alternative -- you'd have to
  separate astrophysical effects from new-physics effects.

- **Gravitational-Wave Strain Data:** Projects like LIGO/Virgo provide
  open time-series data of the detector strain signals. This is \*raw
  interferometer output calibrated to strain (dimensionless)\*\*. The
  calibration here is classical (laser interferometry, mirror motion),
  not based on SM quantum field theory. If your new spacetime geometry
  affects gravitational wave propagation or birefringence, one could
  analyze these time-series for discrepancies from General Relativity
  templates. Again, this tests new gravity more than new QED, but it is
  an example of data you can download and analyze free of SM *particle*
  biases. (Calibration of LIGO data does assume General Relativity for
  the signal shape in some cross-checks, but the raw strain vs time is
  fairly direct.)

**In summary, truly raw, completely SM-uncontaminated data is extremely
rare to find online**. Almost all high-precision experimental data have
some layer of SM-based correction or interpretation. The best you can do
is obtain *the most elementary observables measured* -- time,
frequencies, detector counts, etc. -- **along with documentation of how
they were calibrated** -- and then apply your new theoretical framework
to interpret those observables from scratch. Some concrete steps:

- **Leverage Open Data**: Download LHC open datasets (reconstructed
  level) and focus on observables that your theory might affect. For
  instance, if your new QED predicts a different angular distribution in
  certain decays, you can look at LHC decay data. Keep in mind, you may
  have to *undo* some calibrations or at least be aware of them. The CMS
  open data comes with simulation and software environments to help
  users reproduce analyses. You could rerun parts of reconstruction with
  tweaked parameters if needed.

- **Request or Find Supplemental Data**: Check if collaborations have
  released supplemental data in repositories. Experiments like Muon
  \$g-2\$ sometimes provide **averaged datasets** or plots in papers;
  you might contact the authors to request underlying numbers. Since the
  *g--2* collaboration is now complete (data taking finished), they
  might be more open to distributing some forms of data for novel
  analysis, especially if approached collaboratively.

- **Smaller-scale experiments**: Look for any smaller experiments that
  test QED/QM where data might be more directly accessible. For example,
  the **Proton charge radius measurements** (electron vs muonic hydrogen
  spectroscopy) had some contentious results -- those raw spectra might
  be interesting if a new theory alters energy levels. Similarly,
  **quantum optics experiments** (Bell test data, interferometry
  fringes) sometimes share raw event streams or timestamps for public
  tests of local realism, etc. Those typically only assume basic quantum
  mechanics or classical optics in analysis, and if your new QM
  framework differs, you could reanalyze the raw coincidence counts.

- **Fundamental constants**: Consider data from NIST or national
  metrology labs. They publish results of precision measurements (often
  with raw data in supplementary info or at least detailed error
  budgets). For example, NIST's measurement of the **fine-structure
  constant via atom interferometry** or the **electron \$g-2\$
  measurement from Harvard** might have raw frequency count data or
  interference fringe data available in some form. While these were
  interpreted with standard theory to extract \$\alpha\$, you could
  apply a different interpretation.

**Final Thoughts**

Your observation is astute: **the interconnected nature of modern
experiments and theory means completely model-independent data is
scarce**. The community acknowledges this; for instance, luminosity at
colliders is calibrated using theory (W/Z cross-sections calculated in
SM) *and then cross-checked by independent methods*. This is a way of
*mitigating* model-dependence by not relying on a single assumption. If
you want to validate a revolutionary QM/QED framework, you might
ultimately need to **propose dedicated experiments** or re-run certain
calibrations under new assumptions. In the meantime, you can start with
existing "least processed" data: time-series, frequency ratios,
unfiltered event counts. Use those to see if your theory can reproduce
known phenomena **without** the built-in SM tweaks.

Do keep in mind that any claim of new physics will be scrutinized to
ensure it's not an artifact of misunderstood calibrations. This was one
reason the LHC experiments were initially cautious about open data -- a
naive reanalysis can find "bumps" or deviations that are actually just
unaccounted detector effects. So whatever data you use, be sure to
**understand the experimental setup in detail** (detector response,
selection biases) so you can separate genuine new-physics signals from
calibration residuals.

In summary, **yes, there are data sources you can use online**, but
**no, virtually none of them come totally untouched by Standard
Model-based processing**. The closest alternatives are to look at
**basic measured quantities** (time, counts, ratios) provided by
experiments and then apply your new framework's interpretation on those.
Be prepared for significant work untangling the webs of calibration --
but that's exactly the task if one hopes to confirm a fundamentally new
spacetime geometry or quantum theory using existing data.

**Sources:**

- CERN Open Data Policy (levels of data and raw data limitations)

- ATLAS Open Data documentation (calibration using known processes)

- Muon \$g-2\$ experimental method and use of fundamental constants

- Muon \$g-2\$ field calibration via NMR (probe calibration using CODATA
  values)

- Luminosity calibration at LHC depending on SM cross-sections (W/Z)

- Example of open data release from CMS (100 TB of collision data in
  2011, in analysis format)

- ParticleBites summary of Muon \$g-2\$ (includes wiggle plot and
  discussion of inputs)

- Brookhaven E821 final report (definition of \$a\_\mu\$ and frequency
  ratio \$R\$)
