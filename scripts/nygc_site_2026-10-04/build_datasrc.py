from sg import *
import re

t = read("data-sources.html")
i = t.index("<h2>3. 1000 Genomes Reference Panel</h2>")
j = t.index("<h2>", i + 10)
new = """<h2>3. 1000 Genomes Reference Panel (NYGC 30x high-coverage, current)</h2>
        <p>Used throughout the global ancestry analyses (PCA, F<sub>ST</sub>, PBS, ADMIXTURE) as the external comparison population since 2026-10-02. Source: 1000 Genomes high-coverage GRCh38 call set of the New York Genome Center (3,202 samples; all 22 chromosome files verified against the official md5 manifest), on DRAGEN under <code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x/</code>.</p>
        <table>
            <tr><th>Superpopulation</th><th>Populations included</th><th>N used</th></tr>
            <tr><td>African (AFR)</td><td>YRI, LWK, GWD, MSL, ESN, ACB, ASW</td><td>893</td></tr>
            <tr><td>European (EUR)</td><td>CEU, GBR, FIN, IBS, TSI</td><td>633</td></tr>
            <tr><td>South Asian (SAS)</td><td>GIH, PJL, BEB, STU, ITU</td><td>601</td></tr>
            <tr><td>East Asian (EAS)</td><td>CHB, JPT, CHS, CDX, KHV</td><td>585</td></tr>
            <tr><td>Americas (AMR)</td><td>MXL, PUR, CLM, PEL</td><td>not used in the current runs</td></tr>
            <tr style="font-weight:600;background:#f5f5f5;"><td colspan="2">Total used (EUR + EAS + SAS + AFR)</td><td>2,712</td></tr>
        </table>
        <p style="font-size:0.9em;color:#888;">
            Genome build GRCh38, the same as the Uzbek cohort. Reference samples were not filtered for relatedness. The reference panel has no Central Asian populations; adding some is an open item (see <a href="steps/reference_panels_expansion.html">reference panels</a>).
            The earlier panel (1000G phase-3-era file <code>ALL.chr*.shapeit2_integrated_v1a.GRCh38</code>, 2,548 samples: AFR 671, EUR 522, SAS 492, EAS 515, AMR 348) is deprecated because genotypes at a few dozen sites in that file depend on the sequencing batch; see <a href="steps/deprecated_old_runs.html">deprecated runs</a>.
        </p>

        """
t = t[:i] + new + t[j:]
t = wsub(t, "<td>Samples in final analysis cohort</td><td>1,047 (after the full QC pipeline", "<td>Samples in final analysis cohort</td><td>1,256 in the expanded cohort used for the current population-genetic analyses (original cohort: 1,047; see <a href=\"steps/cohorts_and_sample_sets.html\">cohorts and sample sets</a>) (original pipeline: after the full QC pipeline")
write("data-sources.html", t)
print("datasrc ok")
