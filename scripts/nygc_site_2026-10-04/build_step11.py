from sg import *

t = read("steps/step11.html")
E = R["eth"]; MT = R["mix_test"]; NN = R["nnls_mix"]
reg = {"1.1": "Tashkent (city)", "1.2": "Karakalpakstan", "1.3": "Andijan", "1.4": "Bukhara", "1.5": "Jizzakh", "1.6": "Kashkadarya", "1.7": "Navoiy",
       "1.8": "Namangan", "1.9": "Samarkand", "1.10": "Surxondaryo", "1.11": "Syrdarya", "1.12": "Tashkent Region", "1.13": "Fergana", "1.14": "Khorezm"}


def ms(v):
    return f"{v[0]:.2f} &plusmn; {v[1]:.2f}"


nat_rows = [[r["group"], r["n"], ms(r["EUR"]), ms(r["EAS"]), ms(r["SAS"])] for r in E["nat_table"]]
b_rows = [[reg.get(r["group"], r["group"]), r["n"], ms(r["EUR"]), ms(r["EAS"]), ms(r["SAS"])] for r in E["birth_table"]]
sc = E["scramble"]; pc = E["positive_control"]; kn = E["kw_nat"]; kb = E["kw_birth"]
fmt = lambda p: "&lt; 1e-300" if p == 0 else f"{p:.1e}"
ctl = [["Self-reported nationality (7 groups, n &ge; 5)", fmt(kn["EUR"]), fmt(kn["EAS"]), fmt(kn["SAS"]), "the real test"],
       ["Same test, nationality labels randomly permuted (200 times)", f"median p {sc['EUR']['median_p']:.2f}; {100 * sc['EUR']['frac_below_0p05']:.0f}% below 0.05", f"median p {sc['EAS']['median_p']:.2f}; {100 * sc['EAS']['frac_below_0p05']:.0f}% below 0.05", f"median p {sc['SAS']['median_p']:.2f}; {100 * sc['SAS']['frac_below_0p05']:.0f}% below 0.05", "negative control: should give p near 0.5 and about 5% below 0.05"],
       ["NYGC reference samples by superpopulation label (EUR/EAS/SAS/AFR)", fmt(pc["EUR"]), fmt(pc["EAS"]), fmt(pc["SAS"]), "positive control: must be highly significant"],
       ["Birthplace among self-reported Uzbeks (13 regions with n &ge; 10)", fmt(kb["EUR"]), fmt(kb["EAS"]), fmt(kb["SAS"]), "second real test"]]

blk = f"""
<h2 style="border-bottom:3px solid #16a34a;">ANCESTRY BY SELF-REPORTED GROUP, KEYED BY SAMPLE ID (added 2026-10-04)</h2>
<p>The earlier statement &ldquo;ethnicity is not significant for ADMIXTURE&rdquo; came from joining nationality codes to samples by line position. It is withdrawn. Below, the phenotype sheet (<code>GWAS от 27.08 - For_Plink.csv</code>) is joined to the K=4 ADMIXTURE result (expanded cohort, seed 11, <code>nygc30x_1256/adm/s11/adm_in.4.Q</code>) <strong>by sample ID only</strong>. {E['n_matched']:,} of the {E['n_cohort']:,} cohort samples have a phenotype row with one of the seven listed nationality codes (1,057 match by ID at all; the 2026 batch largely lacks phenotype rows; nationality &ldquo;Other&rdquo; and &ldquo;Missing&rdquo; are not shown). K=4 components were named by the reference means (EUR-like 0.97 in EUR, EAS-like 1.00 in EAS, SAS-like 0.89 in SAS; AFR-like is about 0 in the cohort and not shown).</p>
{table(["Self-reported nationality", "n", "EUR-like (mean &plusmn; SD)", "EAS-like", "SAS-like"], nat_rows)}
<p>Reading: Russians are 85% EUR-like and Tatars 72%, Koreans 87% EAS-like, Kazakhs and Karakalpaks have a large EAS-like share (50% and 42%), and Tajiks (n = 24) are 56% EUR-like, 22% EAS-like and 22% SAS-like, somewhat more EUR-like and less EAS-like than Uzbeks, but the group is small. Self-reported Uzbeks (n = 836) are 48% EUR-like, 31% EAS-like and 21% SAS-like, with large individual spread (see the figure above). Group means describe groups, not individuals.</p>
<h3>Tests, with controls</h3>
{table(["Test (Kruskal&ndash;Wallis, ID-keyed)", "EUR-like", "EAS-like", "SAS-like", "Purpose"], ctl)}
<p>The real test is overwhelmingly significant for all three components; the permuted-label control behaves as it should (median p about 0.5, 4&ndash;7% below 0.05), and the positive control is significant to the limit of floating-point precision. This confirms that the join is correct and that the result is not an artefact of group sizes.</p>
<h3>Birthplace within self-reported Uzbeks</h3>
{table(["Region of birth", "n", "EUR-like", "EAS-like", "SAS-like"], b_rows)}
<p>{E['n_birth_samples']} self-reported Uzbeks with a region of birth in one of the {E['n_birth_groups']} regions with at least 10 people. Differences are significant (table above), but they are small in absolute terms: regional means of the EUR-like share range from {min(r['EUR'][0] for r in E['birth_table']):.2f} to {max(r['EUR'][0] for r in E['birth_table']):.2f} and of the EAS-like share from {min(r['EAS'][0] for r in E['birth_table']):.2f} to {max(r['EAS'][0] for r in E['birth_table']):.2f}. Regions with n &lt; 30 are noisy. The significance comes mostly from the largest regions (Tashkent city n = {E['birth_table'][0]['n']}); the test does not say which regions differ.</p>

<h3>Cross-check of the K=4 proportions without ADMIXTURE</h3>
<p>The ADMIXTURE K=4 means for the Uzbek cohort (EUR-like 0.50, EAS-like 0.30, SAS-like 0.19, AFR-like 0.00) were checked against an independent method: a non-negative least-squares fit of Uzbek allele frequencies as a mixture of the four reference allele frequencies over {R['ue_dist']['n']:,} SNPs (<a href="step9.html">step 9</a>). The fit gives EUR {NN['EUR'] * 100:.0f}%, EAS {NN['EAS'] * 100:.0f}%, SAS {NN['SAS'] * 100:.0f}%, AFR {NN['AFR'] * 100:.0f}% and predicts the observed Uzbek frequency with r = {MT['r_all']:.3f}. The two methods agree on the shares of EUR and EAS within a few points; the fit gives more weight to SAS (23% vs 19%) and the sum is the same by construction. Neither method proves that these sources are the true ancestral populations: the reference panel has no Central Asian groups, which is why K=7 finds a component that no reference carries.</p>
<h3>PCA percentages</h3>
<p>The table above quotes PC1/PC2 as shares of the first 20 PCs (55.0% / 23.4% for the expanded cohort). <a href="step8.html">Step 8</a> and the report use shares of the first 10 PCs (59.0% / 25.1%). Same PCA, different denominators; neither is a share of total variance.</p>
<h3>One sample with a mostly African-like profile</h3>
<p>One cohort sample has an AFR-like component of 0.90 at K=4 and lies inside the African reference cluster in the PCA (<a href="step8.html">step 8</a>). It has not been excluded from any analysis and has been flagged for a provenance check.</p>
"""
anchor = '<h2 style="border-bottom:3px solid #dc2626;color:#991b1b;">DEPRECATED OLD RUNS'
t = put(t, "step11", blk, before=anchor)
t = t.replace("CURRENT RESULTS (NYGC 30x reference, 2026-10-02)", "CURRENT RESULTS (NYGC 30x reference, 2026-10-02; extended 2026-10-04)")
write("steps/step11.html", t)
print("step11 ok")
