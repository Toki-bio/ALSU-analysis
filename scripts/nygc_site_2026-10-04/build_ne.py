from sg import *

t = read("steps/arg_methods_pilot.html")
blk = """
<h2 style="border-bottom:3px solid #16a34a;">IBD-based recent effective population size: attempted 2026-10-09, no usable estimate</h2>
""" + box("key", "<strong>Result.</strong> Two genome-wide attempts to estimate recent Ne(t) from IBD segments (hap-ibd + IBDNe) did not give a usable estimate. The first, on phased <em>imputed</em> data, failed basic sanity checks. The second, on phased <em>typed</em> array SNPs, passed the segment-count check but the IBDNe model does not fit (fit statistic 27.4 for the real data against 0.2&ndash;0.7 in bootstrap resamples) and the curve is implausible in the recent past. <strong>No Ne(t) or divergence-time statement is made from these data.</strong> The item is closed as not feasible with this data set; the earlier chr22 ARG pilot remains a method demonstration only.") + """
<h3>Attempt 1: phased imputed genotypes (1,245 samples, MAF &ge; 0.05)</h3>
<ul>
<li>hap-ibd found 8.4 million IBD segments &ge; 2 cM, about 11 per pair of people; the count per chromosome ranged from 4,700 (chr18) to 1.43 million (chr14), with no relation to chromosome length. For an outbred sample of this size about one segment per pair genome-wide is expected.</li>
<li>IBDNe returned Ne of 6.7&times;10<sup>17</sup> at generation 0 with confidence limits up to 10<sup>100</sup>.</li>
<li>Cause (likely, not tested further): imputed haplotypes are mosaics of reference haplotypes, so people who copied the same reference haplotype appear identical by descent over long stretches; IBD detection on imputed data is known to be unreliable for this reason. This attempt is not used anywhere.</li>
</ul>
<h3>Attempt 2: typed array SNPs, Beagle 5.x phasing (27Feb25 build), hap-ibd, IBDNe</h3>
<ul>
<li><strong>Data:</strong> typed GSA genotypes (<code>FULL_QC_FINAL</code>), 1,243 of the 1,245 final-cohort samples that remain after dropping one of each second-degree pair; SNPs with MAF &ge; 1% and call rate &ge; 95%; phased per chromosome without a reference panel (<code>impute=false</code>); GRCh38 genetic maps (Browning lab); hap-ibd <code>min-seed=2 min-output=2</code>; IBDNe with <code>mincm=2</code>, 80 bootstraps. Script <code>ibdne_typed.sh</code>, output <code>/staging/ALSU-analysis/admixture_analysis/pophistory_typed_full/</code> on DRAGEN.</li>
<li><strong>Segment counts passed the sanity check:</strong> 471,100 segments, 0.61 per pair, roughly proportional to chromosome length (chr21 8,700; chr1 37,600). One outlier chromosome (chr10, 59,800) and 16 high-IBD regions were excluded by IBDNe&rsquo;s own filter (335,582 segments analysed).</li>
<li><strong>The model fit failed:</strong> IBDNe reports a fit of 27.4 for the real data against 0.2&ndash;0.7 for its bootstrap samples; the bootstrap intervals therefore do not describe the uncertainty of the original fit.</li>
<li><strong>The curve is implausible:</strong> Ne falls from 3.6&times;10<sup>7</sup> at generation 0 to 1.3&times;10<sup>5</sup> at generation 13 and then rises to about 5&times;10<sup>5</sup> at generations 28&ndash;36. A drop of two orders of magnitude within ten generations followed by a rise is not a demographic history; it is what the method produces when its assumptions fail.</li>
</ul>
""" + img("nygc_ibdne_typed.png", "IBDNe output for the typed-SNP attempt, shown for transparency only. Do not read it as a population-size history.") + """
<h3>Why the method does not fit this cohort</h3>
<ul>
<li>IBDNe assumes a single randomly mating population. The Uzbek cohort is a mixture of ancestries (<a href="step11.html">step 11</a>; <a href="reference_panels_expansion.html">Kazakh reference</a>); recent admixture creates excess long IBD and distorts the most recent generations.</li>
<li>With about 450,000 typed SNPs and phasing from 1,243 unrelated people, switch errors break long segments, and IBD in the last 5&ndash;10 generations is the most affected.</li>
<li>The cohort is not a random sample of a single population (a clinical cohort drawn from several regions and nationalities).</li>
</ul>
<h3>What would be needed</h3>
<p>A defensible Ne(t) or divergence-time estimate would need high-coverage genomes (for PSMC/MSMC2-type methods, not available here), or dense phased data from a population-specific reference panel, or a method that models admixture explicitly. None of these is available for this project; divergence times between populations cannot be estimated from array-based data. The question is closed for this data set.</p>
"""
t = put(t, "ne", blk, after='<div class="content">')
t = t.replace("Pilot — 2026-09-20", "Pilot 2026-09-20; IBD attempt 2026-10-09") if "Pilot — 2026-09-20" in t else t
write("steps/arg_methods_pilot.html", t)
print("ne ok")
