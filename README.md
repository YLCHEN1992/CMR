![Workflow Diagram](https://github.com/YLCHEN1992/CMR/blob/main/image/background.jpg "Workflow Diagram")
# CMR #
**<span style="font-size:larger;">Circle Mendelian Randomization is a batch Mendelian randomization analysis program which filters genetic instrumental variables using fixed linkage disequilibrium information indicators and ranked significance values and that is different from other methods of fixing significance threshold.</span>**
# Easy to use #
**<span style="font-size:larger;">./CMDnohup.sh</span>**
**<span style="font-size:larger;">UPDATING更新中... ...</span>**

# FQA #
🗣️ Response to Reviewers & Frequently Asked Questions
This document addresses common questions and invaluable feedback from peer reviewers regarding our Mendelian Randomization (MR) analysis methodology and results.

1. On the Limitation to European Ancestry
Reviewer's Comment: "This study is limited to European ancestry, which may restrict the generalizability of the findings. If available, incorporating data from other ancestries could provide a broader perspective and enhance the robustness of the results. If such data are unavailable, the limitation should be explicitly acknowledged, with a discussion of its implications."

<details> <summary><b>Our Response</b></summary>
We sincerely thank the reviewer for this critical comment regarding the ancestral limitations of our study. The genetic data utilized in this analysis was sourced exclusively from the IEU GWAS public database, which contains predominantly European-ancestry populations.

Our analytical approach specifically selected Finnish and European population cohorts for the exposure data to maintain genetic ancestry consistency with the outcome dataset, which was also derived from European populations. We recognize that incorporating genetically mismatched cohorts (e.g., using East Asian exposure data with European outcome data) can introduce significant bias and increased error in MR estimates due to fundamental differences in linkage disequilibrium (LD) structure, allele frequencies, and population-specific genetic architectures.

Furthermore, acquiring well-powered GWAS summary statistics for our specific outcome measures from diverse ancestral backgrounds remains challenging in the current genomic landscape. We have explicitly acknowledged this limitation in the discussion section of our manuscript, noting that our findings should be interpreted within the context of European populations until replicated in more diverse cohorts.

We extend an open invitation to the scientific community: should the reviewer or other researchers have access to appropriate GWAS data for the same outcome trait in non-European populations, we would be exceptionally eager to collaborate and perform complementary analyses to enhance the trans-ancestry robustness and generalizability of our findings.

</details>
2. On the Type of IVW Method Used
Reviewer's Comment: "Please clarify which type of IVW method was used in the analysis. Indicate whether a fixed-effect or random-effect IVW model was applied to ensure transparency and reproducibility."

<details> <summary><b>Our Response</b></summary>
Thank you for prompting this important methodological clarification. In our primary analysis, we employed the fixed-effect inverse-variance weighted (IVW) model as our main analytical approach.

We wish to clarify that the standard implementation of the IVW model in Mendelian Randomization is typically performed without an intercept term, as was done in our study. The random-effects variant of IVW, while sometimes considered in the presence of significant heterogeneity, is conceptually and mathematically more aligned with extended MR models such as MR-Egger, which specifically account for over-dispersion or directional pleiotropy.

Our choice of the fixed-effect IVW model was based on its established performance characteristics and widespread acceptance as the primary method for MR analysis when using multiple genetic instruments.

</details>
3. On Performing Leave-One-Out Analysis
Reviewer's Comment: "I suggest performing a leave-one-out analysis to evaluate the sensitivity of the results. This method will help determine whether any single SNP has a disproportionate influence on the findings."

<details> <summary><b>Our Response</b></summary>
We appreciate this suggestion regarding sensitivity analysis. In our automated MR analysis pipeline, we implemented a specialized form of leave-one-out (LOO) analysis, though we employed it primarily as an iterative filtering mechanism for instrumental variable (SNP) selection rather than as a post-estimation diagnostic tool.

Our algorithm systematically evaluates the influence of individual SNPs by iteratively removing variants that are identified as outliers based on their distance from the IVW model's estimate. Through this process, we typically prune away a substantial portion (approximately 50%) of the most influential SNPs to ensure the robustness of our final results. This approach allows us to maintain a set of genetic instruments that provide consistent causal estimates.

For evaluating the sensitivity of our final results to individual variants, we relied primarily on the Cochran's Q statistic derived from the meta-analysis package integrated within our MR pipeline. This test for heterogeneity effectively serves a similar purpose to a traditional LOO sensitivity analysis by identifying whether the overall causal estimate is disproportionately driven by any single genetic variant or a small subset of variants.

</details>
4. On Multiple Testing Correction (FDR)
Reviewer's Comment: "To control for false positives due to multiple testing, applying False Discovery Rate (FDR) correction is recommended. This adjustment will help ensure the reliability of the reported associations."

<details> <summary><b>Our Response</b></summary>
We thank the reviewer for highlighting this important statistical consideration regarding multiple testing correction. We recognize that conducting numerous statistical tests within a single analysis framework can inflate the family-wise error rate and potentially lead to false positive findings.

In our analytical approach, we made a conscious decision to present uncorrected p-values for several reasons. First, the interpretation of MR results often focuses on specific causal hypotheses for individual exposure-outcome pairs rather than simultaneous inference across all tested relationships. Second, different exposure-outcome combinations represent distinct biological questions with varying prior probabilities of association, making uniform correction methods potentially overly conservative.

We have provided complete results with uncorrected p-values, allowing readers to assess the strength of evidence for each association by examining the magnitude of effects and precision of estimates, in addition to statistical significance. For readers who prefer a more conservative approach, we have included complete results that enable the application of various multiple testing correction methods, including Bonferroni correction (0.05/n) or Benjamini-Hochberg FDR correction, according to their specific analytical preferences.

We acknowledge this approach in the limitations section of our manuscript and encourage readers to interpret our findings with appropriate consideration of the multiple testing issue.

</details>
5. On Adding MR Steiger Filtering
Reviewer's Comment: "Adding MR Steiger filtering is recommended to test the direction of causality for each instrumental variable on exposure and outcome. According to the Steiger filtering assumption, an instrumental variable is valid if it explains more variance in the exposure than in the outcome. Including this step would strengthen the causal inference of the study."

<details> <summary><b>Our Response</b></summary>
We appreciate the expert's suggestion regarding the implementation of MR Steiger directional analysis. We thoroughly considered incorporating this method during our analytical planning phase but encountered specific methodological constraints that limited its applicability to our study.

The MR Steiger method requires estimation of the proportion of variance explained (R²) for each genetic instrument in both the exposure and outcome traits. This calculation typically employs the formula R² = 2 × EAF × (1 - EAF) × β², where EAF is the effect allele frequency and β is the estimated effect size. However, accurate implementation of this approach requires standardized effect estimates and compatible scaling between exposure and outcome measurements.

In our study, the outcome variable's GWAS data was not standardized in a manner that would allow direct comparison of variance explained estimates between exposure and outcome traits. This limitation diminishes the interpretative power and validity of Steiger's Z-test in our specific analytical context.

Additionally, we note that MR Steiger directional analysis is primarily designed for reverse assessment of unidirectional analysis results. Our study employed different sets of SNP loci as instrumental variables for bidirectional MR analyses, which represents a different analytical approach compared to studies that use the same instruments for both directions and then apply Steiger testing to determine causal direction.

Given these methodological considerations, we determined that assessment of causal direction based on the strength and consistency of associations in both directions represented a more appropriate approach for our specific study design and data structure.

</details>
6. On the Management of Linkage Disequilibrium (LD)
Reviewer's Comment: "Regarding the management of linkage disequilibrium (LD), if possible, I suggest implementing LD analysis by selecting SNPs with R² < 0.001 within a 10,000 kb window as part of the SNP selection process. Additionally, performing linkage disequilibrium score regression (LDSC) would help quantify potential genetic confounding due to LD and improve the robustness of the results. If feasible, this analysis and its assumptions should be included in the methodology."

<details> <summary><b>Our Response</b></summary>
Thank you for these sophisticated suggestions regarding LD management. We carefully considered various approaches to LD handling during our methodological development phase and implemented a strategy that balanced statistical rigor with practical analytical constraints.

Regarding the suggestion to select SNPs with R² < 0.001 within a 10,000 kb window: while this approach represents an extremely conservative threshold that would minimize LD contamination, we found that it would render many exposure-outcome combinations analytically intractable by reducing the number of available instruments below the minimum required for reliable MR estimation. Our primary analysis relies on GWAS summary statistics rather than individual-level genetic data, which limits our ability to calculate precise LD measures based on the actual study samples.

Instead, we implemented an LD clumping procedure that retains the most significant SNP within a 10,000 kb window, which represents a conventional approach in MR studies using summary data. This strategy ensures that our instruments are not strongly correlated while maintaining sufficient statistical power for reliable causal estimation.

Concerning LD score regression (LDSC): we acknowledge that this method provides valuable insights into genetic confounding and heritability. However, LDSC requires specific input formats and assumptions that were not compatible with our analytical framework, which primarily operates on processed GWAS summary statistics. The implementation of LDSC would necessitate access to additional reference data and individual-level genetic information that was not available for our study.

We are confident that our current approach to LD management, which combines conservative clumping parameters with sensitivity analyses, provides appropriate protection against false positive results due to LD structure while maintaining analytical feasibility across the wide range of exposure-outcome relationships examined in our study.

</details>
