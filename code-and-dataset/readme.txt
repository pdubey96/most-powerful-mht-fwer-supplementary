Code and Data for: "Necessary Conditions for the Most Powerful Test in Multiple Hypothesis Testing under Exact Family-Wise Error Rate Control — a Path to Fast Computing"

This repository contains the R source code and datasets required to reproduce the simulations, real-data applications, and figures presented in the manuscript.

==============================================================================
FOLDER STRUCTURE
==============================================================================

The "Codes and Datasets" folder contains six primary R scripts, organized by the section of the paper they support.

------------------------------------------------------------------------------
1. Simulation Studies (Section 5)
------------------------------------------------------------------------------
These scripts reproduce the power comparison results presented in the manuscript. Each script implements Algorithm 1 (Coordinate Descent) for a specific distributional setting and compares it against standard FWER-controlling procedures (Bonferroni, Holm, Hochberg, Hommel, Romano-Wolf).

* sim1_truncated_normal.R
  - Manuscript Section: 5.1 (Testing K=3 Independent Normal Means)
  - Description: Simulates K=3 independent tests where the test statistics follow a truncated normal distribution on [-6, 6].
  - Key Output: sim_5.2.pdf (Power curves for Pi_3 and Pi_any).
  - Data: Synthetic data generated within the script.

* sim2_mixture_normal.R
  - Manuscript Section: 5.2 (Testing Two-Sided Mixture Normal Alternatives)
  - Description: Simulates a two-sided testing problem where the alternative hypothesis follows a symmetric Gaussian mixture: 0.5 N(theta, 1) + 0.5 N(-theta, 1).
  - Key Output: sim_7.1.pdf (Power curves).
  - Data: Synthetic data generated within the script.

* sim3_students_t.R
  - Manuscript Section: 5.3 (Testing Heavy-Tailed Alternatives)
  - Description: Evaluates robustness against heavy-tailed alternatives using Student's t-distributions with varying degrees of freedom (df = 2 to 20).
  - Key Output: sim_v8.1.pdf (Power vs. Degrees of Freedom).
  - Data: Synthetic data generated within the script.

* sim4_beta_pvalues.R
  - Manuscript Section: 5.4 (Testing Parameters of the Beta Density Function)
  - Description: Models p-values directly using a Beta distribution under the alternative, strictly satisfying the theoretical boundedness assumptions.
  - Key Output: sim_6.3.pdf (Power curves).
  - Data: Synthetic data generated within the script.

------------------------------------------------------------------------------
2. Real-World Applications (Section 6)
------------------------------------------------------------------------------
These scripts apply the method to empirical datasets to demonstrate practical utility.

* exp1_bcg_subgroup_analysis.R
  - Manuscript Section: 6.1 (BCG Vaccine Subgroup Analysis)
  - Description: Performs fixed-effects meta-analyses on the BCG vaccine dataset across four different subgroup splits (Allocation, Latitude, Risk, Era). It demonstrates that Algorithm 1 detects a significant effect in the "Low Risk" subgroup that standard methods miss.
  - Data Source: 'dat.bcg' from the 'metadat' R package (installed automatically).
  - Key Output: exp_v2.pdf (Dot plot of p-values and rejection decisions) and CSV summary tables.

* exp2_finance_factor_zoo_analysis.R
  - Manuscript Section: 6.2 (Financial Application: Finding Missed Discoveries)
  - Description: Downloads Fama-French factor data, runs time-series regressions for Momentum, Profitability (RMW), and Investment (CMA), and tests the intercepts (alphas). It assumes a heavy-tailed t(4) alternative distribution for the optimal test.
  - Data Source: Automatically downloads data from the Kenneth French Data Library via the 'frenchdata' R package.
  - Key Output: Console output comparing rejection decisions for the three factors.

==============================================================================
SYSTEM REQUIREMENTS
==============================================================================

The code was developed and tested in the following environment:
- R Version: 4.x or higher
- Operating System: macOS / Windows / Linux (Tested on macOS Sequoia 15.0.1)

Required R Packages:
The scripts automatically check for and attempt to install missing dependencies. The primary packages used are:
- Core Logic: metafor, metadat, frenchdata
- Data Manipulation: dplyr, tidyr
- Visualization: ggplot2, patchwork, scales
- Parallel Computing: parallel (optional, for speeding up BCG panel analysis)

==============================================================================
INSTRUCTIONS FOR REPRODUCTION
==============================================================================

1. Set your R working directory to the "Codes and Datasets" folder.
2. Open any script (e.g., sim1_truncated_normal.R) in RStudio or a text editor.
3. Run the script.
   - Simulations: Will generate PDF plots in the same directory. Note that 'nrep' is set to 120,000.
   - Real Data: Will output results to the console and save CSV/PDF files to the working directory.

==============================================================================
LICENSE
==============================================================================
This code is provided for reproducibility purposes associated with the manuscript submission.