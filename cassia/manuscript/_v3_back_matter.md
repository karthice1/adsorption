## Data and code availability

**Location.** All code, the three PHREEQC databases (with checksums and provenance), the virtual corpora, the result tables and the figure scripts are in the `cassia/` directory of the repository `karthice1/adsorption` (CASSIA 0.3.1).

**Fixed snapshot.** The content of this release is identified by `release/MANIFEST.sha256`: SHA-256 checksums of every code, database, data, result and figure file. `release/RELEASE.md` records:
- the code commit `[COMMIT]`;
- the environment (`requirements-lock.txt`, Python 3.13.17);
- the test log and coverage report;
- the exact commands: `python scripts/run_v2.py --stage all` and `python scripts/run_v3.py --stage all`;
- the list of expected output tables.

A tagged GitHub release (`cassia-v0.3.1`) of this commit, archived in Zenodo, will provide the permanent DOI. The tag and archive are created by the author at submission, because they could not be created from the computational environment; the commit hash and manifest above identify the content unambiguously in the meantime.

**Data templates.** The master-dataset template, schema, study registry, search log and PRISMA flow templates are in `cassia/data/`.

**Third-party files.** The PHREEQC databases are distributed with USGS PHREEQC, a U.S. Government work in the public domain, and are included unmodified.

**Software versions.** CASSIA 0.3.1; phreeqpython 1.6.2 (Vitens); PHREEQC databases minteq.v4.dat, minteq.dat and llnl.dat (revisions in Table S2); adsorption 0.1.0 (Amvro23); pyGAPS 4.3.1.dev693; NumPy 2.5.3; SciPy 1.18.1; pandas 3.0.5; scikit-learn 1.9.1; XGBoost 3.4.1; matplotlib 3.11.2; Python 3.13.17.

No experimental data were generated, and no published experimental data were analysed, in this study.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors. [Author to confirm.]

## CRediT authorship contribution statement

**Karthik [Surname]:** Conceptualization, Methodology, Software, Formal analysis, Validation, Visualization, Writing – original draft, Writing – review & editing. [Add co-authors and roles as appropriate.]

## Declaration of competing interest

The author declares no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the author used Claude (Anthropic) to develop and run analysis code and to draft and revise the manuscript text; its use in the research process is described in Section 2.12. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the published article. [The author should confirm this wording against the target journal's current policy.]

## References

Allison, J.D., Brown, D.S., Novo-Gradac, K.J., 1991. MINTEQA2/PRODEFA2, A Geochemical Assessment Model for Environmental Systems: Version 3.0 User's Manual. EPA/600/3-91/021. U.S. Environmental Protection Agency, Athens, GA.

Amvro23, 2022. adsorption: a Python package for isotherm, kinetic and dynamic adsorption models (MIT licence). https://github.com/amvro23/adsorption [accessed 30 October 2026].

Ball, J.W., Nordstrom, D.K., 1998. Critical evaluation and selection of standard state thermodynamic properties for chromium metal and its aqueous ions, hydrolysis species, oxides, and hydroxides. J. Chem. Eng. Data 43, 895–918. https://doi.org/10.1021/je980080a

Barber, R.F., Candès, E.J., Ramdas, A., Tibshirani, R.J., 2021. Predictive inference with the jackknife+. Ann. Stat. 49, 486–507.

Breiman, L., 2001. Random forests. Mach. Learn. 45, 5–32.

Burnham, K.P., Anderson, D.R., 2002. Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer, New York.

Chen, T., Guestrin, C., 2016. XGBoost: A scalable tree boosting system. In: Proc. 22nd ACM SIGKDD Int. Conf. Knowledge Discovery and Data Mining, pp. 785–794.

Dong, X., Ma, L.Q., Li, Y., 2011. Characteristics and mechanisms of hexavalent chromium removal by biochar from sugar beet tailing. J. Hazard. Mater. 190, 909–915.

Duan, N., 1983. Smearing estimate: a nonparametric retransformation method. J. Am. Stat. Assoc. 78, 605–610.

Dunn, R., Wasserman, L., Ramdas, A., 2023. Distribution-free prediction sets for two-layer hierarchical models. J. Am. Stat. Assoc. 118(544), 2491–2502.

Gelman, A., Hill, J., 2007. Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge University Press, Cambridge.

Iacomi, P., Llewellyn, P.L., 2019. pyGAPS: a Python-based framework for adsorption isotherm processing and material characterisation. Adsorption 25, 1533–1542.

Kotaś, J., Stasicka, Z., 2000. Chromium occurrence in the environment and methods of its speciation. Environ. Pollut. 107, 263–283.

Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R.J., Wasserman, L., 2018. Distribution-free predictive inference for regression. J. Am. Stat. Assoc. 113, 1094–1111.

Mohan, D., Pittman, C.U., 2006. Activated carbons and low cost adsorbents for remediation of tri- and hexavalent chromium from water. J. Hazard. Mater. 137, 762–811.

Mohan, D., Rajput, S., Singh, V.K., Steele, P.H., Pittman, C.U., 2011. Modeling and evaluation of chromium remediation from water using low cost bio-char, a green adsorbent. J. Hazard. Mater. 188, 319–333.

Page, M.J., McKenzie, J.E., Bossuyt, P.M., et al., 2021. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ 372, n71.

Park, D., Yun, Y.-S., Park, J.M., 2004. Reduction of hexavalent chromium with the brown seaweed *Ecklonia* biomass. Environ. Sci. Technol. 38, 4860–4864.

Parkhurst, D.L., Appelo, C.A.J., 2013. Description of input and examples for PHREEQC version 3 — A computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations. U.S. Geological Survey Techniques and Methods, book 6, chap. A43.

Pedregosa, F., Varoquaux, G., Gramfort, A., et al., 2011. Scikit-learn: Machine learning in Python. J. Mach. Learn. Res. 12, 2825–2830.

Rai, D., Sass, B.M., Moore, D.A., 1987. Chromium(III) hydrolysis constants and solubility of chromium(III) hydroxide. Inorg. Chem. 26, 345–349.

Ramirez, A., Ocampo, R., Giraldo, S., Padilla, E., Flórez, E., Acelas, N., 2020. Experimental and computational data set on adsorption of Cr(VI) from water using an activated carbon. Data Brief 29, 105292. https://doi.org/10.1016/j.dib.2020.105292

Ran, Y.A., Sharma, S., Balestra, S.R.G., Li, Z., Calero, S., Vlugt, T.J.H., Snurr, R.Q., Dubbeldam, D., 2024. RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids. J. Chem. Phys. 161, 114106.

Rasmussen, C.E., Williams, C.K.I., 2006. Gaussian Processes for Machine Learning. MIT Press, Cambridge, MA.

Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmüller, U., Timmer, J., 2009. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. Bioinformatics 25, 1923–1929.

Roberts, D.R., Bahn, V., Ciuti, S., et al., 2017. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40, 913–929.

Schuirmann, D.J., 1987. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. J. Pharmacokinet. Biopharm. 15, 657–680.

Simonin, J.-P., 2016. On the comparison of pseudo-first order and pseudo-second order rate laws in the modeling of adsorption kinetics. Chem. Eng. J. 300, 254–263.

Tran, H.N., You, S.-J., Hosseini-Bandegharaei, A., Chao, H.-P., 2017. Mistakes and inconsistencies regarding adsorption of contaminants from aqueous solutions: A critical review. Water Res. 120, 88–116.

Wang, J., Guo, X., 2020. Adsorption kinetic models: Physical meanings, applications, and solving methods. J. Hazard. Mater. 390, 122156.

Vovk, V., Gammerman, A., Shafer, G., 2005. Algorithmic Learning in a Random World. Springer, New York.

Wellek, S., 2010. Testing Statistical Hypotheses of Equivalence and Noninferiority, 2nd ed. CRC Press, Boca Raton, FL.

Wilson, E.B., 1927. Probable inference, the law of succession, and statistical inference. J. Am. Stat. Assoc. 22, 209–212.
