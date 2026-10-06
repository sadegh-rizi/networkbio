# References and resources — multi-omics network inference & graph-theoretic controllability

Compiled from the literature searches conducted across this conversation. Organized by the topic sequence discussed.

## 1. Factor-model multi-omics integration (MOFA+)

- Argelaguet et al. **"MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data."** *Genome Biology* (2020). https://link.springer.com/article/10.1186/s13059-020-02015-1
- **MOFA2 FAQ** (developer notes on covariates degrading model fit, factor-count selection). https://biofam.github.io/MOFA2/faq.html
- **MOFA2 documentation/site.** https://biofam.github.io/MOFA2/
- **MOFA GitHub.** https://github.com/bioFAM/MOFA

## 2. Causal network propagation (COSMOS / COSMOS+)

- Dugourd et al. **"Causal integration of multi-omics data with prior knowledge to generate mechanistic hypotheses."** *Molecular Systems Biology* (2021). https://www.embopress.org/doi/full/10.15252/msb.20209730 (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7838823/)
- Rodriguez-Mier et al. **"Modeling causal signal propagation in multi-omic factor space with COSMOS."** *bioRxiv* (2024). https://doi.org/10.1101/2024.07.15.603538
- **cosmosR documentation.** https://saezlab.github.io/cosmosR/

## 3. Optimization-based network inference (CORNETO)

- **"Unifying multi-sample network inference from prior knowledge and omics data with CORNETO."** *Nature Machine Intelligence* (2025). https://www.nature.com/articles/s42256-025-01069-9
- **"Unified knowledge-driven network inference from omics data."** *bioRxiv* (2024). https://www.biorxiv.org/content/10.1101/2024.10.26.620390v1
- **CORNETO GitHub.** https://github.com/saezlab/corneto
- **CORNETO website.** https://corneto.org/

## 4. Multilayer / supra-adjacency spectral methods

- **Frontiers in Molecular Biosciences (2022).** Multilayer omics network review/method. https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.967205/full
- **Frontiers in Genetics (2019).** Multilayer network integration. https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.01381/full
- Zhang et al. **"Spectral clustering of single-cell multi-omics data on multilayer graphs (SCML)."** *Bioinformatics* 38(14):3600 (2022). https://academic.oup.com/bioinformatics/article/38/14/3600/6598796 (preprint: https://www.biorxiv.org/content/10.1101/2022.01.24.477443v1.full)
- **"Characterization of multiple topological scales in multiplex networks through supra-Laplacian eigengaps."** *arXiv:1603.08464.* https://arxiv.org/pdf/1603.08464
- **"Multi-set spectral clustering of time-evolving networks using the supra-Laplacian."** *arXiv:2409.11984.* https://arxiv.org/abs/2409.11984

## 5. Hypergraphs for higher-order regulation

- **"Inferring gene regulatory networks by hypergraph generative model" (HyperG-VAE).** *Cell Reports Methods / ScienceDirect* (2025); preprint *bioRxiv* (2024). https://www.sciencedirect.com/science/article/pii/S2667237525000621 · https://www.biorxiv.org/content/10.1101/2024.04.01.586509v1 · PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12256954/
- **MORE-hypergraph framework** (mRNA + methylation + miRNA integration via self-attention hyperedges). *Briefings in Bioinformatics* 26(1):bbae658. https://academic.oup.com/bib/article/26/1/bbae658/7927591

## 6. Network controllability — foundations

- Liu, Slotine, Barabási. **"Controllability of complex networks."** *Nature* (2011). https://barabasi.com/media/pub_imports/files/329.pdf (doi:10.1038/nature10011)
- **"Network controllability."** *Wikipedia* (background/summary). https://en.wikipedia.org/wiki/Network_controllability
- **"Structural Controllability of Multiplex Networks With the Minimum Number of Driver Nodes."** https://www.researchgate.net/publication/378722151
- **"Structural characteristics in network control of molecular multiplex networks"** (TRN–PPI two-layer multiplex, interlayer coupling changes minimal driver set). PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10062666/
- **"Detection of gene communities in multi-networks reveals cancer drivers."** *arXiv:1507.08415.* https://arxiv.org/pdf/1507.08415
- **"Network controllability: viruses are driver agents in dynamic molecular systems."** *bioRxiv.* https://www.biorxiv.org/content/10.1101/311746.full.pdf

## 7. Single-cell / GRN construction + controllability tools (course-project candidates)

- **CEFCON** — "Deciphering driver regulators of cell fate decisions from single-cell transcriptomics data." PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10733330/
- **Drivergene.net** — "A Cytoscape app for the identification of driver nodes of large-scale complex networks and case studies." PubMed: https://pubmed.ncbi.nlm.nih.gov/39047507
- **scNetViz** — "from single cells to networks using Cytoscape." PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8593621/
- **"Performance assessment of sample-specific network control methods for bulk and single-cell biological data analysis."** *PLOS Computational Biology.* https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1008962
- **Sc-compReg** — "enables the comparison of gene regulatory networks between conditions using single-cell data." *Nature Communications* (2021). https://www.nature.com/articles/s41467-021-25089-2
- **pySCENIC** (GRN inference, GRNBoost2 step). GitHub: https://github.com/aertslab/pySCENIC · Docs: https://pyscenic.readthedocs.io/

## 8. Target control theory / personalized & multi-omics controllability

- **"Target control of linear directed networks based on the path cover problem."** *Scientific Reports* (2024). https://www.nature.com/articles/s41598-024-67442-7 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11266607/)
- **"A novel algorithm for finding optimal driver nodes to target control complex networks and its applications for drug targets."** PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC5780855/
- **"Controllability analysis of the directed human protein interaction network identifies disease genes and drug targets."** *arXiv:1511.07768.* https://arxiv.org/pdf/1511.07768
- **OptiCon** — "Optimal control nodes in disease-perturbed networks as targets for combination therapy." *Nature Communications* (2019). https://www.nature.com/articles/s41467-019-10215-y (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6522545/)
- **CPGD** — "Network controllability-based algorithm to target personalized driver genes for discovering combinatorial drugs of individual patients." *Nucleic Acids Research* (2021). https://academic.oup.com/nar/article/49/7/e37/6090299 (preprint: https://www.biorxiv.org/content/10.1101/571620v2.full)
- **SCS** — "Discovering personalized driver mutation profiles of single samples in cancer by network control strategy." *Bioinformatics* (2018). https://dx.doi.org/10.1093/bioinformatics/bty006
- **PNC** — "A novel network control model for identifying personalized driver genes in cancer." *PLOS Computational Biology.* https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007520
- **"Validation of a multi-omics strategy for prioritizing personalized candidate driver genes."** PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5122402/
- **"Patient-specific driver gene prediction and risk assessment through integrated network analysis of cancer omics profiles."** PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402507/

## 9. Hypergraph controllability (emerging, 2024–2026)

- **"Data-driven Control of Hypergraphs: Leveraging THIS to Damp Noise in Diffusive Hypergraphs."** *arXiv:2511.08647* (2025). https://arxiv.org/html/2511.08647
- Bao, Ma, Ma. **"Controllability of higher-order networks."** *ScienceDirect* (2024). https://www.sciencedirect.com/science/article/pii/S0378437124006174
- Pickard et al. **"Structural Controllability of Large-Scale Hypergraphs."** *arXiv:2603.19955* (2026). https://arxiv.org/abs/2603.19955
- **"On the minimum driver node set of k-uniform linear hypertree networks."** *ScienceDirect* (2024). https://www.sciencedirect.com/science/article/abs/pii/S0096300324001747
- **"Controllability and Observability of Temporal Hypergraphs."** *arXiv:2408.12085.* https://arxiv.org/pdf/2408.12085
- **"High-order Knowledge Based Network Controllability Robustness Prediction: A Hypergraph Neural Network Approach."** *arXiv:2603.02265.* https://arxiv.org/pdf/2603.02265
