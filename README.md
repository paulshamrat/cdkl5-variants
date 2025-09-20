
# CDKL5 Variants Analysis Workflow

![CDKL5 graphical abstract](cdkl5-graphical-abstract-250827.svg)

This repository provides the workflow for the analysis of CDKL5 gene variants, including variant annotation, stability prediction, CDKL5 and it's substrate structural modeling, and reclassification of variants based on $\Delta\Delta G_{\text{Folding}}$ and $\Delta\Delta G_{\text{Binding}}$.

## Methods:

1. **Variant Curation:** CDKL5 variant curation from ClinVar, gnomAD and Literatures.
2. **$\Delta\Delta G_{\text{Folding}}$:** Sequence-based and structure-based $\Delta\Delta G_{\text{Folding}}$ predictions using multiple methods (SAAFEC-SEQ, I-Mutant2.0, INPS, DDGun, mCSM, DDMut, DDGEmb).
3. **Protein-Protein Complex:** Colabfold: CDKL5-Substrate complex modeling using ColabFold.
4. **Protein-Protein Docking:** HADDOCK3: CDKL5-Substrate docking using HADDOCK3.
5. **$\Delta\Delta G_{\text{Binding}}$:** Structure-based $\Delta\Delta G_{\text{Binding}}$ analysis using follwoing tools: DDMutPPI, iSee, mCSM-PPI and SAAMBE-3D.
6. **Pathogenicity prediction:** CDKL5 variants pathogenicity prediction using pathogenicity predctors (PolyPhen-2, MutPred2, ESM-1v, and AlphaMissense).
7. **Variant Reclassification:** Reclassificaiton of CDKL5 variants based on DDG_folding, Binding and Pathogenicity.


## Citation

If this repository is helpful, please cite:  

> Paul, S. K.; Panday, S. K.; Boccuto, L.; Alexov, E. *CDKL5 Deficiency Disorder: Revealing the Molecular Mechanism of Pathogenic Variants*. Preprints 2025, 2025081241.  
> [https://doi.org/10.20944/preprints202508.1241.v1](https://www.preprints.org/manuscript/202508.1241/v1)
---

**Contact:**  
For questions or contributions, please open an issue or contact the repository maintainer.
