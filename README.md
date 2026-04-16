# Cas13IP RNA-Seq Profiles

**Interactive data explorer:** https://farehlab.github.io/Cas13IP-RNA-Seq-Profiles/

---

## Preprint

> **CRISPR-Cas13 Localises Within RNA-Rich Condensates To Drive RNA Recognition and Cleavage** *pending TBC*
> 
> Gurjeet Kaur Gill Jagjeet Singh <sup>1,2</sup>, Wenxin Hu, Carolyn Shembrey, Adrian Hodel, Ilia Viskoboinik, Paul McMillan, Joseph Trapani, Mohamed Fareh
> 
> *bioRxiv* (in preparation)
> 
> DOI: *pending*
> 
> Correspondence: Mohamed.fareh@petermac.org
>
> ---


Syntax: Your text<sup>superscript</sup>
Example: E = mc<sup>2</sup> renders as E = mc²


## Affiliations

¹ Peter MacCallum Cancer Centre, Melbourne 3000, Australia

² Sir Peter MacCallum Department of Oncology, The University of Melbourne, Parkville 3052, Australia

³ Department of Biology, Institute of Molecular Health Sciences, ETH Zürich, Zürich, Switzerland

⁴ Biological Optical Microscopy Platform (BOMP), University of Melbourne, Parkville, Victoria, Australia

---


---

## About

This repository contains differential expression data from RNA immunoprecipitation sequencing (RIP-Seq) experiments characterising the RNA association profile of PspCas13b across subcellular fractions in HEK293T cells. Data are provided as processed differential expression tables and can be explored interactively via the link above. Additionally, the raw count matrix has been provided along with the factor data in COUNTDATA.xlsx

All experiments were performed at the **Peter MacCallum Cancer Centre**, Melbourne, Australia.

---

## Repository Structure

```
Cas13IP-RNA-Seq-Profiles/
├── index.html              # Interactive data explorer
├── README.md               # This file
├── data/
│   ├── COUNTDATA.xlsx      # Raw featureCounts matrix + sample metadata
│   ├── Total-IP_Total.csv
│   ├── CCEF-IP_Total.csv
│   ├── CCEF-IP_CCEF.csv
│   ├── CCEF-IP_Total-IP.csv
│   └── CCEF_Total.csv
└── analysis/
    └── DE_analysis.R       # Full R analysis script (limma-voom pipeline)
```

---

## Experimental Design

12 samples across 4 conditions, 3 biological replicates each:

| Sample Type | Label | Description |
|---|---|---|
| Samples 1–3 | **Total** | Total cytoplasmic RNA |
| Samples 4–6 | **CCEF** | Condensate-enriched cytoplasmic fraction |
| Samples 7–9 | **Total-IP** | Cytoplasmic Cas13 immunoprecipitation |
| Samples 10–12 | **CCEF-IP** | Condensate-enriched Cas13 immunoprecipitation |

*Sample Type numbering reflect data in countmatrix. 

---

## Datasets

Five pairwise comparisons are included:

| File | Comparison | Biological question |
|---|---|---|
| `Total-IP_Total.csv` | Total-IP vs Total | Which RNAs associate with cytoplasmic Cas13? |
| `CCEF-IP_Total.csv` | CCEF-IP vs Total | Which RNAs are enriched in condensate-associated Cas13 vs total? |
| `CCEF-IP_CCEF.csv` | CCEF-IP vs CCEF | Which RNAs are specifically pulled down by Cas13 within condensates? |
| `CCEF-IP_Total-IP.csv` | CCEF-IP vs Total-IP | How does condensate-associated Cas13 differ from cytoplasmic Cas13? |
| `CCEF_Total.csv` | CCEF vs Total | Which RNAs are enriched in the condensate fraction itself? |

---

## Data

Each differential expression file contains:

| Column | Description |
|---|---|
| `labels` | Gene symbol (HGNC) |
| `log2fc` | log₂ fold change |
| `pvalue` | Nominal p-value |
| `fdr` | Adjusted p-value (Benjamini-Hochberg FDR) |

Significance thresholds: **\|log₂FC\| > 1.5**, **FDR < 0.05**

> **Raw sequencing data:** NCBI GEO accession pending. Will be deposited upon preprint submission.

---

## Methods Summary

### Cell culture and fractionation
HEK293T cells transiently expressing dPspCas13b-2xmNeonGreen and cognate crRNA were grown to ~80% confluency. Upon granule formation, cells were lysed (50 mM Tris pH 7.4, 1 mM EDTA, 150 mM NaCl, 0.2% Triton X-100) with RNase inhibitor and protease inhibitor. Fractionation was performed by differential centrifugation: total cytoplasmic lysate (1,000g), and condensate-enriched fraction (CCEF) pellet (18,000g × 2). Immunoprecipitation of FLAG-tagged Cas13 was performed using Pierce Anti-DYKDDDK Magnetic Agarose beads (Invitrogen) with acid elution (0.1M glycine pH 2.8).

### RNA-Seq library preparation and sequencing
RNA was isolated using TRIzol. Library preparation was performed from 100 ng total RNA using the Illumina Ribo-Zero Plus kit. Sequencing was carried out at the Australian Genome Research Facility (AGRF, Melbourne) on the Illumina NextSeq platform, generating paired-end reads.

### Alignment and differential expression
Reads were aligned to the **GRCh38** reference genome using **STAR** and gene-level counts generated using **featureCounts**, both run via Galaxy (June 2024). Count matrices were imported into R for normalisation and differential expression analysis using the **limma-voom pipeline**. Genes with CPM ≥ 0.5 in at least 2 samples were retained. Genes with \|log₂FC\| > 1.5 and FDR < 0.05 (Benjamini-Hochberg) were considered significantly enriched or depleted.

### Software versions

| Tool | Version |
|---|---|
| R | 4.3.3 |
| limma | 3.58.1 |
| edgeR (dependency) | 4.0.16 |
| STAR | run via Galaxy (June 2024) |
| featureCounts | run via Galaxy (June 2024) |
| biomaRt | 2.58.2 |
| ggplot2 | 3.4.4 |
| pheatmap | 1.0.12 |
| VennDiagram | 1.7.3 |

---

## Interactive Explorer

The interactive data explorer at https://farehlab.github.io/Cas13IP-RNA-Seq-Profiles/ includes:

- **Volcano Plot** — visualise differential expression for each comparison with adjustable FDR and log₂FC cutoffs, gene search, and filtered data download
- **Gene Search** — look up any gene to see its profile across all five comparisons, with links to NCBI, UniProt, GeneCards, and Human Protein Atlas
- **Heatmap** — explore top differentially expressed genes across selected datasets, with enriched/depleted split view

All data are embedded client-side — no server required.

---

## Citation

If you use this data, please cite:

> Jagjeet Singh GKG, Hu W, Shembrey C, Hodel A, Viskoboinik I, McMillan P, Trapani J, Fareh M. CRISPR-Cas13 Localises Within RNA-Rich Condensates To Drive RNA Recognition and Cleavage. *bioRxiv* (in preparation).

---

## Issues & Questions

Found a bug or have a question about the data? Please [open an issue](../../issues) on this repository — this way answers are visible to everyone.

---

## License

Code: [MIT License](LICENSE)
Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to reuse with attribution

> ⚠️ This repository is currently private. It will be made public upon preprint deposition.

---

*Repository maintained by the Fareh Lab, Peter MacCallum Cancer Centre. Contact: Mohamed.fareh@petermac.org*
