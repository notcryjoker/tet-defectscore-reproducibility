# TET DefectScore — Analysis Code & Supplementary Data

> **Development source (GitHub):** https://github.com/notcryjoker/tet-defectscore-reproducibility (tag v1.0)
> **Archival snapshot (Zenodo):** https://doi.org/10.5281/zenodo.23643504

**Reproducibility package** for the study *"Thymopoiesis DefectScore is not associated
with myasthenia gravis status in thymic epithelial tumors"* (submitted to *Cancers*, MDPI).

This archive contains the complete analysis pipeline that produces every numerical result
reported in the manuscript, plus the machine-readable supplementary data tables.

> **What is deposited**
> 1. All analysis scripts (`00_config.R` → `10_supplementary.R`) in execution order.
> 2. The pre-locked signature definitions (`signature_definition.R`) and cached
>    `thymo_signatures.rds`.
> 3. Data-acquisition helpers (TCGA / GEO download scripts) so the pipeline is re-runnable
>    from public raw sources.
> 4. Supplementary data tables S1–S2 and the method-concordance values (CSV), as produced
>    by the pipeline.
>
> Raw TCGA-THYM and GEO expression matrices are **public and not included**; download
> instructions are below.

---

## 1. Directory layout

```
zenodo_deposit/
├── README.md
├── LICENSE                      # MIT
├── code/
│   ├── 00_config.R              # global config: ROOT, random seed, library loader
│   ├── 01_tcga_prepare.R       # TCGA-THYM main cohort build (STAR counts + clinical)
│   ├── 02_scores.R             # DefectScore construction (z-score primary, ssGSEA cross-check)
│   ├── 03_associations.R       # DefectScore vs WHO / MG / Masaoka associations
│   ├── 03b_scRNA.R             # single-cell reference validation (GSE147520 stromal scRNA)
│   ├── 04_survival.R           # relapse-free survival (RFS) analysis
│   ├── 05_external.R           # external validation on GSE29695 (Illumina DASL, GPL8432)
│   ├── 06_icb.R                # immune-checkpoint exploration layer (GSE181815)
│   ├── 07_figures.R            # main figures (Fig 4/5/6/7/9/10, 300 dpi)
│   ├── 08_results_fill.R       # compile Results text from .rds outputs
│   ├── 09_fig1_design.R        # Fig 1 study-design flowchart (CN)
│   ├── 09b_fig1_design_EN.R    # Fig 1 study-design flowchart (EN)
│   ├── 10_supplementary.R      # Supplementary Tables S1–S2 + Figure S1
│   ├── signature_definition.R  # pre-locked thymopoiesis cell-state gene sets (constants)
│   ├── thymo_signatures.rds    # cached signature object
│   ├── install_all.R           # install R packages (CRAN + Bioconductor)
│   ├── install_cran.R          # CRAN-only install
│   ├── check_pkgs.R            # dependency check used by the driver
│   ├── dl_tcga_api.sh          # fetch TCGA-THYM STAR counts via GDC API
│   ├── dl_tcga_parallel.sh     # parallel TCGA download
│   ├── dl_gse284479.sh         # design-stage reference (GSE284479, not analysed)
│   ├── dl_sc.py                # fetch GSE147520 single-cell matrix
│   ├── convert_citations.py    # pandoc → numbered-citation conversion (editorial)
│   ├── verify_refs_pubmed.py   # PubMed retraction/erratum sweep (editorial)
│   └── run_pipeline.sh         # driver: 01→08 + 03b + 09 (+ 10 for supplements)
└── data/
    └── supplementary/
        ├── Table_S1_prelocked_gene_sets.csv
        ├── Table_S2_TCGA-THYM_sample_scores.csv
        └── Table_S2_GSE29695_sample_scores.csv
```

## 2. Key method (must match the manuscript)

- **DefectScore = z(CORTICAL) − z(MEDULLARY)**, where each term is the gene-set mean
  *z*-score of the respective pre-locked thymopoiesis signature on the log₂-CPM matrix.
  The *z*-score is the **primary** method; **ssGSEA** (GSVA, `ssgseaParam`, α = 0.25) is a
  cross-check. All manuscript results were computed with the *z*-score method.
- **Random seed = 20261022** is fixed in `00_config.R` for exact reproducibility of scores.

## 3. Dependencies

R (≥ 4.3) with packages: `DESeq2`, `edgeR`, `GSVA`, `TCGAbiolinks`, `GEOquery`,
`org.Hs.eg.db`, `survival`, `survminer`, `clusterProfiler`, `ComplexHeatmap`, `limma`.
Install via `Rscript code/install_all.R`. Python 3 is needed only for the two editorial
helpers (`convert_citations.py`, `verify_refs_pubmed.py`).

## 4. How to run

The deposit root is configurable via the `TET_ROOT` environment variable (defaults to the
current working directory), so the code is portable across machines.

```bash
# from the unzipped deposit root:
export TET_ROOT="$(pwd)"          # optional; driver sets it automatically
Rscript code/install_all.R         # install R dependencies
bash code/run_pipeline.sh          # runs 01→08 + 03b(scRNA) + 09(Fig1)
Rscript code/10_supplementary.R    # regenerates Supplementary Tables S1–S2 + Figure S1
```

To run a single step manually, set `TET_ROOT` (or launch from the deposit root) and call
the script, e.g. `Rscript code/02_scores.R`.

## 5. Data sources (public)

| Cohort | Accession | Use |
|---|---|---|
| TCGA-THYM | dbGaP phs000178 / GDC project TCGA-THYM | main cohort (STAR counts + clinical biotab) |
| GSE29695 | GEO (Badve 2012) | external validation (Illumina HumanRef-8 WG-DASL v3.0, GPL8432) |
| GSE147520 | GEO | single-cell reference validation of pre-locked signatures |
| GSE284479 | GEO | design-stage reference only (not re-analysed) |
| GSE181815 | GEO (SRA raw only) | immune-checkpoint exploration layer (not re-analysed) |

Download helpers are in `code/` (`dl_tcga_api.sh`, `dl_tcga_parallel.sh`, `dl_sc.py`,
`dl_gse284479.sh`).

## 6. Reproduced outputs

Running the pipeline regenerates: the per-sample `DefectScore` / `ActIdx` (`02_scores.R`),
WHO/MG/Masaoka association statistics (`03_associations.R`), the RFS result (`04_survival.R`),
the GSE29695 external-validation gradient (`05_external.R`), all main figures
(`07_figures.R`), the study-design diagram (`09b_fig1_design_EN.R`), and the supplementary
tables/figure (`10_supplementary.R`). The `data/supplementary/*.csv` shipped here are the
exact outputs of `10_supplementary.R`.

## 7. License

Code is **MIT-licensed** (see `LICENSE`). Derived data (DefectScore tables) remain subject
to the original TCGA / GEO data-use terms.
