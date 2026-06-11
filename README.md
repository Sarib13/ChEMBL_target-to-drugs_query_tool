# RNA-seq and Drug Repurposing Toolkit

This repository contains two complementary command-line pipelines:

1. RNA-seq differential expression and pathway/drug discovery pipeline.
2. ChEMBL-based drug repurposing pipeline for one or more target genes.

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Pipeline 1: RNA-seq Analysis](#pipeline-1-rna-seq-analysis)
- [Pipeline 2: ChEMBL Drug Repurposing](#pipeline-2-chembl-drug-repurposing)
- [Typical End-to-End Workflow](#typical-end-to-end-workflow)
- [Limitations](#limitations)
- [License](#license)

## Project Overview

The toolkit helps you go from differential expression output to candidate drug leads.

- Pipeline 1 reads a GEO2R-style TSV, identifies significant genes, performs KEGG pathway enrichment via Enrichr, intersects disease-associated KEGG genes, and retrieves candidate drugs from KEGG.
- Pipeline 2 queries ChEMBL to return approved and investigational drugs for target genes, with mechanism and disease annotations.

## Repository Structure

- `Script/RNA_seq_Pipeline.py`: RNA-seq analysis pipeline (newly added).
- `Script/chembl_drugs_repurposing.py`: ChEMBL drug query and prioritization pipeline.
- `requirements.txt`: Python dependencies.

## Installation

### 1. Clone repository

```bash
git clone https://github.com/Sarib13/ChEMBL_target-to-drugs_query_tool.git
cd ChEMBL_target-to-drugs_query_tool
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
pip install requests
```

Python 3.8+ is recommended.

## Pipeline 1: RNA-seq Analysis

Script: `Script/RNA_seq_Pipeline.py`

### What this script does

Input: differential expression TSV and disease name.

1. Loads TSV and auto-detects gene, p-value, and logFC columns.
2. Filters significant genes using configurable thresholds.
3. Runs Enrichr KEGG 2021 Human pathway enrichment.
4. Searches KEGG DISEASE for disease-associated genes.
5. Finds overlap between significant genes and disease genes.
6. Queries KEGG DRUG for candidate drugs linked to overlapping genes.
7. Writes a multi-sheet Excel report.

### Default thresholds

- P-value threshold: `< 0.05`
- Absolute logFC threshold: `> 2.0`

### Usage

```bash
python Script/RNA_seq_Pipeline.py "Alzheimer disease" --input input.tsv --output rna_seq_analysis_results.xlsx
```

More examples:

```bash
python Script/RNA_seq_Pipeline.py "Parkinson disease" --input GSE203155_AS.top.table.tsv
python Script/RNA_seq_Pipeline.py "breast cancer" --input GSE46517.top.table.tsv --pval 0.01 --logfc 1.5
```

### Command line arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `disease` | Yes | - | Disease name used for KEGG DISEASE search |
| `--input`, `-i` | No | `input.tsv` | Path to differential expression TSV |
| `--output`, `-o` | No | `rna_seq_analysis_results2.xlsx` | Output Excel file |
| `--pval` | No | `0.05` | Significance cutoff for p-value |
| `--logfc` | No | `2.0` | Absolute logFC cutoff |

### Output sheets

The RNA-seq pipeline writes these sheets:

- `Original_Data`: input table after cleaning.
- `Upregulated_Genes`: significant upregulated genes.
- `Downregulated_Genes`: significant downregulated genes.
- `Pathways`: top enriched KEGG pathways from Enrichr.
- `Drugs`: KEGG-derived candidate drugs per disease-related gene.
- `Final_Results`: merged gene-level summary with regulation, pathways, and drugs.

## Pipeline 2: ChEMBL Drug Repurposing

Script: `Script/chembl_drugs_repurposing.py`

### What this script does

- Accepts genes from command line or an Excel sheet.
- Finds corresponding ChEMBL targets and mechanisms.
- Retrieves molecule name, approval year, max phase, mechanism/action type, and disease.
- Supports disease filtering and run logging.
- Exports results to Excel.

### Usage

```bash
python Script/chembl_drugs_repurposing.py --gene EGFR BRCA1 --output chembl_drugs.xlsx
python Script/chembl_drugs_repurposing.py --gene EGFR --disease "colorectal cancer" --logs
```

### Key arguments

| Argument | Description |
|---|---|
| `--gene`, `-g` | One or more target gene symbols |
| `--gene-excel` | RNA-seq output Excel containing genes |
| `--gene-column` | Optional gene column name in `--gene-excel` |
| `--disease`, `-d` | Optional disease filter |
| `--output`, `-o` | Output Excel file name |
| `--delay` | Delay between API calls |
| `--logs` | Enable verbose logs |
| `--log-file` | Custom path for run log file |

## Typical End-to-End Workflow

1. Run `Script/RNA_seq_Pipeline.py` on your DEGs TSV with a disease name.
2. Review `Final_Results` and pathway/drug sheets.
3. Feed selected genes into `Script/chembl_drugs_repurposing.py` for deeper ChEMBL-based drug annotation.

## Limitations

- KEGG and ChEMBL coverage may not include every disease-gene-drug relationship.
- API-driven workflows depend on internet connectivity and endpoint availability.
- Disease text matching can be broad; manually review filtered results.
- Runtime may increase for large gene lists because of rate-limited API queries.

## License

MIT License.
