# Equine80select Remapping Pipeline (EquCab2 to EquCab3)

## Overview
This repository contains a comprehensive computational pipeline designed to remap the **Equine80select genotyping array** from the legacy *EquCab2* reference genome to the modern *EquCab3* assembly.

Standard probe alignment often fails in repetitive regions or inverted loci. This pipeline utilizes a **"Context-Aware" Dual-Alignment strategy**, aligning both the short physical probe (50bp) and the longer genomic context (`TopGenomicSeq`) to ensure high-fidelity coordinate conversion, correct strand assignment, and precise reference/alternative allele determination consistent with VCF standards.

## Key Features
* **Automated Environment Setup:** A master Bash script handles dependency installation and environment creation via Conda/Mamba.
* **Hybrid Re-mapping:** Combines the physical precision of probe alignment with the biological context of `TopGenomicSeq` to resolve paralogs.
* **Strand Authority:** Uses context alignment to correct "complementary allele" errors caused by local inverted repeats.
* **Precise Arithmetic:** Correctly calculates coordinates for Infinium I/II chemistries, including complex edge cases like minus-strand deletions.
* **Quality Control (QC) & Filtration:** Includes post-processing steps to identify and filter markers with array design conflicts, polymorphic site interference, or low mapping quality.

## Prerequisites
* **Linux/macOS**
* **Conda** or **Mamba** package manager installed.
* **EquCab3 Reference Genome** (FASTA format).

## Quick Start

The entire workflow is orchestrated by a master Bash script (`master_pipeline.sh`).

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/equine80-remap.git](https://github.com/your-username/equine80-remap.git)
    cd equine80-remap
    ```

2.  **Run the Pipeline:**
    The script will automatically create the necessary Conda environment, run the remapping tool, and perform quality filtration.

    ```bash
    # Make the script executable
    chmod +x master_pipeline.sh

    # Run the tool
    ./master_pipeline.sh -i original_manifest.csv -r EquCab3.fa -o output_directory
    ```

    * `-i`: Path to the original Illumina manifest (CSV).
    * `-r`: Path to the EquCab3 Reference Genome (FASTA).
    * `-o`: Output directory for remapped files and QC logs.

## Pipeline Methodology

### 1. The Core Remapping Tool (Python)
The internal Python tool (`remap_manifest_final.py`) performs the heavy lifting:
* **Dual Alignment:** Aligns both `AlleleA_ProbeSeq` and `TopGenomicSeq` using `minimap2`.
* **Ref/Alt Determination:** Compares alignments of Allele A and Allele B contexts to computationally determine the Reference allele on EquCab3.
* **Strand Resolution:** Prioritizes the `TopGenomicSeq` strand to prevent allele flipping errors.
* **Deletion Correction:** Automatically adjusts coordinates for deletion events on the negative strand.

### 2. Additional Filtration
Following the remapping, the pipeline executes a QC module to define a high-quality marker set. Additional scripts are designed to identify markers with conflicts with the array design, polymorphic sites, and low-quality probes. Optional filtration of these markers allows selection of a high-quality set of markers for downstream analysis.

Filters include:
* **Mapping Quality (MAPQ):** Removes probes with low confidence alignments (e.g., MAPQ < 30).
* **Array Design Conflicts:** Flags markers where the remapped alleles do not match the expected probe chemistry.
* **Polymorphic Sites:** Identifies and filters probes that align to known polymorphic regions (secondary variants) which may affect binding efficiency.

## Outputs

The pipeline generates the following files in the output directory:

* **`Equine80_EquCab3_Full.csv`**: The raw remapped manifest containing all markers, including strand, coordinates, and MAPQ scores.
* **`Equine80_EquCab3_Filtered.csv`**: A high-quality subset of markers passing all QC filters, ready for downstream analysis (plink/imputation).
* **`QC_Report.txt`**: A summary log detailing how many markers were removed at each filtration step.

## Code Availability
Additional master script is developed in bash to orchestrate the full pipeline from installing dependencies, running the python script, and filtering the low quality marker.

## Citation
If you use this tool in your research, please cite:
> Tamer A. Mansour. "A Context-Aware Computational Pipeline for High-Precision Remapping of Genotyping Arrays: Updating the Equine80select Manifest to EquCab3." https://github.com/drtamermansour/Equine80select_remapper, 2025.

## License
MIT License
