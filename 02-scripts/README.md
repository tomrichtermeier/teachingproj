## 02-scripts
### This directory contains the scripts used in the teachingproj project

-   COMP_Kraken2_Bracken.Snakefile : The Snakemake script for running Kraken2 on all of the samples and evaluating relative abundance using Bracken
-   COMP_Kraken2_Zape.Snakefile : The Snakemake script for running Kraken2 on just the Zape samples
-   COMP_Kraken2_Bracken_SILVA.Snakefile : The Snakemake script for running Kraken2 on all of the samples and evaluating relative abundance using Bracken (with SILVA 16S rRNA databse)
-   COMP_MetaPhlAn4.Snakefile : The Snakemake script for running MetaPhlAn4 on all of the samples
-   COMP_MetaPhlAn_Zape2.Snakefile : The Snakemake script for running MetaPhlAn4 on just the Zape samples
-   PREP_remove_hostDNA.Snakefile : A script for extracting the non-host DNA
-   QUAL_damageprofiles.Snakefile : A script for creating damage profiles for all samples
-   REFG_genome_reconstruction_Pcopri_Smaltophilia.Snakefile : A script for the reconstruction of the genomes of Pcopri and Smaltophilia based (with reference genome)
-   REFG_genome_reconstruction_Tsuccinifaciens.Snakefile : A script for the reconstruction of the genomes of Tsuccinifaciens based (with reference genome)

### types of analysis
1. **PREP** = data preparation
2. **COMP**: = compositional profiling
3. **QUAL**: = quality evaluation
4. **REFG**: = reference genome alignment