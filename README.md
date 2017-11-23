# mk-gatk - pipeline for short Single Nucleotide and small Indel variant calling.

### Abreviations:

**SNVs**: Single Nucleotide Variants
**VQSR**: Variant Quality Score Recalibration
**GATK**: Genome Analysis Toolkit

## About mk-gatk

The mk-gatk pipeline uses the GATK to detect single nucleotide variants and small indels in population data.

GATK offers a wide variety of tools with a primary focus on variant discovery and genotyping.

This pipeline uses the "best practices" from the Broad Institute as recommendations for performing variant discovery analysis in high-throughput sequencing (HTS) data, to go from mapped sequence data all the way to an appropriately filtered variant callset that can be used in downstream analyses.

Using Indel Ralignment, Base Recalibration (of base call quality), Haplotype Caller and Variant Quality Score Recalibration, this pipeline delivers a multisample .vcf file with every genomic variant detected.

## Pipeline configuration.

### Dependencies:

**[GATK](https://software.broadinstitute.org/gatk/)**

The Genome Analysis Toolkit is a structured programming framework designed to ease the development of efficient and robust analysis tools for next-generation DNA sequencers using the functional programming philosophy of MapReduce. The GATK provides a set of data access patterns that encompass the majority of analysis tool needs. Separating specific analysis calculations from common data management infrastructure enables the optimization of the GATK framework for correctness, stability, and CPU and memory efficiency and to enable distributed and shared memory parallelization. \[[1](http://genome.cshlp.org/content/20/9/1297.short)\]

IMPORTANT NOTE: path to GATK must be declared in mk-gatk/analysis/config.mk.

GATK development and installation instructions can be found at: [https://software.broadinstitute.org/gatk/download/](https://software.broadinstitute.org/gatk/download/)

GATK publication can be found at: [McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research, 20(9), 1297-1303.](http://genome.cshlp.org/content/20/9/1297.short)

**[samtools](https://github.com/samtools/bcftools)**

IMPORTANT NOTE: samtools must be called via `samtools` command.

Development and installation instructions can be found at: [https://github.com/samtools](https://github.com/samtools)

**[picard-tools](http://broadinstitute.github.io/picard/)**

IMPORTANT NOTE: Picard tools must be called via `picard-tools` command.

Development and installation instructions can be found at: [https://github.com/broadinstitute/picard](https://github.com/broadinstitute/picard)

**[vcf-validator](https://github.com/EBIvariation/vcf-validator)**

IMPORTANT NOTE: vcf-validator must be called via `vcf-validator` command.

Development and installation instructions can be found at: [https://github.com/EBIvariation/vcf-validator](https://github.com/EBIvariation/vcf-validator)

**[BCFtools](https://samtools.github.io/bcftools/)**

IMPORTANT NOTE: bcftools must be called via `bcftools` command.

development and installation instructions can be found at: [https://github.com/samtools](https://github.com/samtools)

### Input files

mk-gatk requires:
1) Multiple BAM files with marked duplicates, and .bai index. Sample files must located at mk-gatk/data/
1) A common .fasta file for the reference genome, the same version used to generate the sample .bam files. Reference genome files must be located at mk-data/reference/
1) A gold standart INDEL vcf reference file, that must be declared in mk-gatk/analysis/config.mk
1) A dbSNP vcf reference file, that must be declared in mk-gatk/analysis/config.mk
1) A HAPMAP vcf reference file, that must be declared in mk-gatk/analysis/config.mk
1) An OMNI vcf reference file, that must be declared in mk-gatk/analysis/config.mk
1) A one thousand genomes vcf reference file, that must be declared in mk-gatk/analysis/config.mk

### Configuration file

This pipeline includes a config.mk file (located at mk-gatk/analysis/001/config.mk, and propagated to every other module), where you can adjust the following paramters:

PATH=

GATK=

REF=

INDELs=
dbSNP=
HAPMAP=
OMNI=
OTG=

NT=

## Module description.

### 001 -> Indel realigment step

		This step takes an input BAM, locally realign reads such that 
		the number of mismatching bases is minimized across all the reads. 
		This step is important because BWA-like algorithms penalice more 
		indels than mismatches producing bias generating a lot of mismatches where INDELs are supossed to be.

		There are 2 steps to the realignment process:

		1) Determining (small) suspicious intervals which are likely 
		    in need of realignment (RealignerTargetCreator tool).
		2) Running the realigner over those intervals (IndelRealigner)
		     using smith waterman algorithm.

		The indel realigment step, can take a reference of INDEL positions in the 
		1000G-Project in order to detect sites that are more propense to be realigned.
		
	        For more details see: [GATK-forum](https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)

### 002 -> BQSR (Base quality score recalibration) step
		The BQSR apply machine learning to model errors produced by the machine-sequencers.
		Its important because  technical error, leading to over- or under-estimated base quality scores in the data. 
		Some of these errors are due to the physics or the chemistry of how the sequencing reaction 
		works and some are probably due to manufacturing flaws in the equipment.
		
		The base recalibration process involves two key steps: 
		1) The program builds a model of covariation based on the data and a set 
		   of known variants (BaseRecalibrator).
		2) Adjusts the base quality scores in the data based on the model (PrintReads).

		The indel realigment step need a reference for variants:
		The known variants are used to mask out bases at sites of real (expected) variation, to avoid counting real variants as errors. 
		Outside of the masked sites, every mismatch is counted as an error. 

		For more details see: http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr


### 003 -> Individual variant calling step
		This step applies the HaplotypeCaller algorithm of GATK, for stimation of SNPs and short INDELs
		likelihoods for every nucleotide in the genome (using the option: --emitRefConfidence GVCF)

		Note1: All the infered likelihoods still need to be evaluated in step 5, in order to select confident variants. 

                Note2: however that the algorithms used to calculate variant likelihoods is not well 
		suited to extreme allele frequencies (relative to ploidy) so its use is not recommended 
		for somatic (cancer) variant discovery.

		For more details see: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

### 004 -> Joint genotyping step
		In this steps many gVCFs producen in step 3 are used to make a group genotyping,
		this step important to distingh variants in samples that could have low sequencing-depth, 
		taking into account the information in other input samples.

		This tool will combine all spanning records (samples), produce correct genotype likelihoods, re-genotype the newly merged record.

	For more details see: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

	WARNING: In this step many temporary files are requiered to consider, in order to obtain a successful run.:
		1. Change the JAVA tmp directory: export _JAVA_OPTIONS=-Djava.io.tmpdir=/100g/analysis/007a/004/tmp
		2. Extend the limit of open files: ulimit -n 32768

### 005a ->  Variant recalibration (VQSR) steps
		In this step, for each possible variant the quality (likelihood) for each variant is readjusted.
		This is important step to correct multiple testing bias. 

		The approach taken by variant quality score recalibration (VQRS) is to develop a continuous, 
		covarying estimate of the relationship between SNP call annotations (QD, SB, HaplotypeScore, HRun, for example) 
		and the the probability that a SNP is a true genetic variant versus a sequencing 
		or data processing artifact. This model is determined adaptively based on "true sites" provided as input (variant references).

		The VQRS evaluates variants in a two step process:

		1) VariantRecalibrator:
		Create a Gaussian mixture model by looking at the annotations values over a high 
		quality subset of the input call set and then evaluate all input variants. 
		This step produces a recalibration file.

		1) ApplyRecalibration:
		Apply the model parameters to each variant in input VCF files producing 
		a recalibrated VCF file in which each variant is annotated with its VQSLOD value. 
	        In addition, this step will filter the calls based on this new lod score by adding lines to the 
		FILTER column for variants that don't meet the specified lod threshold. And marking with PASS those that pass threshold. 

		The model is run in a first step to recall SNPs and in a second step to recall INDELs.
		
		For more details see: http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr

### 005b -> Hard Filtering
	This module applies hard filters to a variant callset that is too small for VQSR or for which truth/training sets are not available. variants are filtered using parameters recommended by the Broad Institute best practices.

	For more information  go to the following site: [https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)

## mk-gatk directory structure

```
mk-gatk/			##Pipeline main directory.
├── analysis			##Directory for mk modules.
│   ├── 001			##Module for Indel realigment
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data -> ../../data/ ##Symbolic link to data directory at main directory.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   ├── QC			##Directory for evaluation of Quality Control of results.
│   │   │   ├── bin		##Executables directory.
│   │   │   │   └── targets	##Bash script to print required targets to STDOUT.
│   │   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   │   ├── QC-results	##Storage directory for QC files built by mkfile. 
│   │   │   └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   │   ├── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── UT			##Directory for performing Util Test of this module.
│   │       ├── bin		##Executables directory.
│   │       │   └── targets	##Bash script to print required targets to STDOUT.
│   │       ├── Log		##Directory for storage of the logs from running util tests.
│   │       ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │       └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   ├── 002			###Module for BQSR (Base quality score recalibration).
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data -> ../001/results	##Directory for storage of input data for the module.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   ├── QC			##Directory for evaluation of Quality Control of results.
│   │   │   ├── bin		##Executables directory.
│   │   │   │   └── targets	##Bash script to print required targets to STDOUT.
│   │   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   │   ├── QC-results	##Storage directory for QC files built by mkfile. 
│   │   │   └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   │   ├── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── UT			##Directory for performing Util Test of this module.
│   │       ├── bin		##Executables directory.
│   │       │   └── targets	##Bash script to print required targets to STDOUT.
│   │       ├── Log		##Directory for storage of the logs from running util tests.
│   │       ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │       └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   ├── 003			###Module for individual variant calling (g.vcf generation)
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data -> ../002/results	##Directory for storage of input data for the module.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   ├── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── UT			##Directory for performing Util Test of this module.
│   │       ├── bin		##Executables directory.
│   │       │   ├── targets	##Bash script to print required targets to STDOUT.
│   │       │   └── vcf_validator	##tool for validating vcf's in this UT
│   │       ├── Log		##Directory for storage of the logs from running util tests.
│   │       ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │       └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   ├── 004			###Module for joint Genotyping
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data		##Directory for storage of input data for the module.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   ├── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── UT			##Directory for performing Util Test of this module.
│   │       ├── bin		##Executables directory.
│   │       │   ├── targets	##Bash script to print required targets to STDOUT.
│   │       │   └── vcf_validator	##tool for validating vcf's in this UT
│   │       ├── config.mk	##configuration file for UT parameters
│   │       ├── Log		##Directory for storage of the logs from running util tests.
│   │       ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │       └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   ├── 005a			##Module for Variant Quality Score Recalibration
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data -> ../004/results/	##Directory for storage of input data for the module.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   ├── QC			##Directory for evaluation of Quality Control of results.
│   │   │   ├── bin		##Executables directory.
│   │   │   │   ├── Graphicator.R	##Script for QC plots
│   │   │   │   └── targets	##Bash script to print required targets to STDOUT.
│   │   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   │   ├── QC-results	##Storage directory for QC files built by mkfile. 
│   │   │   └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for
│   │   ├── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── UT			##Directory for performing Util Test of this module.
│   │       ├── bin		##Executables directory.
│   │       │   ├── targets	##Bash script to print required targets to STDOUT.
│   │       │   └── vcf_validator	##tool for validating vcf's in this UT
│   │       ├── config.mk	##configuration file for UT parameters
│   │       ├── Log		##Directory for storage of the logs from running util tests.
│   │       ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │       └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
│   ├── 005b			### Module for hard filtering vcfs with low number of variants, that can not be proccessed by VQSR.
│   │   ├── bin			##Executables directory.
│   │   │   └── targets		##Bash script to print required targets to STDOUT.
│   │   ├── data		##Directory for storage of input data for the module.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   │   └── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated by mkfile.
│   │   └── README.txt	##Describes module's objective, software dependencies, input and output data format, and references for module design.
├── data	##Main data directory. Stores BAM files for pipeline processing.
│   ├── sample1.bam	##A BAM file.
│   ├── sample1.bai	##Index for BAM file.
│   ├── sample2.bam	##Another BAM file.
│   ├── sample2.bai	##Index for "another BAM file".
│   ├── sampleN.bam	##More BAM files.
│   └── sampleN.bai	##Indexes for more BAM files.
├── notes/		##Notes about proper execution of modules.
├── README.md		##This document. General workflow description.
├── reference		##Directory for storgae of references used by this pipeline.
└── test-materials	##Directory for data required for running tests in the module.
    ├── bams		##5 different bams from chr22.
    └── reference_vcf	##A vcf file against which to test Specificity and Sensibility of Variant Calling.
```

### References
[1] [McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research, 20(9), 1297-1303.](http://genome.cshlp.org/content/20/9/1297.short)

### Author Info
Developed by Fernando Perez (frpvillatoro@gmail.com) and [Israel Aguilar](https://www.linkedin.com/in/israel-aguilar-ba625949/) (iaguilaror@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/), by request of [INMEGEN](http://www.inmegen.gob.mx/). 2017.

