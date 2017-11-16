Pipeline for GATK short variant calling in the 100GMX project.

	
Folder structure:

##	001: Indel realigment step 
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


		For this MK:

		targets/	BAM files
		results/	BAM realigned from this step and it's .bai index


##	002: BQSR (Base quality score recalibration) step
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
		
		

		For this MK:

		There are perfomed to completed BQRS steps. This is recommended in order to 
		minimize biases in quality scores.

		targets/ BAM files
		
		results/	BAM recalled from this step a recalibration tables
			*.recal-tab1 First recalibration table, result from BaseRecalibrator
			*.recal1.bam	First recalibrated BAM, result from PrintReads
			*.recal-tab2	Second recalibration table, result from BaseRecalibrator
			*.recal2.bam	Second recalibrated BAM, result from PrintReads

		For more details see: http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr


##	003: Individual variant calling step
		This step applies the HaplotypeCaller algorithm of GATK, for stimation of SNPs and short INDELs
		likelihoods for every nucleotide in the genome (using the option: --emitRefConfidence GVCF)

		Note1: All the infered likelihoods still need to be evaluated in step 5, in order to select confident variants. 

                Note2: however that the algorithms used to calculate variant likelihoods is not well 
		suited to extreme allele frequencies (relative to ploidy) so its use is not recommended 
		for somatic (cancer) variant discovery.

		For this MK:

		targets/	BAM files

		results/	gVCFs files generated, those are not recalled


	      For more details see: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

##	004: Joint genotyping step
	   In this steps many gVCFs producen in step 3 are used to make a group genotyping,
	   this step important to distingh variants in samples that could have low sequencing-depth, 
	   taking into account the information in other input samples.

	  This tool will combine all spanning records (samples), produce correct genotype likelihoods
	  re-genotype the newly merged record.
 
		
	 For this MK:

	 Its not necesary to define targets, the script output will be only one VCF for all the input gVCFs.

         For more details see: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php

	WARNING: In this step many temporary files are requiered so consider:
		1. Change the JAVA tmp directory: export _JAVA_OPTIONS=-Djava.io.tmpdir=/100g/analysis/007a/004/tmp
		2. Extend the limit of open files: ulimit -n 32768
		 In order to obtain a successful run.

	Note: 
		- data/95 all amerindians
		- data/94 excluded subject SM-3MGPV due to incorrect FST calculation
		- data/92 equal to data/94 plus 2 aditional subject removed (SM-3MG54 and SM-3MG59) due to missing label in concordance.
##	005:  Variant recalibration steps
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

		2) ApplyRecalibration:
		    Apply the model parameters to each variant in input VCF files producing 
		    a recalibrated VCF file in which each variant is annotated with its VQSLOD value. 
	            In addition, this step will filter the calls based on this new lod score by adding lines to the 
		    FILTER column for variants that don't meet the specified lod threshold. And marking with PASS those that pass threshold. 

	    The model is run in a first step to recall SNPs and in a second step to recall INDELs.
		
	    For this MK:
		targets/   VCF file to apply VQRS
		results/   Recaled VCF

		For more details see: http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr

##	005b:	Hard Filter
	
##	006:	Split VCF between PASS and FAIL VQSR filtered variants

##	007a:	Compressing VCF files
