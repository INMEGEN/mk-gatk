Pipeline for GATK short variant calling in the 100GMX project

	
Folder structure:

	001: Indel realigment step 
		results/	BAM realigned from this step and it's .bai index


	002: BQSR (Base quality score recalibration) step
		results/	BAM recalled from this step a recalibration tables
			*.recal-tab1 First recalibration table, result from BaseRecalibrator
			*.recal1.bam	First recalibrated BAM, result from PrintReads
			*.recal-tab2	Second recalibration table, result from BaseRecalibrator
			*.recal2.bam	Second recalibrated BAM, result from PrintReads

	003: Variant calling step
		results/	gVCFs files generated, those are not recalled
	

	.git
