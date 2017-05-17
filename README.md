Pipeline for GATK short variant calling in the 100GMX project

	
Folder structure:

	##001: Indel realigment step 
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
		
	        For more details see: (https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)


		For this MK:

		targets/	BAM files
		results/	BAM realigned from this step and it's .bai index


	002: BQSR (Base quality score recalibration) step
		results/	BAM recalled from this step a recalibration tables
			*.recal-tab1 First recalibration table, result from BaseRecalibrator
			*.recal1.bam	First recalibrated BAM, result from PrintReads
			*.recal-tab2	Second recalibration table, result from BaseRecalibrator
			*.recal2.bam	Second recalibrated BAM, result from PrintReads

	003: Variant calling step
		results/	gVCFs files generated, those are not recalled
	

	004:
