#Same castillo PATH
PATH=/castle/cfresno/.bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/usr/lib/plan9/bin:/usr/bin/

#GATK directory
GATK="/usr/share/java/GenomeAnalysisTK3.7.jar"

#Archivo de genoma de referencia
##REF="/reference/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta" ##For whole-genome
REF="/remote/reference/mitocondria/NCBI_NC_012920.1_Homo_sapiens_mitochondrion.fasta" ##For mitochondria

#Archivos de variantes de confianza
##INDELs="/reference/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
##dbSNP="/reference/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz"
##HAPMAP="/reference/ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz"
##OMNI="/reference/ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz"
##OTG="/reference/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

##For mitochondria using NC_012920 reference
INDELs="/remote/reference/mitocondria/GATK_bundle_for_NC012920/ALL.chrMT.phase3_callmom.20130502.genotypes.knowindels.vcf.recode.vcf"
OTG="/remote/reference/mitocondria/GATK_bundle_for_NC012920/ALL.chrMT.phase3_callmom.20130502.genotypes.SNVs.recode.vcf"

#Number of threads
NT="4"
