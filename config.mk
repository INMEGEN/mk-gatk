#Same castillo PATH
PATH=/castle/cfresno/.bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/usr/lib/plan9/bin:/usr/bin/

#GATK directory
GATK="/usr/share/java/GenomeAnalysisTK3.7.jar"

#Directorio del genoma de referencia
REF="/reference/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta"

#Directorio con variantes de confianza
INDELs="/reference/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
dbSNP="/reference/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz"
HAPMAP="/reference/ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz" ##Should it be .vcf or vc.gz? originally "$REF/hapmap_3.3.b37.vcf" in 005.
1000G_OMNI="/reference/ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz" ##Should it be .vcf or vc.gz? Originally "$REF/1000G_omni2.5.b37.vcf" in 005.

#Number of threads
NT="4"
