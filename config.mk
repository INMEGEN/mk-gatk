#Same castillo PATH
PATH=/castle/cfresno/.bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/usr/lib/plan9/bin:/usr/bin/

#GATK directory
GATK="/usr/share/java/GenomeAnalysisTK3.7.jar"

#Directorio del genoma de referencia
REF="/reference/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta"

#Directorio con variantes de confianza
INDELs="/reference/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
dbSNP="/reference/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz"
HAPMAP= ##originally "REF/1000G_omni2.5.b37.vcf" in 005; must find correct path in /reference @cluster.inmegen.gob.mx
1000G_OMNI= ##originally "REF/1000G_omni2.5.b37.vcf" in 005; must find correc path in /reference @cluster.inmegen.gob.mx

#Number of threads
NT="4"
