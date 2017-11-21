##################################################
## QC scripts for stage 005 of mk-gatk pipeline ##

For automatic execution of QC test just run:

````
mk QC

````

## variants_by_filter.tsv ##

A tab separated table counting number of variants passing VQSR filters, and variants not passing the same filter in format:

````
Variant_category	Number_of_Variants
PASS			n
NOT_PASS		n

````
## variants_by_filter.tsv.plot.png ##

A pie chart summarizing the table from the .tsv file.
