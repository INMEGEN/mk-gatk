## pipeline mk-GATK for short variant calling
===============

## Module description and objectives

For population studies, it is important to detect genomic variants that occur in a given population.

This pipeline makes SNV and short indel variant calling from .bam files. It uses the Genome Analysis ToolKit (GATK).  
GATK offers a wide variety of tools with a primary focus on variant discovery and genotyping. Its powerful processing engine and high-performance computing features make it capable of taking on projects of any size.

Usage::

Input files should be .bam files with marked duplicates and correctly assigned ReadGroups.

Your @@kind@@ files will be on `results/` when the process ends.

# Options

@@How can you customize the analysis using environment vars or config.mk@@

# Design considerations

@@What was taken into account to build this project?@@

# Requirements

- [`mk`](http://doc.cat-v.org/bell_labs/mk/mk.pdf "A successor for `make`.")

- [`findutils`](https://www.gnu.org/software/findutils/ "Basic directory searching utilities of the GNU operating system.")

- [`coreutils`]( "basic file, shell and text manipulation utilities of the GNU operating system.")

- [@@aditional software@@](@@software URL@@ "@@Description@@")

# References

@@What documents did you used for making this module?@@

@@Where is the documentation for the software used?@@
