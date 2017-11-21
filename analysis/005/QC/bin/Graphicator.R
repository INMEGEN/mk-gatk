args = commandArgs(trailingOnly=TRUE)

wg.variants.df <- read.table(file=args[1], header = T, sep = "\t", stringsAsFactors = F)

png(args[2])

pie(wg.variants.df$Number_of_Variants, labels = wg.variants.df$Variant_category, main="Variants Passing VQSR filter", col=c("green2","red4"))

dev.off()
