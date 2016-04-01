library(cn.mops)

BAMFiles <- list.files(pattern="sort.bam$")

# Setting the window length is important because the default will be extremely large.
bamDataRanges <- getReadCountsFromBAM(BAMFiles, refSeqName = c("Ha4", "Ha16"), mode = "paired", WL = 100)
res <- cn.mops(bamDataRanges, normType = "mean")
res <- calcIntegerCopyNumbers(res)

segm <- as.data.frame(segmentation(res))
CNVs <- as.data.frame(cnvs(res))
CNVRegions <- as.data.frame(cnvr(res))

write.table(segm,file="segmentation.tsv",sep="\t",quote=FALSE)
write.table(CNVs,file="cnvs.tsv",sep="\t",quote=FALSE)	
write.table(CNVRegions,file="cnvr.tsv",quote=FALSE)

# Look in the cnvs file to get the CNV IDs. The code below prints only the first called region
plot(res, which=1)

# By default this prints all chromosomes in the same plot.
segplot(res)