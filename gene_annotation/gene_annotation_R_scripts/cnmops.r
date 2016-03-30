library(cn.mops)

BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, refSeqName = c("Ha4", "Ha16"), mode = "paired", WL = 100)
#bamDataRanges <- getReadCountsFromBAM(BAMFiles, refSeqName = c("Ha4", "Ha16"), mode = "paired")
res <- cn.mops(bamDataRanges, normType = "mean")

segm <- as.data.frame(segmentation(res))
CNVs <- as.data.frame(cnvs(res))
CNVRegions <- as.data.frame(cnvr(res))

write.table(segm,file="segmentation.tsv",sep="\t",quote=FALSE)
write.table(CNVs,file="cnvs.tsv",sep="\t",quote=FALSE)	
write.table(CNVRegions,file="cnvr.tsv",quote=FALSE)

pdf("which_1.pdf")
pdf("which_2.pdf")
pdf("which_3.pdf")