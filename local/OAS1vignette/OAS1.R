#how I created the genePred file
#source /mnt/gluster/home/epantaleo/src/altra/scripts/utils.sh
#genePredView ~/data/annotations/refGene_hg18.txt.gz chr12:111838800-111842895 compressed > hg18.refGene.part.gp
#genePredView ~/data/annotations/refGene_hg19.txt.gz chr12:113354417-113358512 compressed > hg19.OAS1.refGene.part.gp

library(multiseq)

#region      <- "chr12:111838800-111842895" #on hg18
region <- "chr12:113354417-113358512"  #on hg19
samplesheet <- "~/src/multiseq/local/OAS1vignette/samplesheet_OAS1_bwa.txt"  
samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#keep        <- (samples$Type!="NN")
#samples     <- samples[keep,]

x           <- get.counts(samplesheet,region)
#x           <- x[keep,]
transcripts <- get.transcripts(file.path(path.package("multiseq"),"extdata", "hg19.OAS1.refGene.part.gp"), region) 

g           <- factor(samples$Type)
g           <- match(g, levels(g))-1 


OAS1<- NULL   
OAS1$x      <- x
OAS1$region <- region
OAS1$g      <- g
OAS1$read.depth <- samples$ReadDepth
OAS1$assembly   <- "hg19"
OAS1$SNP        <- "rs10774671" 
save(OAS1, file="OAS1.RData")


system.time(res         <- multiseq(x=x, g=g, minobs=1, read.depth=samples$ReadDepth))
system.time(res0         <- multiseq(x=x[g==0,], minobs=1, read.depth=samples$ReadDepth[g==0]))
system.time(res1         <- multiseq(x=x[g==1,], minobs=1, read.depth=samples$ReadDepth[g==1]))
system.time(res2         <- multiseq(x=x[g==2,], minobs=1, read.depth=samples$ReadDepth[g==2]))

pdf("OAS1all.pdf")
par(mfrow=c(5,1), oma = c(0,0,0,0) + 0.1, mar = c(2,2,1,2), mgp=c(0,1,0))

plot(res0, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype AA")
plot(res1, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype AG")
plot(res2, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype GG")
plot(res, z.threshold=2, is.xaxis=FALSE)
plot(transcripts, region)
dev.off()

subset=c(which(g==0)[seq(1,7)],which(g==1)[seq(1,7)],which(g==2)[seq(1,7)])     
OAS1 <- NULL
OAS1$x<-x[subset,]
OAS1$g      <- g[subset]
OAS1$read.depth <- samples$ReadDepth [subset]
OAS1$region <- region
OAS1$assembly   <- "hg19"
OAS1$SNP        <- "rs10774671"
save(OAS1, file="OAS1.RData")


system.time(res0         <- multiseq(x=x[which(g==0)[seq(1,7)],], minobs=1, read.depth=samples$ReadDepth[which(g==0)[seq(1,7)]]))
system.time(res1         <- multiseq(x=x[which(g==1)[seq(1,7)],], minobs=1, read.depth=samples$ReadDepth[which(g==1)[seq(1,7)]]))
system.time(res2         <- multiseq(x=x[which(g==2)[seq(1,7)],], minobs=1, read.depth=samples$ReadDepth[which(g==2)[seq(1,7)]]))
system.time(res          <- multiseq(x=x[subset,], g=g[subset], minobs=1, read.depth=samples$ReadDepth[subset]))


pdf("OAS1.pdf")
par(mfrow=c(5,1), oma = c(0,0,0,0) + 0.1, mar = c(2,2,1,2), mgp=c(0,1,0))

plot(res0, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype AA")
plot(res1, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype AG")
plot(res2, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="genotype GG")
plot(res, z.threshold=2, is.xaxis=FALSE)
plot(transcripts, region)
dev.off()

