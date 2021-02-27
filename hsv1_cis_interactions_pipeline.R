#!/usr/bin/env Rscript

##Mike Mariani UVM 2021

##Let's find the cis 4C peaks for Dave Bloom
##There are 12 samples total: 3 wt and 3 mutant for 2 viewpoints
##Note: there was no "C" because we did not sequence the
##Uninfected cells

##Samples and treatments:
##HSV1_A_1_N701_S1	WT
##HSV1_A_2_N707_S7	WT
##HSV1_B_1_N702_S2	Δctrl2
##HSV1_B_2_N708_S8	Δctrl2
##HSV1_D_1_N703_S3	Δctrl2
##HSV1_D_2_N709_S9	Δctrl2
##HSV1_E_1_N704_S4	Δctrl2
##HSV1_E_2_N710_S10	Δctrl2
##HSV1_F_1_N705_S5	WT
##HSV1_F_2_N711_S11	WT
##HSV1_G_1_N706_S6	WT
##HSV1_G_2_N712_S12	WT

library(peakC)
library(circlize)
library(ggplot2)
library(data.table)
library(dplyr)

##From the peakC paper:
##https://academic.oup.com/nar/article/46/15/e91/5003460
##After splitting the reads into viewpoint specific fastq files, 
##the primer sequence, excluding the restriction site (i.e. GATC for DpnII) 
##is trimmed from the reads. 
##The trimmed reads are subsequently mapped to the hg38 genome with bowtie2 
##using standard settings. Reads with a mapping quality of 1 or higher are retained.

##These appear to have been 85 BP reads, so after removing the 14 bp of the viewpoint
##(excluding the 6bp digestion enzyme as it recommends in the paper) we are left with 
##71bp trimmed reads

####################################### FUNCTIONS ################################################
##################################################################################################

return.score.mm <- function(score.x,score.y)
{
  if(is.na(score.x)){
    return(score.y)
  }
  else if(is.na(score.y)){
    return(score.x)
  }else{
    return(mean(c(score.x,score.y)))
  }
}

mm_bdg_to_wig <- function(bdg, genome.length, frag.width, chrom){
  bdg <- bdg[bdg[,1]==chrom,]
  ##bdg <- bdgs.smc.6a6.rep2
  mean.frag.width <- mean(bdg[,3] - bdg[,2])
  print(paste0("Mean bed width = ",mean.frag.width))
  total.rows <- nrow(bdg)
  print(paste0("Total bed rows = ",total.rows))
  nz.rows <- bdg[bdg[,4]!=0,]
  print(paste0("Number of bed regions with non-zero count = ",nrow(nz.rows)))
  total.bases <- nrow(nz.rows)*mean.frag.width
  print(paste0("Total bases = ", total.bases))
  if(nrow(nz.rows) == total.rows & total.bases < genome.length){
    print("All bed regions have nonzero count and zero interpolation is required")
  }
  chrom.vec <- character(0)
  start.vec <- numeric(0)
  end.vec   <- numeric(0)
  score.vec <- numeric(0)
  for(i in 1:nrow(bdg)){
    pos.start<-numeric(0)
    pos.end<-numeric(0)
    for(j in 1:frag.width){
      start.now <- bdg[i,2] + j
      pos.start <- c(pos.start, start.now)
      end.now <- start.now + 1
      pos.end <- c(pos.end, end.now)
    }
    chrom.vec <- c(chrom.vec, rep(bdg[i,1], each=frag.width))
    start.vec <- c(start.vec, pos.start)
    end.vec   <- c(end.vec, pos.end)
    score.vec <- c(score.vec, rep(bdg[i,4], each=frag.width))
  }
  input.frame <- data.frame(chrom=chrom.vec,
                            chromStart=start.vec,
                            chromEnd=end.vec,
                            score=score.vec,
                            stringsAsFactors = FALSE)
  colnames(input.frame) <- c("chrom", "chromStart", "chromEnd", "score")
  row.num <- nrow(input.frame)
  start.stop <- paste0(input.frame$chromStart,":",input.frame$chromEnd)
  duplicates <- duplicated(start.stop)
  dups <- length(duplicates[duplicates==TRUE])
  duplicated.frame <- input.frame[duplicates,]
  input.frame <- input.frame[!duplicates,]
  row.num.2 <- nrow(input.frame)
  diff <- row.num-(row.num-row.num.2)
  num.dup <- diff/dups + 1
  merged <- dplyr::full_join(duplicated.frame,input.frame, by=c("chrom","chromStart","chromEnd"))
  setDT(merged)[, mean := return.score.mm(score.x,score.y), by = 1:nrow(merged) ]
  merged$score = merged$mean
  merged.correct <- merged[,c(1,2,3,7)]
  ##input.frame <- data.frame(
  ##  pos=floor(bdg[,2] + ((bdg[,3]-bdg[,2])/2)),
  ##  score=bdg[,4],
  ##  stringsAsFactors = FALSE
  zeros.frame <- data.frame(
    chrom="chrHHV6A",
    chromStart=seq(from=0,by=1,to=genome.length-1),
    chromEnd=seq(from=1,by=1,to=genome.length),
    score=0,
    stringsAsFactors = FALSE
  )
  sub.zero.frame <- dplyr::anti_join(zeros.frame, merged.correct, by=c("chromStart","chromEnd"))
  sub.output.frame <- rbind(input.frame, sub.zero.frame)
  if(nrow(sub.output.frame)!=genome.length){
    print(nrow(sub.output.frame))
    print("we have an error!")
  }
  ##Wiggle uses 1-based format (even though it is UCSC)
  ##https://genome.ucsc.edu/goldenPath/help/wiggle.html
  output.frame <- data.frame( 
    pos=sub.output.frame$chromEnd,
    score=sub.output.frame$score,
    stringsAsFactors = FALSE)
  return(output.frame)
}

##PeakC return significant peaks as sequences of single bp positions
##But it is a single vector, so this function breaks continous peak output
##into  list of chunked "regions" for use with other plotting functionality.

mm_peakc_peaks_to_regions <- function(peaks){
  count <- 1
  peaks.list <-list()
  for(i in 1:length(peaks)){
    if(i==1){
      peak.start <- peaks[i]
      peak.now <- peaks[i]
    }else{
      if(peaks[i]-peak.now<=1){
        peak.now <- peaks[i]
        peak.stop <- peaks[i]
      }else if(peaks[i]-peak.now>1){
        peaks.list[[count]] <- c(peak.start,peak.stop)
        count <- count + 1
        peak.start <- peaks[i]
        peak.now <- peak.start
      }
      if(i==length(peaks)){
        peak.stop=peaks[i]
        peaks.list[[count]] <- c(peak.start,peak.stop)
      }
    }
  }
  return(peaks.list)
}

######################################## TRIM FASTQS AND MERGE RUNS ############################
################################################################################################

####Remove first 14 bp from demuxed fastqs (keep the AAGCTT hindiii site for peakC) in bash 
####for both rounds of sequencing:
####cd /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/fastqs
####for i in *.fastq; do awk '{if(NR%4==2 || NR%4==0){print substr($1,15);}else{print;}}' $i > $(basename $i ".fastq")".trimmed.peakc.fastq"; done
####cd /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs
####for i in *.fastq; do awk '{if(NR%4==2 || NR%4==0){print substr($1,15);}else{print;}}' $i > $(basename $i ".fastq")".trimmed.peakc.fastq"; done
####Move the trimmed fastqs to their own folders
##
####Round 1 files = 
##round.1.files <- c("/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_A_1_N701_S10_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_A_2_N707_S16_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_B_1_N702_S11_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_B_2_N708_S17_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_D_1_N703_S12_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_D_2_N709_S18_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_E_1_N704_S13_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_E_2_N710_S19_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_F_1_N705_S14_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_F_2_N711_S20_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_G_1_N706_S15_L002_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs_peakc/HSV1_G_2_N712_S21_L002_R1_001.trimmed.peakc.fastq")
##
####Round 2 files = 
##round.2.files <- c("/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_A_1_N701_S1_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_A_2_N707_S7_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_B_1_N702_S2_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_B_2_N708_S8_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_D_1_N703_S3_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_D_2_N709_S9_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_E_1_N704_S4_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_E_2_N710_S10_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_F_1_N705_S5_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_F_2_N711_S11_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_G_1_N706_S6_L001_R1_001.trimmed.peakc.fastq",
##"/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc/HSV1_G_2_N712_S12_L001_R1_001.trimmed.peakc.fastq")
##
##combined.output.dir <- "/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs_peakc"
##for(i in seq_along(round.1.files)){
##  system(paste0("cat ",round.1.files[[i]]," ",round.2.files[[i]]," > ",combined.output.dir,"/",basename(round.1.files[[i]])))
##}

##Navigate to the combined folder and check multiqc
##cd /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs_peakc
##for i in *.fastq; do fastqc $i ; done
##multiqc ./
##Can move the QC results to its own folder if desired

######################## Initialize Variables ##################################################
################################################################################################

##output.dir <- "/slipstream/home/mmariani/projects/hsv1_4c/output/round_1_and_round_2_combined/cis_peaks"
output.dir <- "/slipstream/home/mmariani/projects/hsv1_4c/output/round_1_and_round_2_combined/cis_peaks_regular_trimming"

setwd(output.dir)
##fastq.dir <- "/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs_peakc"
fastq.dir <- "/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs"
fastq.files <- list.files(fastq.dir,
                          pattern="fastq",
                          full.names=TRUE)

hsv1.genome.length <- 152222
vp.hsv1.1 <- c("HSV-1", 27684,  27990)
vp.hsv1.2 <- c("HSV-1", 133054,  133474)
##NCBI RefSeq human alphaherpesvirus 1 = NC_001806.2  

viewpoint.vp1 <- 27990 - 20
viewpoint.vp2 <- 133474 - 20

hsv1.genome.length <- 152222
##reads.frag.width <- 71 
reads.frag.width <- 65
hsv1.chrom <- "NC_001806.2"

##Run bowtie2 against HSV1 ref genome 
##(Make sure it is indexed for bowtie2)

bowtie2.ref <- "/slipstream/home/mmariani/references/hsv1/hsv1"
bowtie2.path <- "/slipstream/home/mmariani/programs/miniconda3/bin/bowtie2"
num.threads <- 16

for(i in seq_along(fastq.files))
{
  file.out.name <- gsub("_L001_R1_001.trimmed.peakc.fastq|_L001_R1_001.trimmed.fastq","",basename(fastq.files[[i]]))
  print(paste0("Bowtie2 aligning ... \n",fastq.files[[i]]))
  system(paste0(bowtie2.path," ",
                "-p ",num.threads," ",
                "-x ",bowtie2.ref," ", 
                "-U ",fastq.files[[i]]," ",
                "2> ",file.out.name,".bowtie2.log ",
                "| samtools view -@ ",num.threads," -bq 2 -F 4 -o ",
                 paste0(output.dir,"/",file.out.name,".mapqge1.mapped.bam -"))
         )
  ##Sort the bam files:
  system(paste0("samtools sort -@ ",num.threads," -o ",output.dir,"/",file.out.name,".mapqge1.mapped.sorted.bam ",output.dir,"/",file.out.name,".mapqge1.mapped.bam"))
  system(paste0("samtools index ",output.dir,"/",file.out.name,".mapqge1.mapped.sorted.bam"))
  system(paste0("bamCoverage -of bedgraph -b ",output.dir,"/",file.out.name,".mapqge1.mapped.sorted.bam -o ",output.dir,"/",file.out.name,".mapqge1.mapped.sorted.bdg"))
}

################################# load in bedgraphs, convert to wigs, run peakc ###############
###############################################################################################

bdgs.files <- list.files(path=output.dir,pattern=".bdg",full.names = TRUE)
bdgs.frames <- lapply(bdgs.files,read.table,sep="\t",header=FALSE,stringsAsFactors=FALSE)
for(i in 1:length(bdgs.frames))
{
  bdgs.frames[[i]]$sample <- gsub("_L002_R1_001.trimmed.mapped.sorted.bdg","",basename(bdgs.files[i]))
}
bdgs.big.frame <- do.call(rbind,bdgs.frames)

colnames(bdgs.big.frame) <- c("chrom","chromStart","chromEnd","score","sample")

bdgs.big.frame$sample <- gsub(".mapqge1.mapped.sorted.bdg","",bdgs.big.frame$sample)
                       
bdgs.plot <- ggplot(bdgs.big.frame,aes(x=floor((chromStart+chromEnd)/2),y=score,color=sample)) +
  geom_line() +
  facet_wrap(~sample,ncol=3,nrow=4) +
  theme_bw() +
  xlab("bp position") +
  ylab("score") +
  theme(legend.position = "none") +
  ggtitle("Mapped reads to VZV genome ") +
  theme(plot.title = element_text(hjust=0.5))

bdgs.plot.smoothed <- ggplot(bdgs.big.frame,aes(x=floor((chromStart+chromEnd)/2),y=stats::runmed(score,k=11),color=sample)) +
  geom_line() +
  facet_wrap(~sample,ncol=4,nrow=4) +
  theme_bw() +
  xlab("bp position") +
  ylab("score") +
  theme(legend.position = "none") +
  ggtitle("Mapped reads to VZV genome \n(smoothed: running median w/ k=11 bp") +
  theme(plot.title = element_text(hjust=0.5))

ggsave(plot=bdgs.plot,
       filename=paste0(output.dir,"/aligned.vzv.tracks.pdf"),
       width=12,
       height=8,
       device="pdf")

ggsave(plot=bdgs.plot.smoothed,
       filename=paste0(output.dir,"/aligned.vzv.tracks.soothed.median.k11.pdf"),
       width=12,
       height=8,
       device="pdf")

##HSV1_A_1_N701_S1	WT
##HSV1_A_2_N707_S7	WT
##HSV1_B_1_N702_S2	Δctrl2
##HSV1_B_2_N708_S8	Δctrl2
##HSV1_D_1_N703_S3	Δctrl2
##HSV1_D_2_N709_S9	Δctrl2
##HSV1_E_1_N704_S4	Δctrl2
##HSV1_E_2_N710_S10	Δctrl2
##HSV1_F_1_N705_S5	WT
##HSV1_F_2_N711_S11	WT
##HSV1_G_1_N706_S6	WT
##HSV1_G_2_N712_S12	WT

bdgs.vp1.wt.1  <- subset(bdgs.big.frame, sample=="HSV1_A_1_N701_S1")
bdgs.vp1.wt.2  <- subset(bdgs.big.frame, sample=="HSV1_F_1_N705_S5")
bdgs.vp1.wt.3  <- subset(bdgs.big.frame, sample=="HSV1_G_1_N706_S6")
bdgs.vp1.mut.1 <- subset(bdgs.big.frame, sample=="HSV1_B_1_N702_S2")
bdgs.vp1.mut.2 <- subset(bdgs.big.frame, sample=="HSV1_D_1_N703_S3")
bdgs.vp1.mut.3 <- subset(bdgs.big.frame, sample=="HSV1_E_1_N704_S4")

bdgs.vp2.wt.1  <- subset(bdgs.big.frame, sample=="HSV1_A_2_N707_S7")
bdgs.vp2.wt.2  <- subset(bdgs.big.frame, sample=="HSV1_F_2_N711_S11")
bdgs.vp2.wt.3  <- subset(bdgs.big.frame, sample=="HSV1_G_2_N712_S12")
bdgs.vp2.mut.1 <- subset(bdgs.big.frame, sample=="HSV1_B_2_N708_S8")
bdgs.vp2.mut.2 <- subset(bdgs.big.frame, sample=="HSV1_D_2_N709_S9")
bdgs.vp2.mut.3 <- subset(bdgs.big.frame, sample=="HSV1_E_2_N710_S10")

bdgs.vp1.wt <- list(bdgs.vp1.wt.1,
                    bdgs.vp1.wt.2,
                    bdgs.vp1.wt.3)

bdgs.vp1.mut <- list(bdgs.vp1.mut.1,
                     bdgs.vp1.mut.2,
                     bdgs.vp1.mut.3)
                  
bdgs.vp2.wt <- list(bdgs.vp2.wt.1,
                    bdgs.vp2.wt.2,
                    bdgs.vp2.wt.3)

bdgs.vp2.mut <- list(bdgs.vp2.mut.1,
                     bdgs.vp2.mut.2,
                     bdgs.vp2.mut.3)

##Make sure to input into peakC required formatting (
##'pos','score' format similar to wig format)."
##for pos I use the midpoint of the bed positions.

bdg.peakc.in.vp1.wt <- list()
for(i in 1:length(bdgs.vp1.wt)){
  bdg.peakc.in.vp1.wt[[i]] <- mm_bdg_to_wig(bdg=bdgs.vp1.wt[[i]],
                                            genome.length=hsv1.genome.length, 
                                            frag.width=reads.frag.width, 
                                            chrom=hsv1.chrom)
}

bdg.peakc.in.vp1.mut <- list()
for(i in 1:length(bdgs.vp1.mut)){
  bdg.peakc.in.vp1.mut[[i]] <- mm_bdg_to_wig(bdg=bdgs.vp1.mut[[i]],
                                             genome.length=hsv1.genome.length, 
                                             frag.width=reads.frag.width, 
                                             chrom=hsv1.chrom)
}

bdg.peakc.in.vp2.wt <- list()
for(i in 1:length(bdgs.vp2.wt)){
  bdg.peakc.in.vp2.wt[[i]] <- mm_bdg_to_wig(bdg=bdgs.vp2.wt[[i]],
                                            genome.length=hsv1.genome.length, 
                                            frag.width=reads.frag.width, 
                                            chrom=hsv1.chrom)
}

bdg.peakc.in.vp2.mut <- list()
for(i in 1:length(bdgs.vp2.mut)){
  bdg.peakc.in.vp2.mut[[i]] <- mm_bdg_to_wig(bdg=bdgs.vp2.mut[[i]],
                                             genome.length=hsv1.genome.length, 
                                             frag.width=reads.frag.width, 
                                             chrom=hsv1.chrom)
}

#################### Perform peakC analysis ###############################################
###########################################################################################

res.vp1.wt <- peakC::combined.analysis(bdg.peakc.in.vp1.wt, 
                                       num.exp = 3,
                                       vp.pos = viewpoint.vp1, 
                                       wSize=21, 
                                       minDist=0,
                                       alphaFDR = 0.1,
                                       qWr=1)

pdf(file = paste0(output.dir,"/vp1.wt.peakc.pdf"),
    height=8,
    width=8)
  peakC::plot_C(res.vp1.wt,y.max=2e5)
dev.off()

##Check results of first analysis to see if ok:
head(res.vp1.wt$peak)
max(res.vp1.wt$peak)
max.vp1.wt <- max(res.vp1.wt$peak + (0.1*max(res.vp1.wt$peak)))
plot(res.vp1.wt$dbR[,1], 
     res.vp1.wt$dbR[,2], 
     type='h', 
     ylim=c(0,max.vp1.wt), 
     xlab="chromosomal position", 
     ylab="4C signal")
lines(res.vp1.wt$dbR[,1], 
      res.vp1.wt$dbR[,4], 
      lwd=2, 
      col='blue')

res.vp1.mut <- peakC::combined.analysis(bdg.peakc.in.vp1.mut,
                                        num.exp = 3,
                                        vp.pos = viewpoint.vp1, 
                                        wSize=21, 
                                        minDist=0,
                                        alphaFDR = 0.1,
                                        qWr=1)
pdf(file = paste0(output.dir,"/vp1.mut.peakc.pdf"),
    height=8,
    width=8)
  peakC::plot_C(res.vp1.mut,y.max=4e5)
dev.off()

res.vp2.wt <- peakC::combined.analysis(bdg.peakc.in.vp2.wt,
                                        num.exp = 3,
                                        vp.pos = viewpoint.vp2, 
                                        wSize=21, 
                                        minDist=0,
                                        alphaFDR = 0.1,
                                        qWr=1)
pdf(file = paste0(output.dir,"/vp2.wt.peakc.pdf"),
    height=8,
    width=8)
  peakC::plot_C(res.vp2.wt,y.max=4e5)
dev.off()

res.vp2.mut <- peakC::combined.analysis(bdg.peakc.in.vp2.mut,
                                        num.exp = 3,
                                        vp.pos = viewpoint.vp2, 
                                        wSize=21, 
                                        minDist=0,
                                        alphaFDR = 0.1,
                                        qWr=1)
pdf(file = paste0(output.dir,"/vp2.mut.peakc.pdf"),
    height=8,
    width=8)
  peakC::plot_C(res.vp2.mut,y.max=4e5)
dev.off()

####################### Get sig. regions ###########################################################

##Can take a look at regions below if desired:
cat(res.vp1.wt$peak,  sep="\n")
cat(res.vp1.mut$peak, sep="\n")
cat(res.vp2.wt$peak,  sep="\n")
cat(res.vp2.mut$peak, sep="\n")

vp1.wt.regions  <- mm_peakc_peaks_to_regions(res.vp1.wt$peak)
vp1.mut.regions <- mm_peakc_peaks_to_regions(res.vp1.mut$peak)
vp2.wt.regions  <- mm_peakc_peaks_to_regions(res.vp2.wt$peak)
vp2.mut.regions <- mm_peakc_peaks_to_regions(res.vp2.mut$peak)

write.table(x=data.frame(chrom="NC_001806.2",
                           chromStart=do.call(rbind,vp1.wt.regions)[,1],  
                           chromEnd=do.call(rbind,vp1.wt.regions)[,2],
                           stringsAsFactors=FALSE),
            file = paste0(output.dir,"/vp1.wt.regions.bed"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep="\t")

write.table(x=data.frame(chrom="NC_001806.2",
                         chromStart=do.call(rbind,vp1.mut.regions)[,1],  
                         chromEnd=do.call(rbind,vp1.mut.regions)[,2],
                         stringsAsFactors=FALSE),
            file = paste0(output.dir,"/vp1.mut.regions.bed"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep="\t")

write.table(x=data.frame(chrom="NC_001806.2",
                         chromStart=do.call(rbind,vp2.wt.regions)[,1],  
                         chromEnd=do.call(rbind,vp2.wt.regions)[,2],
                         stringsAsFactors=FALSE),
            file = paste0(output.dir,"/vp2.wt.regions.bed"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep="\t")

write.table(x=data.frame(chrom="NC_001806.2",
                         chromStart=do.call(rbind,vp2.mut.regions)[,1],  
                         chromEnd=do.call(rbind,vp2.mut.regions)[,2],
                         stringsAsFactors=FALSE),
            file = paste0(output.dir,"/vp2.mut.regions.bed"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep="\t")

####################### Initialize starting variables for circos with hsv1 ####################
###############################################################################################

vp.hsv1.1 <- list("NC_001806.2", 27684, 27990)
vp.hsv1.2 <- list("NC_001806.2", 133054, 133474)

hsv1.gff3.bed    <- read.table(file = "/slipstream/home/mmariani/references/hsv1/hsv1.genes.bed",
                               sep = "\t",
                               header = FALSE,
                               stringsAsFactors = FALSE)

hsv1.df            <- hsv1.gff3.bed
colnames(hsv1.df)  <- c("chr","start","end","gene","strand")
hsv1.df$gene       <- hsv1.df$gene
hsv1.df$transcript <- hsv1.df$gene
hsv1.df$exon       <- 1
hsv1.df$chr        <- "NC_001806.2"
nrow(hsv1.df)
##79
which(duplicated(hsv1.df$gene))

hsv1.df$gene[62]   <- "LAT_2"
hsv1.df$gene[63]   <- "RL2_2"
hsv1.df$gene[64]   <- "RL1_2"
hsv1.df$gene[79]   <- "RS1_2"

hsv1.df.pos <- hsv1.df[hsv1.df$strand=="+",]
hsv1.df.neg <- hsv1.df[hsv1.df$strand=="-",]

hsv1.vp.df <- data.frame(do.call(rbind,list(vp.hsv1.1,
                                            vp.hsv1.2)),
                         stringsAsFactors = FALSE)

colnames(hsv1.vp.df) <- c("chr","start","end")

hsv1.vp.df$start <- as.numeric(hsv1.vp.df$start)
hsv1.vp.df$end   <- as.numeric(hsv1.vp.df$end)

hsv1.chrom.df <- data.frame(chr="NC_001806.2",
                            start=0,
                            end=152222,
                            stringsAsFactors = FALSE)

col_text="grey25"

cols.selected <- circlize::rand_color(nrow(hsv1.df), 
                                      hue = NULL, 
                                      luminosity = "random", 
                                      transparency = 0)

circos.clear()

##Set par cell padding if desired/necessary:
##circos.par("cell.padding" = c(0.02, 0, 0.02, 0))

##PeakC cis interactions
peakc.files <- list.files(output.dir,
                          pattern=".bed",
                          full.names = TRUE)

peakc.circos.beds <- lapply(peakc.files,read.table,header=FALSE,sep="\t",stringsAsFactors=FALSE)

############################# WT cis circos plot ##########################################
###########################################################################################

vp1.bed.1 <- data.frame(matrix(do.call(rbind,rep(vp.hsv1.1, times=nrow(peakc.circos.beds[[2]]))),ncol=3,byrow=TRUE),stringsAsFactors = FALSE)
vp2.bed.1 <- data.frame(matrix(do.call(rbind,rep(vp.hsv1.2, times=nrow(peakc.circos.beds[[4]]))),ncol=3,byrow=TRUE),stringsAsFactors = FALSE)

vp1.bed.1$X2 <- as.numeric(vp1.bed.1$X2) 
vp1.bed.1$X3 <- as.numeric(vp1.bed.1$X3)
vp2.bed.1$X2 <- as.numeric(vp2.bed.1$X2) 
vp2.bed.1$X3 <- as.numeric(vp2.bed.1$X3)

vp1.bed.2 <- peakc.circos.beds[[2]]
vp2.bed.2 <- peakc.circos.beds[[4]]
vp1.bed.2$V1 <- "NC_001806.2"
vp2.bed.2$V1 <- "NC_001806.2"

colnames(vp1.bed.1) <- c("chr","start","end")
colnames(vp1.bed.2) <- c("chr","start","end")
colnames(vp2.bed.1) <- c("chr","start","end")
colnames(vp2.bed.2) <- c("chr","start","end")

dev.off()
circos.clear()
pdf(file=paste0(output.dir,"/hsv1.wt.circos.plot.pdf"),
    width=4,
    height=4)
circos.genomicInitialize(data=hsv1.chrom.df, 
                         plotType = c("axis")
)
##circos.genomicLabels(hsv1.chrom.df,
##                     labels.column = 1, 
##                     side = "outside", 
##                     niceFacing = TRUE)
##circos.initializeWithIdeogram(plotType = NULL)
##circos.trackPlotRegion(ylim = c(0, 1))

o.cell.padding = circos.par("cell.padding")
arrow.y <- 0.05 
arrow.width <- 0.05
arrow.head.width <- 0.15 
arrow.head.length <- 0.05 
genes.cex <- 0.15
sleepy.time <- 0.01
track.height <- 0.05
current.track <- 2
n.tracks <- 2
pos.high <- 0
for(i in seq_along(unique(hsv1.df$gene))){
  ##print(i)
  if(i!=1){
    ##if(i>5){stop()}
    if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track < n.tracks){
      current.track <- current.track + 1
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##if(direction=="-"){stop()}
      ##print(direction)
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track == n.tracks){
      current.track <- current.track + 1
      circos.trackPlotRegion(ylim = c(0, 1), 
                             track.height = track.height, 
                             bg.border = NA, 
                             ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
                             track.index = current.track
      )
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track,
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      n.tracks <- n.tracks + 1
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else{
      current.track <- 2
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end,hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=2
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      pos.high <- pos2
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }
  }else{
    circos.trackPlotRegion(ylim = c(0, 1), 
                           track.height = track.height, 
                           bg.border = NA, 
                           track.index=2
                           ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
    )
    ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0, 
    ##                track.index = current.track)
    ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0.1, 
    ##                track.index = current.track)
    pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
    pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
    direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
    ##circos.lines(x=c(pos1,pos2),
    ##             y=c(0,0),
    ##             ##clock.wise = TRUE, 
    ##             col = "red", 
    ##             border = "black",
    ##             lwd=3,
    ##             track.index=2)
    circos.arrow(x1 = pos1, 
                 x2 = pos2, 
                 y = arrow.y, 
                 width = arrow.width, 
                 arrow.head.width = arrow.head.width, 
                 arrow.head.length = cm_x(arrow.head.length),
                 arrow.position = ifelse(direction=="+","end","start"),
                 col = ifelse(direction=="+","red","blue"),
                 track.index=2
                 ##tail = tail[CELL_META$sector.numeric.index]
    )
    circos.text(x=(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end+hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end)/2,
                y=-1,
                label=hsv1.df$gene[i],
                ##pos=1,
                cex=genes.cex,
                track.index = 2##,
                ##adj=c(0,0.01)
                ##niceFacing=TRUE,
                ##facing="clockwise"
    )
    pos.high <- pos2
    print(c(i,current.track,unique(hsv1.df$gene)[i]))
    ##stop()
    Sys.sleep(sleepy.time)
    ##stop()
  }
}
##Viewpoints##############################
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = track.height, 
                       bg.border = "black", 
                       track.index=n.tracks+1
                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
)
for(i in 1:2){
  print(i)
  circos.rect(xleft=hsv1.vp.df$start[i], 
              ybottom=0, 
              xright=hsv1.vp.df$end[i], 
              ytop=1,
              track.index=n.tracks+1,
              ##sector.index = get.cell.meta.data("sector.index"),
              ##track.index = get.cell.meta.data("track.index"),
              col="green",
              rot = 0)
  circos.text(x=((hsv1.vp.df$start[i] + hsv1.vp.df$end[i])/2)+1000,
              ##y=-1,
              y=0.5,
              label=paste0("VP ",i),
              ##pos=2,
              cex=0.25,
              track.index = n.tracks+1##,
              ##adj=c(0,0.01)
              ##niceFacing=TRUE,
              ##facing="clockwise"
  )
}
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = 1, 
                       bg.border = NA, 
                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
                       track.index = n.tracks+2
)
col2 <- c("blue","purple")
col2.vec <- Vectorize(adjustcolor)(col2, alpha.f = 0.5)
circos.genomicLink(vp1.bed.1, vp1.bed.2, col = col2.vec[1])
circos.genomicLink(vp2.bed.1, vp2.bed.2, col = col2.vec[2])
##circos.link(sector.index1 = "hsv1",
##            point1=c(1,100),
##            sector.index2 = "hsv1",
##            point2=c(120000,121000)
##            )
##circos.genomicLink(bed1, bed2, col="black", border = "black", track.index=1)
##text(0,0,"hsv1")

title("HSV-1 Wild Type 4C Cis Interactions")

dev.off()

###################################### Mutant Cis Circos Plot ##################################
################################################################################################

vp1.bed.1 <- data.frame(matrix(do.call(rbind,rep(vp.hsv1.1, times=nrow(peakc.circos.beds[[1]]))),ncol=3,byrow=TRUE),stringsAsFactors = FALSE)
vp2.bed.1 <- data.frame(matrix(do.call(rbind,rep(vp.hsv1.2, times=nrow(peakc.circos.beds[[3]]))),ncol=3,byrow=TRUE),stringsAsFactors = FALSE)

vp1.bed.1$X2 <- as.numeric(vp1.bed.1$X2) 
vp1.bed.1$X3 <- as.numeric(vp1.bed.1$X3)
vp2.bed.1$X2 <- as.numeric(vp2.bed.1$X2) 
vp2.bed.1$X3 <- as.numeric(vp2.bed.1$X3)

vp1.bed.2 <- peakc.circos.beds[[1]]
vp2.bed.2 <- peakc.circos.beds[[3]]
vp1.bed.2$V1 <- "NC_001806.2"
vp2.bed.2$V1 <- "NC_001806.2"

colnames(vp1.bed.1) <- c("chr","start","end")
colnames(vp1.bed.2) <- c("chr","start","end")
colnames(vp2.bed.1) <- c("chr","start","end")
colnames(vp2.bed.2) <- c("chr","start","end")

dev.off()
circos.clear()
pdf(file=paste0(output.dir,"/hsv1.mut.circos.plot.pdf"),
    width=4,
    height=4)
circos.genomicInitialize(data=hsv1.chrom.df, 
                         plotType = c("axis")
)
##circos.genomicLabels(hsv1.chrom.df,
##                     labels.column = 1, 
##                     side = "outside", 
##                     niceFacing = TRUE)
##circos.initializeWithIdeogram(plotType = NULL)
##circos.trackPlotRegion(ylim = c(0, 1))

o.cell.padding = circos.par("cell.padding")
arrow.y <- 0.05 
arrow.width <- 0.05
arrow.head.width <- 0.15 
arrow.head.length <- 0.05 
genes.cex <- 0.15
sleepy.time <- 0.01
track.height <- 0.05
current.track <- 2
n.tracks <- 2
pos.high <- 0
for(i in seq_along(unique(hsv1.df$gene))){
  ##print(i)
  if(i!=1){
    ##if(i>5){stop()}
    if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track < n.tracks){
      current.track <- current.track + 1
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##if(direction=="-"){stop()}
      ##print(direction)
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track == n.tracks){
      current.track <- current.track + 1
      circos.trackPlotRegion(ylim = c(0, 1), 
                             track.height = track.height, 
                             bg.border = NA, 
                             ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
                             track.index = current.track
      )
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track,
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      n.tracks <- n.tracks + 1
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else{
      current.track <- 2
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end,hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=2
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=-1,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      pos.high <- pos2
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }
  }else{
    circos.trackPlotRegion(ylim = c(0, 1), 
                           track.height = track.height, 
                           bg.border = NA, 
                           track.index=2
                           ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
    )
    ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0, 
    ##                track.index = current.track)
    ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0.1, 
    ##                track.index = current.track)
    pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
    pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
    direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
    ##circos.lines(x=c(pos1,pos2),
    ##             y=c(0,0),
    ##             ##clock.wise = TRUE, 
    ##             col = "red", 
    ##             border = "black",
    ##             lwd=3,
    ##             track.index=2)
    circos.arrow(x1 = pos1, 
                 x2 = pos2, 
                 y = arrow.y, 
                 width = arrow.width, 
                 arrow.head.width = arrow.head.width, 
                 arrow.head.length = cm_x(arrow.head.length),
                 arrow.position = ifelse(direction=="+","end","start"),
                 col = ifelse(direction=="+","red","blue"),
                 track.index=2
                 ##tail = tail[CELL_META$sector.numeric.index]
    )
    circos.text(x=(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end+hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end)/2,
                y=-1,
                label=hsv1.df$gene[i],
                ##pos=1,
                cex=genes.cex,
                track.index = 2##,
                ##adj=c(0,0.01)
                ##niceFacing=TRUE,
                ##facing="clockwise"
    )
    pos.high <- pos2
    print(c(i,current.track,unique(hsv1.df$gene)[i]))
    ##stop()
    Sys.sleep(sleepy.time)
    ##stop()
  }
}
##Viewpoints##############################
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = track.height, 
                       bg.border = "black", 
                       track.index=n.tracks+1
                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
)
for(i in 1:2){
  print(i)
  circos.rect(xleft=hsv1.vp.df$start[i], 
              ybottom=0, 
              xright=hsv1.vp.df$end[i], 
              ytop=1,
              track.index=n.tracks+1,
              ##sector.index = get.cell.meta.data("sector.index"),
              ##track.index = get.cell.meta.data("track.index"),
              col="green",
              rot = 0)
  circos.text(x=((hsv1.vp.df$start[i] + hsv1.vp.df$end[i])/2)+1000,
              ##y=-1,
              y=0.5,
              label=paste0("VP ",i),
              ##pos=2,
              cex=0.25,
              track.index = n.tracks+1##,
              ##adj=c(0,0.01)
              ##niceFacing=TRUE,
              ##facing="clockwise"
  )
}
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = 1, 
                       bg.border = NA, 
                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
                       track.index = n.tracks+2
)
col2 <- c("red","orange")
col2.vec <- Vectorize(adjustcolor)(col2, alpha.f = 0.5)
circos.genomicLink(vp1.bed.1, vp1.bed.2, col = col2.vec[1])
circos.genomicLink(vp2.bed.1, vp2.bed.2, col = col2.vec[2])
##circos.link(sector.index1 = "hsv1",
##            point1=c(1,100),
##            sector.index2 = "hsv1",
##            point2=c(120000,121000)
##            )
##circos.genomicLink(bed1, bed2, col="black", border = "black", track.index=1)
##text(0,0,"hsv1")

title("HSV-1 Mutant 4C Cis Interactions")

dev.off()

############# Make a circos plot just emphasizing the genes and VP positions ################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

dev.off()
circos.clear()
output.dir <- "/slipstream/home/mmariani/projects/hsv1_4c/output/round_1_and_round_2_combined"
pdf(file=paste0(output.dir,"/virus.features.circos.plot.pdf"),
    width=4,
    height=4)
circos.genomicInitialize(data=hsv1.chrom.df, 
                         plotType = c("axis")
)
##circos.genomicLabels(hsv1.chrom.df,
##                     labels.column = 1, 
##                     side = "outside", 
##                     niceFacing = TRUE)
##circos.initializeWithIdeogram(plotType = NULL)
##circos.trackPlotRegion(ylim = c(0, 1))

o.cell.padding = circos.par("cell.padding")
arrow.y <- 0.1 
arrow.width <- 0.15
arrow.head.width <- 0.25 
arrow.head.length <- 0.1 
genes.cex <- 0.25
sleepy.time <- 0.01
track.height <- 0.1
current.track <- 2
n.tracks <- 2
pos.high <- 0
text.y <- -0.25
adjustment.step <- -0.25
adjustment <- 0
for(i in seq_along(unique(hsv1.df$gene))){
  ##print(i)
  if(i!=1){
    ##if(i>5){stop()}
    if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track < n.tracks){
      adjustment <- 0
      current.track <- current.track + 1
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##if(direction=="-"){stop()}
      ##print(direction)
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=text.y,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else if(subset(hsv1.df, gene==unique(hsv1.df$gene)[i])$start <= pos.high & current.track == n.tracks){
      adjustment <- 0
      current.track <- current.track + 1
      circos.trackPlotRegion(ylim = c(0, 1), 
                             track.height = track.height, 
                             bg.border = NA, 
                             ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
                             track.index = current.track
      )
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=current.track,
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y= text.y,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      n.tracks <- n.tracks + 1
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }else{
      adjustment <- adjustment + adjustment.step
      current.track <- 2
      ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end,hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                0, 
      ##                track.index = current.track)
      ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
      ##                1, 
      ##                track.index = current.track)
      ##draw.sector(pos1[1, "theta"], 
      ##            pos2[1, "theta"], 
      ##            clock.wise = TRUE, 
      ##            col = "red", 
      ##            border = NA)
      pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
      pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
      direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
      ##circos.lines(x=c(pos1,pos2),
      ##             y=c(0,0),
      ##             ##clock.wise = TRUE, 
      ##             col = "red", 
      ##             border = "black",
      ##             lwd=3,
      ##             track.index=current.track)
      circos.arrow(x1 = pos1, 
                   x2 = pos2, 
                   y = arrow.y, 
                   width = arrow.width, 
                   arrow.head.width = arrow.head.width, 
                   arrow.head.length = cm_x(arrow.head.length), 
                   arrow.position = ifelse(direction=="+","end","start"),
                   col = ifelse(direction=="+","red","blue"),
                   track.index=2
                   ##tail = tail[CELL_META$sector.numeric.index]
      )
      circos.text(x=(hsv1.df[i,]$start+hsv1.df[i,]$end)/2,
                  y=text.y + adjustment,
                  label=hsv1.df$gene[i],
                  ##pos=1,
                  cex=genes.cex,
                  track.index = current.track)
      pos.high <- pos2
      print(c(i,current.track,unique(hsv1.df$gene)[i]))
      ##stop()
      Sys.sleep(sleepy.time)
    }
  }else{
    adjustment <- 0
    circos.trackPlotRegion(ylim = c(0, 1), 
                           track.height = track.height, 
                           bg.border = NA, 
                           track.index=2
                           ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
    )
    ##pos1 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0, 
    ##                track.index = current.track)
    ##pos2 = circlize(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end-hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start, 
    ##                0.1, 
    ##                track.index = current.track)
    pos1 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$start
    pos2 = hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end
    direction <- hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$strand
    ##circos.lines(x=c(pos1,pos2),
    ##             y=c(0,0),
    ##             ##clock.wise = TRUE, 
    ##             col = "red", 
    ##             border = "black",
    ##             lwd=3,
    ##             track.index=2)
    circos.arrow(x1 = pos1, 
                 x2 = pos2, 
                 y = arrow.y, 
                 width = arrow.width, 
                 arrow.head.width = arrow.head.width, 
                 arrow.head.length = cm_x(arrow.head.length),
                 arrow.position = ifelse(direction=="+","end","start"),
                 col = ifelse(direction=="+","red","blue"),
                 track.index=2
                 ##tail = tail[CELL_META$sector.numeric.index]
    )
    circos.text(x=(hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end+hsv1.df[hsv1.df$gene==unique(hsv1.df$gene)[i],]$end)/2,
                y=text.y,
                label=hsv1.df$gene[i],
                ##pos=1,
                cex=genes.cex,
                track.index = 2##,
                ##adj=c(0,0.01)
                ##niceFacing=TRUE,
                ##facing="clockwise"
    )
    pos.high <- pos2
    print(c(i,current.track,unique(hsv1.df$gene)[i]))
    ##stop()
    Sys.sleep(sleepy.time)
    ##stop()
  }
}
##CTCF##############################
##circos.trackPlotRegion(ylim = c(0, 1), 
##                       track.height = track.height, 
##                       bg.border = "black", 
##                       track.index=n.tracks+1
##                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
##)
##circos.lines(x=(ctcf.infected.rep1.df.ds$chromStart + ctcf.infected.rep1.df.ds$chromEnd)/2,
##             y=ctcf.infected.rep1.df.ds$score,
##             ##label=paste0("VP ",i),
##             ##pos=1,
##             ##cex=0.3,
##             track.index = n.tracks+1##,
##             ##adj=c(0,0.01)
##             ##niceFacing=TRUE,
##             ##facing="clockwise")
##)
##Viewpoints##############################
##RNA##############################
##circos.trackPlotRegion(ylim = c(0, 1), 
##                       track.height = track.height, 
##                       bg.border = "black", 
##                       track.index=n.tracks+2
##                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
##)
##circos.lines(x=(rna.seq.track.df$chromStart + rna.seq.track.df$chromEnd)/2,
##             y=rna.seq.track.df$score,
##             ##label=paste0("VP ",i),
##             ##pos=1,
##             ##cex=0.3,
##             track.index = n.tracks+2##,
##             ##adj=c(0,0.01)
##             ##niceFacing=TRUE,
##             ##facing="clockwise")
##)
##Viewpoints##############################
circos.trackPlotRegion(ylim = c(0, 1), 
                       track.height = track.height, 
                       ##bg.border = "black", 
                       bg.lty=2,
                       track.index=n.tracks+1
                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
)
for(i in 1:2){
  print(i)
  circos.rect(xleft=hsv1.vp.df$start[i], 
              ybottom=0, 
              xright=hsv1.vp.df$end[i], 
              ytop=1,
              track.index=n.tracks+1,
              ##sector.index = get.cell.meta.data("sector.index"),
              ##track.index = get.cell.meta.data("track.index"),
              col="green",
              rot = 0)
  circos.text(x=(hsv1.vp.df$start[i] + hsv1.vp.df$end[i])/2,
              y=-1,
              label=paste0("VP ",i),
              ##pos=1,
              cex=0.3,
              track.index = n.tracks+1##,
              ##adj=c(0,0.01)
              ##niceFacing=TRUE,
              ##facing="clockwise"
  )
}
##circos.trackPlotRegion(ylim = c(0, 1), 
##                       track.height = 1, 
##                       bg.border = NA, 
##                       ##cell.padding = c(0, o.cell.padding[2], 0, o.cell.padding[4])
##                       track.index = n.tracks+4
##)
##circos.genomicLink(vp1.bed.1, vp1.bed.2, col = "red")
##circos.genomicLink(vp2.bed.1, vp2.bed.2, col = "green")
##circos.genomicLink(vp3.bed.1, vp3.bed.2, col = "blue")
##circos.genomicLink(vp4.bed.1, vp4.bed.2, col = "purple")
##circos.link(sector.index1 = "hsv1",
##            point1=c(1,100),
##            sector.index2 = "hsv1",
##            point2=c(120000,121000)
##            )
##circos.genomicLink(bed1, bed2, col="black", border = "black", track.index=1)

text(0,0,"HSV-1")

title("HSV-1 Genome Characterization")

dev.off()

##circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
##                   border = NA)
##bed1 = generateRandomBed(nr = 100)
##bed1 = bed1[sample(nrow(bed1), 20), ]
##bed2 = generateRandomBed(nr = 100)
##bed2 = bed2[sample(nrow(bed2), 20), ]
##circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
##circos.initializeWithIdeogram()
##circos.genomicLink(bed1, bed2, col = sample(1:5, 20, replace = TRUE), border = NA)
##circos.clear()
##bed1$chr <- "NC_001348.1"
##bed2$chr <- "NC_001348.1"
