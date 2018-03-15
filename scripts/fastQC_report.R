#! /usr/bin/env Rscript

#########################
## check FastQC output ##
#########################
library(fastqcr)
library(ggplot2)
library(reshape2)

#qc.dir <- "fastqs/QC/Test/"
qc.dir <- "fastqs/QC/"
qc <- qc_aggregate(qc.dir)
write.table(summary(qc),
            paste0(qc.dir, "FastQC_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

qc.results <- do.call(cbind.data.frame,
			qc)
write.table(qc.results,
	paste0(qc.dir, "FastQC_results.tsv"),
	sep="\t", quote=FALSE, row.names=FALSE)


# pull out relevant sequencing QC metrics
samples <- list.files(path = qc.dir, pattern = "*.zip")

pbsq <- list()
pbsc <- list()
pbgc <- list()
seq_dup <- list()
rep_seq <- list()
dedup_pc <- list()

for (ii in 1:length(samples)) {
  qcf <- qc_read(paste0(qc.dir, samples[ii]))
  .meanq <- as.data.frame(qcf$per_base_sequence_quality[, c("Base", "Mean")])
  .meanq$Filename <- samples[ii]
  pbsq[[samples[ii]]] <- .meanq
  #.pbsc <- qcf$per_base_sequence_content[, c(2:5)]
  pbsc[[samples[ii]]] <- melt(qcf$per_base_sequence_content, id.vars="Base")
  pbgc[[samples[ii]]] <- qcf$per_sequence_gc_content
  seq_dup[[samples[ii]]] <- qcf$sequence_duplication_levels
  rep_seq[[samples[ii]]] <- qcf$overrepresented_sequences
  dedup_pc[[samples[ii]]] <- qcf$total_deduplicated_percentage
}

print("Binding results together")
seq_qual <- do.call(rbind.data.frame,
                    pbsq)

seq_qual$Base <- factor(seq_qual$Base,
                        levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10-11", "12-13",
                                 "14-15", "16-17", "18-19", "20-21", "22-23", "24-25", "26-27",
                                 "28-29", "30-31", "32-33", "34-35", "36-37", "38-39", "40-41",
                                 "42-43", "44-45", "46-47", "48-49", "50-51", "52-53", "54-55",
                                 "56-57", "58-59", "60-61", "62-63", "64-65", "66-67", "68-69",
                                 "70-71", "72-73", "74-75", "76-77", "78-79", "80-81", "82-83",
                                 "84-85", "86-87", "88-89", "90-91", "92-93", "94-95", "96-97",
                                 "98-99", "100"),
                        labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10-11", "12-13",
                                 "14-15", "16-17", "18-19", "20-21", "22-23", "24-25", "26-27",
                                 "28-29", "30-31", "32-33", "34-35", "36-37", "38-39", "40-41",
                                 "42-43", "44-45", "46-47", "48-49", "50-51", "52-53", "54-55",
                                 "56-57", "58-59", "60-61", "62-63", "64-65", "66-67", "68-69",
                                 "70-71", "72-73", "74-75", "76-77", "78-79", "80-81", "82-83",
                                 "84-85", "86-87", "88-89", "90-91", "92-93", "94-95", "96-97",
                                 "98-99", "100"))

print("Sequence Quality metrics")
write.table(seq_qual,
	gzfile(paste0(qc.dir, "FastQC_BaseQuals.tsv.gz")),
			   sep="\t", quote=FALSE, row.names=FALSE)

seq_content <- do.call(rbind.data.frame,
                       pbsc)

seq_content$FileName <- gsub(rownames(seq_content), pattern="zip.[0-9]+",
                             replacement="zip")

seq_content$Base <- factor(seq_content$Base,
                           levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10-11", "12-13",
                                    "14-15", "16-17", "18-19", "20-21", "22-23", "24-25", "26-27",
                                    "28-29", "30-31", "32-33", "34-35", "36-37", "38-39", "40-41",
                                    "42-43", "44-45", "46-47", "48-49", "50-51", "52-53", "54-55",
                                    "56-57", "58-59", "60-61", "62-63", "64-65", "66-67", "68-69",
                                    "70-71", "72-73", "74-75", "76-77", "78-79", "80-81", "82-83",
                                    "84-85", "86-87", "88-89", "90-91", "92-93", "94-95", "96-97",
                                    "98-99", "100"),
                           labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10-11", "12-13",
                                    "14-15", "16-17", "18-19", "20-21", "22-23", "24-25", "26-27",
                                    "28-29", "30-31", "32-33", "34-35", "36-37", "38-39", "40-41",
                                    "42-43", "44-45", "46-47", "48-49", "50-51", "52-53", "54-55",
                                    "56-57", "58-59", "60-61", "62-63", "64-65", "66-67", "68-69",
                                    "70-71", "72-73", "74-75", "76-77", "78-79", "80-81", "82-83",
                                    "84-85", "86-87", "88-89", "90-91", "92-93", "94-95", "96-97",
                                    "98-99", "100"))

print("Sequence nucleotide content summary")
write.table(seq_content,
	gzfile(paste0(qc.dir, "FastQC_NucContent.tsv.gz")),
				  sep="\t", row.names=FALSE, quote=FALSE)


gc_content <- do.call(cbind.data.frame,
                      pbgc)
gc_content$GC <- qcf$per_sequence_gc_content$`GC Content`
gc_content.melt <- melt(gc_content, id.vars='GC')

print("Sequence GC content summary")
write.table(gc_content,
	gzfile(paste0(qc.dir, "FastQC_GCcontent.tsv.gz")),
				  sep="\t", row.names=FALSE, quote=FALSE)

sequence_dup <- do.call(cbind.data.frame,
                        seq_dup)
sequence_dup$DupLevel <- factor(qcf$sequence_duplication_levels$`Duplication Level`,
                                levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                         ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+"),
                                labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                         ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+"))

sequence_dup.melt <- melt(sequence_dup, id.vars='DupLevel')

print("Sequence duplication levels")
write.table(sequence_dup.melt,
	gzfile(paste0(qc.dir, "FastQC_Duplicates.tsv.gz")),
				  sep="\t", row.names=FALSE, quote=FALSE)

overrep_seq <- do.call(rbind.data.frame,
                       rep_seq)

print("Overrepresented sequences summary")
write.table(overrep_seq, gzfile(paste0(qc.dir, "FastQC_OverRep.tsv.gz")),
			 sep="\t", quote=FALSE, row.names=FALSE)

dedup_seq <- do.call(cbind.data.frame,
                     dedup_pc)
dedup_seq.melt <- melt(dedup_seq)
colnames(dedup_seq.melt) <- c("FileName", "DedupPC")

print("Deduplicated proportions summary")
write.table(dedup_seq, gzfile(paste0(qc.dir, "FastQC_DedupRemain.tsv.gz")),
		       sep="\t", row.names=FALSE, quote=FALSE)