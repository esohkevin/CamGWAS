#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args) < 1) {
   print("",quote=F)
   print("Usage: assoc.R [assoc_result]",quote=F)
   print("assoc_result: should have the following columns; CHR, SNP, BP, P",quote=F)
   print("",quote=F)
   quit(save="no")
} else {
     f <- args[1]
     oMan <- paste0(f,".assoc.png")
     oQq <- paste0(f,".qq.png")
     #?iconv()

     if(!require(data.table))
         install.packages("data.table", repos="https://cloud.r-project.org")
     if(!require(qqman))
         install.packages("qqman", repos="https://cloud.r-project.org")

     #con <- file(description = "qc-camgwas-updated.stats.gz", "r+b")
     #con <- socketConnection(host = "kevine@delgeme.icermali.org", port = 22, server = T, 
     #                 blocking = F, open = "w+", encoding = "", timeout = 100000)
     #con <- gzfile(description = "qc-camgwas-updated.stats.gz", "r+b")

     #bolt_assoc <- readBin(con, numeric(), size = 4)

     #dat <- "qc-camgwas-updated.stats.gz"

     #con <- socketConnection(port = 22, blocking = TRUE)
     #writeLines(paste0(system("whoami", intern = TRUE), "\r"), con)
     #gsub(" *$", "", readLines(con))
     #close(con)

     #assoc <- readBin(dat, n = 20L)
     #head(assoc)
     #View(assoc)

     #assoc <- fread("qc-camgwas-updated.stats", header = T, data.table = F)
     #attach(assoc)
     #assoc$P <- assoc$P_BOLT_LMM
     #png("bolt_assoc.png", height = 480, width = 700, units = "px", res = NA, pointsize = 12)
     #manhattan(assoc)
     #dev.off()
     #qq(assoc$P)

     assoc <- fread(f, h=T, data.table=F, fill=T, nThread = 30)
     attach(assoc)
     png(oMan)
     manhattan(assoc,genomewideline = -log10(5e-08))
     dev.off()
     lamd <- as.numeric(median(qchisq(assoc$P, df=1, lower.tail = F), na.rm = T)/0.456)
     print(paste0("Lambda: ", lamd))
     png(oQq)
     qq(assoc$P)
     text(2, 5, expression(lambda[GC] == lamd))
     dev.off()
}

