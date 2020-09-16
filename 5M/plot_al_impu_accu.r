#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if(length(args)<4) {
   print("usage: ./imput_accuracy.r [AF binsize] [file_base1] [file_base2] [file_base3]")
   quit(save="no")
} else {
     if(!require(data.table))
        install.packages("data.table", repos="https://cloud.r-project.org")
            bs <- as.numeric(args[1])
            bna <- args[2]
            bnb <- args[3]
            bnc <- args[4]
            outpng <- paste0("imputation",bna,bnb,bnc,"accuracy.png")
            bin <- seq(from=0.00, to=0.5, by=bs)
            frqbins <- as.data.frame(bin)
            vc <- c(0)
            chrom <- seq(from=1, to=22, by=1)
            for (chr in 1:22) {
              f <- fread(paste0("chr", chr, bna),h=T, nThread=10)
              colnames(f) <- c("maf","r2")
              for(bindex in 1:(length(bin)-1)) {
                  vc[bindex+1] <- (sum(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]])/length(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]]))
              }
              frqbins[,chrom[chr]] <- vc
            }
            resa <- cbind(bin, frqbins)

            for (chr in 1:22) { 
              f <- fread(paste0("chr", chr, bnb),h=T, nThread=10)
              colnames(f) <- c("maf","r2") 
              for(bindex in 1:(length(bin)-1)) {
                  vc[bindex+1] <- (sum(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]])/length(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]]))
              }
              frqbins[,chrom[chr]] <- vc
            }
            resb <- cbind(bin, frqbins)

            for (chr in 1:22) { 
              f <- fread(paste0("chr", chr, bnc),h=T, nThread=10)
              colnames(f) <- c("maf","r2") 
              for(bindex in 1:(length(bin)-1)) {
                  vc[bindex+1] <- (sum(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]])/length(f$r2[f$maf >= bin[bindex] & f$maf < bin[(bindex+1)]]))
              }
              frqbins[,chrom[chr]] <- vc
            }
            resc <- cbind(bin, frqbins)

            colnames(resa) <- c("MAF","Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr3","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21","Chr22")
            colnames(resb) <- c("MAF","Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr3","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21","Chr22")
            colnames(resc) <- c("MAF","Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr3","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21","Chr22")
            
            png(outpng, width=17, height=16, units="cm", res=100, points=12)
              plot(resa$MAF, resa$Chr1, xlab=paste0("Imputed MAF bin: binsize = ", bs), ylab="Mean imputation accuracy", type="l", lty=4, col=2)
              for (chr in 3:23) { 
            	lines(resa$MAF, resa[,chr], lty=4, col=2)
              }
              for (chr in 2:23) {
                lines(resb$MAF, resb[,chr], lty=4, col=3)
              }
              for (chr in 2:23) {
                lines(resc$MAF, resc[,chr], lty=4, col=4)
              }
            dev.off()
}


