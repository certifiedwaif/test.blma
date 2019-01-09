
source("Load_UScrime.Rs")
time7 <- system.time({blma_result_robust_bayarri2  <- blma(vy, mX, prior="robust_bayarri2"); })[3]
save(time7,blma_result_robust_bayarri2,file=paste(datname,"_blma_result_robust_bayarri2.Rdata",sep=""))