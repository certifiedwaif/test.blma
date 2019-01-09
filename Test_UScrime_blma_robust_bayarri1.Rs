
source("Load_UScrime.Rs")
time6 <- system.time({blma_result_robust_bayarri1  <- blma(vy, mX, prior="robust_bayarri1"); })[3]
save(time6,blma_result_robust_bayarri1,file=paste(datname,"_blma_result_robust_bayarri1.Rdata",sep=""))