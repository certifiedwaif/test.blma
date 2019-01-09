
source("Load_UScrime.Rs")
time4 <- system.time({blma_result_liang_g2         <- blma(vy, mX, prior="liang_g2"); })[3]
save(time4,blma_result_liang_g2,file=paste(datname,"_blma_result_liang_g2.Rdata",sep=""))