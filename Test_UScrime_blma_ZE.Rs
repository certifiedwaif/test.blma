
source("Load_UScrime.Rs")
time2 <- system.time({blma_result_ZE <- blma(vy, mX, prior="ZE"); })[3]
save(time2,blma_result_ZE,file=paste(datname,"_blma_result_ZE.Rdata",sep=""))