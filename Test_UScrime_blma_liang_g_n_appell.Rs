
source("Load_UScrime.Rs")
timeX <- system.time({ blma_result_liang_g_n_appell <- blma(vy, mX, prior="liang_g_n_appell"); })[3]
save(timeX,blma_result_liang_g_n_appell,file=paste(datname,"_blma_result_liang_g_n_appell.Rdata",sep=""))