
source("Load_UScrime.Rs")
time5 <- system.time({blma_result_liang_g_n_quad   <- blma(vy, mX, prior="liang_g_n_quad"); })[3]
save(time5,blma_result_liang_g_n_quad,file=paste(datname,"_blma_result_liang_g_n_quad.Rdata",sep=""))