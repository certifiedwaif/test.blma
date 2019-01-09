
source("Load_UScrime.Rs")
time8 <- system.time({blma_result_ZS_gauss_laguerre  <- blma(vy, mX, prior="zellner_siow_gauss_laguerre"); })[3]
save(time8,blma_result_ZS_gauss_laguerre,file=paste(datname,"_blma_result_ZS_gauss_laguerre.Rdata",sep=""))