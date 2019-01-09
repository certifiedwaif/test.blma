
source("Load_UScrime.Rs")

time1 <- system.time({
	blma_result_BIC <- blma(vy, mX, prior="BIC"); 
})[3]

save(time1,blma_result_BIC,file=paste(datname,"_blma_result_BIC.Rdata",sep=""))