library(shiny)

load("AllEXPProcessedTCGAPDAC.RData")



TCGAsampleTable=read.delim("PublishedSampleTable.txt",sep="\t",header=T,as.is=T)
ADKS=gsub("-","_",TCGAsampleTable$Tumor.Sample.ID)

PFSraw=read.delim("PFS.txt",header=T,as.is=T,sep="\t",dec=",")
rownames(PFSraw)=gsub("-","_",paste0(toupper(PFSraw$patient_barcode),"-01A"))


library(survival)




maindf=data.frame(PFSraw[ADKS,c("PFStime","PFSevent")],OStime= L$samannot[ADKS,"os.time"],OSevent=L$samannot[ADKS,"os.event"],row.names=ADKS,stringsAsFactors = F)



onecut=function(x,cut){
	Q=quantile(x,probs=c(0,cut/100,1))
	if(length(unique(Q))<3)return(NA)

		cutx=cut(x,Q,include.lowest=T)
	levels(cutx)=c("low","high")
	cutx
}

twocuts=function(x,cuts){

	Q=quantile(x,probs=c(0,sort(cuts)/100,1))

	if(length(unique(Q))<4)return(NA)

		cutx=cut(x,Q,include.lowest=T)
	levels(cutx)=c("low",NA,"high")
	cutx
}

string2probe=function(stringinput){
	stringinput=toupper(stringinput)	
	i=which(L$probeannot$GENE_SYMBOL==stringinput)
	if(length(i)==1){
		return(rownames(L$probeannot)[i])	
	}else{
		return(NA)
	}

}


plotSurvOneCut=function(df,cut,probe){
	if(!is.data.frame(df)){
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, paste("No expression values found for gene ",probe),cex = 1.6, col = "black")
		return(NULL)
	}
	gene=strsplit(probe,"\\|")[[1]][1]
	df$geneCut=onecut(df$x,cut)
	if(all(is.na(df$geneCut))){
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, paste("Insufficiant variation of expression for gene ",gene,"at cut",cut,"%"),cex = 1.6, col = "black")
		return(NULL)
	}
	fit=survfit(Surv(time,event)~geneCut,df[ADKS,])
	pv=summary(coxph(Surv(time,event)~geneCut,df[ADKS,]))$sctest[3]
	
	plot(fit,lwd=2,col=c("#00A1D5FF","#DF8F44FF"),mark.time=TRUE,bty='l',xlab="Time in months",ylab="Survival probability",main=NULL,xlim=c(0,75))
	summ=summary(fit)

	legend(45,1, xpd=T,
		legend=paste0(c("low (n: ", "high (n: "),summ$table[,"records"]," ; median: ",signif(summ$table[,"median"],3),")"),
		col=c("#00A1D5FF","#DF8F44FF"), lty=1, cex=1,box.lty=0,lwd=2,title=paste0("log-rank: ",signif(pv,3)))
}

plotSurvTwoCut=function(df,cuts,probe){
	if(!is.data.frame(df)){
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, paste("No expression values found for gene ",probe),cex = 1.6, col = "black")
		return(NULL)
	}
	gene=strsplit(probe,"\\|")[[1]][1]
	df$geneCut=twocuts(df$x,cuts)


	if(all(is.na(df$geneCut))){
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, paste("Insufficiant variation of expression for gene ",gene,"at cut",cut,"%"),cex = 1.6, col = "black")
		return(NULL)
	}
	fit=survfit(Surv(time,event)~geneCut,df[ADKS,])
	pv=summary(coxph(Surv(time,event)~geneCut,df[ADKS,]))$sctest[3]
	
	plot(fit,lwd=2,col=c("#00A1D5FF","#DF8F44FF"),mark.time=TRUE,bty='l',xlab="Time in months",ylab="Survival probability",main=NULL,xlim=c(0,75))
	summ=summary(fit)

	legend(45,1, xpd=T,
		legend=paste0(c("low (n: ", "high (n: "),summ$table[,"records"]," ; median: ",signif(summ$table[,"median"],3),")"),
		col=c("#00A1D5FF","#DF8F44FF"), lty=1, cex=1,box.lty=0,lwd=2,title=paste0("log-rank: ",signif(pv,3)))

}





shinyServer(function(input, output) {

	getprobe<- reactive({
		string2probe(input$gene)
	})


	getdf <-reactive({
		pr=getprobe()
		if(is.na(pr)){return(NA)}
		gene=strsplit(pr,"\\|")[[1]][1]
		data.frame(ID=ADKS,time=maindf[ADKS,paste0(input$type,"time")],event=maindf[ADKS,paste0(input$type,"event")]
			,x=L$exp[pr,ADKS],stringsAsFactors = F,row.names=ADKS)
	})
	output$hvslowplot <- renderPlot({
		plotSurvOneCut( getdf(),cut=input$hvslowcut,getprobe())

	})

	output$intervalplot <- renderPlot({
		plotSurvTwoCut(getdf(), cuts=input$intervalcut,getprobe())

	})


	output$table <-renderTable({

		df=getdf()
		if(!is.data.frame(df))return(data.frame(x="Gene not found"))

			colnames(df)[2:3]=paste0(input$type,"_",colnames(df)[2:3])
		colnames(df)[4]=strsplit(getprobe(),"\\|")[[1]][1]


		df

	})

  output$downloadData <- downloadHandler(
  		

    filename = function() {
      paste0(strsplit(getprobe(),"\\|")[[1]][1],"_",input$type, ".csv")
    },
    content = function(file) {
    	df=getdf()
		if(!is.data.frame(df))return(data.frame(x="Gene not found"))

			colnames(df)[2:3]=paste0(input$type,"_",colnames(df)[2:3])
		colnames(df)[4]=strsplit(getprobe(),"\\|")[[1]][1]

      write.csv(df, file, row.names = FALSE)
    },
    contentType="csv"
  )


	output$HRtxt <- renderTable({
		df=getdf()
		if(!is.data.frame(df))return(data.frame(x="Gene not found"))
			co=coxph(Surv(time,event)~x,data=df[ADKS,])
		summco=summary(co)
		coefs=signif(summco$conf.int,3)


		show=data.frame(matrix(c("HR","log-rank",paste0(coefs[1]," [",coefs[3]," - ",coefs[4],"]"),paste0(signif(summco$sctest[3],3))),ncol=2),stringsAsFactors=F)
		colnames(show)=NULL
		rownames(show)=NULL
		show

	},colnames=F)





})