# Author: PharmCat
# Version = 0.0.1 Alpha


Run <- function(args) {
    #Execute the FCPROCI command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("GROUP", subc="", var="group", ktype="varname"),
        spsspkg.Template("EFFECT", subc="", var="effect", ktype="varname"),
		spsspkg.Template("ALPHA", subc="", var="alpha", ktype="float", vallist=list(0.0, 1.0)),
		spsspkg.Template("CP", subc="CI", var="cicp", ktype="bool"),
		spsspkg.Template("AC", subc="CI", var="ciac", ktype="bool"),
		spsspkg.Template("SC", subc="CI", var="cisc", ktype="bool"),
		spsspkg.Template("SOC", subc="CI", var="cisoc", ktype="bool"),
		spsspkg.Template("BLACKER", subc="CI", var="cibl", ktype="bool"),
		spsspkg.Template("WALD", subc="CI", var="ciwald", ktype="bool"),
		spsspkg.Template("MN", subc="DIFF", var="dmn", ktype="bool"),
		spsspkg.Template("ACW", subc="DIFF", var="dacw", ktype="bool"),
		spsspkg.Template("ASY", subc="DIFF", var="dasy", ktype="bool"),
		spsspkg.Template("ASYY", subc="DIFF", var="dasyy", ktype="bool"),
		spsspkg.Template("NHS", subc="DIFF", var="dnhs", ktype="bool"),
		spsspkg.Template("GNS", subc="RRATIO", var="rrgns", ktype="bool"),
		spsspkg.Template("MHS", subc="RRATIO", var="rrmhs", ktype="bool"),
		spsspkg.Template("MOVER", subc="RRATIO", var="rrmover", ktype="bool"),
		spsspkg.Template("GNC", subc="RRATIO", var="rrgnc", ktype="bool"),
		spsspkg.Template("MNS", subc="ORATIO", var="ormns", ktype="bool"),
		spsspkg.Template("WOOLF", subc="ORATIO", var="orwoolf", ktype="bool"),
		spsspkg.Template("EXACT", subc="ORATIO", var="orexact", ktype="bool")
        ))

    res <- spsspkg.processcmd(oobj, args, "exec")
    
}


exec<-function(group, effect, alpha=0.95, pf=1, cicp=false, ciac=false, cisc=false, cisoc=false, cibl=false, ciwald=false, dmn=false, dacw=false, dasy=false, dasyy=false, dnhs=false, rrgns=false, rrmhs=false, rrmover=false, rrgnc=false, ormns=false, orwoolf=false, orexact=false){
	library(binGroup)
	library(PropCIs)
	library(pairwiseCI)
	
	spsspkg.StartProcedure("Proportion confidence limits")
	
	
	dataset <- spssdata.GetDataFromSPSS(variables=c(group,effect))
	
	
	table <- table(dataset)
	
	
	
	mtt   <- as.data.frame.matrix(table)
	
	print(attributes(table["dim"]))
	
	if (dim(mtt)[1] == 2 && dim(mtt)[2] == 2){
		
	} else {
		print ("Error: this is not 2X2 table")
		return ()
	}
	
	if (pf == 1) {
	mtt <- mtt[,c(2,1)]
	}
	
	
	n1 = mtt[1,1] + mtt [1,2]
	n2 = mtt[2,1] + mtt [2,2]
	
	varlabel <- spssdictionary.GetVariableLabel(group)
	if (varlabel != 0) {
		t1rd = varlabel
	} else {
		t1rd = "Group"
	}
	
	varlabel <- spssdictionary.GetVariableLabel(effect)
	if (varlabel != 0) {
		t1cd = varlabel
	} else {
		t1cd = "Effect"
	}
	
	labels <- data.frame(spssdictionary.GetValueLabels(group))
	
	lab1 <- subset(labels, values==rownames(mtt)[1])[1, "labels"]
	lab2 <- subset(labels, values==rownames(mtt)[2])[1, "labels"]

	rownames(mtt)[1] <- as.character(lab1)
	rownames(mtt)[2] <- as.character(lab2)
	
	spsspivottable.Display(mtt,
title="2 X 2 Table",
rowdim=t1rd,
hiderowdimtitle=FALSE,
coldim = t1cd,
hidecoldimtitle=FALSE,
format=6)
	
	result <- data.frame()
	
	pest = mtt[1,1]/n1
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="CP")
	result <- rbind(result, "Clopper-Pearson (Exact)" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="AC")
	result <- rbind(result, "Agresti-Coull" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="Score")
	result <- rbind(result, "Wilson Score" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="SOC")
	result <- rbind(result, "SOC" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="Blaker")
	result <- rbind(result, "Blaker" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[1,1], conf.level=alpha, method="Wald")
	result <- rbind(result, "Wald" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	colnames(result)[1] <- "Estimate"
	colnames(result)[2] <- "Lower"
	colnames(result)[3] <- "Upper"
	colnames(result)[4] <- "Alpha"
	
	spsspivottable.Display(result,
title=paste0("Proportion Confidence Limits for group: \"", lab1, "\""),
rowdim="Method",
hiderowdimtitle=FALSE,
coldim = "Intrvals",
hidecoldimtitle=FALSE,
format=2)


	result2 <- data.frame()
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="CP")
	result2 <- rbind(result2, "Clopper-Pearson (Exact)" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="AC")
	result2 <- rbind(result2, "Agresti-Coull" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="Score")
	result2 <- rbind(result2, "Wilson Score" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="SOC")
	result2 <- rbind(result2, "SOC" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="Blaker")
	result2 <- rbind(result2, "Blaker" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- binCI(n=n1,y=mtt[2,1], conf.level=alpha, method="Wald")
	result2 <- rbind(result2, "Wald" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	colnames(result2)[1] <- "Estimate"
	colnames(result2)[2] <- "Lower"
	colnames(result2)[3] <- "Upper"
	colnames(result2)[4] <- "Alpha"

	spsspivottable.Display(result2,
title=paste0("Proportion Confidence Limits for group: \"", lab2, "\""),
rowdim="Method",
hiderowdimtitle=FALSE,
coldim = "Intrvals",
hidecoldimtitle=FALSE,
format=2)
	
	#ediff <- data.frame(BinomCI(n1=n1,n2=n2, x=mtt[1,1], y=mtt[2,1],conf.level=alpha,CItype="Two.sided"))
	
	ediff <- data.frame()
	difest = mtt[1,1]/(mtt[1,1]+mtt[1,2])- mtt[2,1]/(mtt[2,1]+mtt[2,2])
	ci <- diffscoreci(x1=mtt[1,1], n1=n1, x2=mtt[2,1], n2=n2, conf.level=alpha)
	ediff <- rbind(ediff, "Miettinen-Nurminen asymptotic score" = c(0, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- wald2ci(x1=mtt[1,1], n1=n1, x2=mtt[2,1], n2=n2, conf.level=alpha, adjust="AC")
	ediff <- rbind(ediff, "Wald (Agresti-Caffo adjusted)" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- prop.test(as.matrix(mtt), correct=FALSE)	
	ediff <- rbind(ediff, "The asymptotic" = c(ci$estimate[1]-ci$estimate[2], ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- prop.test(as.matrix(mtt), correct=TRUE)	
	ediff <- rbind(ediff, "The asymptotic (Yates')" = c(ci$estimate[1]-ci$estimate[2], ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- Prop.diff(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="NHS")
	ediff <- rbind(ediff, "Newcombes Hybrid Score" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	colnames(ediff)[1] <- "Estimate"
	colnames(ediff)[2] <- "Lower"
	colnames(ediff)[3] <- "Upper"
	colnames(ediff)[4] <- "Alpha"
	
	spsspivottable.Display(ediff,
title=paste("Proportion Difference Confidence Limits"),
rowdim="Method",
hiderowdimtitle=FALSE,
coldim = "Intrvals",
hidecoldimtitle=FALSE,
format=2)
	
	rrest = (mtt[1,1]/(mtt[1,1]+mtt[1,2]))/(mtt[2,1]/(mtt[2,1]+mtt[2,2]))
	
	riskr <- data.frame()
	
	ci <- Prop.ratio(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="Score")
	riskr <- rbind(riskr, "Gart-Nam Score" = c(rrest, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- Prop.ratio(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="MNScore")
	riskr <- rbind(riskr, "Miettinen-Nurminen Score" = c(rrest, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- Prop.ratio(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="MOVER")
	riskr <- rbind(riskr, "Method of variance estimates recovery (Donner, Zou, 2012)" = c(rrest, ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- Prop.ratio(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="GNC")
	riskr <- rbind(riskr, "Crude log interval" = c(rrest, ci$conf.int[1], ci$conf.int[2], alpha))
	
	colnames(riskr)[1] <- "Estimate"
	colnames(riskr)[2] <- "Lower"
	colnames(riskr)[3] <- "Upper"
	colnames(riskr)[4] <- "Alpha"
	
	spsspivottable.Display(riskr,
title=paste("Relative Risk Confidence Limits"),
rowdim="Method",
hiderowdimtitle=FALSE,
coldim = "Intrvals",
hidecoldimtitle=FALSE,
format=2)
	
	orest = (mtt[1,1]/mtt[2,1])/(mtt[1,2]/mtt[2,2]) 
	
	oddr <- data.frame()
	ci <- orscoreci(x1=mtt[1,1], n1=n1, x2=mtt[2,1], n2=n2, conf.level=alpha)
	oddr <- rbind(oddr, "MN Score" = c(orest, ci$conf.int[1], ci$conf.int[2], alpha))

	ci <- Prop.or(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="Woolf")
	oddr <- rbind(oddr, "Adjusted Woolf interval" = c(ci$estimate , ci$conf.int[1], ci$conf.int[2], alpha))
	
	ci <- Prop.or(x=c(mtt[1,1],mtt[1,2]),y=c(mtt[2,1],mtt[2,2]), CImethod="Exact")
	oddr <- rbind(oddr, "Exact CI" = c(ci$estimate, ci$conf.int[1], ci$conf.int[2], alpha))
	
	colnames(oddr)[1] <- "Estimate"
	colnames(oddr)[2] <- "Lower"
	colnames(oddr)[3] <- "Upper"
	colnames(oddr)[4] <- "Alpha"
	
	spsspivottable.Display(oddr,
title=paste("Odd Ratio Confidence Limits"),
rowdim="Method",
hiderowdimtitle=FALSE,
coldim = "Intrvals",
hidecoldimtitle=FALSE,
format=2)

	spsspkg.EndProcedure()
    
}