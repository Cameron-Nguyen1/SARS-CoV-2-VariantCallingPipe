library(stringr)
library(ggplot2)
library(data.table)

processGFF = function(x){
    dframe = read.csv(paste0("./GFFs/",x),sep="\t",header = FALSE,comment.char = "#")
    ind_header = grep(">[A-Za-z]",dframe$V1)
    glength = sum(nchar(dframe$V1[(ind_header+1):nrow(dframe)]))
    dframe = dframe[grepl("CDS",dframe$V3)==TRUE,]
    newc = c()
    for (i in 1:length(dframe$V9)){
        temp = str_split(dframe$V9[[i]],";")[[1]][4]
        out = str_match(temp,"=[A-Za-z0-9]*\\[")
        out2 = substr(out,2,(nchar(out)-1))
        newc = c(newc,out2)
    }
    dframe$V9 = newc
    dframe = dframe[,c(1,4,5,9)]
    colnames(dframe)=c("CHR","POS.1","POS.2","Product")
    if (dframe$Product[1] == "ORF1ab"){
        dframe$Product[1] = "ORF1a"
    }
    if (dframe$Product[2] == "ORF1ab"){
        dframe$Product[2] = "ORF1b"
    }
    dframe = cbind(dframe,list("Key"=rep("VD1",nrow(dframe))))
    #names(z) = dframe$Product
    return(dframe)
}

processVCF = function(x){
    key = str_match(x,"VD[0-9]")
    dframe = read.csv(file=paste0("./VCFs/",x),sep="\t",comment.char = "#",header=FALSE)
    decisions = c()
    p = c()
    myl_sb = str_split(dframe$V8,";")
    for (i in 1:length(myl_sb)){
        sb = myl_sb[[i]][myl_sb[[i]] %like% "^SB="]
        sbv = as.numeric(substr(sb,4,nchar(sb)))
        if (sbv <= 50){
            decisions=c(decisions,TRUE)
        }else{decisions=c(decisions,FALSE)}
        if ("INDEL" %in% myl_sb[[i]]){p=c(p,"INDEL")}else{p=c(p,"SNV")}
    }
    frame_filtered = dframe[decisions,]
    frame_filtered = cbind(frame_filtered,list("Product"=p[decisions]))
    frame_filtered = cbind(frame_filtered,list("Key"=rep(key,nrow(frame_filtered))))
    frame_filtered = frame_filtered[,c(1,2,2,9,10)]
    colnames(frame_filtered) = c("CHR","POS.1","POS.2","Product","Key")
    return(frame_filtered)
}


gff_frame = data.frame()
for (file in list.files("./GFFs/")){
    features = processGFF(file)
    gff_frame = rbind(features,gff_frame)
}
vcf_frame = data.frame()
for (file in list.files("./VCFs/")){
    features = processVCF(file)
    vcf_frame = rbind(features,vcf_frame)
}
gff_frame$CHR = factor(gff_frame$CHR,levels=unique(gff_frame$CHR))
vcf_frame$CHR = factor(vcf_frame$CHR,levels=unique(vcf_frame$CHR))
ex = rbind(gff_frame,vcf_frame)

ggplot(data=ex,aes(xmin=POS.1,xmax=POS.2))+
    geom_rect(data=ex[ex$Key=="VD1",],aes(fill = Product, ymin = as.numeric(CHR) - 0.4, ymax = as.numeric(CHR) + 0.4),color="black")+
    geom_rect(data=ex[ex$Key=="VD3",],aes(ymin = as.numeric(CHR) - 0.4, ymax = as.numeric(CHR) + 0.4),color="green")+
    #geom_rect(data=ex[ex$Key=="VD5",],aes(ymin = as.numeric(CHR) - 0.4, ymax = as.numeric(CHR) + 0.4),color="purple")+
    scale_y_continuous(breaks = 1:length(levels(ex$CHR)), labels = levels(ex$CHR)) +
    scale_fill_brewer(palette="Set3") +
    theme_minimal() +
    theme(legend.position="bottom",
          panel.background = element_rect(fill = "darkgrey"),
          panel.grid = element_blank()) +
    labs(x = "Position (BP)", y = "Sample", fill = "Gene")