# PID of current job: 3577287
mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-SaveTransformedData(mSet)
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
