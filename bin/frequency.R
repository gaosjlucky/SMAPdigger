pdf("$outdir/stats_results.pdf", width=8, height=6)
data = read.table("$outdir/R_draw_temp_file_used_once_and_never_seen.data", header=F, as.is=F)
max_num=max(data)
colnames(data) <- c("A->C","A->G","A->T","C->A","C->G","G->C","G->T","T->A","T->C","T->G")
max_num=ceiling(max_num / 1000) * 1000
barplot(as.matrix(data), xlab="Type of Mutations", ylab="No. of Mutations", col="lightblue", ylim=c(0, max_num))
dev.off()
