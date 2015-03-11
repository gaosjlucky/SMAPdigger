args = commandArgs(TRUE)
argc = length(args)
if (argc!=2) stop("\nusage: <ASM.out> <output>\nAuthor:zhouquan\nE-mail:zhouquan@genomics.cn")
data <- read.table(args[1], head=T, stringsAsFactors=F)
options(digits=10)
len = dim(data)[1]
for(index in 1:len){
	n11 = as.numeric(data[index,6]);
	n12 = as.numeric(data[index,7]);
	np1 = as.numeric(data[index,8]);
	np2 = as.numeric(data[index,9]);
	if (n11 + n12 != 0 && np1 + np2 != 0){
		tmp = c(n11,n12,np1,np2);
		ref_avg = n11 / (n11 + n12);
		alt_avg = np1 / (np1 + np2);
		if ((alt_avg != 0 && abs(ref_avg - alt_avg) / alt_avg >= 0.1) || (ref_avg !=0 && abs(ref_avg - alt_avg) / ref_avg >= 0.1)){
#			print("\n")
			p_value = chisq.test(matrix(tmp, ncol=2))$p.value;
			if (p_value < 0.01){
				cat(as.character(data[index,1]),as.character(data[index,2]),as.character(data[index,3]),as.character(data[index,4]),tmp,p_value,"\n",file=args[2],append=T,sep="\t");
			}
		}
	}
}
#dev.flush()
print("all done.")
