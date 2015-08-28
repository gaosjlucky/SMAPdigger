# SMAPdigger
## System requirement
smapDigger works on 64-bit Unix (Linux, Ubuntu, Mac OS, etc) based systems. 

Hardware requirements

	At least one computing node equipped with at least 10 GB Memory
	
Software requirements

	GCC 4.6.0 or higher
	Perl 5.16.3 or higher
	Java 1.6.0 or higher
	zlib 1.2.8 or higher
	

## Step 1 (Setup softwares)

1. cd [smapDigger root directory]
2. Execute command "chmod -R 755 ./"
3. Execute command "./Software_Setup.sh"

## Step 2 (Setup perl modules) 
### If your cannot get access to administrator mode (Not recommanded)
1, Execute command "./Module_Setup.sh [Perl module path (Absolute Full Path)]" in user mode

2, Add the path of your perl module to environmental variable $PARL5LIB in file ~/.bashrc. For most linux system, it is usually something like "export PERL5LIB=$PERL5LIB:<Perl module path>/share/perl5/:<Perl module path>/lib64/perl5/".

3, Execute command "source ~/.bashrc"
### If your can get access to administrator mode (Strongly recommanded)
1, Execute command "./Module_Setup.sh"
## Step 3 (Data prepare)
1, Download example file package "data.tar.gz" from http://gigadb.org/dataset/100143 and save it to the root directory of SMAP. Extract the package with command "tar vxzf data.tar.gz"

2, Put file "hg19.fa" into "./data/common" directory (It could also be downloaded from http://gigadb.org/dataset/100143)

3, Check if file "hg19.chr_len.bed" is in "./data/common" directory

4, Check if file "CpGIsland.bed.seq.example" and "Upstream2k.bed.seq.example" are in "./data/element" directory

5, Check if file adapter file and raw fq files  are in "./data/sample" directory

6, Execute script "./Data_Prepare.sh" when use it first time.
## Step 4 (Configure program) 
1, modify configuration file "Eval.configure" 

2, Execute script "perl ./Monitor.pl -c Eval.configure -o [Output path]"
## Step 5 (Run program) 
1, cd [Output path] 
### If you want to run all the scripts by one step (Recommanded)
2, Execute script “./SMAP_Run.sh"
### If you want to run each script one by one
2, Execute script "perl RRBS_prepare.pl"

3, Execute script "perl RRBS_asm.pl"

4, Execute script "perl RRBS_dmr.pl"

5, Execute script "./RRBS_report.sh"
