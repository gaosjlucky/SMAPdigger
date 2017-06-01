############################################
# Important notice
############################################
The latest version of SMAPdigger has been moved to
https://github.com/hellbelly/SMAPdigger
Several bugs has been fixed and the architecture has been improved in the latest version.


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


## Step 2 (Configure program) 
1, modify configuration file configure 

2, Execute script "perl ./Monitor.pl -c configure -o [Output path]"
## Step 3 (Run program) 
1, cd [Output path] 

    sh RRBS_Run.sh
