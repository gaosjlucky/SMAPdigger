#!/bin/sh
echo "Starting Module Setup...";
if [ $# > 0 ] 
then
	myPath="PREFIX=$1"
	echo "PREFIX=$1";
	if [ ! -d "$1" ] 
	then
		echo "Making directory $1...";
		mkdir "-p" "$1"
	fi
else
	myPath=""
fi
###############################################
# Perl modules
###############################################
echo "Install perl modules..."
cd "./modules"
# Getopt::Long
tar "vxzf" "Getopt-Long-2.42.tar.gz"
cd "Getopt-Long-2.42"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Getopt-Long-2.42"
# Compress-Raw-Bzip2-2.064
tar "vxzf" "Compress-Raw-Bzip2-2.064.tar.gz"
cd "Compress-Raw-Bzip2-2.064"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Compress-Raw-Bzip2-2.064"
# Compress::Zlib
tar "vxzf" "IO-Compress-2.064.tar.gz"
cd "IO-Compress-2.064"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "IO-Compress-2.064"
# Compress::Raw::Zlib
tar "vxzf" "Compress-Raw-Zlib-2.065.tar.gz"
cd "Compress-Raw-Zlib-2.065"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Compress-Raw-Zlib-2.065"
# Statistics::Descriptive
tar "vxzf" "Statistics-Descriptive-3.0607.tar.gz"
cd "Statistics-Descriptive-3.0607"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Statistics-Descriptive-3.0607"
# Statistics::Distributions
tar "vxzf" "Statistics-Distributions-1.02.tar.gz"
cd "Statistics-Distributions-1.02"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Statistics-Distributions-1.02"
# Statistics::TTest
tar "vxzf" "Statistics-TTest-1.1.0.tar.gz"
cd "Statistics-TTest-1.1.0"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Statistics-TTest-1.1.0"
# Math::CDF
tar "vxzf" "Math-CDF-0.1.tar.gz"
cd "Math-CDF-0.1"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Math-CDF-0.1"
# Text::NSP
tar "vxzf" "Text-NSP-1.27.tar.gz"
cd "Text-NSP-1.27"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Text-NSP-1.27"
# Parallel::ForkManager
tar "vxzf" "Parallel-ForkManager-1.06.tar.gz"
cd "Parallel-ForkManager-1.06"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Parallel-ForkManager-1.06"
# PerlIO::gzip
tar "vxzf" "PerlIO-gzip-0.18.tar.gz"
cd "PerlIO-gzip-0.18"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "PerlIO-gzip-0.18"
# List::MoreUtils
tar "vxzf" "List-MoreUtils-0.33.tar.gz"
cd "List-MoreUtils-0.33"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "List-MoreUtils-0.33"
# Test::Simple
tar "vxzf" "Test-Simple-1.001003.tar.gz"
cd "Test-Simple-1.001003"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Test-Simple-1.001003"
# Scalar::List-Utils
tar "vxzf" "Scalar-List-Utils-1.39.tar.gz"
cd "Scalar-List-Utils-1.39"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Scalar-List-Utils-1.39"
# Test::Harness
tar "vxzf" "Test-Harness-3.32.tar.gz"
cd "Test-Harness-3.32"
perl "Makefile.PL" "$myPath"
make
make "install"
cd ".."
rm "-rf" "Test-Harness-3.32"
###############################################