#!/bin/sh
###############################################
# Samtools
###############################################
echo "Install samtools..."
cd "tools"
tar "-jxvf" "./samtools-0.1.19.tar.bz2"
cd "./samtools-0.1.19"
make
cp "-f" "./samtools" "../../bin/samtools"
cp "-f" "./bcftools/bcftools" "../../bin/bcftools"
cp "-f" "./bcftools/vcfutils.pl" "../../bin/vcfutils.pl"
cd ".."
cd ".."
###############################################
# BSMAP
###############################################
echo "Install BSMAP..."
cd "tools"
cd "bsmap-2.74"
make
cp "-f" "./bsmap" "../../bin/bsmap"
cp "-f" "./sam2bam.sh" "../../bin/sam2bam.sh"
cp "-f" "./methratio.py" "../../bin/methratio.py"
cd ".."
cd ".."
###############################################
# Bowtie2
###############################################
echo "Install Bowtie2..."
cd "tools"
unzip "-o" "bowtie2-2.2.2-source.zip"
cd "bowtie2-2.2.2"
make
cp "-f" bowtie2* "../../bin"
chmod "+x" ../../bin/bowtie2*
cd ".."
cd ".."
###############################################
# picard
###############################################
echo "Install picard..."
cd "tools"
unzip "-o" "-d" "../bin" "picard-tools-1.114.zip"
cd ".."
###############################################
# Bs-snper
###############################################
echo "Install Bs-snper"
cd "tools"
unzip "-o" "BS-Snper-master.zip"
cd  "BS-Snper-master"
chmod "+x" BS-Snper.sh
./BS-Snper.sh
cp "-f" "./BS-Snper.pl" "../../bin/"
cp "-f" "./chrLenExtract" "../../bin/"
cp "-f" "./rrbsSnp" "../../bin/"
cd ".."
cd ".."
###############################################
# DMR
###############################################
echo "Install dmr"
cd "tools"
unzip "-o" "dmr.zip"
cd "dmr"
g++ -O2 -o dmr dmr.cpp DivideRange.cpp UnionFile.cpp CoreSearch.cpp cdf.cpp -lm
cp "-f" "./dmr" "../../bin/"
cd ".."
cd ".."
###############################################
# Chmod
###############################################
cd "bin/FastQC/"
chmod "+x" "fastqc"
cd ".."
cd ".."
###############################################
# clear
###############################################
echo "Clear..."
cd "tools"
rm "-rf" "samtools-0.1.19"
rm "-rf" "bowtie2-2.2.2"
rm "-rf" "dmr"
rm "-rf" "Bs-snper"
echo "Complete."
