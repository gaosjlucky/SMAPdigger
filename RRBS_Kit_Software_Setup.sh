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
tar "vxzf" "bsmap-2.74.tgz"
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
# BisSNP
###############################################
echo "Install BisSNP..."
cp "-f" "./tools/BisSNP-0.82.2.jar" "./bin/BisSNP-0.82.2.jar"
###############################################
# clear
###############################################
echo "Clear..."
cd "tools"
rm "-rf" "samtools-0.1.19"
rm "-rf" "bsmap-2.74"
rm "-rf" "bowtie2-2.2.2"
echo "Complete."
