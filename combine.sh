#!/bin/bash

cd output
mkdir -p temp1
mkdir -p temp2
cp ../hadd2.pl ./temp1/
cp ../hadd2.pl ./temp2/

# cp KFParticle* ./temp1/
# cd temp1
# perl hadd2.pl
# rm KFParticle*
# perl hadd2.pl
# cd ../

cp output_* ./temp2/
cd temp2
perl hadd2.pl
rm output_*
perl hadd2.pl
