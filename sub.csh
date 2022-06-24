#!/bin/csh

set MainDir=`pwd`

set XmlDir=./xml
if(! -e $XmlDir) exit
cd $XmlDir
set SubXml=sub.xml
if(-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=40
set nFileTotal=all

# print xml file
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
echo \<command\>$MainDir/run\.csh\</command\> >> $SubXml
echo \<stdout URL=\"file:$MainDir/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml
