#!/bin/csh

set MainDir=`pwd`

set XmlDir=./xml
if (! -e $XmlDir) exit
cd $XmlDir
set SubXml=sub.xml
if (-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=40
set nFileTotal=all

# print xml file
# 19130060 19131037 19135016 19137041 19139063 19140030 19141030 19144012 19144033 19145034 19147021 19147048 19155057 19158020 #Rungroups
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
echo \<command\>$MainDir/run\.csh\</command\> >> $SubXml
echo \<stdout URL=\"file:$MainDir/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
echo \<input URL=\"filelist:/star/data01/pwg/xiatong/git/kfparticle/datalist/run18List27.list\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19147048-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml
