#!/bin/csh
if ( $#argv == 0 ) then
    source clean.csh
endif
set MainDir=`pwd`

# set XmlDir=./xml
# if (! -e $XmlDir) exit
# cd $XmlDir
set SubXml=sub.xml
if (-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=40
set nFileTotal=all

# print xml file
# 19130060 19131037 19135016 19137041 19139063 19140030 19141030 19144012 19144033 19145034 19147021 19147048 19155057 19158020 #Rungroups
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
# add this <shell>singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif</shell>
if ( $HOST =~ *"starsub"* ) then
    echo \<shell\>singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif\</shell\> >> $SubXml
endif
echo \<command\>$MainDir/run\.csh $MainDir\</command\> >> $SubXml
echo \<stdout URL=\"file:$MainDir/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
# echo \<input URL=\"filelist:/star/data01/pwg/xiatong/git/kfparticle/datalist/run18List27.list\" /\> >> $SubXml

### if no argument, use the default list
### if -r option, use unanalyzed list
if ( $#argv == 0) then
    echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19130060-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
else if ( "$argv[1]" == "-r" ) then
    echo \<input URL=\"filelist:$MainDir/unanalyzed_file_list.list\" /\> >> $SubXml
endif
# echo \<input URL=\"catalog:star\.bnl\.gov\?production=P19ib,filetype=daq_reco_PicoDst,trgsetupname~27GeV_production_2018,runnumber\[\]19147048-19268002,sanity=1,tpx=1,storage!=hpss,filename~st_physics\" nFiles=\"$nFileTotal\" /\> >> $SubXml
echo \</job\> >> $SubXml

if ( $HOST =~ *"starsub"* ) star-submit-beta $SubXml
if ( $HOST =~ *"rcas"* ) star-submit $SubXml
sh monitor.sh