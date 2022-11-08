#!/bin/csh

# settings
set nRun=18
set mRun=Run${nRun}
set mEnergy=27.0
set ListDir=/star/data01/pwg/xiatong/git/kfparticle/datalist/ #TODO
set MainDir=/star/data01/pwg/xiatong/git/kfparticle/ #TODO
set TempDir=/home/tmp/maxwoo/ #TODO
# inputs
#set JOBINDEX=$1
##set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/$mEnergy.list.`printf "%.6d" ${JOBINDEX}`
#set FILELIST={$ListDir}/${mEnergy}GeV_${mRun}/test.list
#set JOBID=ScriptTestSandbox

set WorkDir=${TempDir}/$JOBID
mkdir -p $WorkDir
cd $WorkDir
cp -Lr $MainDir/setDEV2.csh .
cp -Lr $MainDir/lMuDst.C .
cp -Lr $MainDir/mix . 
source setDEV2.csh
cp -Lr $MainDir/readPicoDst.C .
cp -Lr $MainDir/.sl73_x8664_gcc485 .

set RootLog=root_${JOBINDEX}.log
if(-e $RootLog) rm $RootLog

root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$JOBINDEX,$nRun,$mEnergy,\"$ListDir\"\) >& $RootLog

set Iter=0
while( `grep -sc '(ret%10)<=kStFatal' $RootLog` )
	@ Iter++
	if( $Iter > 5 ) break
	rm $RootLog
	rm *.root
	root4star -b -q ./readPicoDst.C\(\"$FILELIST\",$JOBINDEX,$nRun,$mEnergy,\"$ListDir\"\) >& $RootLog
end

mv *.log  $MainDir/log/.
mv *.root $MainDir/output/.
