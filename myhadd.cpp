//*****************************************************************************
//
//  Author:
//		Yifei Xu, woshixuyifei@gmail.com
//
//	Description:
//
//		This macro will combine root files in path :
//			$WORKPATH/rootfile_raw/*.root
//		Output files will be in path:
//			$WORKPATH/rootfile_temp/*.root
//
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

int myhadd(TString WorkPath,Int_t TESTMODE=0,TString DestinationDir="Hist",const Int_t Nfilesperjob=10,TString options="")
{
  //bool	TESTMODE	= true;	// 
  //*******************************************************************
  //const int Nfilesperjob	= 10;
  //*******************************************************************
  TString CMD = 0;
  TString ListFileName=0;
  ListFileName.Append("myhadd.list");
  TString ListFileName_random = WorkPath + "/myhadd_random.list";
  CMD = Form("find %s -maxdepth 1 -size +12k -name \"output_*.root\" > %s",WorkPath.Data(),ListFileName.Data());gSystem->Exec(CMD.Data());
  //CMD = Form("find %s -maxdepth 1 -size +12k -name \"DF5A9DA98912C3AB738DAEE21AE85084*\" > %s",WorkPath.Data(),ListFileName.Data());gSystem->Exec(CMD.Data());
  //CMD = Form("find %s -maxdepth 1 -size +12k -name \"F9777B1D60A8A9D6B103EF9B5FFFA6E6*\" > %s",WorkPath.Data(),ListFileName.Data());gSystem->Exec(CMD.Data());
  CMD = Form("sort -R %s > %s",ListFileName.Data(),ListFileName_random.Data());gSystem->Exec(CMD.Data());
  cout<<"================================================================================"<<endl;
  cout<<"             Getting File List   : "<<ListFileName<<endl;
  cout<<"================================================================================"<<endl;
  //*******************************************************************
  TString InputFile  = ListFileName_random;// File list
  TString outFileDir = 0;
  TString errFileDir = 0;
  TString infFileDir = 0;
  TString RootFileDir= 0;

  outFileDir = Form("%s/haddinfo/out",	WorkPath.Data());
  errFileDir = Form("%s/haddinfo/err",	WorkPath.Data());
  infFileDir = Form("%s/haddinfo/location",	WorkPath.Data());
  RootFileDir= Form("%s",	DestinationDir.Data());
  //*******************************************************************
  CMD = Form("mkdir -p %s",infFileDir.Data());gSystem->Exec(CMD.Data());
  CMD = Form("mkdir -p %s",outFileDir.Data());gSystem->Exec(CMD.Data());
  CMD = Form("mkdir -p %s",errFileDir.Data());gSystem->Exec(CMD.Data());
  CMD = Form("mkdir -p %s",RootFileDir.Data());gSystem->Exec(CMD.Data());

  ifstream ifile;ifile.open(InputFile,ios::in);if(!ifile){cerr<<"ifile.open failed!"<<endl;return 0;}
	
  TString ofilename=0;	ofilename.Append("myhadd.xml");
  CMD = Form("touch %s",ofilename.Data());gSystem->Exec(CMD.Data());

  TString CMDcommand = Form("root4star -l -b -q myhadd.C\\\(\\\"output_merged_$JOBID.root\\\",\\\"$FILELIST\\\",$INPUTFILECOUNT,\\\"%s\\\"\\\)",options.Data());

  ofstream ofile;ofile.open(ofilename.Data(),ios::out);if(!ofile){cerr<<"ofile.open failed!"<<endl;return 0;}
  cout<<"================================================================================"<<endl;
  cout<<"             Creating XML File   : "<<ofilename<<endl;
  cout<<"================================================================================"<<endl;

  ofile<<"<?xml version=\"1.0\" encoding=\"utf-8\" ?>"<<endl;
  ofile<<"<job  name=\"hadd\" maxFilesPerProcess=\""<<Nfilesperjob<<"\" filesPerHour=\"100\" fileListSyntax=\"paths\">"<<endl;
  ofile<<"<stdout URL=\"file:"<<outFileDir.Data()<<"/$JOBID.out\"/>"<<endl;
  ofile<<"<stderr URL=\"file:"<<errFileDir.Data()<<"/$JOBID.err\"/>"<<endl;
  ofile<<"<Generator>"<<endl;
  ofile<<"<Location>"<<infFileDir.Data()<<"</Location>"<<endl;
  ofile<<"</Generator>"<<endl;
  ofile<<"<output fromScratch=\"*.root\"    toURL=\""<<RootFileDir.Data()<<"/\" />"<<endl;
  ofile<<"<command>"<<endl;
  ofile<<"starver pro"<<endl;
  ofile<<"ln -s "<<WorkPath.Data()<<"/myhadd.C ."<<endl;
  ofile<<CMDcommand<<endl;
  ofile<<"</command>"<<endl;

  /*for(int j=0;!(ifile.eof());j++)
    {
      TString temp;	ifile>>temp;
      if(ifile.eof())break;
      ofile<<"<input URL=\"file:"<<temp.Data()<<"\"/>"<<endl;
    }*/
  //ofile<<"<input URL=\"file:"<<InputFile.Data()<<"\"/>"<<endl;
ofile<<"<input URL=\"filelist:"<<WorkPath.Data()<<"/&listOfFiles;\"/>"<<endl;
  ofile<<"</job>"<<endl;
  ofile.close();
  //CMD = Form("rm -f %s",ListFileName.Data());gSystem->Exec(CMD.Data());

  if(!TESTMODE) // submit
    {
      cout<<"================================================================================"<<endl;
      cout<<"             Submitting XML File : "<<ofilename<<endl;
      cout<<"================================================================================"<<endl;
      CMD = Form("star-submit-template -template %s -entities listOfFiles=myhadd_random.list",ofilename.Data());gSystem->Exec(CMD.Data());
    }
  return 0;
}
