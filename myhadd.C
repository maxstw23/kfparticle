#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int myhadd(TString OutputFileName,TString Filelist,int NFile,TString options="")
{
  TString CMDhadd=0;CMDhadd.Append("hadd "+options+" "+OutputFileName+" ");

  ifstream ifile;ifile.open(Filelist,ios::in);if(!ifile){cerr<<"ifile.open failed"<<endl;return 0;}

  for(int i=0;i!=NFile;i++)
    {
      TString temp;
      ifile>>temp;
      CMDhadd.Append(temp+" ");
    }
  cout<<CMDhadd.Data()<<endl;
  system(CMDhadd.Data());
  return 0;
}
