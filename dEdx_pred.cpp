#include <iostream>
#include "./StRoot/StBichsel/StdEdxPull.h"

void dEdx_pred()
{	
	gSystem->Load("StBichsel");
	const float mKaon = 0.13957039;
	int p_time_1000 = 0;
	std::vector<float> pvec;
	std::vector<float> dEdxvec;
	while (p_time_1000 < 2000)
	{
		float betagamma = p_time_1000/1000. / mKaon;
		
		pvec.push_back(p_time_1000/1000.);
		dEdxvec.push_back(1e6*StdEdxPull::EvalPred(betagamma, 2, 1));
		p_time_1000 += 10;	
	}

	for (int i = 0; i < pvec.size(); i++) std::cout << pvec[i] << ", ";
	std::cout << std::endl;
	for (int i = 0; i < dEdxvec.size(); i++) std::cout << dEdxvec[i] << ", ";
	std::cout << std::endl;

}
