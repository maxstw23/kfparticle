#include "StLambdaDecayPair.h"


ClassImp(StLambdaDecayPair);


StLambdaDecayPair::StLambdaDecayPair() {
	_p4Lambda.SetPxPyPzE(0,0,0,0);
	_p4Proton.SetPxPyPzE(0,0,0,0);
	//_p4Pion.SetPxPyPzE(0,0,0,0);
	_idxProton = -999;
	_idxPion = -999;
	_dmass = 999;
	_isSelf = true;
}


StLambdaDecayPair::StLambdaDecayPair(TLorentzVector p4Lambda, TLorentzVector p4Proton, int idxProton, int idxPion, bool isSelf) {
	_p4Lambda = p4Lambda;
	_p4Proton = p4Proton;
	_idxProton = idxProton;
	_idxPion = idxPion;
	_isSelf = isSelf;
	_dmass = fabs(p4Lambda.M() - 1.115683);
}


StLambdaDecayPair::StLambdaDecayPair(TLorentzVector p4Lambda, TLorentzVector p4Proton, int idxProton, int idxPion, bool isSelf, double dmass) {
	_p4Lambda = p4Lambda;
	_p4Proton = p4Proton;
	_idxProton = idxProton;
	_idxPion = idxPion;
	_isSelf = isSelf;
	_dmass = dmass;
}


StLambdaDecayPair::~StLambdaDecayPair() {
}
