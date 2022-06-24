#ifndef STLAMBDADECAYPAIR_H
#define STLAMBDADECAYPAIR_H

#include "TLorentzVector.h"


class StLambdaDecayPair
{
private:
	TLorentzVector _p4Lambda;
	TLorentzVector _p4Proton;
	//TLorentzVector p4Pion;
	double _dmass; // fabs(reconstructed_mass - Lambda_mass)
	int _idxProton;
	int _idxPion;
	bool _isSelf; // true: SelfLambda; false: AntiLambda

public:
	StLambdaDecayPair();
	StLambdaDecayPair(TLorentzVector p4Lambda, TLorentzVector p4Proton, int idxProton, int idxPion, bool isSelf);
	StLambdaDecayPair(TLorentzVector p4Lambda, TLorentzVector p4Proton, int idxProton, int idxPion, bool isSelf, double dmass);
	~StLambdaDecayPair();

	TLorentzVector get_p4Lambda() const;
	TLorentzVector get_p4Proton() const;
	TLorentzVector get_p4Pion() const;
	int get_idxProton() const;
	int get_idxPion() const;
	double get_dmass() const;
	bool get_isSelf() const;

	bool operator<(const StLambdaDecayPair &other) const; // for sort

	ClassDef(StLambdaDecayPair, 0);
};

inline TLorentzVector StLambdaDecayPair::get_p4Lambda() const { return _p4Lambda; }
inline TLorentzVector StLambdaDecayPair::get_p4Proton() const { return _p4Proton; }
inline TLorentzVector StLambdaDecayPair::get_p4Pion() const { return _p4Lambda-_p4Proton; }
inline int StLambdaDecayPair::get_idxProton() const { return _idxProton; }
inline int StLambdaDecayPair::get_idxPion() const { return _idxPion; }
inline double StLambdaDecayPair::get_dmass() const { return _dmass; }
inline bool StLambdaDecayPair::get_isSelf() const { return _isSelf; }

inline bool StLambdaDecayPair::operator<(const StLambdaDecayPair &other) const { return _dmass<other.get_dmass(); }


#endif
