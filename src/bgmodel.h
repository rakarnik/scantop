#ifndef _bgmodel
#define _bgmodel

#include "standard.h"
#include "seqset.h"
#include "motif.h"

class BGModel {
	float gc_genome;
	vector<float> model;
	void train_background();                                     // Train 0th order background model

public:
	BGModel(const float bg_gc);
	float gcgenome() const { return gc_genome; }                  // Return overall GC content
	double score_site(const Seqset& seq, const Motif& motif, const int c, const int p, const bool s) const;
	void print_model(ofstream& out);
};

#endif
