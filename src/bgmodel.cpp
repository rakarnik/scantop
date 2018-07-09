#include "bgmodel.h"

BGModel::BGModel(const float bg_gc) : gc_genome(bg_gc), model(4) {
	train_background();
}

double BGModel::score_site(const Seqset& seqset, const Motif& motif, const int c, const int p, const bool s) const {
	const vector<vector <int> >& ss_seq = seqset.seq();
	double L = 0.0;
	int width = motif.get_width();
	vector<int>::const_iterator col_iter = motif.first_column();
	vector<int>::const_iterator last_col = motif.last_column();
	if(s) {
		for(; col_iter != last_col; ++col_iter) {
			L += log2(model[ss_seq[c][p + *col_iter]]);
		}
	} else {
		for(; col_iter != last_col; ++col_iter) {
			L += log2(model[ss_seq[c][p + width - 1 - *col_iter]]);
		}
	}
	return L;
}

void BGModel::train_background() {
	for(int i = 0; i < 4; i++) {
		model[i] = 0;
	}
	
	// Just use genome-wide GC content
	model[0] = (1 - gc_genome)/2;
	model[1] = gc_genome/2;
	model[2] = model[1];
	model[3] = model[0];
}

void BGModel::print_model(ofstream& out) {
	map<int,char> nt;
	nt[0] = 'A';
	nt[1] = 'C';
	nt[2] = 'G';
	nt[3] = 'T';

	out << "Model order:\n";
	for(int i = 0; i < 4; i++) {
		out << nt[i % 4] << '\t' << model[i] << '\n';
	}
	out << "\n\n";
}
