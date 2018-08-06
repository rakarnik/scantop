#include "scantop.h"

int main(int argc, const char* argv[]) {
	if(argc != 6) {
		cerr << "Usage:\n";
		cerr << "scantop <ace_file> <motif_number> <sequence_file> <background_gc> <site_count>\n";
		exit(0);
	}

	string ace_file(argv[1]);
	string motnum(argv[2]);
	string seq_file(argv[3]);
	float bg_gc = atof(argv[4]);
	unsigned int count = atoi(argv[5]);
	vector<char> nucmap(4);

	nucmap[0] = 'A';
	nucmap[1] = 'C';
	nucmap[2] = 'G';
	nucmap[3] = 'T';

	// Read background sequences and create background model
	vector<string> bgseqs, bgnameset;
	BGModel bgm(bg_gc);

	// Read sequences to scan
	vector<string> seqs, nameset;
	get_fasta_fast(seq_file.c_str(), seqs, nameset, 0);
	Seqset scan_sseq(seqs);
	// cerr << scan_sseq.num_seqs() << " sequences read\n";

	// Read motif from ACE file
	vector<string> hits;
	string line;
	ifstream acein(ace_file.c_str());
	bool found = false;
	while(acein.good()) {
		getline(acein, line);
		istringstream iss(line);
		vector<string> fields;
		copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(fields));
		if(found) {
			if(line[0] != '*') {
				hits.push_back(fields[0]);
			} else {
				break;
			}
		} else {
			if(fields.size() > 0 
				&& (fields[0].compare("Motif") == 0)
				&& (fields[1].compare(motnum) == 0)) {
				found = true;
			}
		}
	}
	if(! found) {
		cerr << "Motif number " << motnum << " not found in ACE file '" << ace_file << "'\n";
		exit(1);
	}
	Seqset ace_sseq(hits);
	Motif mot(ace_sseq, 12, bgm.gcgenome());
	mot.clear_cols();
	for(unsigned int i = 0; i < line.length(); i++) {
		if(line[i] == '*') mot.add_col(i);
	}
	for(unsigned int i = 0; i < hits.size(); i++) {
		mot.add_site(i, 0, 1);
	}
	acein.close();

	// Scan sequences
	vector<list <ScoredSite> > topsites(scan_sseq.num_seqs());
	CompareScoredSite comp;
	double* sm = new double[mot.ncols() * 4];
	double* am = new double[mot.ncols()];
	mot.calc_score_matrix(sm, am);
	int width = mot.get_width();
	double scw, scc, msc, bsc;
	for(int c = 0; c < scan_sseq.num_seqs(); c++) {
		for(int p = 0; p < scan_sseq.len_seq(c) - width + 1; p++) {
			msc = mot.score_site(scan_sseq, sm, am, c, p, 1);
			bsc = bgm.score_site(scan_sseq, mot, c, p, 1);
			scw = msc - bsc;
			//cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 1 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scw << '\n';
			if(topsites[c].size() < count || topsites[c].front().score < scw) {
				if(topsites[c].size() == count) {
					topsites[c].pop_front();
				}
				topsites[c].push_front(ScoredSite(scw, c, p, 1));
				topsites[c].sort(comp);
			}
			msc = mot.score_site(scan_sseq, sm, am, c, p, 0);
			bsc = bgm.score_site(scan_sseq, mot, c, p, 0);
			scc = msc - bsc;
			//cerr << c << '\t' << p << '\t' << seqs[c].substr(p, width) << '\t' << 0 << '\t' << setprecision(4) << msc << '\t' << bsc << '\t' << scc << '\n';
			if(topsites[c].size() < count || topsites[c].front().score < scc) {
				if(topsites[c].size() == count) {
					topsites[c].pop_front();
				}
				topsites[c].push_front(ScoredSite(scc, c, p, 1));
				topsites[c].sort(comp);
			}
		}
	}
	delete [] sm;
	delete [] am;

	// Print out the sites
	list<ScoredSite>::iterator it;
	Site s;
	int p;
	for(int c = 0; c < scan_sseq.num_seqs(); c++) {
		for(it = topsites[c].begin(); it != topsites[c].end(); ++it) {
			cout << nameset[c];
			s = it->site;
			p = s.posit();
			cout << '\t' << p;
			cout << '\t' << p + width -1;
			cout << '\t' << s.strand();
			cout << '\t' << it->score;
			cout << '\t';
			if(s.strand() == 1) {
				for(int i = p; i < p + width; i++) {
					cout << nucmap[scan_sseq.seq()[c][i]];
				}
			} else {
				for(int i = p + width - 1; i >= p; i--) {
					cout << nucmap[3 - scan_sseq.seq()[c][i]];
				}
			}
			cout << '\n';
		}
	}
}
