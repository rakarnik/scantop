#include "motif.h"

Motif::Motif() :
seqset(Seqset()) {
}

Motif::Motif(const Seqset& s, int nc, double g) :
seqset(s),
init_nc(nc),
gc(g),
num_seqs(s.num_seqs()),
max_width(3 * init_nc),
columns(init_nc),
num_seqs_with_sites(0),
has_sites(num_seqs)
{
	pseudo[0] = pseudo[3] = (1.0 - gc)/2.0;
	pseudo[1] = pseudo[2] = gc/2.0;
	vector<int>::iterator col_iter = columns.begin();
	for(; col_iter != columns.end(); ++col_iter) {
		*col_iter = distance(columns.begin(), col_iter);
	}
	width = ((int) columns.back()) + 1;
	seq_cutoff = 0.00001;
	expr_cutoff = 0.70;
	iter = "0";
	dejavu = 0;
	spec = 0.0;
}

Motif::Motif(const Motif& m) :
seqset(m.seqset),
init_nc(m.init_nc),
gc(m.gc),
width(m.width),
num_seqs(m.num_seqs),
max_width(m.max_width),
columns(m.columns),
num_seqs_with_sites(m.num_seqs_with_sites),
has_sites(m.has_sites)
{
	for(int i = 0; i < 4; i++)
		pseudo[i] = m.pseudo[i];
  mapsc = m.mapsc;
	spec = m.spec;
  iter = m.iter;
	dejavu = m.dejavu;
	seq_cutoff = m.seq_cutoff;
	expr_cutoff = m.expr_cutoff;
  *this = m;
}

Motif& Motif::operator= (const Motif& m) {
  if(this != &m){
    //assume that the same Seqset is referred to, so ignore some things
		init_nc = m.init_nc;
		gc = m.gc;
		for(int i = 0; i < 4; i++)
			pseudo[i] = m.pseudo[i];
		sitelist.assign(m.sitelist.begin(), m.sitelist.end());
		num_seqs_with_sites = m.num_seqs_with_sites;
		has_sites.assign(m.has_sites.begin(), m.has_sites.end());
    columns.assign(m.columns.begin(), m.columns.end());
		width = ((int) columns.back()) + 1;
		seq_cutoff = m.seq_cutoff;
		expr_cutoff = m.expr_cutoff;
		mapsc = m.mapsc;
		spec = m.spec;
		iter = m.iter;
		dejavu = m.dejavu;
  }
  return *this;
}

void Motif::clear_sites() {
	vector<Site>().swap(sitelist);
	vector<int>(init_nc).swap(columns);
	vector<int>::iterator col_iter = columns.begin();
	for(; col_iter != columns.end(); ++col_iter) {
		*col_iter = distance(columns.begin(), col_iter);
	}
	width = ((int) columns.back()) + 1;
  num_seqs_with_sites = 0;
	has_sites.assign(num_seqs, false);
	mapsc = 0.0;
	spec = 0.0;
}

void Motif::remove_all_sites() {
	sitelist.clear();
	num_seqs_with_sites = 0;
	for(int i = 0; i < num_seqs; i++) {
		has_sites[i] = false;
	}
}

bool Motif::is_open_site(const int c, const int p){
	vector<Site>::iterator site_iter = sitelist.begin();
	for(; site_iter != sitelist.end(); ++site_iter) {
    if(site_iter->chrom() == c) {
      int pp = site_iter->posit();
      if(pp > p - width && pp < p + width) return false;
    }
  }
  return true;
}

void Motif::add_site(const int c, const int p, const bool s){
	assert(p >= 0 && p < seqset.len_seq(c));
	Site st(c, p, s);
	sitelist.push_back(st);
	if(has_sites[c] == 0) num_seqs_with_sites++;
	has_sites[c]++;
}

void Motif::calc_freq_matrix(int *fm) const {
	const vector<vector <int> >& seq = seqset.seq();
	for(int i = 0; i < 4 * ncols(); i++){
		fm[i] = 0;
  }
	
	int g, j, pos, matpos;
	bool s;
  vector<Site>::const_iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
		g = site_iter->chrom();
    j = site_iter->posit();
    s = site_iter->strand();
    if(s) {                              // forward strand
      matpos = 0;
			pos = 0;
			vector<int>::const_iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = j + *col_iter;
				if(pos >= 0 && pos < seqset.len_seq(g)) {
					fm[matpos + seq[g][pos]]++;
				}
				matpos += 4;
      }
    } else {                             // reverse strand
      matpos = 0;
      pos = 0;
			vector<int>::const_iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = j + width - 1 - *col_iter;
				if(pos >= 0 && pos < seqset.len_seq(g)) {
					fm[matpos + 3 - seq[g][pos]]++;
				}
				matpos += 4;
      }
    }
  }
}

bool Motif::column_freq(const int col, int *ret){
  const vector<vector <int> >& seq = seqset.seq();
	for(int i = 0; i < 4; i++) ret[i] = 0;
	int c, p;
	bool s;
	vector<Site>::iterator site_iter = sitelist.begin();
  for(; site_iter != sitelist.end(); ++site_iter) {
    c = site_iter->chrom();
    p = site_iter->posit();
    s = site_iter->strand();
    if(s) {
      if( (p + col > seqset.len_seq(c) - 1) || (p + col < 0) ) return false;
      ret[seq[c][p + col]]++;
    } else {
      if((p + width - 1 - col > seqset.len_seq(c) - 1) || (p + width - 1 - col < 0)) return false;
      ret[3 - seq[c][p + width - 1 - col]]++;
    }
  }
  return true;
}

double Motif::score_site(double* score_matrix, double* avg_matrix, const int c, const int p, const bool s) const {
	return score_site(seqset, score_matrix, avg_matrix, c, p, s);
}

double Motif::score_site(const Seqset& ss, double* score_matrix, double* avg_matrix, const int c, const int p, const bool s) const {
	const vector<vector<int> >& ss_seq = ss.seq();
	double L = 0.0;
	int width = get_width();
	int matpos;
	vector<int>::const_iterator col_iter = columns.begin();
	if(s) {
		matpos = 0;
		for(; col_iter != columns.end() ; ++col_iter) {
			assert(p + *col_iter >= 0);
			assert(p + *col_iter < ss.len_seq(c));
			if(ss_seq[c][p + *col_iter] != -1) {
				L += score_matrix[matpos * 4 + ss_seq[c][p + *col_iter]];
			} else {
				L += avg_matrix[matpos];
			}
			matpos++;
		}
	} else {
		matpos = 0;
		for(; col_iter != columns.end(); ++col_iter) {
			assert(p + width - 1 - *col_iter >= 0);
			assert(p + width - 1 - *col_iter < ss.len_seq(c));
			if(ss_seq[c][p + width - 1 - *col_iter] != -1) {
				L += score_matrix[matpos * 4 + 3 - ss_seq[c][p + width - 1 - *col_iter]];
			} else {
				L += avg_matrix[matpos];
			}
			matpos++;
		}
	}
	return L;
}

int Motif::remove_col(const int c) {
  int ret = 0, ns = 0;         // return number of removed column in new numbering
	bool found = false;
	if(c == 0) {                   // if the column to be removed is the first column, we shift all columns so first is 0
    columns.erase(columns.begin());
		if(columns.size() > 0) {
			ns = columns[0];
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
				*col_iter -= ns;
		}
		ret = -ns;
		shift_sites(ns);
		found = true;
  } else {
		vector<int>::iterator col_iter;
    for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
      if(*col_iter == c) {
				columns.erase(col_iter);
				found = true;
				break;
			}
    }
	}
	if (! found) {
		cerr << "remove_column called for column " << c << " but it was not found!" << endl; 
		abort();
	}
	width = ((int) columns.back()) + 1;
  return ret;
}

void Motif::add_col(const int c) {
  if(c < 0){                   // column to the left of existing, shift existing to the right
		if(columns.size() > 0) {
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter){
				*col_iter -= c;
			}
			shift_sites(c);
		}
		columns.insert(columns.begin(), 0);
  } else {
		bool found = false;
		vector<int>::iterator col_iter;
		for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter){
			if(*col_iter > c) {
				columns.insert(col_iter, c);
				found = true;
				break;
			}
		}
		if(! found)
			columns.push_back(c);
  }
	width = ((int) columns.back()) + 1;
}

bool Motif::has_col(const int c) {
	return binary_search(columns.begin(), columns.end(), c);
}

void Motif::shift_sites(const int shift) {
	vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter)
		site_iter->shift(shift);
}

void Motif::flip_sites() {
  int i;
  for(i = 0; i < number(); i++) {
    sitelist[i].flip();
  }
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
		*col_iter = width - *col_iter - 1;
	}
	vector<int> temp;
	vector<int>::reverse_iterator rev_iter;
	for(rev_iter = columns.rbegin(); rev_iter != columns.rend(); ++rev_iter) {
		temp.push_back(*rev_iter);
	}
	columns.assign(temp.begin(), temp.end());
}

void Motif::orient() {
	int* freq_matrix = new int[4 * ncols()];
  double *freq = new double[4];
	for(int i = 0; i < 4; i++)
		freq[i] = 0.0;
	
  calc_freq_matrix(freq_matrix);
  for(int i = 0; i < 4 * ncols(); i += 4) {
    for(int j = 0; j < 4; j++)
      freq[j] += freq_matrix[i + j];
  }
	for(int i = 0; i < 4; i++)
		freq[i] /= number();
	
  double flip = 1.5 * freq[2] + 1.0 * freq[0] - 1.0 * freq[3] - 1.5 * freq[1];
  if(flip < 0.0) {
		cerr << "\t\t\t\tFlipping motif...\n";
		flip_sites();
	}
	delete [] freq_matrix;
  delete [] freq;
}

int Motif::positions_available() const {
  int ret = 0;
  for(int i = 0; i < num_seqs; i++){
    ret += seqset.len_seq(i) - width + 1;
  }
  return ret;
}

int Motif::positions_available(const vector<bool>& possible) const {
	int ret = 0;
	for(int i = 0; i < num_seqs; i++) {
		if(possible[i])
			ret += seqset.len_seq(i) - width + 1;
	}
	return ret;
}

void Motif::columns_open(int &l, int &r){
	int w = width;
	int c, p, len;
	bool s;
	vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    c = site_iter->chrom();
    p = site_iter->posit();
    s = site_iter->strand();
		len = seqset.len_seq(c);
		l = min(p, l);
		r = min(len - p - w, r);
  }
}

void Motif::freq_matrix_extended(double *fm) const{
  const vector<vector <int> >& seq = seqset.seq();
	int i, col, j;
  int fm_size = (width + 2 * ncols()) * 4;
  for(i = 0; i < fm_size; i++) fm[i] = 0.0;
  if(number() == 0) return;
  for(i = 0; i < number(); i++) {//i = site number
    int c = chrom(i);
    int p = posit(i);
    bool s = strand(i);
    for(j = 0, col = -ncols(); col < width + ncols(); col++, j += 4) {
      if(s) {
				if((p + col <= seqset.len_seq(c) - 1) && (p + col >= 0)) {
					assert(j + seq[c][p + col] < fm_size);
					fm[j + seq[c][p + col]] += 1.0;
				} else {
					for(int k = 0; k < 4; k++) {
						fm[j + k] += 0.25;
					}
				}
      } else {
				if((p + width - 1 - col <= seqset.len_seq(c) - 1) && (p + width - 1 - col >= 0)) {
					assert(j + 3 - seq[c][p + width - 1 - col] < fm_size);
					fm[j + 3 - seq[c][p + width - 1 - col]] += 1.0;
				} else {
					for(int k = 0; k < 4; k++) {
						fm[j + k] += 0.25;
					}
				}
			}
    }
  }
  for(i = 0; i < fm_size; i++) fm[i] /= (double) number();
}

void Motif::calc_score_matrix(double *sm, double *am) const {
	int* fm = new int[4 * ncols()];
	double tot = (double) number();
	for(int j = 0; j < 4; j++)
		tot += pseudo[j];
  calc_freq_matrix(fm);
  for(int i = 0; i < ncols(); i++) {
		for(int j = 0; j < 4; j++) {
			sm[i * 4 + j] = log2((fm[i * 4 + j] + pseudo[j])/tot);
		}
		am[i] = 0.0;
		am[i] += (1.0 - gc) * (sm[i * 4] + sm[i * 4 + 3]);
		am[i] += gc * (sm[i * 4 + 1] + sm[i * 4 + 2]);
	}
	delete [] fm;
}

string Motif::consensus() const {
	int numsites = number();
	if(numsites < 1) return "";
	char nt[] = {'A', 'C', 'G', 'T'};
	const vector<vector <int> >& seq = seqset.seq();
	vector<string> hits;
	hits.reserve(numsites);
	vector<Site>::const_iterator siteit = sitelist.begin();
	int c, p;
	bool s;
	for(; siteit != sitelist.end(); ++siteit) {
		c = siteit->chrom();
    p = siteit->posit();
    s = siteit->strand();
		if(s) {
			string hitseq = "";
			for(int k = p; k < p + width - 1; k++)
				hitseq.append(1, nt[seq[c][k]]);
			hits.push_back(hitseq);
		} else {
			string hitseq = "";
			for(int k = p + width - 1; k > p - 1; k--)
				hitseq.append(1, nt[3 - seq[c][k]]);
			hits.push_back(hitseq);
    }
  }
	
	string cons;
  int wide1 = hits[0].size();
  int num1 = hits.size();
  int num1a, num1c, num1g, num1t;
  for(int i = 0; i < wide1; i++){
    num1a = num1c = num1g = num1t = 0;
    for(int j = 0; j < num1; j++){
      if(hits[j][i] == 'a' || hits[j][i] == 'A') num1a += 1;
      if(hits[j][i] == 'c' || hits[j][i] == 'C') num1c += 1;
      if(hits[j][i] == 'g' || hits[j][i] == 'G') num1g += 1;
      if(hits[j][i] == 't' || hits[j][i] == 'T') num1t += 1;
    }
    if(num1a > num1*0.7) cons += 'A';
    else if(num1c > num1*0.7) cons += 'C';
    else if(num1g > num1*0.7) cons += 'G';
    else if(num1t > num1*0.7) cons += 'T';
    else if((num1a + num1t) > num1 * 0.85) cons +='W';
    else if((num1a + num1c) > num1 * 0.85) cons +='M';
    else if((num1a + num1g) > num1 * 0.85) cons +='R';
    else if((num1t + num1c) > num1 * 0.85) cons +='Y';
    else if((num1t + num1g) > num1 * 0.85) cons +='K';
    else if((num1c + num1g) > num1 * 0.85) cons +='S';
    else cons += '-';
  }
  return cons;
}

void Motif::write(ostream& motout) const {
	char nt[] = {'A', 'C', 'G', 'T'};
	const vector<vector <int> >& seq = seqset.seq();
	vector<Site>::const_iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    int c = site_iter->chrom();
    int p = site_iter->posit();
    bool s = site_iter->strand();
    for(int j = 0; j < width; j++){
      if(s) {
				if(p + j >= 0 && p + j < seqset.len_seq(c))
					motout << nt[seq[c][p + j]];
				else motout << ' ';
      }
      else {
				if(p + width - 1 - j >= 0 && p + width - 1 - j < seqset.len_seq(c))
					motout << nt[3 - seq[c][p + width - 1 - j]];
				else motout << ' ';
      }
    }
    motout << '\t' << c << '\t' << p << '\t' << s << '\n';
  }
	int col = 0, prev_col = 0;
	vector<int>::const_iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {	
		col = *col_iter;
		for(int k = 0; k < (col - prev_col - 1); k++) motout << ' ';
		motout << '*';
		prev_col = col;
	}
  motout << endl;
	
	motout << "MAP Score: " << mapsc << "\n";
	motout << "Specificity Score: " << spec << "\n";
	motout << "Sequence cutoff: " << seq_cutoff << "\n";
	motout << "Expression cutoff: " << expr_cutoff << "\n";
	motout << "Iteration found: " << iter << "\n";
	motout << "Dejavu: " << dejavu << endl << "\n";
}

void Motif::read(istream& motin) {
	char* match;
	char line[200];
	vector<int> c;
	vector<int> p;
	vector<bool> s;
	
	// Read sites
	// (don't add yet, as they will get screwed up by the column changes)
	while(motin.getline(line, 200)) {
		if(line[0] == '*') break;
		match = strtok(line, "\t");
		c.push_back(atoi(strtok(NULL, "\t")));
		p.push_back(atoi(strtok(NULL, "\t")));
		s.push_back(atoi(strtok(NULL, "\0")));
	}
	
	int motwidth = strlen(line);
	columns.clear();
	for(int i = 0; i < motwidth; i++) {
		if(line[i] == '*') add_col(i);
	}
	
	// Add sites
	sitelist.clear();
	int num_sites = c.size();
	for(int i = 0; i < num_sites; i++) {
		assert(p[i] >= 0);
		add_site(c[i], p[i], s[i]);
	}
	
	char* heading;
	// Read MAP score
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_map(atof(strtok(NULL, "\0")));
	
	// Read specificity
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_spec(atof(strtok(NULL, "\0")));
	
	// Read sequence cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_seq_cutoff(atof(strtok(NULL, "\0")));
	
	// Read expression cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_expr_cutoff(atof(strtok(NULL, "\0")));
	
	// Read iteration found
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_iter(strtok(NULL, "\0"));
	
	// Read dejavu
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_dejavu(atoi(strtok(NULL, "\0")));
}

void Motif::destroy() {
  //for when the = operator can fail
  max_width=0;
}

void Motif::print_columns(ostream& out) {
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
		out << " " << *col_iter;
}

bool Motif::check_sites() {
	int c, p;
	bool s;
	vector<Site>::iterator site_iter;
	for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
		c = site_iter->chrom();
		p = site_iter->posit();
		s = site_iter->strand();
		if(p < 0 || p + width - 1 >= seqset.len_seq(c)) {
			cerr << "\t\t\t\t\tc: " << c << " p: " << p << " s: " << s << " w: " << width << " len: " << seqset.len_seq(c) << endl;
			return false;
		}
	}
	return true;
}
