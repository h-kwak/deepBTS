#include "arguments.h"
#include "bedgraph.h"
#include "genelist.h"
#include "fanntrain.h"
#include "datafile.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

// arguments
struct custom_arg
{
  bool do_md, do_tr, do_te;
  string output_file;
  string plus_file, minus_file, tre_file;
  string network_file, tmpnet_file, connet_file;
  string param_file;
  vector<string> train_files;
  vector<int> bin, range;
  int step, win, far;
  float sample_ratio;
  bool sample_near;
  arguments as;
  custom_arg(int argc, char *argv[]);
  vector<string> dummy;
};

// deepPRO functions
void dp_make_data(custom_arg &a);
void dp_train(custom_arg &a);
void dp_test(custom_arg &a);

// main function
int main(int argc, char **argv)
{
  custom_arg a(argc, argv);
  if(a.do_md) {
    dp_make_data(a);
  } else if(a.do_tr) {
    dp_train(a);
  } else {
    dp_test(a);
  }
  return 1;
}

// Functions

void logistic(vector<float> &in, vector<float> &out)
{
  float x = 0.05, y = 0.01;
  float max = *max_element(in.begin(), in.end());
  out.resize(in.size());
  if(max <= 0) {
    fill(out.begin(), out.end(), y);
    return;
  }
  float b = x * max;
  float a = 1 / b * log(1/y - 1);
  for(int i = 0; i < in.size(); ++i)
    out[i] = 1/(1 + exp(-a*(in[i] - b)));
  return;
}

inline void report_count(int i, int e, int t, const char *a, const char *b)
{
  if(i > 0 && i % e == 0) {
    cout << a << " " << i/1000 << "k / " << t/1000 << "k " << b << "\r";
    cout.flush();    
  }
}

// Multires readcount structure
struct mrread
{
  vector<vector<float> > pl, mn;
  vector<int> bin;
  void set_bin(vector<int> &b);
  bool get_from_bg(string chr, bgdata &p, bgdata &m);
  bool gen_eval_vec(int pos, vector<int> &range, vector<float> &output);
  bool check_cutoff(int pos, float co100b, float co1k);
  mrread();
};

mrread::mrread()
{
  bin.push_back(100);
  bin.push_back(50);
  bin.push_back(500);
};

void mrread::set_bin(vector<int> &b)
{
  bin.resize(1);
  bin.insert(bin.end(), b.begin(), b.end());
};

bool mrread::get_from_bg(string chr, bgdata &p, bgdata &m)
{
  pl.clear();
  mn.clear();
  for(int i=0; i<bin.size(); ++i) {
    if(!p.set(chr, bin[i])) return(false);
    if(!m.set(chr, bin[i])) return(false);
    pl.push_back(p.v);
    mn.push_back(m.v);
    for(int j=0;j<mn[i].size();++j) if(mn[i][j]<0) mn[i][j] = -mn[i][j];
  }
  return(true);
};

// Generate evaluation vectors
bool mrread::gen_eval_vec(int pos, vector<int> &range, vector<float> &output)
{
  if(range.size() != bin.size()-1) return(false);
  output.clear();
  for(int i = 0; i < range.size(); ++i) {
    if(range[i] == 0) continue;
    int ss = (pos - range[i])/bin[i+1], se = (pos + range[i])/bin[i+1];
    if(ss < 0 || se > pl[i+1].size() || se > mn[i+1].size()) return(false);
    vector<float> tp(pl[i+1].begin() + ss, pl[i+1].begin() + se);
    vector<float> tm(mn[i+1].begin() + ss, mn[i+1].begin() + se);
    vector<float> lp, lm;
    logistic(tp, lp);
    logistic(tm, lm);
    output.insert(output.end(), lp.begin(), lp.end());
    output.insert(output.end(), lm.begin(), lm.end());
  }
  return(true);
};

// Check if the region contains sufficient reads through heuristic criteria
bool mrread::check_cutoff(int pos, float co100b = 3, float co1k = 1)
{
  int ps = pos/100 - 5, pe = pos/100 + 5;
  if(ps < 0 || pe > pl[0].size() || pe > mn[0].size()) return(false);
  // Check if any of the 100 bp bins in 1 kb region has more than cut-off reads
  if(*max_element(pl[0].begin() + ps, pl[0].begin() + pe) >= co100b ||
     *max_element(mn[0].begin() + ps, mn[0].begin() + pe) >= co100b)
    return(true);
  // Check if sum of the reads in 1 kb regions in both strands more than cut off
  float psum = 0, msum = 0;
  for(int j = ps; j < pe; ++j) {
    psum += pl[0][j];
    msum += mn[0][j];
  }
  if(psum >= co1k && msum >= co1k)
    return(true);
  return(false);
};

custom_arg::custom_arg(int argc, char *argv[]) {
  do_md = do_tr = do_te = false;
  int b[2] = {50, 500};
  int r[2] = {500, 5000};
  bin.insert(bin.begin(), b, b+2);
  range.insert(range.begin(), r, r+2);
  step = 50;
  sample_ratio = 1.5;
  win = 1000, far = 10000;
  param_file = "ANN_parameters.txt";
  tmpnet_file = "temp.net";
  connet_file = "";
  dummy.resize(20);
  as.push("-o", "<Output file>", output_file);
  as.push("-p", "Plus strand bedgraph", plus_file, true);
  as.push("-m", "Minus strand bedgraph", minus_file, true);
  as.push("-n", "ANN network configuration", network_file, true);
  as.push("- - - - - - - - - - - -","- - - - - - - - - - - - - - - - - -", dummy[0], true);
  as.push("make", "Make training data from bedgraph files to the output file", do_md, true);
  as.push("-p ", "Plus strand bedgraph", dummy[1], true);
  as.push("-m ", "Minus strand bedgraph", dummy[2], true);
  as.push("-b", "TRE bed file", tre_file, true);
  as.push("--bin", "Bin sizes (default = 50, 500 bp)", bin, true);
  as.push("--rng", "Scanning ranges (default = 500, 5000 bp)", range, true);
  as.push("--step", "Step size (default = 50 bp, must be a divisor of bins and ranges)", step, true);
  as.push("--sr", "Sample ratio (negative:positive, default = 1.5)", sample_ratio, true);
  as.push("--p", "Parameter file output (default = \"ANN_parameters.txt\")", param_file, true);
  as.push("- - - -","- - - - - - - - - - - - - - - - - - - - - - - - - -", dummy[3], true);
  as.push("train", "Train and save ANN configuration to the output file", do_tr, true);
  as.push("-t", "Input training set files", train_files, true);
  as.push("--p ", "Parameter file input (default = \"ANN_parameters.txt\")", dummy[4], true);
  as.push("--con", "Previous network configuration to continiue training (default = none)", connet_file, true);
  if(!as.get(argc,argv)) exit(0);
};

////////////////
// Train feature

void dp_make_data(custom_arg &a)
{
  // Load PRO-seq bedgraph file
  bgdata pl, mn;
  cout << "LOADING PLUS STRAND BEDGRAPH" << endl;
  pl.load((char *)a.plus_file.c_str());
  cout << "LOADING MINUS STRAND BEDGRAPH" << endl;
  mn.load((char *)a.minus_file.c_str());
	
  // Load PRO-cap derived TREs in bed3 format
  genelist tl;
  cout << "LOADING REFERENCE TRE BED" << endl;
  tl.load((char *)a.tre_file.c_str(), "bed3");

  // Multiple resolution readcounts 
  mrread mr;
  mr.set_bin(a.bin);

  // Evaluation vectors
  vector<vector<float> > tre, ntr;
  
  // Loop through all chromosome
  for(int ch = 0; ch < pl.chr.size(); ++ch) {
    string chr = pl.chr[ch];
    cout << "READING " << chr << "                \r";
    cout.flush();
    if(!mr.get_from_bg(chr, pl, mn)) continue;
    // Initialize TRE region marks
    int arSize = mr.pl[0].size()*mr.bin[0]/a.step;
    if(mr.mn[0].size()*mr.bin[0]/a.step > arSize) arSize = mr.mn[0].size()*mr.bin[0]/a.step;
    vector<bool> tre_mark(arSize, false);
    // Loop through TREs
    for(int i = 0; i < tl.n; ++i) {
      // if the TRE is not in the same chromosome, pass
      if(tl.chrName[i] != chr) continue;
      // TRE window start and end position
      int mp = (tl.txStart[i] + tl.txEnd[i]) / 2;
      int swin = (int)((mp - a.win)/a.step) * a.step;
      int ewin = swin + 2 * a.win;
      // Mark TRE position
      for(int j = swin/a.step; j < ewin/a.step; ++j)
        if(j < arSize && j >=0) tre_mark[j] = true;
      // Check for readcount cutoff
      if(mr.check_cutoff(mp, 4, 2)) {
        // Retrieve evaluation vector and save as tre output
        vector<float> output;
        if(mr.gen_eval_vec(mp, a.range, output)) tre.push_back(output);
      }
    }
	  /* If sampling near, loop through TREs again (DEPRECATED)
	  if(a.sample_near) for(int i = 0; i < tl.n; ++i) {
	    // if the TRE is not in the same chromosome, pass
	    if(tl.chrName[i] != chr) continue;
	    // TRE position
	    int mp = (tl.txStart[i] + tl.txEnd[i]) / 2;
	    // Scan upstream
	    for(int pos = mp - 2000; pos > 0; pos -= 1000) {
	      if(!tre_mark[pos/a.step]) if(mr.check_cutoff(pos, 4, 2)) {
	        vector<float> output;
	        if(mr.gen_eval_vec(pos, a.range, output)) {
	          ntr.push_back(output);
	          break;
	        }
	      }
	    }
	    // Scan downstream
	    for(int pos = mp + 2000; pos < arSize * a.step; pos += 1000) {
	      if(!tre_mark[pos/a.step]) if(mr.check_cutoff(pos, 4, 2)) {
	        vector<float> output;
	        if(mr.gen_eval_vec(pos, a.range, output)) {
	          ntr.push_back(output);
	          break;
	        }
	      }
	    }
	  }
    // Scan through positions for negative training sets
	  if(!a.sample_near) */
    for(int pos = 0; pos < arSize*a.step; pos += a.step) {
      // If near a TRE mark, pass
      if(tre_mark[pos/a.step]) continue;
      // Check for readcount cutoff
      if(!mr.check_cutoff(pos, 4, 2)) continue;
      // Retrieve evaluation vector
      vector<float> output;
      if(mr.gen_eval_vec(pos, a.range, output)) ntr.push_back(output);
    }
  }
  // Finished reading chromosomes
  cout << "FINISHED READING CHROMOSOMES" << endl;
  // Save the evaluation readcounts in a FANN training data format
  ofstream tmp_out(a.output_file);
  int nline;
  float sample_freq = a.sample_ratio * tre.size() / ntr.size();
  if(sample_freq > 1) sample_freq = 1;
  // Sample non-TRE data by the sampling ratio
  srand(time(NULL));
  vector<bool> ntr_sample(ntr.size(), false);
  int data_count = tre.size();
  for(int i = 0; i < ntr_sample.size(); ++i)
    if(rand() < RAND_MAX * sample_freq) {
      ntr_sample[i] = true;
      ++data_count;
    }
  // Write training data header
  int n_input = tre[0].size();
  tmp_out << data_count << "\t" << n_input << "\t" << 2 << endl;
  for(int i = 0; i < tre.size(); ++i) {
    report_count(i, 1000, tre.size()+ntr.size(), "WRITING OUTPUT", "FEATURE VECTORS");
    tmp_out << tre[i][0];
    for(int j = 1; j < n_input; ++j)
      if(j < tre[i].size()) tmp_out << "\t" << tre[i][j];
      else tmp_out<< "\t0";
    tmp_out << endl << "1\t0" <<endl;
  }
  for(int i = 0; i < ntr.size(); ++i) {
    report_count(i + tre.size(), 1000, tre.size()+ntr.size(), "WRITING OUTPUT", "FEATURE VECTORS");
    if(!ntr_sample[i]) continue;
    tmp_out << ntr[i][0];
    for(int j = 1; j < n_input; ++j)
      if(j < ntr[i].size()) tmp_out << "\t" << ntr[i][j];
      else tmp_out<< "\t0";   
    tmp_out << endl << "0\t1" <<endl;
  }
  tmp_out.close();
  // Writing ANN parameter file
  ofstream par_out(a.param_file);
  par_out << "learning_rate=\t0.1" <<endl;
  par_out << "num_layers=\t3" <<endl;
  par_out << "num_input=\t" << tre[0].size() << endl;
  par_out << "num_hidden=\t" << int(sqrt(tre[0].size()*2)) << endl;
  par_out << "num_output=\t2" <<endl;
  par_out << "desired_error=\t0.001" <<endl;
  par_out << "max_iterations=\t2000" << endl;
  par_out << "iterations_between_reports=\t5" << endl;
  par_out << "num_bin=\t" << a.bin.size() << endl;
  par_out << "bins=";
  for(int i = 0; i < a.bin.size(); ++i) par_out << "\t" << a.bin[i];
  par_out << endl;
  par_out << "ranges=";
  for(int i = 0; i < a.range.size(); ++i) par_out << "\t" << a.range[i];
  par_out << endl;
  par_out << "step=\t" << a.step << endl;
  par_out.close();
}

// Read parameter file
void read_param_file(custom_arg &a) {
  ifstream par(a.param_file);
  string buf;
  int num_layers, num_bin;
  getline(par, buf);
  par >> buf >> num_layers;
  for(int i = 0; i < (num_layers + 3) * 2; ++i) par >> buf;
  par >> buf >> num_bin;
  a.bin.resize(num_bin);
  par >> buf;
  for(int i = 0; i < num_bin; ++i) par >> a.bin[i];
  par >> buf;
  a.range.resize(num_bin);
  for(int i = 0; i < num_bin; ++i) par >> a.range[i];
  par >> buf >> a.step;
}

void dp_train(custom_arg &a) {
// Train and save the network using TREs and non-TREs
  FANN::neural_net net;
  FANN::training_data data, data2;
  cout<<"READING TRAINING DATA: "<<a.train_files[0]<<"\r";
  cout.flush();
  data.read_train_from_file(a.train_files[0]);
  
  for(int i = 1; i < a.train_files.size(); ++i) {
    cout<<"READING TRAINING DATA: "<<a.train_files[i]<<"\r";
    cout.flush();
    data2.read_train_from_file(a.train_files[i]);
    data.merge_train_data(data2);
  }
  cout<<endl;
  read_param_file(a);
  stringstream header;
  header << "num_bin= " << a.bin.size() << " bins=";
  for(int i = 0; i < a.bin.size(); ++i) header << " " << a.bin[i];
  header << " ranges=";
  for(int i = 0; i < a.range.size(); ++i) header << " " << a.range[i];
  header << " step= " << a.step;
  simpler_train(net, data, a.output_file, a.param_file, a.connet_file, header.str());
}

void parse_header(string header, custom_arg &a){
  stringstream ss(header);
  string buf;
  int num_bin;
  ss >> buf >> num_bin;
  ss >> buf;
  a.bin.resize(num_bin);
  for(int i = 0; i < num_bin; ++i) ss >> a.bin[i];
  a.range.resize(num_bin);
  ss >> buf;
  for(int i = 0; i < num_bin; ++i) ss >> a.range[i];
  ss >> buf >> a.step;
}

///////////////
// Test feature
void dp_test(custom_arg &a) {
  // Load PRO-seq bedgraph file
  bgdata pl, mn;
  cout << "LOADING PLUS STRAND BEDGRAPH" << endl;
  pl.load((char *)a.plus_file.c_str());
  cout << "LOADING MINUS STRAND BEDGRAPH" << endl;
  mn.load((char *)a.minus_file.c_str());

  // Load network file
  cout << "LOADING NETWORK" << endl;
  FANN::neural_net net;
  string custom_header;
  read_net(net, a.network_file, custom_header);
  if(custom_header != "") parse_header(custom_header, a);

  // Multiple resolution readcounts 
  cout << "READING RESOLUTIONS : " << a.bin[0];
  for(int i = 1; i < a.bin.size(); ++i) cout << ", " << a.bin[i];
  cout << endl;
  cout << "RANGES : " << a.range[0];
  for(int i = 1; i < a.range.size(); ++i) cout << ", " << a.range[i];
  cout << endl;
  cout << "STEP SIZE : " << a.step << endl;

  mrread mr;
  mr.set_bin(a.bin);
  
  // Evaluation vectors
  vector<vector<float> > tr;
  vector<string> tr_chr;
  vector<int> tr_pos;
  
  // Loop through chromosomes
  cout << "READING CHROMOSOMES" << endl;
  for(int ch = 0; ch < pl.chr.size(); ++ch) {
    string chr = pl.chr[ch];
    if(!mr.get_from_bg(chr, pl, mn)) continue;
    cout << "READING " << chr << "\r";
    cout.flush();
    
    // Scan through positions for testing sets
    int chrSize = mr.pl[0].size() * mr.bin[0];
    for(int pos = 0; pos < chrSize; pos += a.step)
    {
      // Check for readcount cutoff
      if(!mr.check_cutoff(pos, 4, 2)) continue;
      // Retrieve evaluation vector
      vector<float> output;
      if(mr.gen_eval_vec(pos, a.range, output)) {
        tr.push_back(output);
        tr_chr.push_back(chr);
        tr_pos.push_back(pos);
      }
    }
  }
  // Finished reading chromosomes
  cout << endl;

  // Sanitize input data
  int n_input = tr[0].size();
  for(int i = 0; i < tr.size(); ++i)
    if(tr[i].size() < n_input) tr[i].resize(n_input, 0.01);

 // Run FANN
  cout << "RUNNING NETWORK" << endl;
  vector<vector<float> > res_score;
  simpler_test(net, tr, res_score);

  ofstream res(a.output_file);
  
  // Loop through all the lines and print as bedgraph fromat
  for(int i = 0; i < tr_chr.size(); ++i) {
    report_count(i, 1000, tr_chr.size(), "WRITING", "SITES"); 
    res << tr_chr[i] << "\t" << tr_pos[i] << "\t" <<
    tr_pos[i]+1 << "\t" << res_score[i][0] << endl;
  }
  res.close();
}
