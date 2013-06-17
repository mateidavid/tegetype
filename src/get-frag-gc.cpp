#include <cstdlib>
#include <iostream>
#include <sstream>
#include <random>

#include "globals.hpp"
#include "Pairing.hpp"
#include "Fasta.hpp"
#include "api/BamReader.h"
#include "igzstream.hpp"
#include "strtk/strtk.hpp"

using namespace std;
using namespace BamTools;


namespace global {
  int progress = 1000;
  int n_gc_bins = 100;
  int gc_window_step_len = 5;
  int max_ns = 10;
  long long seed = -1;
  double fraction = .001;
  vector<SQDict::iterator> bam_to_fa_dict;

  uniform_real_distribution<double> U(0.0, 1.0);
  default_random_engine R;

  vector<vector<int>> frag_gc(global::rg_set.rg_list.size(),
			      vector<int>(global::n_gc_bins, 0));
}


int getTag(const BamAlignment & m) {
  int res = 0;
  if (m.IsPaired()) res += 0x1;
  if (m.IsProperPair()) res += 0x2;
  if (not m.IsMapped()) res += 0x4;
  if (not m.IsMateMapped()) res += 0x8;
  if (m.IsReverseStrand()) res += 0x10;
  if (m.IsMateReverseStrand()) res += 0x20;
  if (m.IsFirstMate()) res += 0x40;
  if (m.IsSecondMate()) res += 0x80;
  if (not m.IsPrimaryAlignment()) res += 0x100;
  if (m.IsFailedQC()) res += 0x200;
  if (m.IsDuplicate()) res += 0x400;
  return res;
}

string getCigar(const BamAlignment & m) {
  ostringstream res;
  for_each(m.CigarData.begin(), m.CigarData.end(), [&] (const CigarOp & op) {
      res << op.Length << op.Type;
    });
  return res.str();
}

string getMapping(const BamAlignment & m, const RefVector & bam_seq) {
  ostringstream res;
  res << m.Name << "\t"
      << getTag(m) << "\t"
      << bam_seq[m.RefID].RefName << "\t"
      << m.Position + 1 << "\t"
      << m.MapQuality << "\t"
      << getCigar(m) << "\t"
      << (m.MateRefID == m.RefID? string("=") : 
	  (m.MateRefID < 0? string("*") : bam_seq[m.MateRefID].RefName)) << "\t"
      << m.MatePosition + 1 << "\t"
      << m.InsertSize;
  return res.str();
}


void
process_file(BamReader & bam_file, const RefVector & bam_seq)
{
  int n_samples = 0;
  BamAlignment m;
  while (bam_file.GetNextAlignmentCore(m)) {
    if (global::U(global::R) > global::fraction) continue;
    if (not m.IsMapped()) continue;
    if (global::bam_to_fa_dict[m.RefID] == global::refDict.end()) continue;
    m.BuildCharData();

    string rg_string;
    m.GetTag(string("RG"), rg_string);
    if (rg_string.size() == 0) {
      cerr << "missing RG in SAM line: " << getMapping(m, bam_seq) << "\n";
      continue;
    }
    ReadGroup * rg_p = global::rg_set.find_by_name(rg_string);
    if (rg_p == NULL) {
      cerr << "unknown RG: " << getMapping(m, bam_seq) << "\n";
      continue;
    }
    size_t rg_idx = rg_p - &global::rg_set.rg_list[0];
    if ((m.IsPaired()
	 and not rg_p->get_pairing()->is_mp_downstream(m.IsFirstMate()? 0 : 1,
						       0,
						       m.IsReverseStrand()? 1 : 0))
	or (not m.IsPaired() and m.IsReverseStrand()))
      continue;
					     
    if (global::verbosity > 1) clog << "sampled: " << getMapping(m, bam_seq);

    // compute gc content at mapped location
    int crt_gc = 0;
    int crt_ns = 0;
    for_each(global::bam_to_fa_dict[m.RefID]->second.seq[0].begin() + m.Position,
	     global::bam_to_fa_dict[m.RefID]->second.seq[0].begin() + m.Position + rg_p->get_pairing()->mean,
	     [&] (char c) {
	       if (c == 'N' or c == 'n') crt_ns++;
	       else if (c == 'G' or c == 'C' or c == 'g' or c == 'c') crt_gc++;
	     });
    if (crt_ns > global::max_ns) {
      if (global::verbosity > 1) clog << ": too many Ns [" << crt_ns << "]\n";
      continue;
    }
    int bin_idx = int((double(crt_gc) / (rg_p->get_pairing()->mean + 1)) * global::n_gc_bins);
    if (global::verbosity > 1) clog << ": bin_idx:[" << bin_idx
				    << "] crt_gc:[" << crt_gc << "]\n";
    if (bin_idx < 10) {
      clog << getMapping(m, bam_seq) << "\n";
      clog << global::bam_to_fa_dict[m.RefID]->second.seq[0].substr(m.Position, rg_p->get_pairing()->mean) << "\n"
	   << "bin_idx:[" << bin_idx << "] crt_gc:[" << crt_gc << "]\n";
    }
    global::frag_gc[rg_idx][bin_idx]++;
    n_samples++;
    if (n_samples % global::progress == 0) clog << n_samples << "\n";
  }
}

void
print_output()
{
  for (size_t rg_idx = 0; rg_idx < global::rg_set.rg_list.size(); ++rg_idx) {
    auto rg_p = global::rg_set.find_by_idx(rg_idx);
    for (int j = 0; j < global::n_gc_bins; ++j) {
      cout << strtk::join(",", rg_p->get_names()) << "\t"
	   << rg_p->get_pairing()->mean << "\t"
	   << j << "\t"
	   << int(ceil((double(j) / global::n_gc_bins) * (rg_p->get_pairing()->mean + 1))) << "\t"
	   << int(double(global::frag_gc[rg_idx][j]) / global::fraction) << "\n";
    }
  }
}


void
get_bam_to_fa_dict(RefVector & bam_seq)
{
  global::bam_to_fa_dict.clear();

  for (size_t i = 0; i < bam_seq.size(); ++i) {
    auto it_first = global::refDict.end();
    auto it_last = global::refDict.end();
    for (auto it = global::refDict.begin(); it != global::refDict.end(); ++it) {
      if (it->second.len != (long long)bam_seq[i].RefLength) continue;
      if (it_first == global::refDict.end()) it_first = it;
      it_last = it;
    }
    if (it_first != it_last) {
      clog << "BAM SQ [" << bam_seq[i].RefName << "] matches more than one fasta seq: ["
	   << it_first->second.name << "," << it_last->second.name << "]; ignoring\n";
      global::bam_to_fa_dict.push_back(global::refDict.end());
    } else if (it_first == global::refDict.end()) {
      clog << "BAM SQ [" << bam_seq[i].RefName << "] not in fasta file; ignoring\n";
      global::bam_to_fa_dict.push_back(it_first);
    } else {
      clog << "BAM SQ [" << bam_seq[i].RefName << "] = fasta SQ [" << it_first->second.name
	   << "]\n";
      global::bam_to_fa_dict.push_back(it_first);
    }
  }
}


int
main(int argc, char * argv[])
{
  string prog_name(argv[0]);
  string fasta_filename;
  string pairing_filename;
  vector<string> mapping_filename;

  char c;
    while ((c = getopt(argc, argv, "vN:s:p:q:f:l:m:")) != -1) {
    switch (c) {
    case 'v':
      ++global::verbosity;
      break;
    case 'N':
      global::num_threads = atoi(optarg);
      break;
    case 's':
      global::seed = atoi(optarg);
      break;
    case 'p':
      global::progress = atoi(optarg);
      break;
    case 'q':
      global::fraction = atof(optarg);
      break;
    case 'f':
      fasta_filename = optarg;
      break;
    case 'l':
      pairing_filename = optarg;
      break;
    case 'm':
      mapping_filename.push_back(string(optarg));
      break;
    default:
      cerr << "unrecognized option: " << c << "\n";
      exit(EXIT_FAILURE);
    }
  }
  if (optind != argc) {
    cerr << "use: " << prog_name << " [options]\n";
    exit(EXIT_FAILURE);
  }

  if (fasta_filename.size() == 0) {
    cerr << "fasta file not given\n";
    exit(EXIT_FAILURE);
  }
  if (pairing_filename.size() == 0) {
    cerr << "pairing file not given\n";
    exit(EXIT_FAILURE);
  }
  if (mapping_filename.size() == 0) {
    cerr << "mapping file(s) not given\n";
    exit(EXIT_FAILURE);
  }

  // init random engine
  if (global::seed < 0) global::seed = time(NULL);
  global::R = default_random_engine(global::seed);

  if (global::verbosity > 0) {
    clog << "using seed [" << global::seed << "]\n";
    clog << "using fraction [" << global::fraction << "]\n";
  }

  // load pairing file
  {
    igzstream pairing_is(pairing_filename);
    if (!pairing_is) {
      cerr << "error opening pairing file: " << pairing_filename << "\n";
      exit(EXIT_FAILURE);
    }
    global::rg_set.load(pairing_is);
    for_each(global::rg_set.rg_list.begin(), global::rg_set.rg_list.end(), [&] (ReadGroup & rg) {
	rg.pairing.mean -= (rg.pairing.mean % global::gc_window_step_len);
      });
  }

  // load reference
  {
    igzstream fasta_is(fasta_filename);
    if (!fasta_is) {
      cerr << "error opening reference fasta file: " << fasta_filename << "\n";
      exit(EXIT_FAILURE);
    }
    readFasta(fasta_is, global::refDict);
  }

  // init global counters
  global::frag_gc = vector<vector<int>>(global::rg_set.rg_list.size(),
					vector<int>(global::n_gc_bins, 0));

  // process mapping files sequentially
  for_each(mapping_filename.begin(), mapping_filename.end(), [&] (const string & fn) {
      BamReader bam_file;
      if (not bam_file.Open(fn)) {
	cerr << "error opening mapping file: " << bam_file.GetErrorString() << "\n";
	exit(EXIT_FAILURE);
      } else if (global::verbosity > 0) {
	clog << "processing file: " << fn << "\n";
      }
      //SamHeader bam_header = bam_file.GetHeader();
      RefVector bam_seq = bam_file.GetReferenceData();

      // for each sq in mappings file, find corresponding sq in fasta file
      get_bam_to_fa_dict(bam_seq);

      process_file(bam_file, bam_seq);
    });

  // output
  print_output();

  return EXIT_SUCCESS;
}
