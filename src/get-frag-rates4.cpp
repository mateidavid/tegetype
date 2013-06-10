#include <cstdlib>
#include <iostream>
#include <sstream>
#include <random>

#include "globals.hpp"
#include "Pairing.hpp"
#include "Fasta.hpp"
#include "api/BamReader.h"
#include "igzstream.hpp"

using namespace std;
using namespace BamTools;


namespace global {
  int progress;
  int n_GC_bins = 100;
  int step_len = 5;
  long long seed = -1;
  double fraction = .001;
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


void
process_file(BamReader & bam_file, SamHeader & bam_header, RefVector & bam_seq)
{
  uniform_real_distribution<double> U(0.0, 1.0);
  default_random_engine R(global::seed >= 0 ? global::seed : time(NULL));

  BamAlignment m;
  while (bam_file.GetNextAlignmentCore(m)) {
    if (U(R) > global::fraction) continue;
    m.BuildCharData();

    if (global::verbosity > 1) {
      clog << m.Name << "\t"
	   << getTag(m) << "\t"
	   << bam_seq[m.RefID].RefName << "\t"
	   << m.Position + 1 << "\t"
	   << m.MapQuality << "\t"
	   << getCigar(m) << "\t"
	   << (m.MateRefID == m.RefID? string("=") : 
	       (m.MateRefID < 0? string("*") : bam_seq[m.MateRefID].RefName)) << "\t"
	   << m.MatePosition + 1 << "\t"
	   << m.InsertSize << "\n";
    }

    string rg_string;
    m.GetTag(string("RG"), rg_string);

  }

}


void
get_bam_to_fa_dict(RefVector & bam_seq, vector<SQDict::iterator> & bam_to_fa_dict)
{
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
      bam_to_fa_dict.push_back(global::refDict.end());
    } else if (it_first == global::refDict.end()) {
      clog << "BAM SQ [" << bam_seq[i].RefName << "] not in fasta file; ignoring\n";
      bam_to_fa_dict.push_back(it_first);
    } else {
      clog << "BAM SQ [" << bam_seq[i].RefName << "] = fasta SQ [" << it_first->second.name
	   << "]\n";
      bam_to_fa_dict.push_back(it_first);
    }
  }
}


int
main(int argc, char * argv[])
{
  string prog_name(argv[0]);
  string fasta_filename;
  string pairing_filename;
  string mapping_filename;

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
      mapping_filename = optarg;
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
    cerr << "mapping file not given\n";
    exit(EXIT_FAILURE);
  }

  {
    igzstream pairing_is(pairing_filename);
    if (!pairing_is) {
      cerr << "error opening pairing file: " << pairing_filename << "\n";
      exit(EXIT_FAILURE);
    }
    load_pairing(pairing_is, global::rg_dict, global::num_rg_dict, global::rg_to_num_rg_dict);
  }

  {
    igzstream fasta_is(fasta_filename);
    if (!fasta_is) {
      cerr << "error opening reference fasta file: " << fasta_filename << "\n";
      exit(EXIT_FAILURE);
    }
    readFasta(fasta_is, global::refDict);
  }

  BamReader bam_file;
  if (not bam_file.Open(mapping_filename)) {
    cerr << "error opening mapping file: " << bam_file.GetErrorString() << "\n";
    exit(EXIT_FAILURE);
  }
  SamHeader bam_header = bam_file.GetHeader();
  RefVector bam_seq = bam_file.GetReferenceData();

  // for each sq in mappings file, find corresponding sq in fasta file
  vector<SQDict::iterator> bam_to_fa_dict;
  get_bam_to_fa_dict(bam_seq, bam_to_fa_dict);

  process_file(bam_file, bam_header, bam_seq);
  
  return EXIT_SUCCESS;
}
