#include <iostream>
#include <cstdlib>
#include <vector>

#include "igzstream.hpp"
#include "strtk/strtk.hpp"
#include "globals.hpp"
#include "Fasta.hpp"

using namespace std;

string prog_name;


namespace global {
  int flank_len = 20;
}


int
get_count_from_frag_list(const string & s)
{
  using namespace strtk;

  int res = 0;
  std_string::token_list_type token_list[3];
  split(";", s, back_inserter(token_list[0]));
  for_each(token_list[0].begin(), token_list[0].end(), [&] (const std_string::token_list_type::value_type & t) {
      token_list[1].clear();
      clog << "t=" << string(t.first, t.second) << "\n";
      split(single_delimiter_predicate<std::string::value_type>(':'),
	    t, back_inserter(token_list[1]));
      clog << "token_list[1][0]=" << string(token_list[1].front().first, token_list[1].front().second) << "\n";
      if (token_list[1].size() != 2) {
	cerr << "error parsing frag list: " << s << "\n";
	exit(EXIT_FAILURE);
      }
      token_list[2].clear();
      split(single_delimiter_predicate<std::string::value_type>(','),
	    token_list[1].back(), back_inserter(token_list[2]));
      res += token_list[2].size();
    });
  return res;
}

void
process_locus(const string & lib_line, const string & ref_evidence_line,
	      const string & alt_evidence_line)
{
  //strtk::ignore_token ignore;

  string locus_name;
  string s[4];
  string ref_chr;
  string alt_chr;
  long long ref_tsd[2][2] = {{0, 0}, {0, 0}};
  long long alt_tsd[2][2] = {{0, 0}, {0, 0}};
  bool is_insertion;

  // parse lib line
  strtk::std_string::token_list_type token_list;
  strtk::split("\t", lib_line, back_inserter(token_list));
  if (token_list.size() < 15) {
    cerr << "could not parse lib line: " << lib_line << "\n";
    exit(EXIT_FAILURE);
  }
  auto itr = token_list.begin();
  locus_name = string(itr->first, itr->second);
  ++itr;
  ref_chr = string(itr->first, itr->second);
  ++itr;
  ++itr;
  ++itr;
  ref_tsd[0][0] = atoll(itr->first);
  ++itr;
  ref_tsd[0][1] = atoll(itr->first);
  ++itr;
  s[0] = string(itr->first, itr->second);
  ++itr;
  s[1] = string(itr->first, itr->second);
  ++itr;
  alt_chr = string(itr->first, itr->second);
  ++itr;
  ++itr;
  ++itr;
  alt_tsd[0][0] = atoll(itr->first);
  ++itr;
  alt_tsd[0][1] = atoll(itr->first);
  ++itr;
  s[2] = string(itr->first, itr->second);
  ++itr;
  s[3] = string(itr->first, itr->second);

  if (s[0] == ".") {
    is_insertion = true;
    alt_tsd[1][0] = atoll(s[2].c_str());
    alt_tsd[1][1] = atoll(s[3].c_str());
  } else {
    is_insertion = false;
    ref_tsd[1][0] = atoll(s[0].c_str());
    ref_tsd[1][1] = atoll(s[1].c_str());
  }

  int count[7];
  // parse ref_evidence & alt_evidence lines
  if (is_insertion) {
    if (not strtk::parse(ref_evidence_line, "\t", count[2], s[0])
	or not strtk::parse(alt_evidence_line, "\t", count[0], count[1],
			    s[1], s[2], s[3])) {
      cerr << "could not parse ref/alt evidence lines:\n" << ref_evidence_line << "\n" << alt_evidence_line << "\n";
      exit(EXIT_FAILURE);
    }
    count[3] = get_count_from_frag_list(s[2]);
    count[4] = get_count_from_frag_list(s[3]);
    count[6] = get_count_from_frag_list(s[1]);
    count[5] = get_count_from_frag_list(s[0]);
  } else {
    if (not strtk::parse(ref_evidence_line, "\t", count[0], count[1], s[0], s[1], s[2])
	or not strtk::parse(alt_evidence_line, "\t", count[2], s[3])) {
      cerr << "could not parse ref/alt evidence lines:\n" << ref_evidence_line << "\n" << alt_evidence_line << "\n";
      exit(EXIT_FAILURE);
    }
    count[3] = get_count_from_frag_list(s[1]);
    count[4] = get_count_from_frag_list(s[2]);
    count[6] = get_count_from_frag_list(s[0]);
    count[5] = get_count_from_frag_list(s[3]);
  }

  clog << locus_name << "\t" << count[0] << "\t" << count[1] << "\t" << count[2] << "\t"
       << count[3] << "\t" << count[4] << "\t" << count[5] << "\t" << count[6] << "\n";

  // check presence of ins and null alleles
  bool have_ins_allele = false;
  bool have_null_allele = false;
  double e_null_cnt;
  double e_ins_cnt[2];
  if (is_insertion)
    {
      // null allele is ref
      e_null_cnt =
	get_expected_complete_span(global::refDict, global::rg_set,
				   ref_chr,
				   ref_tsd[0][0] + 1 - global::flank_len,
				   ref_tsd[0][1] + global::flank_len);

      // ins allele is alt
      e_ins_cnt[0] =
	get_expected_complete_span(global::refDict, global::rg_set,
				   alt_chr,
				   alt_tsd[0][0] + 1 - global::flank_len,
				   alt_tsd[0][1] + global::flank_len);
      e_ins_cnt[1] =
	get_expected_complete_span(global::refDict, global::rg_set,
				   alt_chr,
				   alt_tsd[1][0] + 1 - global::flank_len,
				   alt_tsd[1][1] + global::flank_len);
    }
  else // deletion
    {
      // ins allele is ref
      e_ins_cnt[0] =
	get_expected_complete_span(global::refDict, global::rg_set,
				   ref_chr,
				   ref_tsd[0][0] + 1 - global::flank_len,
				   ref_tsd[0][1] + global::flank_len);
      e_ins_cnt[1] =
	get_expected_complete_span(global::refDict, global::rg_set,
				   ref_chr,
				   ref_tsd[1][0] + 1 - global::flank_len,
				   ref_tsd[1][1] + global::flank_len);

      // null allele is alt
      e_null_cnt =
	get_expected_complete_span(global::refDict, global::rg_set,
				   alt_chr,
				   alt_tsd[0][0] + 1 - global::flank_len,
				   alt_tsd[0][1] + global::flank_len);
    }

  if (global::verbosity > 0)
    clog << "null allele test: have [" << count[5]
	 << "] expect [" << e_null_cnt << "]\n";
  have_null_allele = (double(count[5]) > .5 * e_null_cnt);

  if (global::verbosity > 0)
    clog << "ins allele test: 0: have [" << count[3] << "] expect [" << e_ins_cnt[0]
	 << "] 1: have [" << count[4] << "] expect [" << e_ins_cnt[1] << "]\n";
  have_ins_allele = ((double(count[3]) > .5 * e_ins_cnt[0])
		     or (double(count[4]) > .5 * e_ins_cnt[1]));

}


void
usage(ostream & os)
{
  os << "use: " << prog_name << " -f <ref_fasta_file> -g <alt_fasta_file> -l <pairing_file> -L <lib_file> -r <ref_evidence_file> -a <alt_evidence_file>\n";
}


int
main(int argc, char* argv[])
{
  prog_name = argv[0];
  string ref_fasta_file;
  string alt_fasta_file;
  string pairing_file;
  string lib_file;
  string ref_evidence_file;
  string alt_evidence_file;

  char c;
  while ((c = getopt(argc, argv, "vf:g:l:L:r:a:h")) != -1) {
    switch (c) {
    case 'v':
      global::verbosity++;
      break;
    case 'f':
      ref_fasta_file = optarg;
      break;
    case 'g':
      alt_fasta_file = optarg;
      break;
    case 'l':
      pairing_file = optarg;
      break;
    case 'L':
      lib_file = optarg;
      break;
    case 'r':
      ref_evidence_file = optarg;
      break;
    case 'a':
      alt_evidence_file = optarg;
      break;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
    default:
      cerr << "unrecognized option: " << c << "\n";
      usage(cerr);
      exit(EXIT_FAILURE);
    }
  }
  if (optind != argc) {
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  if (ref_fasta_file == "") { cerr << "missing ref fasta file\n"; exit(EXIT_FAILURE); }
  if (alt_fasta_file == "") { cerr << "missing alt fasta file\n"; exit(EXIT_FAILURE); }
  if (pairing_file == "") { cerr << "missing pairing file\n"; exit(EXIT_FAILURE); }
  if (lib_file == "") { cerr << "missing lib file\n"; exit(EXIT_FAILURE); }
  if (ref_evidence_file == "") { cerr << "missing ref_evidence file\n"; exit(EXIT_FAILURE); }
  if (alt_evidence_file == "") { cerr << "missing alt_evidence file\n"; exit(EXIT_FAILURE); }

  // load pairing file
  {
    igzstream pairing_is(pairing_file);
    global::rg_set.load(pairing_is);
    // check we have fragment rates
    for_each(global::rg_set.rg_list.begin(), global::rg_set.rg_list.end(),
	     [&] (const ReadGroup & rg) {
	       if (rg.get_pairing()->frag_rate.size() == 0) {
		 cerr << "missing fragment rate for read group "
		      << strtk::join(",", rg.get_names()) << "\n";
		 exit(EXIT_FAILURE);
	       }
	     });
  }

  // load ref fasta file
  {
    igzstream fasta_is(ref_fasta_file);
    readFasta(fasta_is, global::refDict);
  }
  // load alt fasta file
  {
    igzstream fasta_is(alt_fasta_file);
    readFasta(fasta_is, global::refDict);
  }

  igzstream lib_is(lib_file);
  igzstream ref_evidence_is(ref_evidence_file);
  igzstream alt_evidence_is(alt_evidence_file);

  int n_lines = 0;
  while (true) {
    string lib_line;
    string ref_evidence_line;
    string alt_evidence_line;
    bool got_lib_line = getline(lib_is, lib_line);
    bool got_ref_evidence_line = getline(ref_evidence_is, ref_evidence_line);
    bool got_alt_evidence_line = getline(alt_evidence_is, alt_evidence_line);

    if (got_lib_line != got_ref_evidence_line or got_lib_line != got_alt_evidence_line) {
      cerr << "error reading line " << n_lines+1 << " from lib/ref_evidence/alt_evidence files\n";
      exit(EXIT_FAILURE);
    }
    if (not got_lib_line) break;

    // got one line from each file
    process_locus(lib_line, ref_evidence_line, alt_evidence_line);
  }

  return EXIT_SUCCESS;
}
