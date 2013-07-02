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
  bool is_male = true;
  bool check_tail_if_head_not_solid = false;
}


int
get_chr_count(const string & chr)
{
  if (chr == "chrX")
    return global::is_male? 1 : 2;
  else if (chr == "chrY")
    return global::is_male? 1 : 0;
  else return 2;
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
      if (global::verbosity >= 3) clog << "t=" << string(t.first, t.second) << "\n";
      split(single_delimiter_predicate<std::string::value_type>(':'),
	    t, back_inserter(token_list[1]));
      if (global::verbosity >= 3) clog << "token_list[1][0]=" << string(token_list[1].front().first, token_list[1].front().second) << "\n";
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
  int ref_strand;
  int tsd_to_check;
  bool is_insertion;
  bool solid_bp[2];

  // parse lib line
  strtk::std_string::token_list_type token_list;
  strtk::split("\t", lib_line, back_inserter(token_list));
  if (token_list.size() < 17) {
    cerr << "could not parse lib line: " << lib_line << "\n";
    exit(EXIT_FAILURE);
  }
  auto itr = token_list.begin(); // name
  locus_name = string(itr->first, itr->second);
  ++itr; // ref_chr
  ref_chr = string(itr->first, itr->second);
  ++itr; // ref_reg_start
  ++itr; // ref_reg_end
  ++itr; // ref_tsd0_start
  ref_tsd[0][0] = atoll(itr->first);
  ++itr; // ref_tsd0_end
  ref_tsd[0][1] = atoll(itr->first);
  ++itr; // ref_tsd1_start
  s[0] = string(itr->first, itr->second);
  ++itr; // ref_tsd1_end
  s[1] = string(itr->first, itr->second);
  ++itr; // strand
  ref_strand = (*(itr->first) == '+'? 0 : 1);
  ++itr; // bp_coverage
  solid_bp[0] = (*(itr->first) == '1');
  solid_bp[1] = (*(itr->first + 1) == '1');
  ++itr; // alt_chr
  alt_chr = string(itr->first, itr->second);
  ++itr; // alt_reg_start
  ++itr; // alt_reg_end
  ++itr; // alt_tsd0_start
  alt_tsd[0][0] = atoll(itr->first);
  ++itr; // alt_tsd0_end
  alt_tsd[0][1] = atoll(itr->first);
  ++itr; // alt_tsd1_start
  s[2] = string(itr->first, itr->second);
  ++itr; // alt_tsd1_end
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
  tsd_to_check = ref_strand;
  if (global::check_tail_if_head_not_solid
      and solid_bp[tsd_to_check] == false and solid_bp[1 - tsd_to_check] == true)
    tsd_to_check = 1 - tsd_to_check;

  int count[7];
  // parse ref_evidence & alt_evidence lines
  //
  // count structure:
  //    --------6--------
  //    ---3---   ---4---
  // ins  -0- ===== -1-
  // ====               ====
  // null      -2-
  //    --------5--------
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

  // check allele presence
  bool ins_allele_present = false;
  bool null_allele_present = false;
  bool ins_allele_absent = false;
  bool null_allele_absent = false;
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

  int chr_count = get_chr_count(ref_chr);
  null_allele_present = (chr_count >= 1
			 and (//count[5] >= 5 or
			      double(count[5]) > max(2.0, .5 * e_null_cnt)
			      )
			 );

  ins_allele_present = (chr_count >= 1
			and ((//count[3 + tsd_to_check] >= 5 or
			      double(count[3 + tsd_to_check])
			      > max(2.0, .5 * e_ins_cnt[0 + tsd_to_check])
			      )
			     )
			);

  if (chr_count == 1 and null_allele_present and ins_allele_present) {
    // cannot both be present, don't make a call
    null_allele_present = false;
    ins_allele_present = false;
  }

  bool more_than_2 = false;
  // check for evidence of more than 2 alleles
  if (null_allele_present and ins_allele_present
      and (double(count[5]) > 1.5 * e_null_cnt
	   or double(count[3 + tsd_to_check]) > 1.5 * e_ins_cnt[0 + tsd_to_check])) {
    more_than_2 = true;
  }

  // ins allele absent?
  if (chr_count == 2 and null_allele_present and not ins_allele_present
      and count[0] == 0
      and count[1] == 0
      and count[3] == 0
      and count[4] == 0
      and double(count[5]) > 1.5 * e_null_cnt)
    ins_allele_absent = true;

  // null allele absent?
  if (chr_count == 2 and ins_allele_present and not null_allele_present
      and count[2] == 0
      and count[5] == 0
      and double(count[3 + tsd_to_check]) > 1.5 * e_ins_cnt[0 + tsd_to_check])
    null_allele_absent = true;

  cout << locus_name << "\t" << (is_insertion? "I" : "D") << "\t";

  char null_allele_name = is_insertion? 'R' : 'A';
  char ins_allele_name = is_insertion? 'A' : 'R';
  string getype_abs = (chr_count == 2? "--": (chr_count == 1? "-X" : "XX"));
  string getype_rel = (chr_count == 2? "--": (chr_count == 1? "-X" : "XX"));

  if (null_allele_present) {
    getype_abs[0] = 'N';
    getype_rel[0] = null_allele_name;
    if (ins_allele_present) {
      getype_abs[1] = 'I';
      getype_rel[1] = ins_allele_name;
    } else if (ins_allele_absent) {
      getype_abs[1] = 'N';
      getype_rel[1] = null_allele_name;
    }
  } else if (ins_allele_present) {
    getype_abs[0] = 'I';
    getype_rel[0] = ins_allele_name;
    if (null_allele_absent) {
      getype_abs[1] = 'I';
      getype_rel[1] = ins_allele_name;
    }
  }
  if (getype_rel == "AR") getype_rel = "RA";

  if (more_than_2) {
    getype_abs += "+";
    getype_rel += "+";
  }

  cout << getype_abs << "\t" << getype_rel << "\n";

  if (global::verbosity >= 2) {
    clog.unsetf(ios_base::floatfield);
    clog << locus_name << "\t" << (is_insertion? "I" : "D")
	 << "\t" << getype_abs << "\t" << getype_rel << "\t"
	 << count[0] << "\t" << count[1] << "\t" << count[2] << "\t"
	 << count[3] << "/" << e_ins_cnt[0] << "\t"
	 << count[4] << "/" << e_ins_cnt[1] << "\t"
	 << count[5] << "/" << e_null_cnt << "\t"
	 << count[6] << "\n";
  }
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
  {
    char * p = getenv("GENDER");
    if (p != NULL and *p == 'f') global::is_male = false;
  }

  string ref_fasta_file;
  string alt_fasta_file;
  string pairing_file;
  string lib_file;
  string ref_evidence_file;
  string alt_evidence_file;

  char c;
  while ((c = getopt(argc, argv, "vf:g:l:L:r:a:x:th")) != -1) {
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
    case 't':
      global::check_tail_if_head_not_solid = true;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
    case 'x':
      if (optarg[0] == 'm') {
	global::is_male = true;
      } else if (optarg[0] == 'f') {
	global::is_male = false;
      } else {
	cerr << "unrecognized gender: " << optarg << "\n";
	exit(EXIT_FAILURE);
      }
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

  if (global::verbosity >= 1) {
    clog << "gender: " << (global::is_male? "m" : "f") << "\n";
  }

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
