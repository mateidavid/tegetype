#include <iostream>
#include <cstdlib>
#include <vector>

#include "igzstream.hpp"
#include "strtk/strtk.hpp"
#include "globals.hpp"
#include "Clone.hpp"
#include "CloneGen.hpp"
#include "SamMapping.hpp"
#include "SamMappingSetGen.hpp"
#include "Pairing.hpp"
#include "common.hpp"
#include "Fasta.hpp"

using namespace std;


class TSD {
public:
  int start;
  int end;
  int count;
  TSD(int start_, int end_) : start(start_), end(end_), count(0) {}
};


string prog_name;
vector<TSD> tsd;
vector<vector<vector<int>>> cluster;
string (*cnp)(const string &);
void (*fnp)(const string &, Clone &, int &);
int flank_len = 30;
int min_non_repeat_bp = 20;
int min_read_len = 20;
int min_read_len_left = 20; // after removing NM from either side
int expected_insert_size = 320; // used to place limit on allowable fragment sizes
int min_mqv = 0;
int max_nm = 10;
bool is_alt = false;
int reg_start;
int reg_end;
int total_bp;
int total_bp_left;
int total_bp_right;
int total_bp_mid;


void
process_mapping_set(const string & clone_name, vector<SamMapping> & v_sm)
{
  if ((global::rg_set.rg_list.size() == 0 and v_sm.size() != 1)
      or (global::rg_set.rg_list.size() > 0 and v_sm.size() != 2)) {
    cerr << "incorrect number of mappings for clone [" << clone_name << "]\n";
    exit(EXIT_FAILURE);
  }

  // convert to Mapping structures
  vector<Mapping> v_m(2);
  for (size_t j = 0; j < v_sm.size(); ++j) {
    if (v_sm[j].mapped) v_m[j] = convert_SamMapping_to_Mapping(v_sm[j]);
  }

  // trim mapping positions by NM value on either side;
  //if (tsd.size() == 2)
  {
    long long min_pos = -1;
    long long max_pos = -1;
    // count non-repeat bp
    int non_repeat_bp = 0;

    for (size_t j = 0; j < v_m.size(); ++j) {
      if (!v_sm[j].mapped) continue;
      int edit_dist = 0;
      for (size_t i = 0; i < v_sm[j].rest.size(); ++i) {
	if (v_sm[j].rest[i].key == "NM") {
	  edit_dist = atoi(v_sm[j].rest[i].value.c_str());
	  break;
	}
      }
      // discard fragment if NM too large
      if (edit_dist > max_nm) {
	LOG(1) << "[" << v_sm[0].name << "]: discarding; large NM\n";
	return;
      }
      if (edit_dist > 0) {
	v_m[j].dbPos[0] += edit_dist;
	v_m[j].dbPos[1] -= edit_dist;
      }
      if (v_m[j].dbPos[1] - v_m[j].dbPos[0] + 1 < min_read_len_left) {
	LOG(1) << "[" << v_sm[0].name << "]: discarding; small read len after NM trim\n";
	return;
      }
      if (min_pos < 0 or v_m[j].dbPos[0] < min_pos) min_pos = v_m[j].dbPos[0];
      if (max_pos < 0 or v_m[j].dbPos[1] > max_pos) max_pos = v_m[j].dbPos[1];

      if (is_alt)
	for_each(v_m[j].db->seq[0].begin() + v_m[j].dbPos[0],
		 v_m[j].db->seq[0].begin() + v_m[j].dbPos[1],
		 [&] (char c) {
		   non_repeat_bp += (c >= 'A' and c <= 'Z'? 1 : 0);
		 });
    }

    /*
    // if fragment does not capture 20bp outside of the inner repeat, ignore
    if (tsd.size() == 2
	and (min_pos < 0
	     or max_pos < 0
	     or (min_pos > tsd[0].end - min_non_repeat
		 and max_pos < tsd[1].start + min_non_repeat))) {
    */
    if (is_alt and non_repeat_bp < min_non_repeat_bp) {
      LOG(1) << "[" << v_sm[0].name << "]: discarding: not enough non-repeat bp\n";
      return;
    }
  }

  // count bp mapped left/right/between TSDs
  for (size_t j = 0; j < v_sm.size(); ++j) {
    if (not v_sm[j].mapped) continue;
    total_bp += int(v_m[j].dbPos[1] - v_m[j].dbPos[0] + 1);
    if (v_m[j].dbPos[1] < tsd[0].start - flank_len)
      total_bp_left += int(v_m[j].dbPos[1] - v_m[j].dbPos[0] + 1);
    else if (v_m[j].dbPos[0] > tsd[tsd.size() - 1].end + flank_len)
      total_bp_right += int(v_m[j].dbPos[1] - v_m[j].dbPos[0] + 1);
    else if (tsd.size() == 2 and v_m[j].dbPos[0] > tsd[0].end + flank_len
	and v_m[j].dbPos[1] < tsd[1].start - flank_len)
      total_bp_mid += int(v_m[j].dbPos[1] - v_m[j].dbPos[0] + 1);
  }

  // check if fragment completely captures either TSD
  for (size_t i = 0; i < tsd.size(); ++i) {
    for (size_t j = 0; j < v_sm.size(); ++j) {
      if (v_sm[j].mapped
	  and v_sm[j].mqv >= min_mqv
	  and v_m[j].dbPos[0] <= tsd[i].start - flank_len
	  and v_m[j].dbPos[1] >= tsd[i].end + flank_len) {
	// captures this TSD!
	tsd[i].count++;
	LOG(1) << "[" << v_sm[0].name << "]: captures tsd [" << i + 1 << "]\n";
	//return;
      }
    }
  }

  if (v_sm.size() != 2) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; unpaired\n";
    return;
  }

  if (v_m[0].dbPos[1] - v_m[0].dbPos[0] + 1 < min_read_len
      or v_m[1].dbPos[1] - v_m[1].dbPos[0] + 1 < min_read_len) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; one read too small\n";
    return;
  }

  if (not v_sm[0].mapped or not v_sm[1].mapped) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; not both mapped\n";
    return;
  }

  if (v_sm[0].mqv < min_mqv and v_sm[1].mqv < min_mqv) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; neither read has min mqv\n";
    return;
  }

  if (not v_sm[0].flags[1] or not v_sm[1].flags[1]) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; not proper pair\n";
    return;
  }

  if (v_sm[0].mqv < min_mqv or v_sm[1].mqv < min_mqv) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; one read doesn't have min mqv\n";
    return;
  }

  const Pairing * pairing = NULL;
  int rg_idx = -1;
  if (fnp == NULL) {
    // get read group info from SAM tags
    size_t k = 0;
    while (k < v_sm[0].rest.size() and v_sm[0].rest[k].key != "RG") ++k;
    if (k >= v_sm[0].rest.size()) {
      cerr << "could not determine read group for clone: " << clone_name << "\n";
      exit(EXIT_FAILURE);
    }
    ReadGroup * rg_p = global::rg_set.find_by_name(v_sm[0].rest[k].value);
    if (rg_p == NULL) {
      cerr << "error: missing read group [" << v_sm[0].rest[k].value
	   << "] of clone [" << clone_name << "]\n";
      exit(EXIT_FAILURE);
    }
    pairing = rg_p->get_pairing();
    rg_idx = global::rg_set.get_idx(rg_p);
  } else {
    // use full name parser to get read group info
    Clone c;
    int nip;
    fnp(v_sm[0].name, c, nip);
    pairing = c.pairing;

    // HACK: clone should store read group, not pairing pointer
    for (size_t i = 0; i < global::rg_set.rg_list.size(); ++i) {
      if (pairing == &global::rg_set.rg_list[i].pairing) {
	rg_idx = i;
	break;
      }
    }
    if (rg_idx < 0) {
      cerr << "hack failed; pairing pointer not found\n";
      exit(EXIT_FAILURE);
    }
  }

  // if this is the null allele and fragment length is too small, ignore
  int frag_len = pairing->get_t_len(v_m[0], 0, v_m[1], 0);
  if (tsd.size() == 1 and frag_len < pairing->mean - expected_insert_size/2) {
    LOG(1) << "[" << v_sm[0].name << "]: discarding; frag_len too small\n";
    return;
  }

  // read pair, both mapped, neither read captures a TSD
  long long left_end = min(v_m[0].dbPos[0], v_m[1].dbPos[0]);
  long long right_end = max(v_m[0].dbPos[1], v_m[1].dbPos[1]);

  if (tsd.size() == 1) {
    if (left_end <= tsd[0].start - flank_len
	and right_end >= tsd[0].end + flank_len) {
      LOG(1) << "[" << v_sm[0].name << "]: straddles single tsd\n";
      cluster[0][rg_idx].push_back(frag_len);
    }
  } else {
    if (left_end <= tsd[0].start - flank_len and
	right_end >= tsd[1].end + flank_len) {
      LOG(1) << "[" << v_sm[0].name << "]: straddles both tsds\n";

      cluster[0][rg_idx].push_back(pairing->get_t_len(v_m[0], 0, v_m[1], 0));
    } else if (left_end <= tsd[0].start - flank_len and
	       right_end >= tsd[0].end + flank_len and
	       right_end <= tsd[1].start - flank_len) {
      LOG(1) << "[" << v_sm[0].name << "]: straddles left tsd\n";
      cluster[1][rg_idx].push_back(frag_len);
    } else if (left_end >= tsd[0].end + flank_len and
	       left_end <= tsd[1].start - flank_len and
	       right_end >= tsd[1].end + flank_len) {
      LOG(1) << "[" << v_sm[0].name << "]: straddles right tsd\n";
      cluster[2][rg_idx].push_back(frag_len);
    }
  }
}

void
print_cluster_evidence(const vector<vector<int>> & v, ostream & os)
{
  for (size_t i = 0; i < v.size(); ++i) {
    if (i > 0) os << ";";
    os << global::rg_set.rg_list[i].get_num_id() << ":"
       << strtk::join(",", v[i]);
  }
}

string
default_cnp(const string& s)
{
  return string(s);
}

void
usage(ostream& os)
{
  os << "use: " << prog_name << " [ -l <pairing_file> ] -t <l_start>,<l_end> [ -t <r_start>,<r_end> ] [ <file> ]\n";
}

int
main(int argc, char* argv[])
{
  prog_name = string(argv[0]);
  {
    char * p;
    if ((p = getenv("MIN_MQV")) != NULL) min_mqv = atoi(p);
    if ((p = getenv("MAX_NM")) != NULL) max_nm = atoi(p);
    if ((p = getenv("MIN_READ_LEN")) != NULL) min_read_len = atoi(p);
    if ((p = getenv("MIN_READ_LEN_LEFT")) != NULL) min_read_len_left = atoi(p);
    if ((p = getenv("FLANK_LEN")) != NULL) flank_len = atoi(p);
    if ((p = getenv("MIN_NON_REPEAT_BP")) != NULL) min_non_repeat_bp = atoi(p);
  }

  cnp = default_cnp;
  string pairing_file;
  string fasta_file;

  char c;
  while ((c = getopt(argc, argv, "af:l:t:s:PN:g:vh")) != -1) {
    switch (c) {
    case 'a':
      is_alt = true;
      break;
    case 'f':
      fasta_file = optarg;
      break;
    case 'l':
      pairing_file = optarg;
      break;
    case 't':
      if (optarg[0] != '.') {
	if (tsd.size() >= 2) {
	  cerr << "wrong number of tsds\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	string tmp(optarg);
	size_t i = tmp.find(',');
	if (i == string::npos) {
	  cerr << "error parsing tsd specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	int start = atoi(tmp.substr(0, i).c_str());
	int end = atoi(tmp.substr(i+1).c_str());
	if (start < 0 or end < start) {
	  cerr << "error parsing tsd specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	if (tsd.size() == 1 and start < tsd[0].end) {
	  cerr << "tsds in wrong order\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	tsd.push_back(TSD(start, end));
      }
      break;
    case 's':
      {
	string tmp(optarg);
	size_t i = tmp.find(',');
	if (i == string::npos) {
	  cerr << "error parsing region start/end specification: " << tmp << "\n";
	  usage(cerr);
	  exit(EXIT_FAILURE);
	}
	reg_start = atoi(tmp.substr(0, i).c_str());
	reg_end = atoi(tmp.substr(i+1).c_str());
      }
      break;
    case 'P':
      cnp = cloneNameParser;
      fnp = fullNameParser;
      break;
    case 'g':
      global::default_rg_name = optarg;
      break;
    case 'N':
      global::num_threads = atoi(optarg);
      break;
    case 'v':
      global::verbosity++;
      break;
    case 'h':
      usage(cout);
      exit(EXIT_SUCCESS);
    default:
      cerr << "unrecognized option: " << c << endl;
      usage(cerr);
      exit(EXIT_FAILURE);
    }
  }

  if (optind + 1 < argc) {
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  if (is_alt) {
    // retrieve actual sequence from fasta file
    igzstream fasta_is(fasta_file);
    readFasta(fasta_is, global::refDict);
  }

  // load pairing file
  {
    igzstream pairing_is(pairing_file);
    if (!pairing_is) {
      cerr << "error opening pairing file: " << pairing_file << "\n";
      exit(EXIT_FAILURE);
    }
    global::rg_set.load(pairing_is);
  }

  //if (global::verbosity > 0) clog << "number of threads: " << num_threads << '\n';
  if (tsd.size() == 1) {
    cluster = vector<vector<vector<int>>>
      (1, vector<vector<int>>(global::rg_set.rg_list.size()));
  } else if (tsd.size() == 2) {
    cluster = vector<vector<vector<int>>>
      (3, vector<vector<int>>(global::rg_set.rg_list.size()));
  } else {
    cerr << "wrong number of tsds\n";
    usage(cerr);
    exit(EXIT_FAILURE);
  }

  LOG(1) << "pairing file: [" << pairing_file << "]\n";
  LOG(1) << "tsd 1: [" << tsd[0].start << "," << tsd[0].end << ")\n";
  if (tsd.size() > 1) {
    LOG(1) << "tsd 2: [" << tsd[1].start << "," << tsd[1].end << ")\n";
  }
  LOG(1) << "region limits: [" << reg_start << "," << reg_end << "]\n";
  LOG(1) << "min_mqv: [" << min_mqv << "]\n";
  LOG(1) << "max_nm: [" << max_nm << "]\n";
  LOG(1) << "min_read_len: [" << min_read_len << "]\n";
  LOG(1) << "min_read_len_left: [" << min_read_len_left << "]\n";
  LOG(1) << "flank_len: [" << flank_len << "]\n";
  LOG(1) << "internal naming: [" << (cnp == default_cnp? "no" : "yes") << "]\n";

  {
    igzstream mapIn(optind < argc? argv[optind] : "-");
    if (!mapIn) {
      cerr << "error opening mappings file: " << argv[optind] << endl;
      exit(EXIT_FAILURE);
    }

    SamMappingSetGen map_gen(&mapIn, cnp, NULL, &global::refDict, not is_alt);
    pair<string,vector<SamMapping> >* m = map_gen.get_next();
    int n_fragments = 0;
    while (m != NULL) {
      ++n_fragments;
      process_mapping_set(m->first, m->second);
      delete m;
      m = map_gen.get_next();
    }
  }

  if (tsd.size() == 1) {
    cout << tsd[0].count << "\t";
    print_cluster_evidence(cluster[0], cout);
  } else {
    cout << tsd[0].count << "\t";
    cout << tsd[1].count << "\t";
    print_cluster_evidence(cluster[0], cout);
    cout << "\t";
    print_cluster_evidence(cluster[1], cout);
    cout << "\t";
    print_cluster_evidence(cluster[2], cout);
  }
  cout << "\t"
       << (reg_start < tsd[0].start - flank_len ?
	   double(total_bp_left) / double (tsd[0].start - flank_len - reg_start + 1)
	   : 0)
       << "\t"
       << (reg_end > tsd[tsd.size() - 1].end + flank_len ?
	   double(total_bp_right)
	   / double (reg_end - tsd[tsd.size() - 1].end - flank_len + 1)
	   : 0);
  if (tsd.size() == 2) {
    cout << "\t"
	 << (tsd[0].end + flank_len < tsd[1].end - flank_len ?
	     double(total_bp_mid)
	     / double (tsd[1].end - flank_len - tsd[0].end - flank_len + 1)
	     : 0);
  }
  cout << "\n";

  return EXIT_SUCCESS;
}
