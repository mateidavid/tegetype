#include "Pairing.hpp"

#include <iomanip>

#include "strtk/strtk.hpp"
#include "Read.hpp"
#include "globals.hpp"

void
Pairing::parse_token(const string& t)
{
  strtk::std_string::token_list_type token_list;
  strtk::split("=", t, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();
  if (itr == token_list.end()) {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }
  string key(itr->first, itr->second);
  ++itr;
  if (itr == token_list.end()) {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }
  int i_val = atoi(itr->first);
  if (!key.compare("paired")) {
    paired = i_val;
  } else if (!key.compare("st_diff")) {
    st_diff = i_val;
  } else if (!key.compare("min")) {
    min = i_val;
  } else if (!key.compare("max")) {
    max = i_val;
  } else if (!key.compare("mean")) {
    mean = i_val;
  } else if (!key.compare("stddev")) {
    stddev = i_val;
  } else if (!key.compare("r1_len")) {
    r_len[0] = i_val;
  } else if (!key.compare("r2_len")) {
    r_len[1] = i_val;
  } else if (!key.compare(0, 2, "gc")) {
    if (frag_rate.size() == 0) frag_rate = vector<double>(100);
    int bin_idx = atoi(key.substr(2).c_str());
    frag_rate[bin_idx] = atof(itr->first);
  } else {
    cerr << "could not parse pairing token: " << t << endl;
    exit(1);
  }    
}

Pairing::Pairing(const string& s)
{
  strtk::std_string::token_list_type token_list;
  strtk::split(",", s, back_inserter(token_list));
  strtk::std_string::token_list_type::iterator itr = token_list.begin();

  while (itr != token_list.end()) {
    parse_token(string(itr->first, itr->second));
    ++itr;
  }
}

int
Pairing::get_mp_st(const Mapping& mapping, int r_st) const
{
  int st = (r_st + mapping.st) % 2;
  return (st + st_diff) % 2;
}

int
Pairing::get_mp_st(int read_st, int mapping_st) const
{
  int st = (read_st + mapping_st) % 2;
  return (st + st_diff) % 2;
}

bool
Pairing::is_mp_downstream(const Mapping& mapping, int r_st) const
{
  assert(mapping.qr != NULL);
  int st = (r_st + mapping.st) % 2;
  int nip = mapping.qr->nip;
  int mp_st = get_mp_st(mapping, r_st);
  int sign_5p_diff;

  if (nip == 0) {
    sign_5p_diff = (st == 0? 1 : -1);
  } else {
    sign_5p_diff = (mp_st == 0? -1 : 1);
  }

  return mean * sign_5p_diff > 0;
}

bool
Pairing::is_mp_downstream(int read_nip, int read_st, int mapping_st) const
{
  int st = (read_st + mapping_st) % 2;
  int mp_st = get_mp_st(read_st, mapping_st);
  int sign_5p_diff;

  if (read_nip == 0) {
    sign_5p_diff = (st == 0? 1 : -1);
  } else {
    sign_5p_diff = (mp_st == 0? -1 : 1);
  }

  return mean * sign_5p_diff > 0;
}

vector<Interval<long long int> >
Pairing::get_mp_pos(const Mapping& mapping, int r_st) const
{
  vector<Interval<long long int> > result;

  assert(mapping.qr != NULL);
  assert(mapping.qr->mp != NULL);
  int st = (r_st + mapping.st) % 2;
  int nip = mapping.qr->nip;
  int mp_st = get_mp_st(mapping, r_st);
  int mp_rlen = mapping.qr->mp->len;

  long long int pos_5p = mapping.dbPos[st];
  int sign_5p_diff;
  if (nip == 0) {
    sign_5p_diff = (st == 0? 1 : -1);
  } else {
    sign_5p_diff = (mp_st == 0? -1 : 1);
  }

  Interval<long long int> mp_pos_5p;
  mp_pos_5p[0] = pos_5p + sign_5p_diff * min;
  mp_pos_5p[1] = pos_5p + sign_5p_diff * max;
  mp_pos_5p.sort();

  int sign_mp_st = (mp_st == 0? 1 : -1);
  Interval<long long int> tmp;
  tmp[0] = mp_pos_5p[0];
  tmp[1] = mp_pos_5p[0] + sign_mp_st * (mp_rlen - 1);
  tmp.sort();
  result.push_back(tmp);
  tmp[0] = mp_pos_5p[1];
  tmp[1] = mp_pos_5p[1] + sign_mp_st * (mp_rlen - 1);
  tmp.sort();
  result.push_back(tmp);
    
  return result;
}

int
Pairing::get_t_len(const Mapping& mapping0, int r0_st,
		   const Mapping& mapping1, int r1_st) const
{
  int st0 = (r0_st + mapping0.st) % 2;
  int st1 = (r1_st + mapping1.st) % 2;
  long long pos_5p_0 = mapping0.dbPos[st0];
  long long pos_5p_1 = mapping1.dbPos[st1];
  if (st0 == 0) {
    return int(pos_5p_1 - pos_5p_0);
  } else {
    return int(pos_5p_0 - pos_5p_1);
  }
}


bool
Pairing::pair_concordant(const Mapping& mapping0, int r0_st,
			 const Mapping& mapping1, int r1_st) const
{
  //int st = (r0_st + mapping0.st) % 2;
  int mp_st = get_mp_st(mapping0, r0_st);
  vector<Interval<long long int> > mp_pos = get_mp_pos(mapping0, r0_st);

  return (r1_st + mapping1.st) % 2 == mp_st and
    mapping1.dbPos[0] >= mp_pos[0][0] and mapping1.dbPos[0] <= mp_pos[1][0];
}


ostream &
operator <<(ostream & os, const Pairing& pairing)
{
  if (!pairing.paired) {
    os << "paired=0";
  } else {
    os << "paired=1"
       << ",st_diff=" << pairing.st_diff
       << ",min=" << pairing.min
       << ",max=" << pairing.max
       << ",mean=" << pairing.mean
       << ",stddev=" << pairing.stddev;
    if (pairing.frag_rate.size() > 0) {
      for (int i = 0; i < 100; ++i) {
	os << ",gc" << i << "=" << scientific << setprecision(3) << pairing.frag_rate[i];
      }
    }
  }
  return os;
}


ReadGroup::ReadGroup(const string & s)
{
  strtk::std_string::token_list_type token_list;
  strtk::split("\t", s, back_inserter(token_list));
  if (token_list.size() < 3) {
    cerr << "cannot parse read group line: " << s << endl;
    exit(1);
  }
  strtk::std_string::token_list_type::iterator itr = token_list.begin();
  strtk::split(",", string(itr->first, itr->second),
	       strtk::range_to_type_back_inserter(name));
  ++itr;
  num_id = string(itr->first, itr->second);
  ++itr;
  pairing = Pairing(string(itr->first, itr->second));
}


ostream &
operator <<(ostream & os, const ReadGroup & rg)
{
  os << strtk::join(",", rg.get_names()) << "\t"
     << rg.get_num_id() << "\t"
     << *rg.get_pairing();
  return os;
}


void
ReadGroupSet::add(const ReadGroup & rg)
{
  auto names = rg.get_names();
  for_each(names.begin(), names.end(), [&] (const string & name) {
      if (rg_name_dict.count(name) > 0) {
	cerr << "error: duplicate read group name: " << name << "\n";
	exit(EXIT_FAILURE);
      }
    });
  string num_id = rg.get_num_id();
  if (rg_num_id_dict.count(num_id) > 0) {
    cerr << "error:duplicate read group num id: " << num_id << "\n";
    exit(EXIT_FAILURE);
  }
  if (rg_list.size() > 0) {
    if (num_id.size() != rg_num_id_len) {
      cerr << "error: rg num ids of different sizes: " << rg_num_id_len
	   << " vs " << num_id.size() << "\n";
      exit(EXIT_FAILURE);
    }
  } else {
    rg_num_id_len = num_id.size();
  }

  int rg_idx = rg_list.size();
  rg_list.push_back(rg);
  for_each(names.begin(), names.end(), [&] (const string & name) {
      rg_name_dict.insert(pair<string,int>(name, rg_idx));
    });
  rg_num_id_dict.insert(pair<string,int>(num_id, rg_idx));
  if (global::verbosity > 0) clog << "added rg [" << rg << "]\n";
}

void
ReadGroupSet::load(istream & is)
{
  string s;
  while (getline(is, s))
    add(s);
  if (is.bad()) {
    cerr << "error reading pairing file\n";
    exit(EXIT_FAILURE);
  }
}

ReadGroup *
ReadGroupSet::find_by_name(const string & s)
{
  if (rg_name_dict.count(s) == 0)
    return NULL;
  return find_by_idx(rg_name_dict[s]);
}

ReadGroup *
ReadGroupSet::find_by_num_id(const string & s) {
  if (rg_num_id_dict.count(s) == 0)
    return NULL;
  return find_by_idx(rg_num_id_dict[s]);
}

int
ReadGroupSet::get_idx(const ReadGroup * rg_p) {
  if (rg_list.size() == 0 or rg_p < &rg_list[0] or rg_p > &rg_list[rg_list.size() - 1]) {
    cerr << "get_idx error: read group [" << *rg_p << "] not in current set:" << *this;
    exit(EXIT_FAILURE);
  }
  return int(rg_p - &rg_list[0]);
}


ostream &
operator << (ostream & os, const ReadGroupSet & rgs) {
  for_each(rgs.rg_list.begin(), rgs.rg_list.end(), [&os] (const ReadGroup & rg) {
      os << rg << "\n";
    });
  return os;
}


double
get_expected_complete_span(const SQDict & sq_dict, const ReadGroupSet & rg_set,
			   const string & chr, long long start_1, long long end_1)
{
  if (start_1 < end_1) {
    cerr << "error: start=" << start_1 << " < end=" << end_1 << "\n";
    exit(EXIT_FAILURE);
  }
  auto it = sq_dict.find(chr);
  if (it == sq_dict.end()) {
    cerr << "contig " << chr << " not found\n";
    exit(EXIT_FAILURE);
  }
  const Contig& ctg = it->second;
  if (ctg.name.length() == 0) {
    cerr << "error: couldn't find chr=" << chr << "\n";
    exit(EXIT_FAILURE);
  } else if (start_1 > ctg.len) {
    cerr << "error: start=" << start_1 << " outside chr=" << chr
	 << " len=" << ctg.len << "\n";
    exit(EXIT_FAILURE);
  } else if (end_1 > ctg.len) {
    cerr << "error: end=" << end_1 << " outside chr=" << chr
	 << " len=" << ctg.len << "\n";
    exit(EXIT_FAILURE);
  }

  double res = 0;

  for (size_t i = 0; i < rg_set.rg_list.size(); ++i) {
    const ReadGroup & rg = rg_set.rg_list[i];
    int rounded_mean = rg.get_pairing()->mean - (rg.get_pairing()->mean % 5);

    // check if frags from this rg can fully span interval
    if (end_1 - start_1 + 1 > rounded_mean) continue;

    long long reg_start_1 = max(1ll, end_1 - rounded_mean + 1);
    long long reg_end_1 = min(start_1 + rounded_mean - 1, ctg.len);

    int count_n = 0;
    int count_gc = 0;
    long long first_1 = reg_start_1;
    long long last_1 = reg_start_1 - 1;
    while (last_1 < reg_end_1) {
      ++last_1;
      char c = ctg.seq[0][last_1 - 1];
      if (c == 'N' or c == 'n') ++count_n;
      if (c == 'G' or c == 'g' or c == 'C' or c == 'c') ++count_gc;
      if (last_1 - first_1 + 1 > rounded_mean) {
	c = ctg.seq[0][first_1 - 1];
	if (c == 'N' or c == 'n') --count_n;
	if (c == 'G' or c == 'g' or c == 'C' or c == 'c') --count_gc;
	++first_1;
      }
      if (last_1 - first_1 + 1 == rounded_mean) {
	// process current region
	int bin_idx = int((double(count_gc) / (rounded_mean + 1)) * 100);
	res += rg.get_pairing()->frag_rate[bin_idx];
      }
    }
  }

  return res;
}
