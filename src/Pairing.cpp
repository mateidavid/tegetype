#include "Pairing.hpp"

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
  int val=atoi(itr->first);
  if (!key.compare("paired")) {
    paired = val;
  } else if (!key.compare("st_diff")) {
    st_diff = val;
  } else if (!key.compare("min")) {
    min = val;
  } else if (!key.compare("max")) {
    max = val;
  } else if (!key.compare("mean")) {
    mean = val;
  } else if (!key.compare("stddev")) {
    stddev = val;
  } else if (!key.compare("r1_len")) {
    r_len[0] = val;
  } else if (!key.compare("r2_len")) {
    r_len[1] = val;
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

ostream&
operator <<(ostream& ostr, const Pairing& pairing)
{
  if (!pairing.paired) 
    ostr << "paired=0";
  else
    ostr << "paired=1"
	 << ",st_diff=" << pairing.st_diff
	 << ",min=" << pairing.min
	 << ",max=" << pairing.max
	 << ",mean=" << pairing.mean
	 << ",stddev=" << pairing.stddev;
  return ostr;
}

void
load_pairing(istream& istr, RGDict& rg_dict, RGDict& num_rg_dict, RGRGDict& rg_to_num_rg_dict)
{
  string s;
  while (getline(istr, s)) {
    ReadGroup rg(s);
    for_each(rg.name.begin(), rg.name.end(), [&] (const string & name) {
	if (rg_dict.count(name) != 0) {
	  cerr << "error: duplicate read group: " << name << endl;
	  exit(1);
	}
      });

    if (num_rg_dict.count(rg.num_id) != 0) {
      cerr << "error: duplicate numeric read group: " << rg.num_id << endl;
      exit(1);
    }
    if (num_rg_dict.size() > 0) {
      if ((int)rg.num_id.size() != global::num_rg_len) {
	cerr << "error: numeric RG ids of different sizes: " << rg.num_id.size()
	     << " vs " << global::num_rg_len << endl;
	exit(1);
      }
    } else {
      global::num_rg_len = rg.num_id.size();
    }

    rg.idx = num_rg_dict.size();
    for_each(rg.name.begin(), rg.name.end(), [&] (const string & name) {
	rg_dict.insert(pair<string,ReadGroup>(name, rg));
	rg_to_num_rg_dict.insert(pair<string,string>(name, rg.num_id));
      });
    num_rg_dict.insert(pair<string,ReadGroup>(rg.num_id, rg));
    if (global::verbosity > 0) clog << "added RG [" << strtk::join(",", rg.name)
				    << "] with pairing [" << rg << "]\n";
  }
  if (istr.bad()) {
    cerr << "error reading pairing file" << endl;
    exit(1);
  }
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
  (Pairing &)(*this) = Pairing(string(itr->first, itr->second));
}