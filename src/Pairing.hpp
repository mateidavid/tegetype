#ifndef Pairing_hpp_
#define Pairing_hpp_

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "Interval.hpp"
#include "Mapping.hpp"
#include "DNASequence.hpp"


class Pairing
{
public:
  bool paired;
  int st_diff;
  int min;
  int max;
  int mean;
  int stddev;
  int r_len[2];
  vector<double> frag_rate;

  Pairing() { paired = false; }
  Pairing(const string&);

  void parse_token(const string&);
  int get_mp_st(const Mapping&, int) const;
  int get_mp_st(int, int) const;
  bool is_mp_downstream(const Mapping&, int) const;
  bool is_mp_downstream(int, int, int) const;
  vector<Interval<long long int> > get_mp_pos(const Mapping&, int) const;
  int get_t_len(const Mapping&, int, const Mapping&, int) const;
  bool pair_concordant(const Mapping&, int, const Mapping&, int) const;
};

ostream & operator <<(ostream &, const Pairing &);


class ReadGroup {
public:
  vector<string> name;
  string num_id;
  Pairing pairing;

  ReadGroup() {}
  ReadGroup(const string &);

  vector<string> get_names() const { return name; }
  string get_num_id() const { return num_id; }
  const Pairing * get_pairing() const { return &pairing; }
};

ostream & operator <<(ostream &, const ReadGroup &);


class ReadGroupSet {
public:
  vector<ReadGroup> rg_list;
  map<string,int> rg_name_dict;
  map<string,int> rg_num_id_dict;

  size_t rg_num_id_len;

  ReadGroupSet() : rg_num_id_len(0) {}

  void add(const ReadGroup &);
  void add(const string & s) { add(ReadGroup(s)); }
  void load(istream &);
  ReadGroup * find_by_name(const string &);
  ReadGroup * find_by_num_id(const string &);
  ReadGroup * find_by_idx(int i) { return &rg_list[i]; }
  int get_idx(const ReadGroup *);
};

ostream & operator <<(ostream &, const ReadGroupSet &);

// return expected number of fragments fully spanning the given closed 1-based interval,
// assuming a single allelic region
double get_expected_complete_span(const SQDict &, const ReadGroupSet &,
				  const string &, long long, long long);

#endif
