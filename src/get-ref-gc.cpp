#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <boost/regex.hpp>

#include "strtk/strtk.hpp"
#include "igzstream.hpp"

using namespace std;
using strtk::parse;
using strtk::for_each_line;


namespace global {
  int verbosity;
  int progress;
  int n_gc_bins = 100;
  int step_len = 5;
}


int
main(int argc, char* argv[])
{
  string prog_name(argv[0]);
  vector<int> win_size;
  vector<string> chr_group;

  char c;
  int w;
  while ((c = getopt(argc, argv, "vp:w:c:")) != -1) {
    switch (c) {
    case 'v':
      ++global::verbosity;
      break;
    case 'p':
      global::progress = atoi(optarg);
      break;
    case 'w':
      w = atoi(optarg);
      if (w > 0 and w % 5 == 0)
	win_size.push_back(atoi(optarg));
      else
	clog << "ignoring window size [" << w << "]\n";
      break;
    case 'c':
      chr_group.push_back(optarg);
      break;
    default:
      cerr << "unrecognized option: " << c << "\n";
      exit(EXIT_FAILURE);
    }
  }

  if (win_size.size() == 0) {
    cerr << "no window sizes\n";
    exit(EXIT_FAILURE);
  }

  if (optind < argc - 1) {
    cerr << "use: " << prog_name << " [options] [<gc5_file>]\n";
    exit(EXIT_FAILURE);
  }

  igzstream in_file(optind < argc? argv[optind] : "-");

  if (global::verbosity > 0) {
    clog << "using chr_group:";
    for_each(chr_group.begin(), chr_group.end(), [] (const string& s) {
	clog << " " << s;
      });
    clog << " default\n";
    clog << "using win_size:";
    for_each(win_size.begin(), win_size.end(), [] (int w) {
	clog << " " << w;
      });
    clog << "\n";
  }

  // create structures
  vector<deque<int> > q(win_size.size());
  vector<int> gc_cnt(win_size.size(), 0);
  vector<pair<long long,long long> > q_pos(win_size.size());
  vector<vector<vector<int> > > cnt(chr_group.size() + 1,
				    vector<vector<int> >(win_size.size(),
							 vector<int>(global::n_gc_bins, 0)));
  vector<boost::regex> p;
  for_each(chr_group.begin(), chr_group.end(), [&p] (const string& s) {
      p.push_back(boost::regex(s));
    });

  size_t crt_chr_group;

  for_each_line(in_file, [&] (const string& line) {
      strtk::ignore_token ignore;
      string s;
      string chr;

      if (parse(line, " \t", ignore, s, ignore) and parse(s, "=", ignore, chr))
	{
	  if (global::verbosity > 0) clog << "found chr [" << chr << "]: ";

	  // find group
	  int first_match = -1;
	  int last_match = -1;
	  for (size_t i = 0; i < p.size(); ++i) {
	    if (boost::regex_match(chr, p[i])) {
	      if (first_match < 0) first_match = i;
	      last_match = i;
	    }
	  }
	  if (last_match == -1) {
	    crt_chr_group = chr_group.size();
	    if (global::verbosity > 0) clog << "using leftover group\n";
	  } else if (last_match == first_match) {
	    crt_chr_group = size_t(last_match);
	    if (global::verbosity > 0)
	      clog << "using group [" << chr_group[crt_chr_group] << "]\n";
	  } else {
	    clog << "matches more than one group: [" << chr_group[first_match] << "], [" << chr_group[last_match] << "]\n";
	    exit(EXIT_FAILURE);
	  }

	  // clear all queues
	  for (size_t i = 0; i < q.size(); ++i) {
	    q[i].clear();
	    gc_cnt[i] = 0;
	    q_pos[i].first = -1;
	    q_pos[i].second = -1;
	  }
	}
      else
	{
	  long long crt_pos;
	  int crt_cnt;
	  int a;
	  double b;
	  if (parse(line, " \t", crt_pos, a)) {
	    crt_cnt = (a * 5) / 100;
	  } else if (parse(line, " \t", crt_pos, b)) {
	    crt_cnt = int((b * 5) / 100 + 0.5);
	  } else {
	    cerr << "error parsing line: [" << line << "]\n";
	    exit(EXIT_FAILURE);
	  }

	  if (global::verbosity > 1) clog << "crt_pos=[" << crt_pos << "] crt_cnt=[" << crt_cnt << "]\n";

	  for (size_t i = 0; i < q.size(); ++i) {
	    // if range is not contiguous, clear queue first
	    if (q_pos[i].second >= 0 and crt_pos != q_pos[i].second + 1) {
	      if (global::verbosity > 1) clog << "skipping to pos=[" << crt_pos << "]\n";
	      q[i].clear();
	      gc_cnt[i] = 0;
	      q_pos[i].first = -1;
	      q_pos[i].second = -1;
	    }

	    // add current region to queue
	    q[i].push_back(crt_cnt);
	    gc_cnt[i] += crt_cnt;
	    if (q_pos[i].first < 0) q_pos[i].first = crt_pos;
	    q_pos[i].second = crt_pos + (global::step_len - 1);

	    // remove excess from queue
	    while (q_pos[i].second - q_pos[i].first + 1 > win_size[i]) {
	      gc_cnt[i] -= q[i].front();
	      q[i].pop_front();
	      q_pos[i].first += global::step_len;
	    }

	    // if full region captured, record gc count
	    if (q_pos[i].second - q_pos[i].first + 1 == win_size[i]) {
	      int bin_idx = int((double(gc_cnt[i]) / (win_size[i] + 1)) * global::n_gc_bins);
	      cnt[crt_chr_group][i][bin_idx]++;
	      if (global::verbosity > 1) clog << "region:[" << q_pos[i].first << "," << q_pos[i].second << "] gc:[" << gc_cnt[i] << "] bin_idx:[" << bin_idx << "]\n";
	    }
	  }

	}
    });

  // print by win_size, then bin number
  for (size_t i = 0; i < win_size.size(); ++i) {
    for (int j = 0; j < global::n_gc_bins; ++j) {
      cout << win_size[i] << "\t"
	   << j << "\t"
	   << int(ceil((double(j) / global::n_gc_bins) * (win_size[i] + 1)));
      for (size_t k = 0; k <= chr_group.size(); ++k) {
	cout << "\t" << cnt[k][i][j] * global::step_len;
      }
      cout << "\n";
    }
  }

  return EXIT_SUCCESS;
}
