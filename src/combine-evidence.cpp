#include <iostream>
#include <cstdlib>
#include <vector>

#include "igzstream.hpp"
#include "strtk/strtk.hpp"
#include "globals.hpp"
#include "Fasta.hpp"

using namespace std;

string prog_name;


void
usage(ostream & os)
{
  os << "use: " << prog_name << " -f <fasta_file> -l <pairing_file> -L <lib_file> -r <ref_evidence_file> -a <alt_evidence_file>\n";
}


int
main(int argc, char* argv[])
{
  prog_name = argv[0];
  string fasta_file;
  string pairing_file;
  string lib_file;
  string ref_evidence_file;
  string alt_evidence_file;

  char c;
  while ((c = getopt(argc, argv, "vf:l:L:r:a:h")) != -1) {
    switch (c) {
    case 'v':
      global::verbosity++;
      break;
    case 'f':
      fasta_file = optarg;
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

  if (fasta_file == "") { cerr << "missing fasta file\n"; exit(EXIT_FAILURE); }
  if (pairing_file == "") { cerr << "missing pairing file\n"; exit(EXIT_FAILURE); }
  if (lib_file == "") { cerr << "missing lib file\n"; exit(EXIT_FAILURE); }
  if (ref_evidence_file == "") { cerr << "missing ref_evidence file\n"; exit(EXIT_FAILURE); }
  if (alt_evidence_file == "") { cerr << "missing alt_evidence file\n"; exit(EXIT_FAILURE); }

  // load pairing file
  {
    igzstream pairing_is(pairing_file);
    global::rg_set.load(pairing_is);
  }

  // load fasta file
  {
    igzstream fasta_is(fasta_file);
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
  }

  return EXIT_SUCCESS;
}
