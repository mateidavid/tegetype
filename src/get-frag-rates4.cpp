#include <cstdlib>
#include <iostream>
#include <random>
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;


int
main(int argc, char * argv[])
{
  uniform_real_distribution<double> unif(0.0, 1.0);
  default_random_engine re;

  BamReader reader;
  reader.Open(argv[1]);
  SamHeader header = reader.GetHeader();

  BamAlignment m;
  while (reader.GetNextAlignmentCore(m)) {
    if (unif(re) > .001) continue;

    m.BuildCharData();
    string rg;
    char t;
    m.GetTagType(string("RG"), t);
    m.GetTag(string("RG"), rg);
    cout << m.Name << "\t" << reader.GetReferenceData()[m.RefID].RefName << ":" << m.Position << "\t" << reader.GetReferenceData()[m.MateRefID].RefName << ":" << m.MatePosition << "\t" << m.InsertSize << "\tRG:" << t << ":" << rg << "\n";

  }

  /*
  string rg;
  char t;
  m.GetTagType(string("RG"), t);
  m.GetTag(string("RG"), rg);
  cout << m.Name << "\t" << reader.GetReferenceData()[m.RefID].RefName << ":" << m.Position << "\t" << reader.GetReferenceData()[m.MateRefID].RefName << ":" << m.MatePosition << "\t" << m.InsertSize << "\tRG:" << t << ":" << rg << "\n";

  cout << "HasIndex: " << reader.HasIndex() << "\n";
  cout << "LocateIndex: " << reader.LocateIndex() << "\n";
  if (!reader.LocateIndex()) {
    cout << "CreateIndex: ... ";
    reader.CreateIndex();
    cout << "done\n";
  }
  */

  return EXIT_SUCCESS;
}
