#include <string>
#include <list>
#include <map>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


typedef struct exon {
  pair<int, int> range;
  string strand;
  string chr;
  double pgc;
  int len;
  int gcc;

  exon(string c, string s, pair<int, int> r)
    : chr(c), strand(s), range(r), len(0), gcc(0), pgc(0.0) {}
} exon;

typedef list<exon> transcript;

typedef struct gene {
  map<string, transcript> t;
  pair<int, int> range;
  string strand;
  list<exon> e;
  string chr;
  string id;
  int len;

  gene() {}
  gene(string i, string c, string s)
    : id(i), chr(c), strand(s), len(0), range(make_pair(0,0)) {}
} gene;

