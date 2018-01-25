#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <list>
#include <map>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <stdlib.h>
#include <stdio.h>
#include <Rcpp.h>
#include "gtf_parser.h"

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;


////////////////////////////////////////////////////////////////////////////////
// itoa: Returns an integer's string representation.
string itoa(int value, int base = 10) {

  // check that the base if valid
  if (base < 2 || base > 36) return '\0';

  char result[32] = {0};
  char* ptr = result, *ptr1 = result, tmp_char;
  int tmp_value, max = 32;
  char d[72] = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz";

  do {
    tmp_value = value;
    value /= base;
    *ptr++ = d[35 + (tmp_value - value * base)];
  } while ( value && max--);

  // Apply negative sign
  if (tmp_value < 0) *ptr++ = '-';
  *ptr-- = '\0';
  while(ptr1 < ptr) {
    tmp_char = *ptr;
    *ptr--= *ptr1;
    *ptr1++ = tmp_char;
  }

  return string(result);
}

////////////////////////////////////////////////////////////////////////////////
// unquote: Removes quotes and semicolon from GTF attributes.
string unquote(const string& s) {
  int n = s.find('"');

  return s.substr(n + 1, n + s.rfind('"') - 1);
}

////////////////////////////////////////////////////////////////////////////////
// defineContiguousRegion: Determines if two exons do define a contiguous region
// Overlaps version sould remove the '+ 1' from the greater side.
inline bool defineContiguousRegion(exon e, exon e2) {

  return e2.range.first <= e.range.second + 1;
  // and e2.range.second + 1 >= e.range.first; // BUT e < e2 in our setting
}

////////////////////////////////////////////////////////////////////////////////
// isLower: Comparator to order exons based in their genomic coordinates.
inline bool isLower(exon e, exon e2) { return e.range.first < e2.range.first; }

////////////////////////////////////////////////////////////////////////////////
// extend: Fusion two contiguous/overlapping exons into one. Return value needed
// to know if we should remove e2.
bool extend(exon& e, exon& e2) {

  if (not defineContiguousRegion(e, e2)) return false;
  e.range.first = min(e.range.first, e2.range.first);
  e.range.second = max(e.range.second, e2.range.second);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// reduce: Reduce contiguous/overlapping exons while adding exon lenghts to gene
// model.
void reduce(gene& g) {

  g.e.sort(isLower);
  for (list<exon>::iterator it = g.e.begin(); it != g.e.end(); ) {
    exon& e = *it++;
    while (it != g.e.end() and extend(e, *it)) it = g.e.erase(it);
    e.len = e.range.second - e.range.first + 1;
    g.len += e.len;
    g.range.first = g.range.first? min(g.range.first, e.range.first)
      : e.range.first;
    g.range.second = max(g.range.second, e.range.second);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Parallelized computation of lengths
struct SummarizeGene : public RcppParallel::Worker {
  map<string, gene>& annotation;

  // Basic constructor
  SummarizeGene(map<string, gene>& ann) : annotation(ann) {}

  // Worker job
  void operator()(size_t begin, size_t end) {
    map<string, gene>::iterator it_end = annotation.begin();
    map<string, gene>::iterator it = annotation.begin();
    advance(it_end, end);
    advance(it, begin);

    while (it != it_end) reduce((it++)->second);
  }
};

////////////////////////////////////////////////////////////////////////////////
// buildAnnotDF: Returns formatted DataFrame with annotation from genes.
// SAF version (with per-row exonic ranges) needed for Rsubread.
// Compressed version with per-row genes.
DataFrame buildAnnotDF(map<string, gene>& ann, string level = "exon") {
  vector<string> gid, strand, start, end, chr;
  map<string, gene>::iterator it;
  vector<int> len;

  if(level == "exon") {
    for (it = ann.begin(); it != ann.end(); ++it) {
      gene& g = it->second;
      list<exon>::iterator e;

      for(e = g.e.begin(); e != g.e.end(); ++e) {
        gid.push_back(g.id);
        chr.push_back(e->chr);
        start.push_back(itoa((e->range).first));
        end.push_back(itoa((e->range).second));
        strand.push_back(e->strand);
        len.push_back(e->len);
      }
    }
  } else if (level == "gene") {
    for (it = ann.begin(); it != ann.end(); ++it) {
      gene& g = it->second;

      gid.push_back(g.id);
      chr.push_back(g.chr);
      start.push_back(itoa(g.range.first));
      end.push_back(itoa(g.range.second));
      strand.push_back(g.strand);
      len.push_back(g.len);
    }
  }

  return DataFrame::create(_["GeneID"] = gid,
                           _["Chr"] = chr,
                           _["Start"] = start,
                           _["End"] = end,
                           _["Strand"] = strand,
                           _["Length"] = len,
                           _["stringsAsFactors"] = false);
}


////////////////////////////////////////////////////////////////////////////////
// parseGtf: Builds Annotation DataFrame upon gene length and GC value
// computation from a GTF file.
// [[Rcpp::export]]
List gtf2SAF(std::string fileName) {
  ifstream gtfFile(fileName.c_str());
  map<string, gene> ann;
  string line;

  // Load features
  while (getline(gtfFile, line)) {
    if (line[0] == '#') continue; // Ignore comment lines
    string chr, s, featType, strand, gid, tid;
    istringstream iss(line);
    int start, end;

    iss >> chr >> s >> featType >> start >> end >> s >> strand >> s;
    if (featType != "exon") continue; // Only exonic features
    while (iss >> s) {
      if (s == "gene_id") {
        iss >> gid;
        gid = unquote(gid);
      } else if (s == "transcript_id") {
        iss >> tid;
        tid = unquote(tid);
      }
    }
    if (!ann.count(gid)) ann[gid] = gene(gid, chr, strand);
    ann[gid].e.push_back(exon(chr, strand, make_pair(start, end)));
    ann[gid].t[tid].push_back(exon(ann[gid].e.back()));
  }

  // Compute lengths ang %GC values parallelized by genes
  SummarizeGene ret(ann);
  parallelFor(0, ann.size(), ret);

  // Returns both versions, SAF and compressed formatted annotations
  return List::create(_["exon"] = buildAnnotDF(ann),
                      _["gene"] = buildAnnotDF(ann, "gene"));
}

