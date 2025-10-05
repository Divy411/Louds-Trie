#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "louds-trie.hpp"

using namespace std;
using namespace std::chrono;

static vector<string> read_lines(const char* path) {
  ifstream f(path);
  if (!f) { cerr << "Error opening " << path << "\n"; std::exit(1); }
  vector<string> out;
  string line;
  while (getline(f, line)) out.push_back(line);
  return out;
}

// Enforce add()’s precondition (strictly increasing, no dups)
static void sort_and_dedupe(vector<string>& v) {
  sort(v.begin(), v.end());

  v.erase(unique(v.begin(), v.end()), v.end());
}

static louds::Trie build_trie_from(const vector<string>& keys) {
  louds::Trie t;
  for (const auto& k : keys) t.add(k);
  t.build();

  return t;
}

static vector<string> set_union_vec(vector<string> A, vector<string> B) {
  A.insert(A.end(), B.begin(), B.end());
  sort_and_dedupe(A);
  return A;
}

// Check presence and that ids form a permutation of 0..n-1
static void verify_permutation_ids(const louds::Trie& t, const vector<string>& U) {
  vector<int64_t> ids; ids.reserve(U.size());

  for (const auto& s : U) {

    int64_t id = t.lookup(s);
    if (id < 0) {
      cerr << "[FAIL] missing key: \"" << s << "\"\n";
      assert(false);
    }
    ids.push_back(id);
  }
  sort(ids.begin(), ids.end());
  for (size_t i = 0; i < ids.size(); ++i) {

    if (ids[i] != static_cast<int64_t>(i)) {

      cerr << "[FAIL] ids not a permutation at i=" << i

           << " expected " << i << " got " << ids[i] << "\n";
      assert(false);
      
    }
  }
}

int main(int argc, char* argv[]) {
  ios::sync_with_stdio(false);

  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <keys1.txt> <keys2.txt>\n";
    return 1;
  }

  // Read
  vector<string> keys1 = read_lines(argv[1]);
  vector<string> keys2 = read_lines(argv[2]);

  // Ensure add() precondition even if files are already sorted
  sort_and_dedupe(keys1);
  sort_and_dedupe(keys2);

  // Build trie1
  auto t0 = high_resolution_clock::now();
  auto trie1 = build_trie_from(keys1);
  auto t1 = high_resolution_clock::now();

  // Build trie2
  auto trie2 = build_trie_from(keys2);
  auto t2 = high_resolution_clock::now();

  // Merge
  auto tm0 = high_resolution_clock::now();
  louds::Trie merged = louds::merge_trie(trie1, trie2);  // if merge_trie doesn't build, call merged.build() here
  auto tm1 = high_resolution_clock::now();

  const auto U = set_union_vec(keys1, keys2);

  // Print stats
  auto ns = [](auto a, auto b){ return (double)duration_cast<nanoseconds>(b-a).count(); };
  cout << "[Trie1] #keys=" << trie1.n_keys() << "  #nodes=" << trie1.n_nodes()<< "  size=" << trie1.size()
   << "  build=" << (keys1.empty()?0.0:ns(t0,t1)/keys1.size()) << " ns/key\n";

  cout << "[Trie2] #keys=" << trie2.n_keys() << "  #nodes=" << trie2.n_nodes() << "  size=" << trie2.size()
       << "  build=" << (keys2.empty()?0.0:ns(t1,t2)/keys2.size()) << " ns/key\n";

  cout << "[Merge] #keys=" << merged.n_keys() << "  #nodes=" << merged.n_nodes() << "  size=" << merged.size()
       << "  merge=" << (U.empty()?0.0:ns(tm0,tm1)/U.size()) << " ns/key\n";

  // Verify
  verify_permutation_ids(merged, U);
  cout << "Verification: \n";

  // Quick lookup timings
  auto s0 = high_resolution_clock::now();

  for (const auto& k : U) (void)merged.lookup(k);

  auto s1 = high_resolution_clock::now();
  cout << "[Merge] seq.lookup " << (U.empty()?0.0:ns(s0,s1)/U.size()) << " ns/key\n";

  vector<string> shuffled = U;
  std::mt19937 rng(std::random_device{}());
  std::shuffle(shuffled.begin(), shuffled.end(), rng);
  auto r0 = high_resolution_clock::now();
  for (const auto& k : shuffled) (void)merged.lookup(k);
  auto r1 = high_resolution_clock::now();
  cout << "[Merge] rnd.lookup " << (U.empty()?0.0:ns(r0,r1)/U.size()) << " ns/key\n";

  return 0;
}
