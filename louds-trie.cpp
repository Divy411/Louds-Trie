#include "louds-trie.hpp"

#ifdef _MSC_VER
 #include <intrin.h>
 #include <immintrin.h>
#else  // _MSC_VER
 #include <x86intrin.h>
#endif  // _MSC_VER

#include <cassert>
#include <vector>

namespace louds {
namespace {

uint64_t Popcnt(uint64_t x) {
#ifdef _MSC_VER
  return __popcnt64(x);
#else  // _MSC_VER
  return __builtin_popcountll(x);
#endif  // _MSC_VER
}

uint64_t Ctz(uint64_t x) {
#ifdef _MSC_VER
  return _tzcnt_u64(x);
#else  // _MSC_VER
  return __builtin_ctzll(x);
#endif  // _MSC_VER
}

struct BitVector {
  struct Rank {
    uint32_t abs_hi;
    uint8_t abs_lo;
    uint8_t rels[3];

    uint64_t abs() const {
      return ((uint64_t)abs_hi << 8) | abs_lo;
    }
    void set_abs(uint64_t abs) {
      abs_hi = (uint32_t)(abs >> 8);
      abs_lo = (uint8_t)abs;
    }
  };

  vector<uint64_t> words;
  vector<Rank> ranks;
  vector<uint32_t> selects;
  uint64_t n_bits;

  BitVector() : words(), ranks(), selects(), n_bits(0) {}

  uint64_t get(uint64_t i) const {
    return (words[i / 64] >> (i % 64)) & 1UL;
  }
  void set(uint64_t i, uint64_t bit) {
    if (bit) {
      words[i / 64] |= (1UL << (i % 64));
    } else {
      words[i / 64] &= ~(1UL << (i % 64));
    }
  }

  void add(uint64_t bit) {
    if (n_bits % 256 == 0) {
      words.resize((n_bits + 256) / 64);
    }
    set(n_bits, bit);
    ++n_bits;
  }
  // build builds indexes for rank and select.
  void build() {
    uint64_t n_blocks = words.size() / 4;
    uint64_t n_ones = 0;
    ranks.resize(n_blocks + 1);
    for (uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      ranks[block_id].set_abs(n_ones);
      for (uint64_t j = 0; j < 4; ++j) {
        if (j != 0) {
          uint64_t rel = n_ones - ranks[block_id].abs();
          ranks[block_id].rels[j - 1] = rel;
        }

        uint64_t word_id = (block_id * 4) + j;
        uint64_t word = words[word_id];
        uint64_t n_pops = Popcnt(word);
        uint64_t new_n_ones = n_ones + n_pops;
        if (((n_ones + 255) / 256) != ((new_n_ones + 255) / 256)) {
          uint64_t count = n_ones;
          while (word != 0) {
            uint64_t pos = Ctz(word);
            if (count % 256 == 0) {
              selects.push_back(((word_id * 64) + pos) / 256);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_ones = new_n_ones;
      }
    }
    ranks.back().set_abs(n_ones);
    selects.push_back(words.size() * 64 / 256);
  }

  // rank returns the number of 1-bits in the range [0, i).
  uint64_t rank(uint64_t i) const {
    uint64_t word_id = i / 64;
    uint64_t bit_id = i % 64;
    uint64_t rank_id = word_id / 4;
    uint64_t rel_id = word_id % 4;
    uint64_t n = ranks[rank_id].abs();
    if (rel_id != 0) {
      n += ranks[rank_id].rels[rel_id - 1];
    }
    n += Popcnt(words[word_id] & ((1UL << bit_id) - 1));
    return n;
  }
  // select returns the position of the (i+1)-th 1-bit.
  uint64_t select(uint64_t i) const {
    const uint64_t block_id = i / 256;
    uint64_t begin = selects[block_id];
    uint64_t end = selects[block_id + 1] + 1UL;
    if (begin + 10 >= end) {
      while (i >= ranks[begin + 1].abs()) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks[middle].abs()) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks[rank_id].abs();

    uint64_t word_id = rank_id * 4;
    if (i < ranks[rank_id].rels[1]) {
      if (i >= ranks[rank_id].rels[0]) {
        word_id += 1;
        i -= ranks[rank_id].rels[0];
      }
    } else if (i < ranks[rank_id].rels[2]) {
      word_id += 2;
      i -= ranks[rank_id].rels[1];
    } else {
      word_id += 3;
      i -= ranks[rank_id].rels[2];
    }
    return (word_id * 64) + Ctz(_pdep_u64(1UL << i, words[word_id]));
  }

  uint64_t size() const {
    return sizeof(uint64_t) * words.size()
      + sizeof(Rank) * ranks.size()
      + sizeof(uint32_t) * selects.size();
  }
};

struct Level {
  BitVector louds;
  BitVector outs;
  vector<uint8_t> labels;
  uint64_t offset;

  Level() : louds(), outs(), labels(), offset(0) {}

  uint64_t size() const;
};

uint64_t Level::size() const {
  return louds.size() + outs.size() + labels.size();
}

}  // namespace

class TrieImpl {
 public:
  TrieImpl();
  ~TrieImpl() {}

  void add(const string &key);
  void build();
  int64_t lookup(const string &query) const;

  uint64_t n_keys() const {
    return n_keys_;
  }
  uint64_t n_nodes() const {
    return n_nodes_;
  }
  uint64_t size() const {
    return size_;
  }
  // Gives read only access to the internal levels of the trie
  const std::vector<Level>& levels() const { 
    return levels_; 
  }

  // Checks if the trie includes the empty string as a key
  bool has_empty_key() const { 
    return levels_[0].outs.get(0) != 0; 
  }

 private:
  vector<Level> levels_;
  uint64_t n_keys_;
  uint64_t n_nodes_;
  uint64_t size_;
  string last_key_;
};

TrieImpl::TrieImpl()
  : levels_(2), n_keys_(0), n_nodes_(1), size_(0), last_key_() {
  levels_[0].louds.add(0);
  levels_[0].louds.add(1);
  levels_[1].louds.add(1);
  levels_[0].outs.add(0);
  levels_[0].labels.push_back(' ');
}

void TrieImpl::add(const string &key) {
  if (key.empty()) {
    levels_[0].outs.set(0, 1);
    ++levels_[1].offset;
    ++n_keys_;
    return;
  }
  assert(key > last_key_); //moved it downward so that the trie could first handle the special case of adding  
                          //the empty string before enforcing the lexicographic order for all subsequent (non-empty) keys.
  if (key.length() + 1 >= levels_.size()) {
    levels_.resize(key.length() + 2);
  }
  uint64_t i = 0;
  for ( ; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    uint8_t byte = key[i];
    if ((i == last_key_.length()) || (byte != level.labels.back())) {
      level.louds.set(levels_[i + 1].louds.n_bits - 1, 0);
      level.louds.add(1);
      level.outs.add(0);
      level.labels.push_back(key[i]);
      ++n_nodes_;
      break;
    }
  }
  for (++i; i < key.length(); ++i) {
    Level &level = levels_[i + 1];
    level.louds.add(0);
    level.louds.add(1);
    level.outs.add(0);
    level.labels.push_back(key[i]);
    ++n_nodes_;
  }
  levels_[key.length() + 1].louds.add(1);
  ++levels_[key.length() + 1].offset;
  levels_[key.length()].outs.set(levels_[key.length()].outs.n_bits - 1, 1);
  ++n_keys_;
  last_key_ = key;
}

void TrieImpl::build() {
  uint64_t offset = 0;
  for (uint64_t i = 0; i < levels_.size(); ++i) {
    Level &level = levels_[i];
    level.louds.build();
    level.outs.build();
    offset += levels_[i].offset;
    level.offset = offset;
    size_ += level.size();
  }
}

int64_t TrieImpl::lookup(const string &query) const {
  if (query.length() >= levels_.size()) {
    return false;
  }
  uint64_t node_id = 0;
  for (uint64_t i = 0; i < query.length(); ++i) {
    const Level &level = levels_[i + 1];
    uint64_t node_pos;
    if (node_id != 0) {
      node_pos = level.louds.select(node_id - 1) + 1;
      node_id = node_pos - node_id;
    } else {
      node_pos = 0;
    }

    // Linear search implementation
    // for (uint8_t byte = query[i]; ; ++node_pos, ++node_id) {
    //   if (level.louds.get(node_pos) || level.labels[node_id] > byte) {
    //     return -1;
    //   }
    //   if (level.labels[node_id] == byte) {
    //     break;
    //   }
    // }

    // Binary search implementation
    uint64_t end = node_pos;
    uint64_t word = level.louds.words[end / 64] >> (end % 64);
    if (word == 0) {
      end += 64 - (end % 64);
      word = level.louds.words[end / 64];
      while (word == 0) {
        end += 64;
        word = level.louds.words[end / 64];
      }
    }
    end += Ctz(word);
    uint64_t begin = node_id;
    end = begin + end - node_pos;

    uint8_t byte = query[i];
    while (begin < end) {
      node_id = (begin + end) / 2;
      if (byte < level.labels[node_id]) {
        end = node_id;
      } else if (byte > level.labels[node_id]) {
        begin = node_id + 1;
      } else {
        break;
      }
    }
    if (begin >= end) {
      return -1;
    }
  }
  const Level &level = levels_[query.length()];
  if (!level.outs.get(node_id)) {
    return false;
  }
  return level.offset + level.outs.rank(node_id);
}

Trie::Trie() : impl_(new TrieImpl) {}

Trie::~Trie() {
  delete impl_;
}

void Trie::add(const string &key) {
  return impl_->add(key);
}

void Trie::build() {
  impl_->build();
}

int64_t Trie::lookup(const string &query) const {
  return impl_->lookup(query);
}

uint64_t Trie::n_keys() const {
  return impl_->n_keys();
}

uint64_t Trie::n_nodes() const {
  return impl_->n_nodes();
}

uint64_t Trie::size() const {
  return impl_->size();
}




// Determines the range of child nodes in the next level for a given node.
// reused from node traversal logic in the lookup() function.
inline void ChildRange(const Level& level, uint64_t& nodeId, uint64_t& begin, uint64_t& endBitsBegin) {
    uint64_t nodePos = 0;

    if (nodeId != 0) {
        // Select locates the position of the (nodeId - 1)th 1-bit and move to the next bit
        nodePos = level.louds.select(nodeId - 1) + 1;

        // // Adjust nodeId to index within labels
        nodeId = nodePos - nodeId;
    }
  
    uint64_t end = nodePos;

    // Extract the bits starting at the position of this node
    uint64_t word = level.louds.words[end / 64] >> (end % 64);

    // If the current word is all zero scan forward until we find a 1-bit
    if (word == 0) {
        end += 64 - (end % 64);
        word = level.louds.words[end / 64];

        while (word == 0) {
            end += 64;
            word = level.louds.words[end / 64];
        }
    }

    // Count trailing zeros to locate the next 1-bit
    end += Ctz(word);

    // Begin index is updated nodeId. endBitsBegin is relative position in labels
    begin = nodeId;
    endBitsBegin = begin + end - nodePos;
}
// A depth-first iterator that walks through all the keys stored in a LOUDS trie
class KeyIterator {
 public:
 // Constructor - it needs a TrieImpl to walk over.
  explicit KeyIterator(const TrieImpl& T)
      : levels_(T.levels()), finished_(false) {

    //  If the trie has no key, we can return.
    if (levels_.size() <= 1) { 
      finished_ = true; 
      return; 
    }

    // Start with the root frame. Actual keys start from level 1.
    Frame root;
    root.depth = 0;
    root.node_id = 0;
    // Determine the child range for root from level 1
    {
      uint64_t nid = 0, b = 0, e = 0;
      ChildRange(levels_[1], nid, b, e);
      root.begin = b; // Start index of root’s children
      root.end = e;  // End index of root’s children
      root.next = b; // Next child to explore is at begin
    }
    root.yielded = true;  // root itself is not a key 
    frames_.push_back(root); // Pushing root frame onto DFS stack
  }

  // Returns true and writes the next key into out if available.
  // Otherwise returns false to mean that it is end of iteration.
  bool next(std::string& out) {
    if (finished_) return false;

    while (!frames_.empty()) {
      Frame& f = frames_.back();

      // Yield the key if this node is terminal and not already returned
      if (f.depth >= 1 && !f.yielded) {
        f.yielded = true;

        if (levels_[f.depth].outs.get(f.node_id)) {
          out.assign(key_.begin(), key_.end());
          return true;
        }
      }

      // Try to descend into the next child node if one is available
      if (f.next < f.end) {
        const uint64_t child_id = f.next++;
        const uint64_t child_depth = f.depth + 1;

        // Append this child's label to the key being built
        key_.push_back(static_cast<char>(levels_[child_depth].labels[child_id]));

        Frame child;
        child.depth = child_depth;
        child.node_id = child_id;
        child.yielded = false;

        // Check if this node has children by looking into the next level
        if (child_depth + 1 < levels_.size()) {
          uint64_t nid = child.node_id, b = 0, e = 0;
          ChildRange(levels_[child_depth + 1], nid, b, e); 
          child.begin = b; 
          child.end = e; 
          child.next = b;
        } else {
          // Reached a leaf node
          child.begin = child.end = child.next = 0;
        }
        frames_.push_back(child);
        continue;
      }

      // If there are no more children to explore at this node we will backtrack
      if (f.depth >= 1) 
      {
        key_.pop_back(); // Remove the label we added when descending
      }
      frames_.pop_back();  // Remove this frame from the stack
    }

    // No more keys left
    finished_ = true;
    return false;
  }

 private:
  struct Frame {
    uint64_t depth = 0;       
    uint64_t node_id = 0;     // global child index at this level
    uint64_t begin = 0, end = 0, next = 0;  // children range in labels[]
    bool yielded = false;     // whether terminal check was done
  };
  const std::vector<Level>& levels_; 
  std::vector<Frame> frames_; 
  std::string key_; //current key being built during traversal
                   // whem we descend through trie levels. It is modified when we push/pop frames.
  bool finished_;
};


// Merge two built tries into a new LOUDS trie (set union, duplicates removed).
// Complexity: Θ(total length of unique output keys); Extra space: O(max depth).
Trie merge_trie(const Trie& a, const Trie& b) {
  Trie out;

  // // If either trie contains the empty string we add it.
    if (a.impl_->has_empty_key() || b.impl_->has_empty_key()) {

    out.add(std::string());
  }
  
  // Set up depth-first iterators for both input tries.
  KeyIterator iter_a(*a.impl_);
  KeyIterator iter_ab(*b.impl_);
  
  std::string key_a, key_b;

  // Initialize both iterators to fetch their first key (if any).
  bool has_a = iter_a.next(key_a);
  bool has_b = iter_ab.next(key_b);
  
  // Walk through both tries in sorted order and merge keys.
  while (has_a || has_b) {
    if (has_a && has_b) {
      if (key_a < key_b) {
        out.add(key_a);     // key_a is smaller, add it to the result

        has_a = iter_a.next(key_a); 

      } else if (key_b < key_a) {
        out.add(key_b);  

        has_b = iter_ab.next(key_b);

      } else {
        // If both keys are equal, add it only once
        out.add(key_a);

        has_a = iter_a.next(key_a);  // advance both to skip duplicates

        has_b = iter_ab.next(key_b);
      }
    } else if (has_a) {
      // Only iterator A has keys left. dump them all
      out.add(key_a);
      has_a = iter_a.next(key_a);
    } else {
      // Only iterator B has keys left. dump them all
      out.add(key_b);
      has_b = iter_ab.next(key_b);
    }
  }

  out.build();
  return out;
}


}  // namespace louds