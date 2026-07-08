#ifndef JMOTIF_h
#define JMOTIF_h
//
//#include <RcppArmadillo.h>
#include <vector>
#include <deque>
#include <random>
#include <Rcpp.h>
using namespace Rcpp;
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp14")]]

//
// Define the letters array
//
const char LETTERS[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
                        'q', 'r', 's', 't', 'u',  'v', 'w', 'x', 'y', 'z'};

//
// SAX stack
//
// NumericVector znorm(NumericVector ts, double threshold);
//
NumericVector paa(NumericVector ts, int paa_num);
//
NumericVector alphabet_to_cuts(int a_size);
//
CharacterVector series_to_chars(NumericVector ts, int a_size);
//
CharacterVector series_to_string(NumericVector ts, int a_size);
//
std::map<int, std::string> sax_via_window(NumericVector ts, int w_size, int paa_size, int a_size,
    CharacterVector nr_strategy, double n_threshold);
//
std::map<int, CharacterVector> sax_by_chunking(NumericVector ts, int paa_size,
                                               int a_size, double n_threshold);

//
// Strings operations
//
char idx_to_letter(int idx);
//
int letter_to_idx(char letter);
//
IntegerVector letters_to_idx(CharacterVector str);
//
bool is_equal_str(CharacterVector a, CharacterVector b);
//
bool is_equal_mindist(CharacterVector a, CharacterVector b);

//
// Distance
//
double euclidean_dist(NumericVector seq1, NumericVector seq2);
//
double early_abandoned_dist(NumericVector seq1, NumericVector seq2, double upper_limit);

//
// SAX-VSM
//
Rcpp::DataFrame series_to_wordbag(NumericVector ts, int w_size, int paa_size, int a_size,
    CharacterVector nr_strategy, double n_threshold);
//
Rcpp::DataFrame manyseries_to_wordbag(NumericMatrix data, int w_size, int paa_size, int a_size,
    CharacterVector nr_strategy, double n_threshold);
//
Rcpp::DataFrame bags_to_tfidf(Rcpp::List data);
//
Rcpp::DataFrame cosine_sim(Rcpp::List data);

//
// Discord
//
// Defines the discord record which consists of its index and the distance to NN
struct discord_record {
  int index;
  double nn_distance;
  int dist_calls;
};
//
class VisitRegistry {
public:
  int size;
  std::vector<bool> registry;
  std::vector<int> indexes;
  int unvisited_count;
  VisitRegistry( int capacity, int seed = -1 );
  int getNextUnvisited();
  void markVisited(int idx);
  void markVisited(int start, int end);
  bool isVisited(int idx);
  ~VisitRegistry();
};
//
Rcpp::DataFrame find_discords_brute_force(NumericVector ts, int w_size, int discords_num,
                                          double n_threshold, int seed);
//

//
// HOT-SAX
//
Rcpp::DataFrame find_discords_hotsax(NumericVector ts, int w_size, int paa_size,
                                      int a_size, double n_threshold, int discords_num,
                                      int seed);

//
// REPAIR
//

// the basic work string token
//
class repair_symbol {
public:
  std::string payload;
  int str_index;
  repair_symbol() { // a generic constructor
    str_index = -1;
  };
  repair_symbol(const std::string str, int index); // a more advance constructor
  virtual bool is_guard(){ // if this is a guard? no... overriden in Guard
    return false;
  }
  std::string to_string(){ // get the payload
    return std::string(payload);
  }
};

// the symbol (token) wrapper for the string data structure0
//
class repair_symbol_record {
public:
  repair_symbol* payload;
  repair_symbol_record* prev;
  repair_symbol_record* next;
  repair_symbol_record( repair_symbol* symbol );
};

// the grammar rule
//
class repair_rule {
public:
  int id;
  int rule_use;
  repair_symbol* first;
  repair_symbol* second;
  std::string expanded_rule_string;
  std::vector<int> occurrences;
  repair_rule(){
    id = -1; rule_use = 0; first = nullptr; second = nullptr;
  };
  std::string get_rule_string();
};

// the guard symbol: a special symbol that holds the rule
//
class repair_guard: public repair_symbol {
public:
  repair_rule* r;
  repair_guard();
  repair_guard(repair_rule* rule, int idx){
    r = rule;
    payload = r->get_rule_string();
    str_index = idx;
  }
  bool is_guard(){
    return true;
  }
  std::string get_expanded_string(){
    return r->expanded_rule_string;
  }
};

//
// the priority queue elements...
//
class repair_digram {
public:
  std::string digram;
  int freq;
  repair_digram(const std::string str, int index);
};

//
// the priority queue (a dobly-linked list) node
//
class repair_pqueue_node {
public:
  repair_pqueue_node* prev;
  repair_pqueue_node* next;
  repair_digram* payload;
  repair_pqueue_node() {
    payload = nullptr;
    prev = nullptr;
    next = nullptr;
  }
  repair_pqueue_node(repair_digram* d) {
    payload = d;
    prev = nullptr;
    next = nullptr;
  }
};

//
// Arena/object pool for the RePair working objects.
//
// Every symbol/record/guard/rule/digram/pqueue-node the grammar inference allocates was
// previously `new`'d and NEVER freed -- a leak that is benign for a single one-shot call
// but accumulates (20x big_none peaked ~397 MB RSS). The pool owns all of them in
// std::deque containers (NOT std::vector: the algorithm holds raw pointers --
// r0[i]->next, rule->first/second, pqueue links -- across many later allocations, and
// deque guarantees stable element addresses where vector realloc would dangle them).
// A single pool, declared as the FIRST local of _str_to_repair_grammar, is destroyed
// LAST (after the result is built), freeing everything in one shot on every exit path
// (normal or exception). Object lifetimes within the call are unchanged.
//
class repair_pool {
public:
  std::deque<repair_symbol> symbols;
  std::deque<repair_symbol_record> records;
  std::deque<repair_guard> guards;
  std::deque<repair_rule> rules;
  std::deque<repair_digram> digrams;
  std::deque<repair_pqueue_node> nodes;
  repair_symbol* new_symbol(const std::string& str, int index) {
    symbols.emplace_back(str, index);
    return &symbols.back();
  }
  repair_symbol_record* new_record(repair_symbol* sym) {
    records.emplace_back(sym);
    return &records.back();
  }
  repair_guard* new_guard(repair_rule* rule, int idx) {
    guards.emplace_back(rule, idx);
    return &guards.back();
  }
  repair_rule* new_rule() {
    rules.emplace_back();
    return &rules.back();
  }
  repair_digram* new_digram(const std::string& str, int freq) {
    digrams.emplace_back(str, freq);
    return &digrams.back();
  }
  repair_pqueue_node* new_node(repair_digram* d) {
    nodes.emplace_back(d);
    return &nodes.back();
  }
};

//
// the priority queue taking care about repair digrams ordering.
//
// Bucketed (count-indexed) design (Larsson-Moffat): buckets[f] heads a doubly-linked
// list of every digram whose frequency is exactly f; max_count tracks the highest
// non-empty bucket. Max-select (dequeue) is an O(1) array index + downward scan past
// empties; every +/-1 frequency change is an O(1) unlink + push-front. This replaces
// the previous frequency-sorted single linked list whose enqueue / update did O(Q)
// linear band-walks. The `nodes` map (digram -> node) is retained for O(1) locate.
//
// Order discipline (preserves rule numbering): push-FRONT within a bucket reproduces
// the previous list's LIFO-within-band order exactly for the seed pass (pure enqueues),
// so the dequeue order -- hence rule-id assignment -- is unchanged.
//
class repair_priority_queue {
public:
  std::vector<repair_pqueue_node*> buckets; // buckets[freq] -> head of that freq's list
  int max_count;                            // highest non-empty bucket index
  std::unordered_map<std::string, repair_pqueue_node*> nodes; // fastmap <digram> -> <node>
  repair_pool* pool;                        // arena that owns the queue's nodes (not owned here)
  repair_priority_queue() {
    max_count = 0;
    pool = nullptr;
  }
  repair_digram* enqueue(repair_digram* digram);
  repair_digram* dequeue();
  repair_digram* update_digram_frequency(std::string *digram_string, int new_value);
  bool contains_digram(std::string *digram_string);
  void remove_node(repair_pqueue_node* node);
  std::string to_string();
};

//
// RRA
//
class rule_record {
public:
  int rule_id;
  std::string rule_string;
  std::string expanded_rule_string;
  std::vector<int> rule_occurrences;
  std::vector<std::pair<int, int>> rule_intervals;
  // int rule_use;
  rule_record() {
    rule_id = -1;
    // rule_use = 0;
  }
};

struct rra_discord_record {
  int rule;
  int start;
  int end;
  double nn_distance;
  int distance_calls;
};

std::unordered_map<int, rule_record> _str_to_repair_grammar(std::string s);
Rcpp::List str_to_repair_grammar(CharacterVector str);

Rcpp::DataFrame find_discords_rra(NumericVector series, int w_size, int paa_size,
  int a_size, CharacterVector nr_strategy, double n_threshold, int discords_num,
  int seed);

//
// Utilities
//
NumericVector col_means(NumericMatrix m);
//
NumericVector subseries(NumericVector ts, int start, int end);
//
int armaRand();

// internal high performance computing
//
std::vector<double> _alphabet_to_cuts(int a_size);
int _count_spaces(std::string *s);
double _mean(std::vector<double> *ts, int *start, int *end);
std::vector<double> _znorm(const std::vector<double>& ts, double threshold);
void _znorm_slice(const std::vector<double>& ts, int start, int end,
                  double threshold, std::vector<double>& out);
std::vector<double> _paa2(const std::vector<double>& ts, int paa_num);
void _paa2_range(const std::vector<double>& ts, int start, int end, int paa_num,
                 std::vector<double>& out);
double _euclidean_dist(std::vector<double>* seq1, std::vector<double>* seq2);
double _early_abandoned_dist(std::vector<double>* seq1, std::vector<double>* seq2,
                             double upper_limit);
std::vector<double> _subseries(std::vector<double>* ts, int start, int end);
//
std::string _series_to_string(std::vector<double> ts, int a_size);
bool _is_equal_mindist(std::string a, std::string b);
//
std::unordered_map<int, std::string> _sax_via_window(
    const std::vector<double>& ts, int w_size, int paa_size, int a_size,
    std::string nr_strategy, double n_threshold);

#endif
