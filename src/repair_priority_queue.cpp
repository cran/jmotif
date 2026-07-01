#include <jmotif.h>

// Bucketed (count-indexed) priority queue for RePair digrams (Larsson-Moffat).
//
// buckets[f] heads a doubly-linked list of all digrams with frequency exactly f.
// max_count is the highest non-empty bucket. The `nodes` map gives O(1) locate by
// digram string. Compared with the previous frequency-sorted single list, max-select
// and every frequency change are O(1) amortized instead of O(Q) linear walks.
//
// Order discipline: nodes are pushed to the FRONT of their bucket. During the initial
// seed pass (pure enqueues) this reproduces the previous list's LIFO-within-frequency
// order, so dequeue order -- and hence the rule-id assignment that depends on it --
// is byte-for-byte unchanged.

// ensure buckets[] can be indexed at `freq`
static inline void ensure_capacity(std::vector<repair_pqueue_node*>& buckets, int freq) {
  if (freq >= (int) buckets.size()) {
    buckets.resize(freq + 1, nullptr);
  }
}

// unlink a node from its current bucket list (does NOT touch the `nodes` map)
void repair_priority_queue::remove_node(repair_pqueue_node* node) {
  int f = node->payload->freq;
  if (nullptr != node->prev) {
    node->prev->next = node->next;
  } else if (f >= 0 && f < (int) buckets.size() && buckets[f] == node) {
    // node was the bucket head
    buckets[f] = node->next;
  }
  if (nullptr != node->next) {
    node->next->prev = node->prev;
  }
  node->prev = nullptr;
  node->next = nullptr;
  nodes.erase(node->payload->digram);
}

// enqueue a new digram; throws if the key already exists (matches prior contract)
repair_digram* repair_priority_queue::enqueue(repair_digram* digram) {

  if (nodes.find(digram->digram) != nodes.end()) {
    throw std::range_error("Inadmissible value, key exists...");
    return nullptr;
  }

  // allocate the node from the arena when one is attached, else fall back to new
  repair_pqueue_node* nn = (nullptr != pool) ? pool->new_node(digram)
                                             : new repair_pqueue_node(digram);

  int f = digram->freq;
  ensure_capacity(buckets, f);

  // push to the FRONT of bucket f (LIFO within the frequency band)
  nn->next = buckets[f];
  nn->prev = nullptr;
  if (nullptr != buckets[f]) {
    buckets[f]->prev = nn;
  }
  buckets[f] = nn;

  if (f > max_count) {
    max_count = f;
  }

  nodes.emplace(digram->digram, nn);
  return nn->payload;
}

// update a digram's frequency, relocating it to the proper bucket; evict if < 2
repair_digram* repair_priority_queue::update_digram_frequency(
    std::string *digram_string, int new_value) {

  std::unordered_map<std::string, repair_pqueue_node*>::iterator it =
    nodes.find(*digram_string);
  if (it == nodes.end()) {
    return nullptr;
  }
  repair_pqueue_node* node = it->second;

  // trivial case
  if (new_value == node->payload->freq) {
    return node->payload;
  }

  // a digram occurring < 2 times is no longer a replacement candidate -- evict it
  if (2 > new_value) {
    remove_node(node);
    return nullptr;
  }

  // relocate: unlink from the old bucket, set new freq, push-front into the new bucket.
  // (remove_node erases the nodes-map entry; re-add after the move.)
  remove_node(node);
  node->payload->freq = new_value;
  ensure_capacity(buckets, new_value);
  node->next = buckets[new_value];
  node->prev = nullptr;
  if (nullptr != buckets[new_value]) {
    buckets[new_value]->prev = node;
  }
  buckets[new_value] = node;
  if (new_value > max_count) {
    max_count = new_value;
  }
  nodes.emplace(node->payload->digram, node);

  return node->payload;
}

// pop the most-frequent digram: scan down from max_count to the first non-empty bucket
repair_digram* repair_priority_queue::dequeue() {
  while (max_count > 0 &&
         (max_count >= (int) buckets.size() || nullptr == buckets[max_count])) {
    max_count--;
  }
  if (max_count <= 0 || nullptr == buckets[max_count]) {
    return nullptr;
  }
  repair_pqueue_node* head = buckets[max_count];
  repair_digram* res = head->payload;
  // pop the bucket head
  buckets[max_count] = head->next;
  if (nullptr != head->next) {
    head->next->prev = nullptr;
  }
  nodes.erase(res->digram);
  return res;
}

// is this digram currently in the queue?
bool repair_priority_queue::contains_digram(std::string *digram_string) {
  return nodes.find(*digram_string) != nodes.end();
}

// print the whole queue, highest frequency first (debug)
std::string repair_priority_queue::to_string() {
  std::stringstream res;
  for (int f = max_count; f >= 0; f--) {
    if (f >= (int) buckets.size()) continue;
    for (repair_pqueue_node* p = buckets[f]; nullptr != p; p = p->next) {
      res << p->payload->digram << " : " << p->payload->freq << std::endl;
    }
  }
  return res.str();
}
