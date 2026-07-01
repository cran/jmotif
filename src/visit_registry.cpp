#include <jmotif.h>

// constructor
//
// ``seed`` makes the shuffled visit order -- and hence the early-abandoning
// search trajectory and its distance-call count -- reproducible. A negative
// seed (the default) keeps the historical std::random_device behavior (a fresh,
// nondeterministic order each run); seed >= 0 seeds the std::mt19937 engine so
// repeated runs are byte-identical. The discord *result* is order-independent
// either way.
//
VisitRegistry::VisitRegistry( int capacity, int seed ) {

  indexes = std::vector<int>(capacity);
  registry = std::vector<bool>(capacity);

  for( int i = 0; i < capacity; i++ ) {
    indexes[i] = i;
  }

  std::mt19937 g;
  if (seed < 0) {
    std::random_device rd;
    g.seed(rd());
  } else {
    g.seed((unsigned) seed);
  }
  std::shuffle(std::begin(indexes), std::end(indexes), g);

  unvisited_count = capacity;
  size = capacity;

}

// destructor
//
VisitRegistry::~VisitRegistry() {
  registry = std::vector<bool>();
  indexes = std::vector<int>();
}

// gets next unvisited position... following the HOTSAX heuristics this need to be a random
// position, so we can abandon the search earlier, if possible
//
int VisitRegistry::getNextUnvisited(){
  int random_index = -1;
    do {
      if(indexes.size()<=0){
        return -1;
      }
      random_index = indexes.back();
      indexes.pop_back();
    } while ( registry[random_index] );
    return random_index;
}

// marks a position visited, takes care about the counter
//
void VisitRegistry::markVisited(int idx){
  if(idx>=0 && idx<registry.size()){
    if(registry[idx]){
      return;
    }else{
      unvisited_count = unvisited_count - 1;
      registry[idx] = true;
    }
  }
}

// marks an interval as visited
//
void VisitRegistry::markVisited(int start, int end){
  for(int i=start; i<end; i++){
    if(i>=0 && i<registry.size()){
      if(registry[i]){
        continue;
      }else{
        unvisited_count = unvisited_count - 1;
        registry[i] = true;
      }
    }
  }
}

// check if the position has been marked as visited
//
bool VisitRegistry::isVisited(int idx){
  if(idx>=0 && idx<registry.size()){
    return(registry[idx]);
  }
  return false;
}
