// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include "maxclique_pattabiraman_heuristic.h"

#include <boost/foreach.hpp>

MaxCliquePattabiramanHeuristic::MaxCliquePattabiramanHeuristic()  {

}

std::vector<int> MaxCliquePattabiramanHeuristic::max_clique(const MaxClique::Graph &g, const size_t N, const size_t lower_bound) {
  return static_max_clique(g, N, lower_bound);
}

std::vector<int> MaxCliquePattabiramanHeuristic::static_max_clique(const MaxClique::Graph &g, const size_t N, const size_t lower_bound) {
        // Size of max clique 
        size_t global_max = lower_bound;
        std::vector<int> global_best_clique;

        // Go through every vertex
//#pragma omp parallel for schedule(dynamic) default(none) shared(g, N, global_max, global_best_clique)
        for (size_t i = 0; i < N; i++) {
    // Check for abort
    R_CheckUserInterrupt();

    // if at least as many degrees as the current top candidate, continue

    if (boost::out_degree(i, g) >= global_max) {
      std::set<int> U; // clique of this vertex, initially empty

      // Go through outer edges of vertex i
      for(out_edge_iterator_range_t outEdges = boost::out_edges(i, g); outEdges.first != outEdges.second; ++outEdges.first) {
        // Get number of this vertex (slighty confusing to do it this way)
        int j = (int) outEdges.first->m_target;

        // If j has enough out degrees, add j to the clique
        if (boost::out_degree(j,g) >= global_max) {
                U.insert(j);
        }
      }
      // Calculate clique closure of vertex i
      bool found = false;
      std::vector<int> clique = clique_heuristic(g, U, 1, global_max, found);
      if (found) {
        // *ding* *ding* *ding*
        clique.push_back(i);
//#pragma omp critical
//{
        global_max = clique.size(); 
        global_best_clique = clique;
//}
      }
    }
  }
  return global_best_clique;
}

std::vector<int> MaxCliquePattabiramanHeuristic::clique_heuristic(const MaxClique::Graph &g, std::set<int> &U, const size_t size, const size_t max, bool &found) {
  std::vector<int> clique;

  if (U.empty()) {
    if (size > max) {
            /* *ding* *ding* *ding* */
            // Tail end of recursion
            found = true;
    }
    return clique;
  }

  // Find vertex with max degree
  int vertex_u = -1; // Will actually be initialised in the next foreach loop
  int max_degree = -1;

  BOOST_FOREACH( int test_vertex, U ) {
    int d_test_vertex = boost::out_degree( test_vertex, g );

    if (d_test_vertex > max_degree) {
      max_degree = d_test_vertex;
      vertex_u = test_vertex;
    }
  }

  // remove u from U
  U.erase(vertex_u);

  // find all out edges of u
  std::set<int> Nprime;
  for(out_edge_iterator_range_t outEdges = boost::out_edges(vertex_u, g); outEdges.first != outEdges.second; ++outEdges.first) {
    int vertex_w = (int) outEdges.first->m_target;
    size_t d_w = boost::out_degree(vertex_w, g);

    if (d_w >= max) {
      Nprime.insert(vertex_w);
    }
  }

  // intersect out vertices of u with U
  std::set<int> Unew;
  set_intersection(U.begin(),U.end(),Nprime.begin(),Nprime.end(), std::inserter( Unew, Unew.begin() ) );

  // call recursively
  clique = clique_heuristic(g, Unew, size+1, max, found);
  if (found) {
    /* *ding* *ding* *ding* */
    clique.push_back(vertex_u);
  }
  return clique;
}

