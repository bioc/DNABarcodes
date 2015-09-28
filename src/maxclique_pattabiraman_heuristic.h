// maxclique_pattabiraman_heuristic.h

#ifndef __MAXCLIQUE_PATTABIRAMAN_HEURISTIC_H_INCLUDED__ 
#define __MAXCLIQUE_PATTABIRAMAN_HEURISTIC_H_INCLUDED__

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "maxclique.h"

#include <iostream>
#include <utility>
#include <memory>

#include <boost/graph/adjacency_list.hpp>

class MaxCliquePattabiramanHeuristic : public MaxClique { 

        public:
                MaxCliquePattabiramanHeuristic(); // default constructor

                std::vector<int> max_clique(const MaxClique::Graph &g, const size_t N, const size_t lower_bound);
                static std::vector<int> static_max_clique(const MaxClique::Graph &g, const size_t N, const size_t lower_bound);

        private:
                typedef boost::graph_traits< MaxClique::Graph >::out_edge_iterator out_edge_iterator_t;
                typedef std::pair<out_edge_iterator_t, out_edge_iterator_t> out_edge_iterator_range_t;
                typedef MaxClique::Graph::vertex_descriptor vertex_t;

                static std::vector<int> clique_heuristic(const MaxClique::Graph &g, std::set<int> &U, const size_t size, const size_t max, bool &found);

};

#endif 
