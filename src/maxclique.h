// maxclique.h

#ifndef __MAXCLIQUE_H_INCLUDED__ 
#define __MAXCLIQUE_H_INCLUDED__

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include <iostream>
#include <utility>

#include <boost/graph/adjacency_list.hpp>

class MaxClique { 

        public:
                typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS > Graph;

                virtual std::vector<int> max_clique(const Graph &, const size_t, const size_t) = 0;

};

#endif 
