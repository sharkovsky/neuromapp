#include <boost/range/algorithm_ext/iota.hpp>
#include <iostream>
#include <time.h>
#include <ctime>
#include <numeric>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/range/algorithm_ext/iota.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>

#include "coreneuron_1.0/event_passing/environment/presyn_maker.h"

namespace environment {
namespace connectionRules {

void fixedindegree::connect(gid2in_type& inputs, gid2out_type& outputs, neurondistribution* neuron_dist, int rank) {
    //create local presyns with empty vectors
    for(int i = 0; i < neuron_dist->getlocalcells(); ++i){
        const int gid = neuron_dist->local2global(i);
        outputs[gid];
    }

    //used for random presyn and netcon selection
    boost::mt19937 rng(time(NULL) + rank);
    boost::random::uniform_int_distribution<> uni_d(0, neuron_dist->getglobalcells()-1);

    //foreach gid, select srcs
    for(int i = 0; i < neuron_dist->getlocalcells(); ++i){
	for(int j = 0; j < fanin_; ++j){
	    const int cur = uni_d(rng);
	    //local GID
	    if(neuron_dist->isLocal(cur)){
		//add self to src gid
		const int g_i = neuron_dist->local2global(i);
		outputs[cur].push_back(g_i);
	    }
	    //remote GID
	    else{
		//add self to input presyn for gid
		const int g_i = neuron_dist->local2global(i);
		inputs[cur].push_back(g_i);
	    }
	}
    }
}



void fixedoutdegree::connect(gid2in_type& inputs, gid2out_type& outputs, neurondistribution* neuron_dist, int rank) {
    //create local presyns with empty vectors
    for(int i = 0; i < neuron_dist->getlocalcells(); ++i){
        const int gid = neuron_dist->local2global(i);
        outputs[gid];
    }

    //used for random presyn and netcon selection
    //seed has to be fixed across the ranks to generate proper fixed out degree
    boost::mt19937 rng(15623);
    boost::random::uniform_int_distribution<> uni_d(0, neuron_dist->getglobalcells()-1);

    for(int cur = 0; cur < neuron_dist->getglobalcells(); ++cur){
	for(int j = 0; j < fanout_; ++j){
	    int picked = uni_d(rng);
	    //connect if target neuron is local
	    //only one rank stores the connection
	    if(neuron_dist->isLocal(picked)) {
		if(neuron_dist->isLocal(cur)){
		    //add self to src gid
		    outputs[cur].push_back(picked);
		}
		//remote GID
		else{
		    //add self to input presyn for gid
		    inputs[cur].push_back(picked);
		}
	    }
	}
    }
}


void fromfile::connect (gid2in_type& inputs, gid2out_type& outputs, neurondistribution* neuron_dist, int rank) {
    //create local presyns with empty vectors
    for(int i = 0; i < neuron_dist->getlocalcells(); ++i){
        const int gid = neuron_dist->local2global(i);
        outputs[gid];
    }
    std::stringstream fname;
    std::fstream gid2src_f;
    size_t target;

    for(int i = 0; i < neuron_dist->getlocalcells(); ++i){
	target = neuron_dist->local2global(i);
	fname << filepath_root_ << "/gid2src_gid" << target << ".dat";
	gid2src_f.open(fname.str().c_str(), std::ios::in);
	size_t src;
	while ( gid2src_f >> src ) {
	    if(neuron_dist->isLocal(src)){
		//add self to src gid
		//if ( rank == 0 ) { std::cerr << "loc target: " << target << " src: " << src << std::endl; }
		outputs[src].push_back(target);
	    }
	    //remote GID
	    else{
		//if ( rank == 0 ) { std::cerr << "rem target: " << target << " src: " << src << std::endl; }
		//add self to input presyn for gid
		inputs[src].push_back(target);
	    }
	}
	gid2src_f.close();
	fname.str(std::string());
    }
}

} // namespace connectionRules

const presyn* presyn_maker::find_input(int key) const{
    std::map<int, std::vector<int> >::const_iterator it = inputs_.begin();
    const presyn* input_ptr = NULL;
    it = inputs_.find(key);
    if(it != inputs_.end()){
        input_ptr = &(it->second);
    }
    return input_ptr;
}

const presyn* presyn_maker::find_output(int key) const{
    std::map<int, std::vector<int> >::const_iterator it = outputs_.begin();
    const presyn* output_ptr = NULL;
    it = outputs_.find(key);
    if(it != outputs_.end()){
        output_ptr = &(it->second);
    }
    return output_ptr;
}

} //end of namespace
