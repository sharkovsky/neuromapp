/*
 * neurondistribution.cpp
 *
 *  Created on: Jul 11, 2016
 *      Author: schumann
 */
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "coreneuron_1.0/event_passing/environment/neurondistribution.h"

environment::continousdistribution::continousdistribution(size_t groups, size_t me, size_t cells):
        global_number(cells)
{
    //neuron distribution
    const int offset = cells % groups;
    const bool hasonemore = offset > me;
    local_number = cells / groups + hasonemore;

    start = me * local_number;
    if (!hasonemore)
        start += offset;
}

environment::continousdistribution::continousdistribution(size_t groups, size_t me, environment::continousdistribution* parent_distr):
    global_number(parent_distr->getglobalcells())
{
    //neuron distribution
    const int offset = parent_distr->getlocalcells() % groups;
    const bool hasonemore = offset > me;
    local_number = parent_distr->getlocalcells() / groups + hasonemore;

    start = me * local_number + parent_distr->start;
    if (!hasonemore)
        start += offset;
}

environment::distribution_from_file::distribution_from_file(const size_t me, const std::string gid_file_dir_path) {
    std::stringstream filename;
    filename << gid_file_dir_path << "/gid" << me << ".dat";
    std::fstream gid_file;
    gid_file.open( filename.str().c_str() , std::ios::in );

    size_t this_gid;
    while ( gid_file >> this_gid ) {
	local_gids_.push_back(this_gid);
    }

    gid_file.close();

    size_t local_num_of_cells = local_gids_.size();
    MPI_Allreduce(&local_num_of_cells, &global_number_, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

}

bool environment::distribution_from_file::isLocal(size_t id) const {
    return ( std::find( local_gids_.begin(), local_gids_.end(), id ) != local_gids_.end() );
}

size_t environment::distribution_from_file::global2local(size_t glo) const {
    const_iterator it = std::find( local_gids_.begin(), local_gids_.end(), glo );
    assert( it != local_gids_.end()  );
    return ( std::distance( local_gids_.begin(), it) );
}


