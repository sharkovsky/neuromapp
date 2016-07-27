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

environment::distribution_from_file::distribution_from_file(const size_t me, const size_t glob_ncells, const std::string gid_file_dir_path) : global_number_(glob_ncells) {
    std::stringstream filename;
    filename << gid_file_dir_path << "/gid" << me << ".dat";
    std::fstream gid_file;
    gid_file.open( filename.str().c_str() , std::ios::in );

    size_t this_gid;
    while ( gid_file >> this_gid ) {
	local_gids_.push_back(this_gid);
    }

    gid_file.close();

}


