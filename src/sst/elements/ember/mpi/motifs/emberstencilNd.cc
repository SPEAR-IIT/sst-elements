// Copyright 2009-2021 NTESS. Under the terms
// of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.
//
// Copyright (c) 2009-2021, NTESS
// All rights reserved.
//
// Portions are copyright of other developers:
// See the file CONTRIBUTORS.TXT in the top level directory
// the distribution for more information.
//
// This file is part of the SST software package. For license
// information, see the LICENSE file in the top level directory of the
// distribution.

#include <sst_config.h>
#include "emberstencilNd.h"
using namespace SST::Ember;

EmberStencilNdGenerator::EmberStencilNdGenerator(SST::ComponentId_t id, Params& params) :
	EmberMessagePassingGenerator(id, params, "StencilNd"),
	m_loopIndex(0),
	m_startTime(0),
	m_stopTime(0)
{
	params.find_array<int>("arg.gridsize", gridsize);
	stencildim = gridsize.size();

	motifname = "emberStencil" + std::to_string(stencildim) + "d";

	mycoord.resize(stencildim);
	posneighbors.resize(stencildim);
	negneighbors.resize(stencildim);

	items_per_cell = params.find<uint32_t>("arg.fields_per_cell", 1);
	performReduction = (params.find<int>("arg.doreduce", 1) == 1);
	sizeof_cell = params.find<uint32_t>("arg.datatype_width", 8);

	nsCompute  = params.find<uint64_t>("arg.computetime", (uint64_t) 0);

	iterations = params.find<uint32_t>("arg.iterations", 1);

    jobId = params.find<int>("_jobId"); //NetworkSim

	m_wait2start = params.find<uint64_t>("arg.wait2start", 1);

	//get my cooridnate and calculate neighbors
	configure(params);

	mkdir("./ember_stats/" ,0755);
	outfile = "./ember_stats/" + motifname + "_rank" + std::to_string(rank()) + "_" + std::to_string(size()) + ".stats";
	FILE *file = fopen(outfile.c_str(), "w");
	fprintf(file, "rank,ite,start(ns),stop(ns),comm(ns)\n");
	fclose(file);
}

void EmberStencilNdGenerator::configure(Params& params)
{
	int worldSize = 1;
	for (int axlen : gridsize){
		worldSize *= axlen;
	}
	assert( worldSize == size() );

	lex_coords(mycoord, gridsize, rank() );

	//determine rank of neighbors
    std::vector<int> tmp_coords(stencildim);
    for (int d = 0; d < stencildim; d++){
        neighbor_coords(mycoord, tmp_coords, d, 1);
        posneighbors[d] = lex_rank(tmp_coords, gridsize);
		neighbor_coords(mycoord, tmp_coords, d, -1);
		negneighbors[d] = lex_rank(tmp_coords, gridsize);
    }

    if(0 == rank()) {
		output("\n");
		output("%s processor grid: ", motifname.c_str());
		for (int axlen : gridsize){
			output("%dx", axlen);
		}
		output("\n");
		output("%s compute time: %" PRIu32 " ns\n", motifname.c_str(), nsCompute);
		output("%s iterations: %" PRIu32 "\n", motifname.c_str(), iterations);
		output("%s iterms/cell: %" PRIu32 "\n", motifname.c_str(), items_per_cell);
		//
		output("%s size of cell: %" PRIu32 "\n", motifname.c_str(), sizeof_cell);
		output("%s do reduction: %" PRIu32 "\n", motifname.c_str(), performReduction);
		//
		output("%s wait2start: %lu ns\n", motifname.c_str(), m_wait2start);	

		output("%s rank 0, neighbors:", motifname.c_str() );
		for (int nei : posneighbors){
			output(" %d ", nei);
		}
		for (int nei : negneighbors){
			output(" %d ", nei);
		}
		output("\n");
	}
	// int newrank = lex_rank( mycoord, stencildim, gridsize);
	// output("%s, rank %d, newrank %d, coor: ", motifname.c_str(), rank(), newrank);
	// for (int idx : mycoord){
	// 	output("%d ", idx);
	// }
	// output("\n");
}

/*------------------------------------------------------------------*/
/* Convert rank to machine coordinates */
void EmberStencilNdGenerator::lex_coords( std::vector<int>& coords, std::vector<int>& gridsize, const uint32_t rank){
  int d;
  uint32_t r = rank;

  for(d = 0; d < stencildim; d++){
    //
    if(gridsize[d] == 0) printf("rank %d: gridsize[%d] = 0\n", r, d);
    
    coords[d] = r % gridsize[d];
    r /= gridsize[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert machine coordinate to linear rank (inverse of
   lex_coords) */
//added two checks to return -1 if we are trying to access a neighbor that
//is beyond the scope of the lattice.

int EmberStencilNdGenerator::lex_rank(const std::vector<int>& coords, const std::vector<int>& gridsize)
{
    int rank = coords[stencildim-1];
	for(int axdim=0; axdim < stencildim; axdim++){
		if ( coords[axdim] == -1 || coords[axdim]>=gridsize[axdim]){
			return -1;
		}
	}
    for(int d = stencildim-2; d >= 0; d--){
        rank = rank * gridsize[d] + coords[d];
    }
	assert( rank < size());
    return rank;
}

void EmberStencilNdGenerator::neighbor_coords(const std::vector<int>& coords, std::vector<int>& n_coords, const int dim, const int offset){
	assert(dim < stencildim);
	for(int d=0; d<stencildim; d++){
		n_coords[d]=coords[d];
	}
    n_coords[dim] = n_coords[dim]+offset;
}

bool EmberStencilNdGenerator::generate( std::queue<EmberEvent*>& evQ )
{	
	uint64_t comm_time = m_stopTime - m_startTime - nsCompute;
	if(rank()==0){
	// progress msg
		output("%s rank%d, ite %d, start %lu, end %lu, comm.time %lu @ %lu\n", motifname.c_str(), rank(), m_loopIndex, m_startTime, m_stopTime, comm_time, getCurrentSimTimeNano() );
	}

	// IO
	if(m_loopIndex>0){
		std::ofstream iofile;
		iofile.open(outfile.c_str(), std::ios::app);
		if (iofile.is_open())
		{	
			iofile <<rank()<<","<<m_loopIndex<<","<<m_startTime<<","<<m_stopTime<<","<<comm_time<<"\n";
			iofile.close();
		}
		else assert(0);
	}

	//change place of exit assertion so that last iteration is record
	if ( m_loopIndex == iterations ) {
		//
		if(rank()==0){
			output("%s rank 0 finished %d iterations\n", motifname.c_str(), iterations);
		}
		return true;
	}

	// wait 2 start
	if ( 0 == m_loopIndex ){
		enQ_compute( evQ, m_wait2start );
	}

	// ite begins
	enQ_getTime( evQ, &m_startTime );
	enQ_compute( evQ, nsCompute);

	std::vector<MessageRequest*> requests;

	// positive loop
	for (int d = 0; d < stencildim; d++){
		if (posneighbors[d] != -1){
			MessageRequest*  req  = new MessageRequest();
			requests.push_back(req);
			enQ_irecv( evQ, posneighbors[d],
					items_per_cell * sizeof_cell, 0, GroupWorld, req);
		}
	}
	for (int d = 0; d < stencildim; d++){
		if ( posneighbors[d] != -1 ){
			enQ_send( evQ, posneighbors[d], items_per_cell * sizeof_cell, 0, GroupWorld);
		}
	}

	// negative loop
	for (int d = 0; d < stencildim; d++){
		if (negneighbors[d] != -1){
			MessageRequest*  req  = new MessageRequest();
			requests.push_back(req);
			enQ_irecv( evQ, negneighbors[d],
					items_per_cell * sizeof_cell, 0, GroupWorld, req);
		}
	}
	for (int d = 0; d < stencildim; d++){
		if ( negneighbors[d] != -1 ){
			enQ_send( evQ, negneighbors[d], items_per_cell * sizeof_cell, 0, GroupWorld);
		}
	}

	for(uint32_t i = 0; i < requests.size(); ++i) {
		enQ_wait( evQ, requests[i]);
	}

	requests.clear();

	if(performReduction) {
		enQ_allreduce( evQ, NULL, NULL, 1, DOUBLE, MP::SUM, GroupWorld);
	}

	// comm. ends
	enQ_getTime( evQ, &m_stopTime );

	// change exit place to begine
	m_loopIndex++;
	return false;
}

