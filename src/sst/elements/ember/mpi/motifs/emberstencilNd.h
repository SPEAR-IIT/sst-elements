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

// yao
// updated from Halo3D Motif
#ifndef _H_EMBER_STENCIL_ND_BLOCKING
#define _H_EMBER_STENCIL_ND_BLOCKING

#include "mpi/embermpigen.h"

//yao
#include <sys/stat.h>
#include <fstream>

namespace SST {
namespace Ember {

class EmberStencilNdGenerator : public EmberMessagePassingGenerator {

public:
    SST_ELI_REGISTER_SUBCOMPONENT_DERIVED(
        EmberStencilNdGenerator,
        "ember",
        "StencilNdMotif",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "Stencil Nd motif",
        SST::Ember::EmberGenerator
    )

    SST_ELI_DOCUMENT_PARAMS(
        {   "arg.gridsize",  "performs a N d stencil pattern", "[1,1,1,1,1]"},
        {   "arg.fields_per_cell",  "Specify how many variables are being computed per cell (this is one of the dimensions in message size. Default is 1", "1"},
        {   "arg.doreduce",     "How often to do reduces, 1 = each iteration",      "1"},
        {   "arg.datatype_width",   "Specify the size of a single variable, single grid point, typically 8 for double, 4 for float, default is 8 (double). This scales message size to ensure byte count is correct.", "8"},
        {   "arg.computetime",      "Sets the number of nanoseconds to compute for",    "10"},
        {   "arg.iterations",       "Sets the number of halo3d operations to perform",  "10"},
        //yao
        {   "arg.wait2start",      "wait how long to start?",        "1"}, 
    )

    SST_ELI_DOCUMENT_STATISTICS(
        { "time-Init", "Time spent in Init event",          "ns",  0},
        { "time-Finalize", "Time spent in Finalize event",  "ns", 0},
        { "time-Rank", "Time spent in Rank event",          "ns", 0},
        { "time-Size", "Time spent in Size event",          "ns", 0},
        { "time-Send", "Time spent in Recv event",          "ns", 0},
        { "time-Recv", "Time spent in Recv event",          "ns", 0},
        { "time-Irecv", "Time spent in Irecv event",        "ns", 0},
        { "time-Isend", "Time spent in Isend event",        "ns", 0},
        { "time-Wait", "Time spent in Wait event",          "ns", 0},
        { "time-Waitall", "Time spent in Waitall event",    "ns", 0},
        { "time-Waitany", "Time spent in Waitany event",    "ns", 0},
        { "time-Compute", "Time spent in Compute event",    "ns", 0},
        { "time-Barrier", "Time spent in Barrier event",    "ns", 0},
        { "time-Alltoallv", "Time spent in Alltoallv event", "ns", 0},
        { "time-Alltoall", "Time spent in Alltoall event",  "ns", 0},
        { "time-Allreduce", "Time spent in Allreduce event", "ns", 0},
        { "time-Reduce", "Time spent in Reduce event",      "ns", 0},
        { "time-Bcast", "Time spent in Bcast event",        "ns", 0},
        { "time-Gettime", "Time spent in Gettime event",    "ns", 0},
        { "time-Commsplit", "Time spent in Commsplit event", "ns", 0},
        { "time-Commcreate", "Time spent in Commcreate event", "ns", 0},
    )


public:
	EmberStencilNdGenerator(SST::ComponentId_t, Params& params);
	~EmberStencilNdGenerator() {}
	//yao
    // void configure();
    void configure(Params& params);
	bool generate( std::queue<EmberEvent*>& evQ );

private:

    void lex_coords( std::vector<int>& coords, std::vector<int>& size, const uint32_t rank);
    int lex_rank(const std::vector<int>& coords, const std::vector<int>& size);
    void neighbor_coords(const std::vector<int>& coords, std::vector<int>& n_coords, const int dim, const int offset);

	uint32_t m_loopIndex;
	bool performReduction;
    uint32_t iterations;

    //yao
    int stencildim;
    std::vector<int> gridsize;
    std::vector<int> mycoord;
	uint32_t nsCompute;
	uint32_t items_per_cell;
	uint32_t sizeof_cell;
    std::vector<int> posneighbors;
    std::vector<int> negneighbors;

    int jobId; //NetworkSim
	uint64_t m_startTime;  //NetworkSim
    uint64_t m_stopTime;  //NetworkSim

    //yao
    std::string motifname;
    uint64_t m_wait2start;
    std::string outfile; //IO file 

};

}
}

#endif
