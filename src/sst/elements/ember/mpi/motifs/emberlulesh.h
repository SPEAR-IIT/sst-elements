// Ember LULESH
// Roth, Philip C. et.al "Automated characterization of parallel application communication patterns." HPDC'15
// Developled based on Halo3d26 and sweep3d
// yao kang ykang17@hawk.iit.edu

#ifndef _H_EMBER_LULESH_MOTIF
#define _H_EMBER_LULESH_MOTIF

#include "mpi/embermpigen.h"

#include <sys/stat.h>
#include <fstream>

namespace SST {
namespace Ember {

class EmberLuleshGenerator : public EmberMessagePassingGenerator {

public:
    SST_ELI_REGISTER_SUBCOMPONENT_DERIVED(
        EmberLuleshGenerator,
        "ember",
        "LuleshMotif",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "LULESH benchmark",
        SST::Ember::EmberGenerator
    )

    SST_ELI_DOCUMENT_PARAMS(
        {   "arg.iterations",       "Sets the number of bcast operations to perform",   "1"},
        {   "arg.compute",      "Sets the time spent computing",        "1"},

        {   "arg.msgscale",      "total msg quantity scale (0 --1]",        "1"},

        // for 3d grid
        {   "arg.pex",          "Sets the processors in X-dimension (0=auto)",      "0"},
        {   "arg.pey",          "Sets the processors in Y-dimension (0=auto)",      "0"},
        {   "arg.pez",          "Sets the processors in Z-dimension (0=auto)",      "0"},

        // for 2d grid
        {   "arg.pep",          "Sets the processors in P-dimension (0=auto)",      "0"},
        {   "arg.peq",          "Sets the processors in Q-dimension (0=auto)",      "0"},  
        {   "arg.kba",          "Sweep3d number of slices per Z-axis",      "0"},       

        {   "arg.wait2start",      "wait how long to start?",        "1"}, 

        {   "arg.dobcast",      "dobcast?",        "1"}, 
        {   "arg.doreduce",     "doreduce?",       "1"}, 
        {   "arg.do3dnn26",     "do3dnn26?",       "1"}, 
        {   "arg.do3dsweep",    "do3dsweep?",      "1"}, 
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
	EmberLuleshGenerator(SST::ComponentId_t, Params& params);
    bool generate( std::queue<EmberEvent*>& evQ);

private:

    enum comm_phase {
        BCAST = 0,   
        REDUCE,     
        NN3D26,
        SWEEP3D,     
        DONE,
    };

    int curr_phase;
    double msgscale;

    void config3d26(Params& params);
    void issue3D26stencil( std::queue<EmberEvent*>& evQ);
    void configsweep();
    void issue3dsweep( std::queue<EmberEvent*>& evQ);

    bool dobcast;
    bool doreduce;
    bool do3dnn26;
    bool do3dsweep;

    std::string motifname;
    uint64_t m_wait2start;
    std::string outfile; //IO file 

    uint32_t bcast_vol;
    uint32_t reduce_vol;
    uint32_t nn3d26_vol;
    uint32_t sweep3d_vol;

    uint64_t m_startTime;
    uint64_t m_stopTime;
    uint64_t m_compute;
    uint32_t m_iterations;
    uint32_t m_count;
    void*    m_sendBuf;
    void*    m_recvBuf;
    int      m_root;
    uint32_t m_loopIndex;

    uint64_t m_bcastdoneTime;
    uint64_t m_reducedoneTime;
    uint64_t m_nn3d26doneTime;


    //halo3d26
    std::vector<MessageRequest> requests;
    
    uint32_t peX;
	uint32_t peY;
	uint32_t peZ;

	int32_t  xface_down;
	int32_t  xface_up;

	int32_t  yface_down;
	int32_t  yface_up;

	int32_t  zface_down;
	int32_t  zface_up;

	int32_t  line_a;
	int32_t  line_b;
	int32_t  line_c;
	int32_t  line_d;
	int32_t  line_e;
	int32_t  line_f;
	int32_t  line_g;
	int32_t  line_h;
	int32_t  line_i;
	int32_t  line_j;
	int32_t  line_k;
	int32_t  line_l;

	int32_t  corner_a;
	int32_t  corner_b;
	int32_t  corner_c;
	int32_t  corner_d;
	int32_t  corner_e;
	int32_t  corner_f;
	int32_t  corner_g;
	int32_t  corner_h;

        //sweep3d
    uint32_t px;
	uint32_t py;
    uint32_t kba;

    int32_t  x_up;
	int32_t  x_down;
	int32_t  y_up;
	int32_t  y_down;


};

}
}

#endif
