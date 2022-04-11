// -*- mode: c++ -*-

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


#ifndef COMPONENTS_MERLIN_GENERATORS_TRAFFICEGEN_H
#define COMPONENTS_MERLIN_GENERATORS_TRAFFICEGEN_H

#include <cstdlib>

#include <sst/core/rng/mersenne.h>
#include <sst/core/rng/gaussian.h>
#include <sst/core/rng/discrete.h>
#include <sst/core/rng/expon.h>
#include <sst/core/rng/uniform.h>

#include <sst/core/component.h>
#include <sst/core/event.h>
#include <sst/core/link.h>
#include <sst/core/timeConverter.h>
#include <sst/core/output.h>
#include <sst/core/interfaces/simpleNetwork.h>



#include "sst/elements/merlin/merlin.h"

#define ENABLE_FINISH_HACK 0

namespace SST {
namespace Merlin {

class TrafficGen : public Component {

public:

    SST_ELI_REGISTER_COMPONENT(
        TrafficGen,
        "merlin",
        "trafficgen",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "Pattern-based traffic generator.",
        COMPONENT_CATEGORY_NETWORK)

    SST_ELI_DOCUMENT_PARAMS(
        {"job_id",                                "job id if not the only job on the system."},

        {"id",                                    "Network ID of endpoint."},
        {"num_peers",                             "Total number of endpoints in network."},
        {"num_vns",                               "Number of requested virtual networks."},
        {"link_bw",                               "Bandwidth of the router link specified in either b/s or B/s (can include SI prefix)."},
        {"topology",                              "Name of the topology subcomponent that should be loaded to control routing."},
        {"buffer_length",                         "Length of input and output buffers.","1kB"},
        {"packets_to_send",                       "Number of packets to send in the test.","1000"},
        {"packet_size",                           "Packet size specified in either b or B (can include SI prefix).","5"},
        {"delay_between_packets",                 "","0"},
        {"message_rate",                          "","1GHz"},
        //Tornado, Stencil_3D, all2all-fft, random neighbor added
        {"PacketDest:pattern",                    "Address pattern to be used (NearestNeighbor, Uniform, HotSpot, Normal, Binomial)",NULL},
        {"PacketDest:Seed",                       "Sets the seed of the RNG", "11" },
        {"PacketDest:RangeMax",                   "Minumum address to send packets.","0"},
        {"PacketDest:RangeMin",                   "Maximum address to send packets.","INT_MAX"},
        {"PacketDest:NearestNeighbor:3DSize",     "For Nearest Neighbors, the 3D size \"x y z\" of the mesh", ""},
        {"PacketDest:HotSpot:target",             "For HotSpot, which node is the target", ""},
        {"PacketDest:HotSpot:targetProbability",  "For HotSpot, with what probability is the target targeted", ""},
        {"PacketDest:Normal:Mean",                "In a normal distribution, the mean", ""},
        {"PacketDest:Normal:Sigma",               "In a normal distribution, the mean variance", ""},
        {"PacketDest:Binomial:Mean",              "In a binomial distribution, the mean", ""},
        {"PacketDest:Binomial:Sigma",             "In a binomial distribution, the variance", ""},
        {"PacketSize:pattern",                    "Address pattern to be used (Uniform, HotSpot, Normal, Binomial)",NULL},
        {"PacketSize:Seed",                       "Sets the seed of the RNG", "11" },
        {"PacketSize:RangeMax",                   "Minumum size of packets.","0"},
        {"PacketSize:RangeMin",                   "Maximum size of packets.","INT_MAX"},
        {"PacketSize:HotSpot:target",             "For HotSpot, the target packet size", ""},
        {"PacketSize:HotSpot:targetProbability",  "For HotSpot, with what probability is the target targeted", ""},
        {"PacketSize:Normal:Mean",                "In a normal distribution, the mean", ""},
        {"PacketSize:Normal:Sigma",               "In a normal distribution, the mean variance", "1.0"},
        {"PacketSize:Binomial:Mean",              "In a binomial distribution, the mean", ""},
        {"PacketSize:Binomial:Sigma",             "In a binomial distribution, the variance", "0.5"},
        
        //add step delay pattern
        {"PacketDelay:pattern",                   "Address pattern to be used (Uniform, HotSpot, Normal, Binomial, Step)",NULL},
        {"PacketDelay:Seed",                      "Sets the seed of the RNG", "11" },
        {"PacketDelay:RangeMax",                  "Minumum delay between packets.","0"},
        {"PacketDelay:RangeMin",                  "Maximum delay between packets.","INT_MAX"},
        {"PacketDelay:HotSpot:target",            "For HotSpot, the target packet delay", ""},
        {"PacketDelay:HotSpot:targetProbability", "For HotSpot, with what probability is the target targeted", ""},
        {"PacketDelay:Normal:Mean",               "In a normal distribution, the mean", ""},
        {"PacketDelay:Normal:Sigma",              "In a normal distribution, the mean variance", "1.0"},
        {"PacketDelay:Binomial:Mean",             "In a binomial distribution, the mean", ""},
        {"PacketDelay:Binomial:Sigma",            "In a binomial distribution, the variance", "0.5"},
        
        {"tornadoUrMixThreshold",                  "mix threshold for tornado and ur mix", "0.5"},

        {"PacketDelay:packet_delay_list",   "array of packet delays", ""},
        {"PacketDelay:packet_num_list",   "array of packet nums", ""},

        {"Tornado:shift", "ADV + i pattern", "1"},
    )

    SST_ELI_DOCUMENT_PORTS(
        {"rtr",  "Port that hooks up to router.", { "merlin.RtrEvent", "merlin.credit_event" } }
    )

    SST_ELI_DOCUMENT_SUBCOMPONENT_SLOTS(
        {"networkIF", "Network interface", "SST::Interfaces::SimpleNetwork" }
    )


private:

#if ENABLE_FINISH_HACK
    static int count;
    static int received;
    static int min_lat;
    static int max_lat;
    static int mean_sum;
#endif

    class Generator {
    public:
        virtual int getNextValue(void) = 0;
        virtual void seed(uint32_t val) = 0;

        virtual void dump_status() {};
        virtual int get_num_targets(){ return 1; };

        Generator(int id) : id(id) { current_target_idx=0; };
        Generator(){ id=0; current_target_idx=0; };

        int id;
        int current_target_idx; //used for packetDest generator in case of one-to-many traffic pattern
    };

    //step pattern, only used for packet delay
    class Step : public Generator{
    public:
        Step(int id, std::vector<int> delay_list, std::vector<int> num_list) : Generator(id) 
        {   
            assert(delay_list.size() == num_list.size());

            msg_count = 0;
            delays = delay_list;
            num_msgs = num_list;

            int step_value = 0;
            for( int x :  num_msgs){
                step_value += x;
                num_steps.push_back(step_value);
            }
        }

        ~Step(){}

        int getNextValue(void)
        {   

            assert(msg_count < num_steps.back());
            int cur_idx = delays.size()-1;
            for(cur_idx; cur_idx >= 0; cur_idx--){
                if(msg_count >= num_steps[cur_idx])
                    break;
            }

            assert(cur_idx+1 < delays.size());

            int new_delay = delays[cur_idx+1];
            msg_count++;
            return new_delay;
        }

        void seed(uint32_t val) {}

        void dump_status(){

            printf("TrafficGen Endpoint: Step PacketDelay Generator\n");
            printf("Msg interval list(ns): ");
            for(auto x : delays ){
                printf("%d ", x );
            }
            printf("\nMsg num list: ");
            for(auto x : num_msgs ){
                printf("%d ", x );
            }
            printf("\nMsg num step: ");
            for(auto x : num_steps ){
                printf("%d ", x );
            }
            printf("\n");
        }
            
    private:
        int msg_count;
        std::vector<int> delays;
        std::vector<int> num_msgs;
        std::vector<int> num_steps;
    };


    //tornado pattern, only used for dragonfly topo
    class Tornado : public Generator{
    public:
        Tornado(int id, int hpr, int rpg, int num_g, int shift) : Generator(id) 
        {
            num_groups = num_g;
            group_size = hpr * rpg;
            num_peers = group_size * num_groups;
            g_shift = shift;

            dist_size = group_size - 0;
            gen = new MersenneRNG();
            dist = new SSTUniformDistribution(dist_size, gen);
        }

        ~Tornado(){
            delete dist;
            delete gen;
        }

        int getNextValue(void)
        {
            int group_id = id / group_size;
            int local_id = id % group_size;

            int dest_g = (group_id + g_shift ) % num_groups;
            int tmpid = ( (int) dist->getNextDouble() ) + 0 - 1;
            int dest_id = dest_g * group_size + tmpid;
            assert(dest_id < num_peers); 
            return dest_id;
        }

        void seed(uint32_t val) {
            delete dist;
            delete gen;
            gen = new MersenneRNG((unsigned int) val);
            dist = new SSTUniformDistribution(dist_size, gen);
        }

    private:
        int num_groups;
        int group_size;
        int num_peers;
        MersenneRNG* gen;
        SSTUniformDistribution* dist;
        int dist_size;
        int g_shift;
    };

    // mix of tornado and ur
    class TornadoURmix : public Generator{
    public:
        TornadoURmix(int id, int hpr, int rpg, int num_g, int min, int max, double threshold) : 
        Generator(id),
        min(min) 
        {
            // tornado part
            num_groups = num_g;
            group_size = hpr * rpg;
            num_peers = group_size * num_groups;

            // ur part
            gen = new MersenneRNG();
            dist_size = std::max(1, max-min);
            dist = new SSTUniformDistribution(dist_size, gen);

            mixthreshold = threshold;
        }

        ~TornadoURmix(){
            delete dist;
            delete gen;
            // delete rng;
        }

        int getNextValue(void)
        {
            int group_id = id / group_size;
            int local_id = id % group_size;
            int dest_g = (group_id+1)%num_groups;
            int dest_id_tornado = dest_g * group_size + local_id;
            assert(dest_id_tornado < num_peers); 
            int dest_id_ur = ((int) dist->getNextDouble()) +min-1;
            double prob = gen->nextUniform();

            int final_dest = -1;
            if(prob < mixthreshold){
                final_dest = dest_id_tornado;
            }
            else{
                final_dest = dest_id_ur;
            }
            assert(final_dest>=0 && final_dest< num_peers);
            return final_dest;
        }

        void seed(uint32_t val) {
            delete dist;
            delete gen;
            gen = new MersenneRNG((unsigned int) val);
            dist = new SSTUniformDistribution(dist_size,gen);
        }

    private:
        // for tornado
        int num_groups;
        int group_size;
        int num_peers;
        // for  ur
        MersenneRNG* gen;
        SSTUniformDistribution* dist;
        int dist_size;
        int min;

        double mixthreshold;  // percentage of tornado packets
    };
  
    //derived from NearestNeighbor
    class Stencil_3D : public Generator {
        int* neighbors;
        int numNeighbors = 6;
    public:
        Stencil_3D(int id, int maxX, int maxY, int maxZ) :
        Generator(id) 
        {  
            int myX = id % maxX;
            int myY = (id / maxX) % maxY;
            int myZ = id / (maxX*maxY);

            neighbors = new int[numNeighbors];
            neighbors[0] = (((myX-1) + maxX) % maxX) + myY*maxX                         + myZ*(maxX*maxY);
            neighbors[1] = ((myX+1) % maxX)          + myY*maxX                         + myZ*(maxX*maxY);
            neighbors[2] = myX                       + (((myY-1) + maxY) % maxY) * maxX + myZ*(maxX*maxY);
            neighbors[3] = myX                       + (((myY+1)) % maxY) * maxX        + myZ*(maxX*maxY);
            neighbors[4] = myX                       + myY*maxX                         + (((myZ-1) + maxZ) % maxZ) * (maxX*maxY);
            neighbors[5] = myX                       + myY*maxX                         + (((myZ+1)) % maxZ) * (maxX*maxY);
        }

        void seed(uint32_t val)
        {
            //pass;
        }

        int getNextValue(void)
        {   
            assert(current_target_idx < numNeighbors);
            int nbr = neighbors[current_target_idx];
            current_target_idx++;
            return nbr;
        }

        int get_num_targets(){
            return numNeighbors;
        }
    };
  
    //3d fft, nodes arrange in a 3d-grid, group along z-axis.
    //all-to-all within communicator group
    class FFT3D_all2all : public Generator {
        int myX, myY, myZ;
        int sizeX, sizeY, sizeZ;

    public:
        FFT3D_all2all(int id, int maxX, int maxY, int maxZ) :
        Generator(id),
        sizeX(maxX), sizeY(maxY), sizeZ(maxZ) 
        {  
            myX = id % maxX;
            myY = (id / maxX) % maxY;
            myZ = id / (maxX*maxY);
        }

        void seed(uint32_t val)
        {
            //pass;
        }

        int getNextValue(void)
        {   

            assert(current_target_idx < sizeZ);
            int next_Z = current_target_idx >= myZ? (current_target_idx+1):current_target_idx; 
            assert(next_Z < sizeZ);
            int nbr = myX + myY * sizeX + next_Z*(sizeX*sizeY);        
            current_target_idx++; 
            return nbr;
        }

        int get_num_targets(){
            return sizeZ-1;
        }
    };

    class RandomNeighbors: public Generator {
        MersenneRNG* gen_num_nbrs;
        SSTUniformDistribution* dist_num_nbrs;
        int dist_size_num_nbrs; 

        MersenneRNG* gen_target;
        SSTUniformDistribution* dist_target;

        int num_neighbors;
        int min;
        int job_size;
    
    public:
        RandomNeighbors(int id, int min, int max, int max_size) :
        Generator(id), min(min), job_size(max_size)
        {
    		gen_num_nbrs = new MersenneRNG();
    		dist_size_num_nbrs = std::max(1, max-min);
    		dist_num_nbrs = new SSTUniformDistribution(dist_size_num_nbrs, gen_num_nbrs);

            gen_target = new MersenneRNG();
    		dist_target = new SSTUniformDistribution(max_size, gen_target);
        }

    	~RandomNeighbors() {
    		delete dist_num_nbrs;
    		delete gen_num_nbrs;

            delete dist_target;
    		delete gen_target;
    	}

        int get_num_targets(){
            //begining of a new round of sending packets
            if (current_target_idx == 0){
                num_neighbors = (int) dist_num_nbrs->getNextDouble() + min - 1;
            }
            // else, this is resumed from the previous round due to nic buffer size overflow
            // thus number of neighbors unchanged
            return  num_neighbors;
        }

        int getNextValue(void)
        {      
            assert(num_neighbors > 0);
            assert(current_target_idx < num_neighbors);
            int nbr = (int) dist_target->getNextDouble() + 0 - 1;
            current_target_idx++;
            return nbr;
        }

        void seed(uint32_t val)
        {
            //pass
        }
    };

    class NearestNeighbor : public Generator {
        Generator *dist;
        int *neighbors;
    public:
        NearestNeighbor(Generator *dist, int id, int maxX, int maxY, int maxZ, int numNeighbors) :
            dist(dist)
        {
            int myX = id % maxX;
            int myY = (id / maxX) % maxY;
            int myZ = id / (maxX*maxY);

            neighbors = new int[numNeighbors];
            switch (numNeighbors) {
            case 6: {
                neighbors[0] = (((myX-1) + maxX) % maxX) + myY*maxX                         + myZ*(maxX*maxY);
                neighbors[1] = ((myX+1) % maxX)          + myY*maxX                         + myZ*(maxX*maxY);
                neighbors[2] = myX                       + (((myY-1) + maxY) % maxY) * maxX + myZ*(maxX*maxY);
                neighbors[3] = myX                       + (((myY+1)) % maxY) * maxX        + myZ*(maxX*maxY);
                neighbors[4] = myX                       + myY*maxX                         + (((myZ-1) + maxZ) % maxZ) * (maxX*maxY);
                neighbors[5] = myX                       + myY*maxX                         + (((myZ+1)) % maxZ) * (maxX*maxY);
                break;
            }
            default:
                Simulation::getSimulation()->getSimulationOutput().fatal(CALL_INFO, -1,
                    "Unsure how to deal with %d neighbors\n", numNeighbors);
            }
        }

        int getNextValue(void)
        {
            int neighbor = dist->getNextValue();
            return neighbors[neighbor];
        }

        void seed(uint32_t val)
        {
            dist->seed(val);
        }
    };

    class ExponentialDist : public Generator {
        MersenneRNG* gen;
        SSTExponentialDistribution* dist;

    public:
        ExponentialDist(int lambda)
        {
            gen = new MersenneRNG();
	    dist = new SSTExponentialDistribution((double) lambda);
        }

        ~ExponentialDist() {
            delete dist;
            delete gen;
        }

        int getNextValue(void)
        {
            return (int) dist->getNextDouble();
        }

        void seed(uint32_t val)
        {
            delete gen;
            gen = new MersenneRNG((unsigned int) val);
        }
    };


    class UniformDist : public Generator {
        MersenneRNG* gen;
        SSTUniformDistribution* dist;

        int dist_size;
        int min;

    public:
        UniformDist(int min, int max): min(min)
        {
    		gen = new MersenneRNG();

    		dist_size = std::max(1, max-min);
    		dist = new SSTUniformDistribution(dist_size, gen);
        }

    	~UniformDist() {
    		delete dist;
    		delete gen;
    	}

        int getNextValue(void)
        {
            return (int) dist->getNextDouble() + min - 1;
        }

        void seed(uint32_t val)
        {
            delete dist;
            delete gen;
            gen = new MersenneRNG((unsigned int) val);
            dist = new SSTUniformDistribution(dist_size,gen);
        }
    };


    class DiscreteDist : public Generator {
	MersenneRNG* gen;
	SSTDiscreteDistribution* dist;
        int minValue;
    public:
        DiscreteDist(int min, int max, int target, double targetProb) : minValue(min)
        {
            int size = std::max(max - min, 1);
            double dfltP = (1.0 - targetProb) / (size - 1);
            std::vector<double> probs(size);
            for ( int i = 0 ; i < size ; i++ ) {
                probs[i] = dfltP;
            }
            probs[target] = targetProb;

	    gen = new MersenneRNG();
	    dist = new SSTDiscreteDistribution(&probs[0], size, gen);
        }

	~DiscreteDist() {
		delete dist;
		delete gen;
	}

        int getNextValue(void)
        {
            return ((int) dist->getNextDouble()) + minValue;
        }

        void seed(uint32_t val)
        {
            gen = new MersenneRNG((unsigned int) val);
        }
    };

    class NormalDist : public Generator {
	SSTGaussianDistribution* dist;
	MersenneRNG* gen;

        int minValue;
        int maxValue;
    public:
        NormalDist(int min, int max, double mean, double stddev) : minValue(min), maxValue(max)
        {
            gen = new MersenneRNG();
            dist = new SSTGaussianDistribution(mean, stddev);
        }

	~NormalDist() {
		delete dist;
		delete gen;
	}

        int getNextValue(void) {
            double val = -1.0;
            while ((int)val >= maxValue || (int)val < minValue || val < 0){
                val = dist->getNextDouble();
            }
            return (int) val;
        }

        void seed(uint32_t val)
        {
            gen = new MersenneRNG((unsigned int) val);
        }
    };

    class BinomialDist : public Generator {
        int minValue;
    public:
        BinomialDist(int min, int max, int trials, float probability) : minValue(min)
        {
            merlin_abort.fatal(CALL_INFO, -1, "BinomialDist is not currently supported\n");
        }
        virtual int getNextValue(void)
        {
            // return dist(gen) + minValue;
            return 0;
        }
        virtual void seed(uint32_t val)
        {
            // gen.seed(val);
        }
    };



    enum AddressMode { SEQUENTIAL, FATTREE_IP };

    AddressMode addressMode;

    Output out;
    int id;
    int ft_loading;
    int ft_radix;
    int num_peers;
    int num_vns;
//    int last_vc;
    
    int job_id;
    std::string traffic_pattern;

    uint64_t packets_sent;
    uint64_t packets_recd;

    bool done;

    SST::Interfaces::SimpleNetwork* link_control;
    SST::Interfaces::SimpleNetwork::Handler<TrafficGen>* send_notify_functor;
    Clock::Handler<TrafficGen>* clock_functor;
    TimeConverter* clock_tc;

    int base_packet_size;
    uint64_t packets_to_send;

    int base_packet_delay;
    int packet_delay;

    Generator *packetDestGen;
    Generator *packetSizeGen;
    Generator *packetDelayGen;

    double torUrmixth;

    SimTime_t next_send_time;

    int tor_group_shift=1;

public:
    TrafficGen(ComponentId_t cid, Params& params);
    ~TrafficGen();

    void init(unsigned int phase);
    void setup();
    void finish();


private:
    Generator* buildGenerator(const std::string &prefix, Params& params);
    bool clock_handler(Cycle_t cycle);
    int fattree_ID_to_IP(int id);
    int IP_to_fattree_ID(int id);
    bool handle_receives(int vn);
    bool send_notify(int vn);

protected:
    int getPacketDest(void);
    int getPacketSize(void);
    int getDelayNextPacket(void);

};

} //namespace Merlin
} //namespace SST

#endif
