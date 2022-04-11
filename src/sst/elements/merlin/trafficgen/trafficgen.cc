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
#include "trafficgen/trafficgen.h"
#include <unistd.h>
#include <climits>
#include <signal.h>

#include <sst/core/params.h>
#include <sst/core/simulation.h>
#include <sst/core/timeLord.h>

using namespace SST::Merlin;
using namespace SST::Interfaces;

#if ENABLE_FINISH_HACK
int TrafficGen::count = 0;
int TrafficGen::received = 0;
int TrafficGen::min_lat = 0xffffffff;
int TrafficGen::max_lat = 0;
int TrafficGen::mean_sum = 0;
#endif

TrafficGen::TrafficGen(ComponentId_t cid, Params& params) :
    Component(cid),
//    last_vc(0),
    packets_sent(0),
    packets_recd(0),
    done(false),
    packet_delay(0),
    packetDestGen(NULL),
    packetSizeGen(NULL),
    packetDelayGen(NULL),
    next_send_time(0)
{
    out.init(getName() + ": ", 0, 0, Output::STDOUT);

    job_id = params.find<int>("job_id",0);

    id = params.find<int>("id",-1);
    if ( id == -1 ) {
        out.fatal(CALL_INFO, -1, "id must be set!\n");
    }

    num_peers = params.find<int>("num_peers",-1);
    if ( num_peers == -1 ) {
        out.fatal(CALL_INFO, -1, "num_peers must be set!\n");
    }

    // num_vns = params.find_integer("num_vns");
    // if ( num_vns == -1 ) {
    //     out.fatal(CALL_INFO, -1, "num_vns must be set!\n");
    // }
    num_vns = 1;

    addressMode = SEQUENTIAL;

    // Create a LinkControl object
    // First see if it is defined in the python
    link_control = loadUserSubComponent<SST::Interfaces::SimpleNetwork>
        ("networkIF", ComponentInfo::SHARE_NONE, num_vns );

    if ( !link_control ) {

        // Just load the default
        Params if_params;

        if_params.insert("link_bw",params.find<std::string>("link_bw"));
        if_params.insert("input_buf_size",params.find<std::string>("buffer_length","1kB"));
        if_params.insert("output_buf_size",params.find<std::string>("buffer_length","1kB"));
        if_params.insert("port_name","rtr");

        link_control = loadAnonymousSubComponent<SST::Interfaces::SimpleNetwork>
            ("merlin.linkcontrol", "networkIF", 0,
             ComponentInfo::SHARE_PORTS | ComponentInfo::INSERT_STATS, if_params, num_vns /* vns */);
    }


    packets_to_send = params.find<uint64_t>("packets_to_send", 1000);

    /* Distribution selection */
    packetDestGen = buildGenerator("PacketDest", params);
    assert(packetDestGen);
    packetDestGen->seed(id);

    /* Packet size */
    // base_packet_size = params.find_integer("packet_size", 64); // In Bits
    packetSizeGen = buildGenerator("PacketSize", params);
    if ( packetSizeGen ) {
        printf("EP trafficgen, packetSizeGen enabled\n");
        packetSizeGen->seed(id);
    } 

    std::string packet_size_s = params.find<std::string>("packet_size", "8B");
    UnitAlgebra packet_size(packet_size_s);
    if ( packet_size.hasUnits("B") ) {
        packet_size *= UnitAlgebra("8b/B");
    }

    if ( !packet_size.hasUnits("b") ) {
        out.fatal(CALL_INFO, -1, "packet_size must be specified in units of either B or b!\n");
    }

    base_packet_size = packet_size.getRoundedValue();


    // base_packet_delay = params.find_integer("delay_between_packets", 0);
    packetDelayGen = buildGenerator("PacketDelay", params);
    if ( packetDelayGen ) {
        packetDelayGen->seed(id);
    } 

    std::string packet_delay_s = params.find<std::string>("delay_between_packets", "0s");

    UnitAlgebra packet_delay(packet_delay_s);

    if ( !packet_delay.hasUnits("s") ) {
        out.fatal(CALL_INFO, -1, "packet_delay must be specified in units of s!\n");
    }

    // assume 1Ghz message rate
    base_packet_delay = (packet_delay / UnitAlgebra("1ns")).getRoundedValue();

    registerAsPrimaryComponent();
    primaryComponentDoNotEndSim();
    clock_functor = new Clock::Handler<TrafficGen>(this,&TrafficGen::clock_handler);
    clock_tc = registerClock( params.find<std::string>("message_rate", "1GHz"), clock_functor, false);

    std::string tmp_s = params.find<std::string>("message_rate", "1GHz");
    if(tmp_s.compare("1GHz")) out.fatal(CALL_INFO, -1, "delay_between_packets use ns as unit, here message rate has to be 1GHz");


    // Register a receive handler which will simply strip the events as they arrive
    link_control->setNotifyOnReceive(new SST::Interfaces::SimpleNetwork::Handler<TrafficGen>(this,&TrafficGen::handle_receives));
    send_notify_functor = new SST::Interfaces::SimpleNetwork::Handler<TrafficGen>(this,&TrafficGen::send_notify);

    traffic_pattern = params.find<std::string>("PacketDest:pattern");
    if(id == 0){

        printf("\n------- EP trafficgen---------\n");
        printf("  traffic pattern\t%s\n", traffic_pattern.c_str());

        if (!traffic_pattern.compare("Tornado")){
            printf("  \tADV+%d\n", tor_group_shift);
        } 

        else if (!traffic_pattern.compare("Stencil_3D")) {
            std::string shape = params.find<std::string>("PacketDest:Stencil_3D:3DSize");
            int maxX, maxY, maxZ;
            assert (sscanf(shape.c_str(), "%d %d %d", &maxX, &maxY, &maxZ) == 3);
            printf("  \tCube size %dx%dx%d\n", maxX, maxY, maxZ);
        }

        else if (!traffic_pattern.compare("FFT3D_all2all")) {
            std::string shape = params.find<std::string>("PacketDest:FFT3D_all2all:3DSize");
            int maxX, maxY, maxZ;
            assert (sscanf(shape.c_str(), "%d %d %d", &maxX, &maxY, &maxZ) == 3);
            printf("  \tCube size %dx%dx%d\n", maxX, maxY, maxZ);
        }

        else if (!traffic_pattern.compare("RandomNeighbors")) {
            int rn_min = params.find<int>("PacketDest:RandomNeighbors:range_min");
            int rn_max = params.find<int>("PacketDest:RandomNeighbors:range_max");
                printf("  \tNeighbors number range [%d, %d)\n", rn_min, rn_max);
        }

        printf("  job id(size)\t%d(%d)\n", job_id, num_peers);
       
        printf("  num packets\t%lu\n", packets_to_send); 
        printf("  base pkt size\t%d bit\n", base_packet_size);

        printf("  base delay\t%d\n", base_packet_delay);
        printf("  pkt size gen");
        if(packetSizeGen) printf("\tYes\n");
        else printf("\tNo\n");
        printf("  pkt delay gen");
        if(packetDelayGen){
            printf("\tEnabled\n");
            packetDelayGen->dump_status();
        } 
        else printf("\tNo\n");
        printf("\n------------------------------\n\n");
    }
    // printf("EP constructor, cid %lu, cur_t %lu\n", cid, getCurrentSimTime(clock_tc));
    
}


TrafficGen::~TrafficGen()
{
    delete link_control;
}


TrafficGen::Generator* TrafficGen::buildGenerator(const std::string &prefix, Params &params)
{
    

    Generator* gen = NULL;
    std::string pattern = params.find<std::string>(prefix + ":pattern");
    std::pair<int, int> range = std::make_pair(
        params.find<int>(prefix + ":RangeMin", 0),
        params.find<int>(prefix + ":RangeMax", INT_MAX));

    uint32_t rng_seed = params.find<uint32_t>(prefix + ":Seed", 1010101);

    if ( !pattern.compare("NearestNeighbor") ) {
        std::string shape = params.find<std::string>(prefix + ":NearestNeighbor:3DSize");
        int maxX, maxY, maxZ;
        assert (sscanf(shape.c_str(), "%d %d %d", &maxX, &maxY, &maxZ) == 3);
        gen = new NearestNeighbor(new UniformDist(0, 6), id, maxX, maxY, maxZ, 6);

    } else if (!pattern.compare("Stencil_3D")) {
        std::string shape = params.find<std::string>(prefix + ":Stencil_3D:3DSize");
        int maxX, maxY, maxZ;
        assert (sscanf(shape.c_str(), "%d %d %d", &maxX, &maxY, &maxZ) == 3);
        gen = new Stencil_3D( id, maxX, maxY, maxZ);
    
    } else if (!pattern.compare("FFT3D_all2all")) {
        std::string shape = params.find<std::string>(prefix + ":FFT3D_all2all:3DSize");
        int maxX, maxY, maxZ;
        assert (sscanf(shape.c_str(), "%d %d %d", &maxX, &maxY, &maxZ) == 3);
        gen = new FFT3D_all2all( id, maxX, maxY, maxZ);
    
    } else if (!pattern.compare("RandomNeighbors")) {
        int min_nodes = params.find<int>(prefix + ":RandomNeighbors:range_min");
        int max_nodes = params.find<int>(prefix + ":RandomNeighbors:range_max");

        gen = new RandomNeighbors( id, min_nodes, max_nodes, num_peers);

    // Tornado pattern, only for dragonfly
    } else if (!pattern.compare("Tornado")) {
        std::string tmp_topo = params.find<std::string>("topology");
        if (tmp_topo.find("dragonfly") == std::string::npos) {
            out.fatal(CALL_INFO, -1, "Tornado pattern only designed for dragonfly topology for now, not for '%s'\n", tmp_topo.c_str());
        }
        bool found1, found2, found3;
        int hpr = params.find<int>("dragonfly:hosts_per_router",found1);
        if(!found1){
            out.fatal(CALL_INFO, -1, "Missing values for dragonfly:hosts_per_router");
        }
        
        int rpg = params.find<int>("dragonfly:routers_per_group",found2);
        if(!found2){
            out.fatal(CALL_INFO, -1, "Missing values for dragonfly:routers_per_group");
        }
        int num_g = params.find<int>("dragonfly:num_groups",found3);
        if(!found3){
            out.fatal(CALL_INFO, -1, "Missing values for dragonfly:num_groups");
        }
        if(!found1 || !found2 || !found3){
            out.fatal(CALL_INFO, -1, "Missing values");
        }
        tor_group_shift = params.find<int>("Tornado:shift", 1, found3);
        gen = new Tornado(id, hpr, rpg, num_g, tor_group_shift);

    } else if (!pattern.compare("TornadoURmix")){
        std::string tmp_topo = params.find<std::string>("topology");
        if (tmp_topo.find("dragonfly") == std::string::npos) {
            out.fatal(CALL_INFO, -1, "Tornado pattern only designed for dragonfly topology for now, not for '%s'\n", tmp_topo.c_str());
        }

        bool found;
        double th = params.find<double>("tornadoUrMixThreshold",found);
        assert(found);

        bool found1, found2, found3;
        int hpr = params.find<int>("dragonfly:hosts_per_router",found1);
        int rpg = params.find<int>("dragonfly:routers_per_group",found2);
        int num_g = params.find<int>("dragonfly:num_groups",found3);
        if(!found1 || !found2 || !found3){
            out.fatal(CALL_INFO, -1, "Missing values");
        }

        gen = new TornadoURmix(id, hpr, rpg, num_g, range.first, range.second, th);

    
    } else if ( !pattern.compare("Step") ) {

        std::vector<std::string> delay_step;
        params.find_array<std::string>(prefix + ":packet_delay_list", delay_step);

        if(delay_step.size() == 0){
            out.fatal(CALL_INFO, -1, "packet_delay_list not found for Step pattern");
        }

        std::vector<int> delay_step_int;

        for(std::string delay_str : delay_step){
            UnitAlgebra packet_delay(delay_str);
            if ( !packet_delay.hasUnits("s") ) {
                out.fatal(CALL_INFO, -1, "packet_delay must be specified in units of s!\n");
            }
            // assume 1Ghz message rate
            int delay_int = (packet_delay / UnitAlgebra("1ns")).getRoundedValue();
            delay_step_int.push_back(delay_int);
        }

        std::vector<int> num_step;
        params.find_array<int>(prefix + ":packet_num_list", num_step);

        if(num_step.size() == 0){
            out.fatal(CALL_INFO, -1, "packet_num_list not found for Step pattern");
        }
        gen = new Step(id, delay_step_int, num_step);  

        packets_to_send = 0;
        for(int num : num_step ){
            packets_to_send += num;
        }

    } else if ( !pattern.compare("Uniform") ) {
        gen = new UniformDist(range.first, range.second);
    } else if ( !pattern.compare("HotSpot") ) {
        int target = params.find<int>(prefix + ":HotSpot:target");
        float targetProb = params.find<float>(prefix + ":HotSpot:targetProbability");
        gen = new DiscreteDist(range.first, range.second, target, targetProb);
    } else if ( !pattern.compare("Normal") ) {
        float mean = params.find<float>(prefix + ":Normal:Mean", range.second/2.0f);
        float sigma = params.find<float>(prefix + ":Normal:Sigma", 1.0f);
        gen = new NormalDist(range.first, range.second, mean, sigma);
    } else if ( !pattern.compare("Exponential") ) {
        float lambda = params.find<float>(prefix + ":Exponential:Lambda", range.first);
        gen = new ExponentialDist(lambda);
    } else if ( !pattern.compare("Binomial") ) {
        int trials = params.find<int>(prefix + ":Binomial:Mean", range.second);
        float probability = params.find<float>(prefix + ":Binomial:Sigma", 0.5f);
        gen = new BinomialDist(range.first, range.second, trials, probability);
    } else if ( pattern.compare("") ) { // Allow none - non-pattern
        out.fatal(CALL_INFO, -1, "Unknown pattern '%s'\n", pattern.c_str());
    }

    if ( gen ) gen->seed(rng_seed);

    return gen;
}

void TrafficGen::finish()
{
    link_control->finish();
}

void TrafficGen::setup()
{

    link_control->setup();
#if ENABLE_FINISH_HACK
    count++;
#endif
}

void
TrafficGen::init(unsigned int phase) {
    link_control->init(phase);
    return;
}


bool
TrafficGen::clock_handler(Cycle_t cycle)
{   
    if ( done ) return true;
    else if (packets_sent >= packets_to_send) {
        primaryComponentOKToEndSim();
        done = true;
    }

    uint64_t curr_time = getCurrentSimTimeNano();
    if(curr_time >= next_send_time){
    
        // packets_sent
        // 1-to-1 pattern:     this is the actual number of packets
        // 1-to-many pattern:  this is the round of a iteration, 
        //                     that multiple packets are sent
        if ( packets_sent < packets_to_send ) {

            int num_targets = packetDestGen->get_num_targets();
            int curr_target_idx = packetDestGen->current_target_idx;

            int packet_size = getPacketSize();

            //print progress
            if( (curr_target_idx == 0) && 
            (packets_to_send/5 !=0) &&  id == 0 && ( packets_sent %  (packets_to_send/5) == 0)) {
                printf("%lu/%lu message sent\n", packets_sent, packets_to_send);
            }
                
            for(; curr_target_idx<num_targets; curr_target_idx++) {

                if ( link_control->spaceToSend(0,packet_size) ) {

                    int target = getPacketDest();

                    if(target < 0 || target >= num_peers){
                        printf("EP %d, job %d, cal dest for %d\n", id, job_id, target);
                        assert ( target >= 0 && target < num_peers);
                    }

                    SimpleNetwork::Request* req = new SimpleNetwork::Request();
                    req->head = true;
                    req->tail = true;

                    switch ( addressMode ) {
                    case SEQUENTIAL:
                        req->dest = target;
                        req->src = id;
                        break;
                    case FATTREE_IP:
                        req->dest = fattree_ID_to_IP(target);
                        req->src = fattree_ID_to_IP(id);
                        break;
                    }
                    req->vn = 0;
                    // ev->size_in_flits = packet_size;
                    req->size_in_bits = packet_size;

                    req->setTraceID(packets_sent);
                    req->special_index = (3*packets_to_send)/5;
                    
                    bool sent = link_control->send(req,0);
                    assert( sent );
                }
                else {      
                    link_control->setNotifyOnSend(send_notify_functor);
                    return true;
                }
            }

            // This round of sending packets is finished. This round may send 1 packet or multiple packets
            ++packets_sent;
            next_send_time = getCurrentSimTimeNano() + (uint64_t) getDelayNextPacket();
            packetDestGen->current_target_idx = 0;

        }
    }
    return false;
}

int TrafficGen::fattree_ID_to_IP(int id)
{
    union Addr {
        uint8_t x[4];
        int32_t s;
    };

    Addr addr;

    int edge_switch = (id / ft_loading);
    int pod = edge_switch / (ft_radix/2);
    int subnet = edge_switch % (ft_radix/2);

    addr.x[0] = 10;
    addr.x[1] = pod;
    addr.x[2] = subnet;
    addr.x[3] = 2 + (id % ft_loading);

#if 0
    out.output("Converted NIC id %d to %u.%u.%u.%u.\n", id, addr.x[0], addr.x[1], addr.x[2], addr.x[3]\n);
#endif

    return addr.s;
}


int TrafficGen::IP_to_fattree_ID(int ip)
{
    union Addr {
        uint8_t x[4];
        int32_t s;
    };

    Addr addr;
    addr.s = ip;

    int id = 0;
    id += addr.x[1] * (ft_radix/2) * ft_loading;
    id += addr.x[2] * ft_loading;
    id += addr.x[3] -2;

    return id;
}

bool
TrafficGen::handle_receives(int vn)
{   
    int tmp_msg_id = -1;

    SimpleNetwork::Request* req = link_control->recv(vn);
    if ( req != NULL ) {
        packets_recd++;
        tmp_msg_id = req->getTraceID();
        delete req;
    }
    return true;
}


bool
TrafficGen::send_notify(int vn)
{   

    reregisterClock(clock_tc, clock_functor);
    return false;
}


int TrafficGen::getPacketDest(void)
{   
    int dest = packetDestGen->getNextValue();
    return dest;
}


int TrafficGen::getPacketSize(void)
{
    if ( packetSizeGen ) {
        return packetSizeGen->getNextValue();
    } else {
        return base_packet_size;
    }
}


int TrafficGen::getDelayNextPacket(void)
{
    if ( packetDelayGen ) {
        return packetDelayGen->getNextValue();
    } else {
        return base_packet_delay;
    }
}
