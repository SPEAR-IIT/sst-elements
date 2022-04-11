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
//

#include <sst_config.h>
#include <sst/core/shared/sharedArray.h>
#include "sst/core/rng/xorshift.h"

#include "merlin.h"
#include "dragonfly.h"

#include <stdlib.h>
#include <sstream>

using namespace SST::Merlin;

// const uint8_t bit_array::masks[8] = { 0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f };


void
RouteToGroup::init_write(const std::string& basename, int group_id, global_route_mode_t route_mode,
                         const dgnflyParams& params, const std::vector<int64_t>& global_link_map,
                         bool config_failed_links, const std::vector<FailedLink>& failed_links_vec)
{   
    if( (params.g-1) * params.n != global_link_map.size() ){
        printf("MerlinDF:  (params.g-1) * params.n : %d,  global_link_map.size() %ld\n", (params.g-1) * params.n, global_link_map.size());
    }
    assert( (params.g-1) * params.n == global_link_map.size() );

    // Get a shared region
    data.initialize(basename+"group_to_global_port",
                    ((params.g-1) * params.n) * sizeof(RouterPortPair));


    groups = params.g;
    routers = params.a;
    slices = params.n;
    links = params.h;
    gid = group_id;
    global_start = params.p + params.a - 1;
    mode = route_mode;
    consider_failed_links = config_failed_links;

    // Fill in the data
    for ( int i = 0; i < global_link_map.size(); i++ ) {
        // Figure out all the mappings
        int64_t value = global_link_map[i];
        if ( value == -1 ) continue;

        int group = value % (params.g - 1);
        int route_num = value / (params.g - 1);
        int router = i / params.h;
        int port = (i % params.h) + params.p + params.a - 1;

        RouterPortPair rpp;
        rpp.router = router;
        rpp.port = port;
        // sr->modifyArray(group * params.n + route_num, rpp);
        data.write(group * params.n + route_num, rpp);
    }
    data.publish();
    // sr->publish();
    // data = sr->getPtr<const RouterPortPair*>();

    // If we're not doing failed links, we're done
    if ( !config_failed_links ) return;

    // For now, we only support the same number of links between
    // each pair of groups, so we'll start with that many and
    // subtract for failed links.  We'll build it all in a private
    // array and then copy it becuase multiple modification within
    // the shared region is inefficient at best and breaks at
    // worst.
    link_counts.initialize(basename + "group_link_counts",groups * groups,
                           params.n,Shared::SharedObject::NO_VERIFY);

    // Need to have a bit field where each element tells us whether or
    // not that global link is active.  Start at group 0, global link
    // 0 and count from there.
    size_t bf_size = groups * params.a * params.h;
    failed_links.initialize(basename + "downed_links", bf_size);


    // uint8_t* link_counts_ld = new uint8_t[groups * groups];
    // std::fill_n(link_counts_ld, groups * groups, params.n);

    // uint8_t* bf_down = static_cast<uint8_t*>(sr_dl->getRawPtr());
    // uint8_t* bf_down = new uint8_t[bf_size];
    // std::fill_n(bf_down, bf_size, 0);

    for ( auto x : failed_links_vec ) {
        if ( x.low_group >= groups || x.high_group >= groups || x.slice >= params.n ) {
            merlin_abort.fatal(CALL_INFO,1,"Illegal link specification: %d:%d:%d\n",x.low_group, x.high_group, x.slice);
        }
        link_counts.write(x.low_group * groups + x.high_group, link_counts[x.low_group * groups + x.high_group] - 1);
        link_counts.write(x.high_group * groups + x.low_group, link_counts[x.high_group * groups + x.low_group] - 1);

        const RouterPortPair& rp = getRouterPortPairForGroup(x.low_group, x.high_group, x.slice);
        int bit_index = (x.low_group * params.a * params.h) + (rp.router * params.h + rp.port - global_start);
        failed_links.write(bit_index,true);

        const RouterPortPair& rp2 = getRouterPortPairForGroup(x.high_group, x.low_group, x.slice);
        bit_index = (x.high_group * params.a * params.h) + (rp2.router * params.h + rp2.port - global_start);
        failed_links.write(bit_index,true);
    }

    // Publish the link counts
    link_counts.publish();

    // Published the down links data.  Already have a pointer to the
    // data.
    failed_links.publish();
}

void
RouteToGroup::init(const std::string& basename, int group_id, global_route_mode_t route_mode,
                   const dgnflyParams& params, bool config_failed_links)
{
    data.initialize(basename+"group_to_global_port");
    data.publish();

    groups = params.g;
    routers = params.a;
    slices = params.n;
    links = params.h;
    gid = group_id;
    global_start = params.p + params.a - 1;
    mode = route_mode;
    consider_failed_links = config_failed_links;

    if ( !config_failed_links ) return;

    // Get the shared regions and publish them.
    link_counts.initialize(basename + "group_link_counts");
    link_counts.publish();

    size_t bf_size = (groups * params.a * params.h + 7) / 8;
    failed_links.initialize(basename + "downed_links");
    failed_links.publish();
}



const RouterPortPair&
RouteToGroup::getRouterPortPair(int group, int route_number) const
{
    return getRouterPortPairForGroup(gid,group,route_number);
}

const RouterPortPair&
RouteToGroup::getRouterPortPairForGroup(uint32_t src_group, uint32_t dest_group, uint32_t slice) const
{
    // Look up global port to use
    switch ( mode ) {
    case ABSOLUTE:
        if ( dest_group >= src_group ) dest_group--;
        break;
    case RELATIVE:
        if ( dest_group > src_group ) {
            dest_group = dest_group - src_group - 1;
        }
        else {
            dest_group = groups - src_group + dest_group - 1;
        }
        break;
    default:
        break;
    }

    return data[dest_group*slices + slice];
}

int
RouteToGroup::getValiantGroup(int dest_group, RNG::SSTRandom* rng) const
{
    if ( !consider_failed_links )  {
        int group;
        do {
            group = rng->generateNextUInt32() % groups;
        } while ( group == gid || group == dest_group );
        return group;
    }

    // Need to worry about failed links
    uint32_t possible_mid_groups = groups - 2;

    uint8_t min_links;
    uint8_t req_links;
    uint32_t mid;
    do {
        // We need to generate two random numbers: mid_group plus a
        // random numbers that will determin how many links are needed
        // between groups in order to be considered a valid route.
        // This will weight the valiant groups based on how many
        // downed links there are along the path.  We compare against
        // the minimum link count between src->mid and mid->dest.
        uint32_t state = rng->generateNextUInt32() % (possible_mid_groups * slices);
        mid = state /  slices;
        req_links = state % slices + 1;

        // Need to adjust mid.  We don't consider the src or dest
        // group.
        if ( mid >= gid || mid >= dest_group ) mid++;
        if ( mid >= gid && mid >= dest_group ) mid++;

        // Need to get the lowest number of links between src->mid and
        // mid->dest and compare against the req_links
        min_links = getLinkCount(gid,mid);
        uint8_t links = getLinkCount(mid,dest_group);
        if ( links < min_links ) min_links = links;
    } while ( min_links < req_links );
    // printf("src_group = %d, dest_group = %d, mid_group = %d\n",gid,dest_group,mid);
    return mid;
}


// void
// RouteToGroup::setRouterPortPair(int group, int route_number, const RouterPortPair& pair) {
//     region->modifyArray(group*routes+route_number,pair);
// }


topo_dragonfly::topo_dragonfly(ComponentId_t cid, Params &p, int num_ports, int rtr_id, int num_vns) :
    Topology(cid),
    num_vns(num_vns),
    rtr_id(rtr_id),
    router_id_global(rtr_id), // router_id_global == rtr_id
    useQrouting(0)
{

    qtablefile = qtableFileDir+"rtr"+ std::to_string(router_id_global) +"_qtable"; 
    out2file.init("", 0, 0, Output::FILE, qtablefile);

    bool found;
    // use_VC = p.find<bool>("use_VC", true);

    UnitAlgebra link_latency_ua = p.find<UnitAlgebra>("link_lat_global","100ns");
    if ( !link_latency_ua.hasUnits("s") ) {
        output.fatal(CALL_INFO,-1,"link_lat must specified in seconds");
    }
    link_latency_global = (link_latency_ua / UnitAlgebra("1ns")).getRoundedValue();

    link_latency_ua = p.find<UnitAlgebra>("link_lat_local","100ns");
    if ( !link_latency_ua.hasUnits("s") ) {
        output.fatal(CALL_INFO,-1,"link_lat must specified in seconds");
    }
    link_latency_local = (link_latency_ua / UnitAlgebra("1ns")).getRoundedValue();

    learning_rate = p.find<float>("learning_rate", 0.5, found);

    learning_rate2 = p.find<float>("learning_rate2", 0.5, found);

    if(!found){
        learning_rate2 = learning_rate;
    }

    epsilon = p.find<float>("epsilon", 0.1, found);
    save_qtable = p.find<bool>("save_qtable", false, found);
    save_qtable_time = p.find<int>("save_qtable_time", 1000, found);  //us
    save_qtable_time *= 1000; // convert to ns

    // q_vc = p.find<int>("q_vc", 4, found);

    load_qtable = p.find<bool>("load_qtable", false, found);
    max_hops = p.find<int>("max_hops", 6, found);

    assert(max_hops >= 0);

    src_group_q = p.find<bool>("src_group_q", false, found);
    src_mid_group_q = p.find<bool>("src_mid_group_q", false, found);

    // printf("Rtr df %lu construct, use_VC? %d, link latency %llu, learning_rate %f\n", (unsigned long)router_id_global, use_VC, link_latency, learning_rate);
    // printf("Rtr df %d, epsilon %f\n", router_id_global, epsilon);
    pathToQtableFile = p.find<std::string>("pathToQtableFile","",found);
    if(load_qtable && pathToQtableFile.empty()){
        output.fatal(CALL_INFO, -1, "Need to specify the qtable file for loading\n");

    }

    // printf("Rtr df %d, save qtable? %d, file path %s\n", router_id_global, save_qtable, pathToQtableFile.c_str());

    // printf("Rtr df %d, max_hops allowed %d\n", router_id_global, max_hops);
    //===>|

    params.p = p.find<uint32_t>("hosts_per_router");
    params.a = p.find<uint32_t>("routers_per_group");
    params.k = num_ports;
    params.h = p.find<uint32_t>("intergroup_per_router");
    params.g = p.find<uint32_t>("num_groups");
    params.n = p.find<uint32_t>("intergroup_links");

    group_id = rtr_id / params.a;
    router_id = rtr_id % params.a;

    std::string global_route_mode_s = p.find<std::string>("global_route_mode","absolute");
    if ( global_route_mode_s == "absolute" ) global_route_mode = ABSOLUTE;
    else if ( global_route_mode_s == "relative" ) global_route_mode = RELATIVE;
    else {
        output.fatal(CALL_INFO, -1, "Invalid global_route_mode specified: %s.\n",global_route_mode_s.c_str());
    }

    vns = new vn_info[num_vns];
    std::vector<std::string> vn_route_algos;
    if ( p.is_value_array("algorithm") ) {
        p.find_array<std::string>("algorithm", vn_route_algos);
        if ( vn_route_algos.size() != num_vns ) {
            fatal(CALL_INFO, -1, "ERROR: When specifying routing algorithms per VN, algorithm list length must match number of VNs (%d VNs, %lu algorithms).\n",num_vns,vn_route_algos.size());
        }
    }
    else {
        std::string route_algo = p.find<std::string>("algorithm", "minimal");
        for ( int i = 0; i < num_vns; ++i ) vn_route_algos.push_back(route_algo);
    }

    adaptive_threshold = p.find<double>("adaptive_threshold",2.0);

    bool config_failed_links = p.find<bool>("config_failed_links","false");

    //TODO: under assumption of no failed links
    assert(!config_failed_links);

    // Set up the RouteToGroup object

    if ( rtr_id == 0 ) {
        // Get the global link map
        std::vector<int64_t> global_link_map;
        p.find_array<int64_t>("global_link_map", global_link_map);

        std::vector<FailedLink> failed_links;
        p.find_array<FailedLink>("failed_links", failed_links);
        group_to_global_port.init_write("network_", group_id, global_route_mode, params, global_link_map,
                                        config_failed_links, failed_links);
    }
    else {
        group_to_global_port.init("network_", group_id, global_route_mode, params, config_failed_links);
    }

    bool show_adprouting = false;
    
    // Setup the routing algorithms
    int curr_vc = 0;
    for ( int i = 0; i < num_vns; ++i ) {
        vns[i].start_vc = curr_vc;
        vns[i].bias = 50;
        if ( !vn_route_algos[i].compare("valiant") ) {
            if ( params.g <= 2 ) {
                /* 2 or less groups... no point in valiant */
                vns[i].algorithm = MINIMAL;
                vns[i].num_vcs = 2;
            } else {
                vns[i].algorithm = VALIANT;
                // vns[i].num_vcs = 3;

                // 4vc with intermediate router instead of group
                vns[i].num_vcs = 4;
            }
        }
        // else if ( !vn_route_algos[i].compare("adaptive-local") ) {
        else if ( !vn_route_algos[i].compare("ugal-3vc") ) {    
            vns[i].algorithm = UGAL_3VC;
            vns[i].num_vcs = 3;
            show_adprouting = true;
        }
        else if ( !vn_route_algos[i].compare("ugal-4vc") ) {    
            vns[i].algorithm = UGAL_4VC;
            vns[i].num_vcs = 4;
            show_adprouting = true;
        }

        else if ( !vn_route_algos[i].compare("par") ) {
            vns[i].algorithm = PAR;
            vns[i].num_vcs = 5;
            show_adprouting = true;
        }

        else if ( !vn_route_algos[i].compare("minimal") ) {
            vns[i].algorithm = MINIMAL;
            vns[i].num_vcs = 2;
        }
        else if ( !vn_route_algos[i].compare("ugal") ) {
            vns[i].algorithm = UGAL;
            vns[i].num_vcs = 3;
        }
        else if ( !vn_route_algos[i].compare("min-a") ) {
            vns[i].algorithm = MIN_A;
            vns[i].num_vcs = 2;
        }

        else if ( !vn_route_algos[i].compare("q-adaptive") ) {
            vns[i].algorithm = Q1;
            if(src_group_q){
                assert(!src_mid_group_q);
                vns[i].num_vcs = max_hops + 3;
            }
            else if (src_mid_group_q) {
                assert(!src_group_q);
                vns[i].num_vcs = 2 + 3;
            }
            else{
                vns[i].num_vcs = max_hops + 3;
            }
            useQrouting = 1;
        }
        else {
            fatal(CALL_INFO_LONG,1,"ERROR: Unknown routing algorithm specified: %s\n",vn_route_algos[i].c_str());
        }
        curr_vc += vns[i].num_vcs;
    }

    rng = new RNG::XORShiftRNG(rtr_id+1);
    rng_q = new RNG::XORShiftRNG(rtr_id+1010);

    output.verbose(CALL_INFO, 1, 1, "%u:%u:  ID: %u   Params:  p = %u  a = %u  k = %u  h = %u  g = %u\n",
            group_id, router_id, rtr_id, params.p, params.a, params.k, params.h, params.g);

    // printf("RtrDF, G %u:rtr %u:  ID: %u   Params:  p = %u  a = %u  k = %u  h = %u  g = %u\n", group_id, router_id, rtr_id, params.p, params.a, params.k, params.h, params.g);

    q_threshold1 = p.find<double>("q_threshold1", 0.0, found);
    assert(q_threshold1>=0);

    q_threshold2 = p.find<double>("q_threshold2", 0.0, found);
    assert(q_threshold2>=0);


    qtable_row_type = p.find<std::string>("qtable_row_type", "g");
    qtable_rows = 0;
    if ( qtable_row_type == "g" ) {
        qtable_rows = params.g -1;
    }
    else if ( qtable_row_type == "r" ) {
        qtable_rows = (params.g * params.a ) -1;
    }
    else if ( qtable_row_type == "n" ) {
        qtable_rows = (params.g * params.a * params.p ); // no need to -1, store every destination Est.
    }
    else if ( qtable_row_type == "destG_srcN" ){
        qtable_rows = (params.g -1) * params.p;
    }
    else{
        output.fatal(CALL_INFO, -1, "Unknown qtable row type: %s\n", qtable_row_type.c_str());
    }

    //set the correct tables
    if (useQrouting > 0){
        // yao 2021-10-8 setQtable is called by hr_router
        // setQtable();

        for(int i = params.a-1; i<params.k - params.p; i++){
            //global port index - num_host_port
            global_port_idx.push_back(i);
        }
 
    }

    // a periodic clock handler for general debug purpose 
    prid_func = p.find<bool>("perid_func", false, found);
    if(prid_func && rtr_id == 2){

        FILE *file = fopen("rtr2vc.csv", "w");
        fclose(file);

        UnitAlgebra periodic_interval("100ns");
        if ( !periodic_interval.hasUnits("s") ) {
            output.fatal(CALL_INFO,-1,"Dragonfly?? something is wrong");
        }

        peridClock_handler = new Clock::Handler<topo_dragonfly>(this,&topo_dragonfly::perid_funct_handler);
        perid_tc = registerClock( periodic_interval, peridClock_handler, false);

    }

    qtable_bcast = NOBCAST;

    if(rtr_id == 0 ){
        
        output.output("\n---------");
        output.output("Dragonfly Topology: %d nodes", params.p * params.a * params.g);
        output.output("---------\n");
        output.output("config_failed_links: (bool) %d\n", config_failed_links);
        output.output("host/router %d\trouter/group %d\tnum_group %d\n", 
            params.p, params.a, params.g);

        output.output("num_ports %d\tglobal_link %d\tglobal_link/router %d\n", 
            params.k, params.n, params.h);
        
        output.output("link latency:local %ld\tglobal %ld\n", link_latency_local, link_latency_global);

        output.output("total VC num: %d\n", curr_vc);
        output.output("number of VN : %d\n", num_vns);
        output.output("periodic func %d\n", prid_func);

        output.output("routing algorithm:\n\t");
        for(auto algo : vn_route_algos){
            output.output("%s ", algo.c_str());
        }
        output.output("\n");

        if(show_adprouting){
            output.output("adp threhold %.2f\n", adaptive_threshold);
        }
        if( useQrouting ){ 

            output.output("lr %.4f\tlr2 %.4f\tq_threshold1 %.4f, q_threshold2 %.4f\n", learning_rate, learning_rate2, q_threshold1, q_threshold2);

            output.output("exp %.4f\tsrc_group_q %d\tsrc_mid_group_q %d\tmax hops %d\n", 
                epsilon, src_group_q, src_mid_group_q, max_hops);
            
            output.output("qtable type %s\tqtable rows: %d\tloadqtable %d\tpathToQFile %s\n", qtable_row_type.c_str(), qtable_rows, load_qtable, pathToQtableFile.c_str());

            output.output("saveqtable %d\ttsaveqtable_time(us) %d\n", save_qtable, save_qtable_time/1000);
            
            output.output("qtable Bcast %d\n", qtable_bcast);

            output.output("global ports are ( - num_hosts): ");
            for(auto x : global_port_idx){
                output.output("%d ",x);
            }
            output.output("\n");
        }
        output.output("-------------------------------\n");

    }
}

topo_dragonfly::~topo_dragonfly()
{
    delete[] vns;
    if (useQrouting > 0){
        delete[] qtable;
    }
}

void topo_dragonfly::route_nonadaptive(int port, int vc, internal_router_event* ev)
{
    int msgid = ev->getTraceID();
    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    // Break this up by port type
    uint32_t next_port = 0;
    if ( (uint32_t)port < params.p ) {
        // Host ports
        if ( td_ev->dest.group == td_ev->src_group ) {
            // Packets stays within the group
            if ( td_ev->dest.router == router_id ) {
                // Stays within the router
                next_port = td_ev->dest.host;
            }
            else {
                // Route to the router specified by mid_group.  If
                // this is a direct route then mid_group will equal
                // router and packet will go direct.
                next_port = port_for_router(td_ev->dest.mid_group);
            }
        }
        else {
            // Packet is leaving group.  Simply route to group
            // specified by mid_group.  If this is a direct route then
            // mid_group will be set to group.
            next_port = port_for_group(td_ev->dest.mid_group, td_ev->global_slice);
        }
    }
    else if ( (uint32_t)port < ( params.p + params.a - 1) ) {
        // Intragroup links
        if ( td_ev->dest.group == group_id ) {
            // In final group
            if ( td_ev->dest.router == router_id ) {
                // In final router, route to host port
                next_port = td_ev->dest.host;
            }
            else {
                // This is a valiantly routed packet within a group.
                // Need to increment the VC and route to correct
                // router.
                // td_ev->setVC(vc+1);

                // yao as q routing implement change vc differently, only happens in qrouting() funciton
                // if (algorithm != Q)
                td_ev->setVC(vc+1);

                next_port = port_for_router(td_ev->dest.router);
            }
        }
        else {
            // Not in correct group, should route out one of the
            // global links
            if ( td_ev->dest.mid_group != group_id ) {
                next_port = port_for_group(td_ev->dest.mid_group, td_ev->global_slice);
            } else {
                next_port = port_for_group(td_ev->dest.group, td_ev->global_slice);
            }
        }
    }
    else { // global
        /* Came in from another group.  Increment VC */
        // td_ev->setVC(vc+1);
        td_ev->setVC(vc+1);
        // if (algorithm != Q)
        //     td_ev->setVC(vc+1);

        if ( td_ev->dest.group == group_id ) {
            if ( td_ev->dest.router == router_id ) {
                // In final router, route to host port
                next_port = td_ev->dest.host;
            }
            else {
                // Go to final router
                next_port = port_for_router(td_ev->dest.router);
            }
        }
        else {
            // Just passing through on a valiant route.  Route
            // directly to final group
            next_port = port_for_group(td_ev->dest.group, td_ev->global_slice);
        }
    }

    output.verbose(CALL_INFO, 1, 1, "%u:%u, Recv: %d/%d  Setting Next Port/VC:  %u/%u\n", group_id, router_id, port, vc, next_port, td_ev->getVC());
    td_ev->setNextPort(next_port);
}


void topo_dragonfly::route_ugal(int port, int vc, internal_router_event* ev)
{
    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);
    int vn = ev->getVN();

    // Need to determine where we are in the routing.  We can check
    // this based on input port

    // Input port
    if ( port < params.p ) {
        // Packet stays in group.
        if ( td_ev->dest.group == group_id ) {
            // Check to see if the dest is in the same router
            if ( td_ev->dest.router == router_id ) {
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }

            // Adaptive routing when packet stays in group.
            // Check to see if we take the valiant route or direct route.

            int direct_route_port = port_for_router(td_ev->dest.router);
            int direct_route_weight = output_queue_lengths[direct_route_port * num_vcs + vc];

            int valiant_route_port = port_for_router(td_ev->dest.mid_group);
            int valiant_route_weight = output_queue_lengths[valiant_route_port * num_vcs + vc];

            if ( direct_route_weight <= 2 * valiant_route_weight + vns[vn].bias ) {
                td_ev->setNextPort(direct_route_port);
            }
            else {
                td_ev->setNextPort(valiant_route_port);
            }

            return;
        }

        // Packet leaves the group
        else {
            // Need to find the lowest weighted route.  Loop over all
            // the slices.
            int min_weight = std::numeric_limits<int>::max();
            std::vector<std::pair<int,int> > min_ports;
            for ( int i = 0; i < params.n; ++i ) {
                // Direct routes
                int weight;
                int port = port_for_group(td_ev->dest.group, i);
                if ( port != -1 ) {
                    weight = output_queue_lengths[port * num_vcs + vc];

                    if ( weight == min_weight ) {
                        min_ports.emplace_back(port,i);
                    }
                    else if ( weight < min_weight ) {
                        min_weight = weight;
                        min_ports.clear();
                        min_ports.emplace_back(port,i);
                    }
                }

                // Valiant routes
                port = port_for_group(td_ev->dest.mid_group, i);
                if ( port != -1 ) {
                    weight = 2 * output_queue_lengths[port * num_vcs + vc] + vns[vn].bias;

                    if ( weight == min_weight ) {
                        min_ports.emplace_back(port,i);
                    }
                    else if ( weight < min_weight ) {
                        min_weight = weight;
                        min_ports.clear();
                        min_ports.emplace_back(port,i);
                    }
                }
            }

            auto& route = min_ports[rng->generateNextUInt32() % min_ports.size()];
            td_ev->setNextPort(route.first);
            td_ev->global_slice = route.second;
            return;
        }
    }

    // Intragroup links
    else if ( port < ( params.p + params.a -1 ) ) {
        // In final group
        if ( td_ev->dest.group == group_id ) {
            if ( td_ev->dest.router == router_id ) {
                // In final router, route to host port
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }
            else {
                // This is a valiantly routed packet within a group.
                // Need to increment the VC and route to correct
                // router.  We'll use VC2 to avoid interfering with
                // valiant traffic routing through the group.
                td_ev->setVC(vc+2);
                td_ev->setNextPort(port_for_router(td_ev->dest.router));
                return;
            }
        }
        // Need to route it over the global link
        else {
            if ( td_ev->dest.mid_group == group_id ) {
                // In valiant group, just route out to next group
                td_ev->setNextPort( port_for_group(td_ev->dest.group, td_ev->global_slice) );
                return;
            }

            // It's possible that there are links to both the valiant
            // and direct group on the same slice (we don't detect the
            // case where there are routes to both but on a different
            // slice).  We'll need to check both and if they both
            // exist take the port with the least weight.

            // Check direct route first
            int direct_port = port_for_group(td_ev->dest.group, td_ev->global_slice);
            int valiant_port = port_for_group(td_ev->dest.mid_group, td_ev->global_slice);

            // Need to see if these are global ports on this router.
            // If not, then we won't consider them.  At least one of
            // these will be a global port from this router.
            int min_port = std::numeric_limits<int>::max();
            if ( direct_port != -1 && is_port_global(direct_port) ) {
                min_port = direct_port;
            }

            if ( valiant_port != -1 && is_port_global(valiant_port) ) {
                if ( min_port == std::numeric_limits<int>::max() ) {
                    min_port = valiant_port;
                }
                else {
                    // Need to check weights
                    int direct_weight = output_queue_lengths[direct_port * num_vcs + vc];
                    int valiant_weight = 2 * output_queue_lengths[valiant_port * num_vcs + vc] + vn[vns].bias;
                    if ( direct_weight > valiant_weight ) {
                        min_port = valiant_port;
                    }
                }
            }
            td_ev->setNextPort(min_port);
            return;
        }
    }
    // Came in from global routes
    else {
        // Need to increment the VC
        vc++;
        td_ev->setVC(vc);

        // See if we are in the target group
        if ( td_ev->dest.group == group_id ) {
            if ( td_ev->dest.router == router_id ) {
                // Deliver to output port
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }
            else {
                // Need to route to the dest router
                td_ev->setNextPort(port_for_router(td_ev->dest.router));
                return;
            }
        }

        // Just routing through.  Need to look at all possible routes
        // to the dest group and pick the lowest weighted route
        int min_weight = std::numeric_limits<int>::max();
        std::vector<std::pair<int,int> > min_ports;

        // Look through all routes.  If the port is in current router,
        // weight with 1, other weight with 2
        for ( int i = 0; i < params.n; ++i ) {
            int port = port_for_group(td_ev->dest.group, i);
            if ( port == -1 ) continue;
            int weight = output_queue_lengths[port * num_vcs + vc];

            if ( !is_port_global(port) ) weight *= 2;

            if ( weight == min_weight ) {
                min_ports.emplace_back(port,i);
            }
            else if ( weight < min_weight ) {
                min_weight = weight;
                min_ports.clear();
                min_ports.emplace_back(port,i);
            }
        }
        auto& route = min_ports[rng->generateNextUInt32() % min_ports.size()];
        td_ev->setNextPort(route.first);
        td_ev->global_slice = route.second;
        return;
    }

}

void topo_dragonfly::route_mina(int port, int vc, internal_router_event* ev)
{
    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);
    int vn = ev->getVN();

    // Need to determine where we are in the routing.  We can check
    // this based on input port

    // Input port
    if ( port < params.p ) {
        // Packet stays in group.
        if ( td_ev->dest.group == group_id ) {

            // Check to see if the dest is in the same router
            if ( td_ev->dest.router == router_id ) {
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }

            // No adaptive routing, just route to destination router
            td_ev->setNextPort(port_for_router(td_ev->dest.router));
            return;
        }

        // Packet leaves the group
        else {
            // Need to find the lowest weighted route.  Loop over all
            // the slices, looking only at minimal routes.  For now,
            // just weight all paths equally.
            int min_weight = std::numeric_limits<int>::max();
            std::vector<std::pair<int,int> > min_ports;
            for ( int i = 0; i < params.n; ++i ) {

                // Direct routes
                int weight;
                int port = port_for_group(td_ev->dest.group, i);
                if ( port != -1 ) {
                    int hops = hops_to_router(td_ev->dest.group, td_ev->dest.router, i);
                    // Weight by hop count, thus favoring shorter
                    // paths.  The "+ hops" on the end is to make
                    // shorter paths win ties.
                    weight = hops * output_queue_lengths[port * num_vcs + vc] + hops;

                    if ( weight == min_weight ) {
                        min_ports.emplace_back(port,i);
                    }
                    else if ( weight < min_weight ) {
                        min_weight = weight;
                        min_ports.clear();
                        min_ports.emplace_back(port,i);
                    }
                }
            }

            auto& route = min_ports[rng->generateNextUInt32() % min_ports.size()];
            td_ev->setNextPort(route.first);
            td_ev->global_slice = route.second;
            return;
        }
    }

    // Intragroup links
    else if ( port < ( params.p + params.a -1 ) ) {
        // In final group
        if ( td_ev->dest.group == group_id ) {
            if ( td_ev->dest.router == router_id ) {
                // In final router, route to host port
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }
            else {
                // Shouldn't happen with minimal traffic
                merlin_abort.fatal(CALL_INFO,1,"INTERNAL ERROR: routing error with min-a routing.");
            }
        }
        // Need to route it over the global link
        else {
            // Find the route.  We stored the global slice when
            // initially routing, so this should be a global link.
            td_ev->setNextPort(port_for_group(td_ev->dest.group, td_ev->global_slice));
            return;
        }
    }
    // Came in from global routes
    else {
        // Need to increment the VC
        vc++;
        td_ev->setVC(vc);

        // See if we are in the target group
        if ( td_ev->dest.group == group_id ) {
            if ( td_ev->dest.router == router_id ) {
                // Deliver to output port
                td_ev->setNextPort(td_ev->dest.host);
                return;
            }
            else {
                // Need to route to the dest router
                td_ev->setNextPort(port_for_router(td_ev->dest.router));
                return;
            }
        }
        else {
            // Shouldn't happen for minimal routing
            merlin_abort.fatal(CALL_INFO,1,"INTERNAL ERROR: routing error with min-a routing.");
        }
    }
    return;

}

void topo_dragonfly::route_adaptive_local(int port, int vc, internal_router_event* ev)
{   

    int vn = ev->getVN();
    assert( vns[vn].algorithm == ADAPTIVE_LOCAL );

    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();

    // For now, we make the adaptive routing decision only at the
    // input to the network and at the input to a group for adaptively
    // routed packets
    if ( port >= params.p && port < (params.p + params.a-1) ) return;


    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    // Adaptive routing when packet stays in group
    if ( port < params.p && td_ev->dest.group == group_id ) {

        // If we're at the correct router, no adaptive needed
        if ( td_ev->dest.router == router_id) return;


        int direct_route_port = port_for_router(td_ev->dest.router);
        int direct_route_credits = output_credits[direct_route_port * num_vcs + vc];

        int valiant_route_port = port_for_router(td_ev->dest.mid_group);
        int valiant_route_credits = output_credits[valiant_route_port * num_vcs + vc];

        if ( valiant_route_credits > (int)((double)direct_route_credits * adaptive_threshold) ) {
            
            td_ev->setNextPort(valiant_route_port);
            rtr_ev->setAdpRouted();
        }
        else {
            td_ev->setNextPort(direct_route_port);
        }        
        return;
    }

    // If the dest is in the same group, no need to adaptively route
    if ( td_ev->dest.group == group_id ) return;


    // Based on the output queue depths, choose minimal or valiant
    // route.  We'll chose the direct route unless the remaining
    // output credits for the direct route is half of the valiant
    // route.  We'll look at two slice for each direct and indirect,
    // giving us a total of 4 routes we are looking at.  For packets
    // which came in adaptively on global links, look at two direct
    // routes and chose between the.
    int direct_slice1 = td_ev->global_slice_shadow;
    // int direct_slice2 = td_ev->global_slice;
    int direct_slice2 = (td_ev->global_slice_shadow + 1) % params.n;
    int direct_route_port1 = port_for_group(td_ev->dest.group, direct_slice1, 0 );
    int direct_route_port2 = port_for_group(td_ev->dest.group, direct_slice2, 1 );
    int direct_route_credits1 = output_credits[direct_route_port1 * num_vcs + vc];
    int direct_route_credits2 = output_credits[direct_route_port2 * num_vcs + vc];
    int direct_slice;
    int direct_route_port;
    int direct_route_credits;
    if ( direct_route_credits1 > direct_route_credits2 ) {
        direct_slice = direct_slice1;
        direct_route_port = direct_route_port1;
        direct_route_credits = direct_route_credits1;
    }
    else {
        direct_slice = direct_slice2;
        direct_route_port = direct_route_port2;
        direct_route_credits = direct_route_credits2;
    }

    int valiant_slice = 0;
    int valiant_route_port = 0;
    int valiant_route_credits = 0;

    if ( port >= (params.p + params.a-1) ) {
        // Global port, no indirect routes.  Set credits negative so
        // it will never get chosen
        valiant_route_credits = -1;
    }
    else {
        int valiant_slice1 = td_ev->global_slice;
        // int valiant_slice2 = td_ev->global_slice;
        int valiant_slice2 = (td_ev->global_slice + 1) % params.n;
        int valiant_route_port1 = port_for_group(td_ev->dest.mid_group_shadow, valiant_slice1, 2 );
        int valiant_route_port2 = port_for_group(td_ev->dest.mid_group_shadow, valiant_slice2, 3 );
        int valiant_route_credits1 = output_credits[valiant_route_port1 * num_vcs + vc];
        int valiant_route_credits2 = output_credits[valiant_route_port2 * num_vcs + vc];
        if ( valiant_route_credits1 > valiant_route_credits2 ) {
            valiant_slice = valiant_slice1;
            valiant_route_port = valiant_route_port1;
            valiant_route_credits = valiant_route_credits1;
        }
        else {
            valiant_slice = valiant_slice2;
            valiant_route_port = valiant_route_port2;
            valiant_route_credits = valiant_route_credits2;
        }
    }
    
    if ( valiant_route_credits > (int)((double)direct_route_credits * adaptive_threshold) ) { // Use valiant route
        td_ev->dest.mid_group = td_ev->dest.mid_group_shadow;
        td_ev->setNextPort(valiant_route_port);
        td_ev->global_slice = valiant_slice;
        rtr_ev->setAdpRouted();

    }
    else { // Use direct route
        td_ev->dest.mid_group = td_ev->dest.group;
        td_ev->setNextPort(direct_route_port);
        td_ev->global_slice = direct_slice;
    }
}

void topo_dragonfly::route_packet(int port, int vc, internal_router_event* ev) {
    int vn = ev->getVN();
    if( vns[vn].algorithm == UGAL_3VC ) {
        route_ugal_3vc(port,vc,ev);
    }
    else if (vns[vn].algorithm == UGAL_4VC){
        route_ugal_4vc(port,vc,ev);
    }
    else if ( vns[vn].algorithm == Q1 )
    {
        q_adaptive(port,vc,ev);
    }

    else if (vns[vn].algorithm == MINIMAL){
        route_minimal(port,vc,ev);
    }
    else if (vns[vn].algorithm == VALIANT){
        route_valiant(port,vc,ev);
    }

    else if (vns[vn].algorithm == PAR){
        route_PAR(port,vc,ev);
    }
    else{
        output.fatal(CALL_INFO, -1, "Dragonfly Unknown rting: %d \n", vns[vn].algorithm);
    }
}

internal_router_event* topo_dragonfly::process_input(RtrEvent* ev)
{
    int msgid = ev->getTraceID();
    dgnflyAddr dstAddr = {0, 0, 0, 0, 0, 0, 0}; 
    idToLocation(ev->getDest(), &dstAddr); 

    int vn = ev->getRouteVN();

    switch (vns[vn].algorithm) {
    case MINIMAL:
    case Q1:
        if ( dstAddr.group == group_id ) {
            dstAddr.mid_group = dstAddr.router;
        }
        else {
            dstAddr.mid_group = dstAddr.group;
        }
        break;
    case VALIANT:
    case ADAPTIVE_LOCAL:
    case UGAL:
    case UGAL_3VC:
    case UGAL_4VC:
    case PAR:        
        dstAddr.mid_router = params.a+10; 
        dstAddr.mid_group = params.g+10;
        break;
    default:
        output.fatal(CALL_INFO, -1, "Unknown rting: %d \n", vns[vn].algorithm);
    }
    dstAddr.mid_group_shadow = dstAddr.mid_group;

    topo_dragonfly_event *td_ev = new topo_dragonfly_event(dstAddr);
    td_ev->src_group = group_id;
    td_ev->setEncapsulatedEvent(ev);
    td_ev->setVC(vns[vn].start_vc);
    td_ev->global_slice = ev->getTrustedSrc() % params.n;
    td_ev->global_slice_shadow = ev->getTrustedSrc() % params.n;

    //TODO: so far, consider system at its maximum size
    assert(params.n == 1);

    if ( td_ev->getTraceType() != SST::Interfaces::SimpleNetwork::Request::NONE ) {
        output.output("TRACE(%d): process_input():"
                      " mid_group_shadow = %u\n",
                      td_ev->getTraceID(),
                      td_ev->dest.mid_group_shadow);
    }
    return td_ev;
}

std::pair<int,int>
topo_dragonfly::getDeliveryPortForEndpointID(int ep_id)
{
    return std::make_pair<int,int>(ep_id / params.p, ep_id % params.p);
}


int
topo_dragonfly::routeControlPacket(CtrlRtrEvent* ev)
{
    const auto& dest = ev->getDest();

    int dest_id = dest.addr;

    // Check to see if this event is for the current router. If so,
    // return -1/
    if ( dest.addr_is_router ) {
        if ( dest_id == rtr_id ) return -1;
        // If this is a router ID, multiply by number of hosts per
        // router (params.p) in order to get the id for an endpoint in
        // this router.  Then we can use the idToLocation() function.
        else dest_id *= params.p;
    }
    // Addr is an endpoint id, check to see if the event is actually
    // for a router and, if so, if it is for the current router
    else if ( dest.addr_for_router && ((dest_id / params.p) == rtr_id) ) {
        return -1;
    }

    dgnflyAddr addr;
    idToLocation(dest_id,&addr);

    // Just get next port on minimal route
    int next_port;
    if ( addr.group != group_id ) {
        next_port = port_for_group(addr.group, 0 /* global slice */);
    }
    else if ( addr.router != router_id ) {
        next_port = port_for_router(addr.router);
    }
    else {
        next_port = addr.host;
    }
    return next_port;


}


void topo_dragonfly::routeInitData(int port, internal_router_event* ev, std::vector<int> &outPorts)
{
    RtrEvent* tmp_ev =  ev->getEncapsulatedEvent();
    int msgid = tmp_ev->getTraceID();
    bool broadcast_to_groups = false;
    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);
    if ( td_ev->dest.host == (uint32_t)INIT_BROADCAST_ADDR ) {
        if ( (uint32_t)port >= (params.p + params.a-1) ) {
            /* Came in from another group.
             * Send to locals, and other routers in group
             */
            for ( uint32_t p  = 0 ; p < (params.p + params.a-1) ; p++ ) {
                outPorts.push_back((int)p);
            }
        } else if ( (uint32_t)port >= params.p ) {
            /* Came in from another router in group.
             * send to hosts
             * if this is the source group, send to other groups
             */
            for ( uint32_t p = 0 ; p < params.p ; p++ ) {
                outPorts.push_back((int)p);
            }
            if ( td_ev->src_group == group_id ) {
                broadcast_to_groups = true;
            }
        } else {
            /* Came in from a host
             * Send to all other hosts and routers in group, and all groups
             */
            // for ( int p = 0 ; p < (int)params.k ; p++ ) {
            for ( int p = 0 ; p < (int)(params.p + params.a - 1) ; p++ ) {
                if ( p != port )
                    outPorts.push_back((int)p);
            }
            broadcast_to_groups = true;
        }

        if ( broadcast_to_groups ) {
            for ( int p = 0; p < (int)(params.g - 1); p++ ) {
                const RouterPortPair& pair = group_to_global_port.getRouterPortPair(p,0);
                if ( pair.router == router_id ) outPorts.push_back((int)(pair.port));
            }
        }
    } else {
        // Not all data structures used for routing during run are
        // initialized yet, so we need to just do a quick minimal
        // routing scheme for init.
        // route(port, 0, ev);

        // Minimal Route
        int next_port;
        if ( td_ev->dest.group != group_id ) {
            next_port = port_for_group_init(td_ev->dest.group, td_ev->global_slice);
        }
        else if ( td_ev->dest.router != router_id ) {
            next_port = port_for_router(td_ev->dest.router);
        }
        else {
            next_port = td_ev->dest.host;
        }
        outPorts.push_back(next_port);
    }

}


internal_router_event* topo_dragonfly::process_InitData_input(RtrEvent* ev)
{
    int msgid = ev->getTraceID();
    dgnflyAddr dstAddr;
    idToLocation(ev->getDest(), &dstAddr);
    topo_dragonfly_event *td_ev = new topo_dragonfly_event(dstAddr);
    td_ev->src_group = group_id;
    td_ev->setEncapsulatedEvent(ev);

    return td_ev;
}


Topology::PortState topo_dragonfly::getPortState(int port) const
{
    if ( is_port_endpoint(port) ) return R2N;
    else if ( is_port_local_group(port) ) return R2R;
    else {
        // These are global ports, see if they are failed
        if ( group_to_global_port.isFailedPort(RouterPortPair(router_id,port)) ) return FAILED;
        else return R2R;
    }
}

std::string topo_dragonfly::getPortLogicalGroup(int port) const
{
    if ( (uint32_t)port < params.p ) return "host";
    if ( (uint32_t)port >= params.p && (uint32_t)port < (params.p + params.a - 1) ) return "group";
    else return "global";
}

int
topo_dragonfly::getEndpointID(int port)
{
    return (group_id * (params.a /*rtr_per_group*/ * params.p /*hosts_per_rtr*/)) +
        (router_id * params.p /*hosts_per_rtr*/) + port;
}

void
topo_dragonfly::setOutputBufferCreditArray(int const* array, int vcs)
{
    output_credits = array;
    num_vcs = vcs;
}


void
topo_dragonfly::setOutputQueueLengthsArray(int const* array, int vcs)
{
    output_queue_lengths = array;
    num_vcs = vcs;
}

void topo_dragonfly::idToLocation(int id, dgnflyAddr *location)
{
    if ( id == INIT_BROADCAST_ADDR) {
        location->group = (uint32_t)INIT_BROADCAST_ADDR;
        location->mid_group = (uint32_t)INIT_BROADCAST_ADDR;
        location->router = (uint32_t)INIT_BROADCAST_ADDR;
        location->host = (uint32_t)INIT_BROADCAST_ADDR;
    } else {
        uint32_t hosts_per_group = params.p * params.a;
        location->group = id / hosts_per_group;
        location->router = (id % hosts_per_group) / params.p;
        location->host = id % params.p;
        location->router_global = location->group * params.a + location->router;
    }
}


int32_t topo_dragonfly::router_to_group(uint32_t group)
{

    if ( group < group_id ) {
        return group / params.h;
    } else if ( group > group_id ) {
        return (group-1) / params.h;
    } else {
        output.fatal(CALL_INFO, -1, "Trying to find router to own group.\n");
        return 0;
    }
}


int32_t topo_dragonfly::hops_to_router(uint32_t group, uint32_t router, uint32_t slice)
{
    int hops = 1;
    const RouterPortPair& pair = group_to_global_port.getRouterPortPair(group,slice);
    if ( pair.router != router_id ) hops++;
    const RouterPortPair& pair2 = group_to_global_port.getRouterPortPairForGroup(group, group_id, slice);
    if ( pair2.router != router ) hops++;
    return hops;
}

/* returns local router port if group can't be reached from this router */
int32_t topo_dragonfly::port_for_group(uint32_t group, uint32_t slice, int id)
{
    assert(group >=0 && group < params.g);

    const RouterPortPair& pair = group_to_global_port.getRouterPortPair(group,slice);
    if ( group_to_global_port.isFailedPort(pair) ) {
        // printf("******** Skipping failed port ********\n");
        return -1;
    }

    if ( pair.router == router_id ) {
        return pair.port;
    } else {
        return port_for_router(pair.router);
    }
}

// Always ignore failed links during init
int32_t topo_dragonfly::port_for_group_init(uint32_t group, uint32_t slice)
{
    // TraceFunction trace(CALL_INFO_LONG);
    // trace.output("Routing from group %d to group %u over slice %u\n",group_id,group,slice);
    const RouterPortPair& pair = group_to_global_port.getRouterPortPair(group,slice);
    // trace.output("This maps to router %d and port %d\n",pair.router,pair.port);

    if ( pair.router == router_id ) {
        return pair.port;
    } else {
        return port_for_router(pair.router);
    }
}


int32_t topo_dragonfly::port_for_router(uint32_t router)
{
    assert(router >= 0 && router < params.a);
    
    uint32_t tgt = params.p + router;
    if ( router > router_id ) tgt--;
    return tgt;
}


// additional functions

int topo_dragonfly::choose_port_for_group(topo_dragonfly_event* td_ev, int dest_g){

    assert(dest_g>=0);

    int direct_slice1 = td_ev->global_slice_shadow;
    int direct_slice2 = (td_ev->global_slice_shadow + 1) % params.n;
    int direct_route_port1 = port_for_group(dest_g, direct_slice1, 0 );
    int direct_route_port2 = port_for_group(dest_g, direct_slice2, 1 );

    int direct_crt1 = 0;
    int direct_crt2 = 0;

    int port_to_g = -1;

    const int vn = td_ev->getVN();
    const int start_vc = vns[vn].start_vc;
    const int end_vc = start_vc + vns[vn].num_vcs;

    assert(start_vc>=0);
    assert(end_vc <= num_vcs);
    
    // direct_crt1 = output_queue_lengths[direct_route_port1 * num_vcs + vc_direct];
    // direct_crt2 = output_queue_lengths[direct_route_port2 * num_vcs + vc_direct];

    for(int i=start_vc; i<end_vc; i++){
        direct_crt1 += output_queue_lengths[direct_route_port1 * num_vcs + i];
        direct_crt2 += output_queue_lengths[direct_route_port2 * num_vcs + i];
    }

    if ( direct_crt1 < direct_crt2 ){
        td_ev->global_slice  = direct_slice1;
        port_to_g = direct_route_port1;
    }
    else {
        td_ev->global_slice = direct_slice2;
        port_to_g = direct_route_port2;       
    }

    assert(port_to_g >= params.p);

    return port_to_g;

}


void topo_dragonfly::route_minimal(int port, int vc, internal_router_event* ev){

    int next_port = -1;
    int next_vc = -1;

    int vn = ev->getVN();
    int start_vc = vns[vn].start_vc;

    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();

    int link_hops = rtr_ev->getNumHops();

    assert(link_hops < 4);

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    if(group_id == td_ev->dest.group){
        next_vc = start_vc + 1;
        if(router_id == td_ev->dest.router){
            next_port = td_ev->dest.host;
        }
        else {
            next_port = port_for_router(td_ev->dest.router);
        }
    }
    else{
        assert(group_id == td_ev->src_group);
        next_vc = start_vc;
        next_port = choose_port_for_group(td_ev, td_ev->dest.group);
    }

    assert(next_vc >=0);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

//same group message min routing
//no local valiant 
void topo_dragonfly::route_valiant(int port, int vc, internal_router_event* ev){
    int next_port = -1;
    int next_vc = -1;

    const int vn = ev->getVN();
    assert( vns[vn].algorithm == VALIANT );

    const int start_vc = vns[vn].start_vc;

    RtrEvent* rtr_ev = ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();

    int link_hops = rtr_ev->getNumHops();

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    int current_vc = td_ev->getVC();
    assert(current_vc == vc);

    assert(link_hops < 7);

    if(group_id == td_ev->dest.group)
    {
        next_vc = start_vc + 3;
        if(router_id == td_ev->dest.router){
            next_port = td_ev->dest.host;
        }
        else{
            next_port = port_for_router(td_ev->dest.router);
        }
    }
    else if(group_id == td_ev->src_group){
        assert(td_ev->src_group != td_ev->dest.group);
        assert(current_vc == start_vc);
        if(port<params.p){
            //local port
            assert(td_ev->dest.mid_group == params.g + 10);
            assert(td_ev->dest.mid_router == params.a + 10);

            uint32_t tmp_router = params.g * params.a + 10;
            uint32_t tmp_group = params.g+ 10;

            do {
                tmp_router = rng->generateNextUInt32() % (params.g * params.a);
                tmp_group = tmp_router / params.a;
            } while ( tmp_group == group_id || tmp_group == td_ev->dest.group );

            assert(tmp_group>=0 && tmp_group<params.g);

            td_ev->dest.mid_group = tmp_group;
            td_ev->dest.mid_router = (tmp_router % params.a);
            rtr_ev->midgroup = tmp_group;
        }

        assert(td_ev->dest.mid_group >= 0 && td_ev->dest.mid_group < params.g );

        next_vc = start_vc;
        next_port = choose_port_for_group(td_ev, td_ev->dest.mid_group);
    }
    else {
        assert(group_id == td_ev->dest.mid_group);
        assert(link_hops<5);
        assert(port >= params.p);
        assert(td_ev->dest.mid_router>=0 && td_ev->dest.mid_router<params.a);

        if(router_id == td_ev->dest.mid_router){
            assert(current_vc == start_vc + 1 || current_vc == start_vc + 0);
            next_vc = start_vc + 2;
            next_port = choose_port_for_group(td_ev, td_ev->dest.group);
        }
        else{
            if(current_vc == start_vc + 0){
                assert(port > params.p + params.a - 2); // global port
                next_vc = start_vc + 1;
                next_port = port_for_router(td_ev->dest.mid_router);
            }
            else{
                assert(current_vc == start_vc + 2);
                assert(port >= params.p && port < params.p + params.a - 1); //local port
                next_vc = start_vc + 2;
                next_port = choose_port_for_group(td_ev, td_ev->dest.group);
            }
        }
    }

    assert(next_vc >= start_vc);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

//TODO: slice is used to select differnt global links in case multiple global links exist in current group connect to dest group

std::pair<int,int> topo_dragonfly::adp_select_port(topo_dragonfly_event *td_ev, int port, int vc_start, int vc_end, std::vector<int> vc_min, std::vector<int> vc_nonmin, std::string debuginfo)
{
    assert(vc_min.size()>0);
    assert(vc_nonmin.size()>0);

    const int vn = td_ev->getVN();

    //support for multi VNs, where each VN may have difference VCs
    assert(vc_start == vns[vn].start_vc);
    assert(vc_end == vc_start + vns[vn].num_vcs );

    RtrEvent* rtr_ev =  td_ev->getEncapsulatedEvent();

    int next_port = -1;
    int take_non_min = -1;

    int direct_slice1 = td_ev->global_slice_shadow;
    int direct_slice2 = (td_ev->global_slice_shadow + 1) % params.n;
    int direct_route_port1 = port_for_group(td_ev->dest.group, direct_slice1, 0 );
    int direct_route_port2 = port_for_group(td_ev->dest.group, direct_slice2, 1 );

    int direct_route_port;
    int direct_queue;
    int direct_slice;

    int direct_queue1 = 0;
    int direct_queue2 = 0;
    
    // direct_crt1 = output_queue_lengths[direct_route_port1 * num_vcs + vc_direct];
    // direct_crt2 = output_queue_lengths[direct_route_port2 * num_vcs + vc_direct];

    for(int i=vc_start; i<vc_end; i++){
        direct_queue1 += output_used_credits[direct_route_port1 * num_vcs + i];
        direct_queue2 += output_used_credits[direct_route_port2 * num_vcs + i];

        direct_queue1 += output_queue_lengths[direct_route_port1 * num_vcs + i];
        direct_queue2 += output_queue_lengths[direct_route_port2 * num_vcs + i];
    }

    if ( direct_queue1 < direct_queue2 ){
        direct_slice = direct_slice1;
        direct_route_port = direct_route_port1;
        direct_queue = direct_queue1;
    }
    else {
        direct_slice = direct_slice2;
        direct_route_port = direct_route_port2;
        direct_queue = direct_queue2;        
    }

    int valiant_slice; 
    int valiant_route_port; 
    int valiant_queue;

    int valiant_slice1 = td_ev->global_slice;
    int valiant_slice2 = (td_ev->global_slice + 1) % params.n;
    int valiant_route_port1 = port_for_group(td_ev->dest.mid_group, valiant_slice1, 2 );
    int valiant_route_port2 = port_for_group(td_ev->dest.mid_group, valiant_slice2, 3 );            
    int valiant_queue1 = 0;
    int valiant_queue2 = 0;

    // valiant_crt1 = output_queue_lengths[valiant_route_port1 * num_vcs + vc_valiant];
    // valiant_crt2 = output_queue_lengths[valiant_route_port2 * num_vcs + vc_valiant];

    for(int i=vc_start; i<vc_end; i++){
        valiant_queue1 += output_used_credits[valiant_route_port1 * num_vcs + i];
        valiant_queue2 += output_used_credits[valiant_route_port2 * num_vcs + i];

        valiant_queue1 += output_queue_lengths[valiant_route_port1 * num_vcs + i];
        valiant_queue2 += output_queue_lengths[valiant_route_port2 * num_vcs + i];
    }

    if ( valiant_queue1 < valiant_queue2 ) {
        valiant_slice = valiant_slice1;
        valiant_route_port = valiant_route_port1;
        valiant_queue = valiant_queue1;
    }
    else {
        valiant_slice = valiant_slice2;
        valiant_route_port = valiant_route_port2;
        valiant_queue = valiant_queue2;        
    }

    if(direct_route_port == valiant_route_port){
        direct_queue = 0;
        for(auto vc_m : vc_min){
            direct_queue += output_used_credits[direct_route_port * num_vcs + vc_m];
            direct_queue += output_queue_lengths[direct_route_port * num_vcs + vc_m];
        }

        valiant_queue = 0;
        for(auto vc_nm : vc_nonmin){
            valiant_queue += output_used_credits[valiant_route_port * num_vcs + vc_nm];
            valiant_queue += output_queue_lengths[valiant_route_port * num_vcs + vc_nm];
        }
    }

    // if ( valiant_crt > (int)((double)direct_crt * adaptive_threshold) ){
    if(  direct_queue > (int)( (double) valiant_queue * adaptive_threshold) ) {
        // non min path
        next_port = valiant_route_port;
        td_ev->global_slice = valiant_slice;
        rtr_ev->setAdpRouted();

        take_non_min = 1;
    }
    else{
        // min path
        next_port = direct_route_port;
        td_ev->global_slice = direct_slice;

        take_non_min = 0;
    }
    assert(next_port >=0);
    assert(take_non_min >=0);

    std::pair <int,int> adp_result (take_non_min, next_port);

    return adp_result;
}

void topo_dragonfly::route_ugal_3vc(int port, int vc, internal_router_event* ev)
{   
    const int vn = ev->getVN();
    assert( vns[vn].algorithm == UGAL_3VC );
    const int start_vc = vns[vn].start_vc;
    const int ugal_vcs = vns[vn].num_vcs;
    assert(ugal_vcs == 3);
    const int end_vc = start_vc + ugal_vcs;

    // only using 1 vn so far, not true with multiple VNs
    // assert(ugal_vcs == num_vcs); 

    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();

    int link_hops = rtr_ev->getNumHops();
    
    assert(link_hops < 6);

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    const int current_vc = td_ev->getVC();

    int next_port = -1;
    int next_vc = -1;

    //no local adp routing
    if(group_id == td_ev->dest.group)
    {
        next_vc = start_vc + 2;
        if(router_id == td_ev->dest.router){
            next_port = td_ev->dest.host;
        }
        else{
            next_port = port_for_router(td_ev->dest.router);
        }
    }
    else if(group_id == td_ev->src_group)
    {
        assert(td_ev->src_group != td_ev->dest.group);
        if(port<params.p){
            //host port, do adp routing
            assert(!rtr_ev->getAdpRouted());
            assert(current_vc == start_vc);
            assert(td_ev->dest.mid_group == params.g + 10);
            assert(td_ev->dest.mid_router == params.a + 10);

            uint32_t tmp_router = params.g * params.a + 10;
            uint32_t tmp_group = params.g+ 10;

            do {
                tmp_router = rng->generateNextUInt32() % (params.g * params.a);
                tmp_group = tmp_router / params.a;
            } while ( tmp_group == group_id || tmp_group == td_ev->dest.group );

            assert(tmp_group>=0 && tmp_group<params.g);

            td_ev->dest.mid_group = tmp_group;
            td_ev->dest.mid_router = (tmp_router % params.a);
            rtr_ev->midgroup = tmp_group;

            std::vector<int> vc_min {start_vc+1};
            std::vector<int> vc_nonmin {start_vc+0};

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, end_vc, vc_min, vc_nonmin);

            if( adp_result.first ) {
                // non min path
                assert(rtr_ev->getAdpRouted());
                next_vc = start_vc + 0;
            }
            else{
                // min path
                assert(!rtr_ev->getAdpRouted());
                next_vc = start_vc + 1;
            }
            next_port = adp_result.second;
        }
        
        else {
            assert(port < params.p + params.a + 1); //must be local
            assert(td_ev->dest.mid_group>=0 && td_ev->dest.mid_group<params.g);
            assert(td_ev->dest.mid_router>=0 && td_ev->dest.mid_router<params.a);

            int to_group = -1;
            if(rtr_ev->getAdpRouted()){
                assert(current_vc == start_vc + 0);
                to_group = td_ev->dest.mid_group;
                next_vc = start_vc + 0;
            }
            else{
                assert(current_vc == start_vc + 1);
                to_group = td_ev->dest.group;
                next_vc = start_vc + 1;
            }
            next_port = choose_port_for_group(td_ev, to_group);
            assert(next_port >= params.p + params.a - 1);
        }
    }
    else {
        assert(group_id == td_ev->dest.mid_group);
        assert(link_hops<4);
        assert(port >= params.p);
        assert(current_vc == start_vc + 1 || current_vc == start_vc + 0);
        
        next_port = choose_port_for_group(td_ev, td_ev->dest.group);
        next_vc = start_vc + 1;
    }

    assert(next_vc >= start_vc);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

void topo_dragonfly::route_ugal_4vc(int port, int vc, internal_router_event* ev)
{   
    const int vn = ev->getVN();
    assert( vns[vn].algorithm == UGAL_4VC );
    const int start_vc = vns[vn].start_vc;
    const int ugal_vcs = vns[vn].num_vcs;
    assert(ugal_vcs == 4);

    // this is no longer true since different routing has differnt number of VCs. and each VN may use different routing
    // may have multi vn for QoS
    // assert(ugal_vcs * num_vns == num_vcs ); 
    // assert(start_vc == ugal_vcs*vn);
    const int end_vc = start_vc + ugal_vcs;

    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();

    int link_hops = rtr_ev->getNumHops();
    
    assert(link_hops < 7);

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    const int current_vc = td_ev->getVC();

    int next_port = -1;
    int next_vc = -1;

    //no local adp routing
    if(group_id == td_ev->dest.group)
    {
        next_vc = start_vc + 3;
        if(router_id == td_ev->dest.router){
            next_port = td_ev->dest.host;
        }
        else{
            next_port = port_for_router(td_ev->dest.router);
        }
    }
    else if(group_id == td_ev->src_group)
    {
        assert(td_ev->src_group != td_ev->dest.group);
        if(port<params.p){
            //host port, do adp routing
            assert(!rtr_ev->getAdpRouted());
            assert(current_vc == start_vc);
            assert(td_ev->dest.mid_group == params.g + 10);
            assert(td_ev->dest.mid_router == params.a + 10);

            uint32_t tmp_router = params.g * params.a + 10;
            uint32_t tmp_group = params.g+ 10;

            do {
                tmp_router = rng->generateNextUInt32() % (params.g * params.a);
                tmp_group = tmp_router / params.a;
            } while ( tmp_group == group_id || tmp_group == td_ev->dest.group );

            assert(tmp_group>=0 && tmp_group<params.g);

            td_ev->dest.mid_group = tmp_group;
            td_ev->dest.mid_router = (tmp_router % params.a);
            rtr_ev->midgroup = tmp_group;

            std::vector<int> vc_min {start_vc+2};
            // std::vector<int> vc_min {1,2};
            std::vector<int> vc_nonmin {start_vc+0};

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, end_vc, vc_min, vc_nonmin);

            if( adp_result.first ) {
                // non min path
                assert(rtr_ev->getAdpRouted());
                next_vc = start_vc + 0;

                assert(rtr_ev->val_route_pos == 0);
                rtr_ev->val_route_pos = 1;

            }
            else{
                // min path
                assert(!rtr_ev->getAdpRouted());
                next_vc = start_vc + 2;
            }
            next_port = adp_result.second;
        }
        
        else {
            assert(port < params.p + params.a + 1); //must be local
            assert(td_ev->dest.mid_group>=0 && td_ev->dest.mid_group<params.g);
            assert(td_ev->dest.mid_router>=0 && td_ev->dest.mid_router<params.a);

            int to_group = -1;
            if(rtr_ev->getAdpRouted()){
                assert(current_vc == start_vc + 0);
                to_group = td_ev->dest.mid_group;
                next_vc = start_vc + 0;
            }
            else{
                assert(current_vc == start_vc + 2);
                to_group = td_ev->dest.group;
                next_vc = start_vc + 2;
            }
            next_port = choose_port_for_group(td_ev, to_group);
            assert(next_port >= params.p + params.a - 1);
        }
    }
    else {
        assert(group_id == td_ev->dest.mid_group);
        assert(link_hops<5);
        assert(port >= params.p);
        
        if(router_id == td_ev->dest.mid_router){
            assert(current_vc == start_vc + 1 || current_vc == start_vc + 0);

            next_vc = start_vc + 2;
            next_port = choose_port_for_group(td_ev, td_ev->dest.group);
        }
        else{
            if(current_vc == start_vc + 0){
                assert(port > params.p + params.a - 2); // global port
                next_vc = start_vc + 1;
                next_port = port_for_router(td_ev->dest.mid_router);

                // for debug purpose
                assert(rtr_ev->val_route_pos ==1);
                int minpathport = choose_port_for_group(td_ev, td_ev->dest.group);
                if(minpathport != next_port){
                    rtr_ev->val_route_pos = 3;

                }
            }
            else{
                assert(current_vc == start_vc + 2);
                assert(port >= params.p && port < params.p + params.a - 1); //local port
                next_vc = start_vc + 2;
                next_port = choose_port_for_group(td_ev, td_ev->dest.group);
            }
        }
    }

    assert(next_vc >= start_vc);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

// yao Progressive Adaptive Routing
void topo_dragonfly::route_PAR(int port, int vc, internal_router_event* ev)
{   
    const int vn = ev->getVN();
    assert( vns[vn].algorithm == PAR );
    const int start_vc = vns[vn].start_vc;
    const int par_vcs = vns[vn].num_vcs;
    const int end_vc = start_vc + par_vcs;
    
    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();
    int link_hops = rtr_ev->getNumHops();
    
    assert(link_hops < 8);
    assert(par_vcs == 5);
    // no longer true with multiple VNs
    // assert(par_vcs == num_vcs); // only using 1 vn so far

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    const int current_vc = td_ev->getVC();
    assert(current_vc == vc);

    int next_port = -1;
    int next_vc = -1;

    std::vector<int> vc_min {start_vc+0};
    std::vector<int> vc_nonmin {start_vc+1};

    //no local adp routing
    if(group_id == td_ev->dest.group)
    {
        next_vc = start_vc + 4;
        if(router_id == td_ev->dest.router){
            next_port = td_ev->dest.host;
        }
        else{
            next_port = port_for_router(td_ev->dest.router);
        }
    }
    else if(group_id == td_ev->src_group){
        assert(td_ev->src_group != td_ev->dest.group);

        if(port<params.p){
            //host port, do adp routing
            assert(td_ev->dest.mid_group == params.g + 10);
            assert(td_ev->dest.mid_router == params.a + 10);
            assert(!rtr_ev->getAdpRouted());
            assert(current_vc == start_vc);

            uint32_t tmp_router = params.g * params.a + 10;
            uint32_t tmp_group = params.g+ 10;

            do {
                tmp_router = rng->generateNextUInt32() % (params.g * params.a);
                tmp_group = tmp_router / params.a;
            } while ( tmp_group == group_id || tmp_group == td_ev->dest.group );

            assert(tmp_group>=0 && tmp_group<params.g);

            td_ev->dest.mid_group = tmp_group;
            td_ev->dest.mid_router = (tmp_router % params.a);
            assert(td_ev->dest.mid_group>=0 && td_ev->dest.mid_group<params.g);
            rtr_ev->midgroup = tmp_group;

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, end_vc, vc_min, vc_nonmin, "1st");

            if(adp_result.first){
                //non min routing
                assert(rtr_ev->getAdpRouted());
                assert(rtr_ev->val_route_pos == 0);
                next_vc = start_vc + 1;
                rtr_ev->val_route_pos = 1;
            }
            else{
                // min path
                assert(!rtr_ev->getAdpRouted());
                next_vc = start_vc + 0;
            }

            next_port = adp_result.second;
        }
        else{
            assert(port < params.p + params.a + 1); //must be local
            assert(td_ev->dest.mid_group>=0 && td_ev->dest.mid_group<params.g);

            if(rtr_ev->getAdpRouted()){
                assert(current_vc == start_vc + 1);
                assert(rtr_ev->val_route_pos == 1 || rtr_ev->val_route_pos == 2);
                next_port = choose_port_for_group(td_ev, td_ev->dest.mid_group);
                assert(next_port >= params.p + params.a - 1);
                next_vc = start_vc + 1;
            }
            else{
                assert(current_vc == start_vc + 0);
                assert(rtr_ev->val_route_pos == 0);
                //make a second adp routing decision
                std::pair <int,int> adp_result = adp_select_port(
                    td_ev, port, start_vc, end_vc, vc_min, vc_nonmin, "2nd");

                if(adp_result.first){
                    //take non min
                    assert(rtr_ev->getAdpRouted());
                    next_vc = start_vc + 1;
                    next_port = adp_result.second;
                    rtr_ev->val_route_pos = 2;
                }
                else{
                    assert(!rtr_ev->getAdpRouted());
                    next_vc = start_vc + 3;
                    next_port = adp_result.second;
                    assert(next_port >= params.p + params.a -1); //must be a global port
                }
            }
        }
    }
    else {
        assert(group_id == td_ev->dest.mid_group);
        assert(link_hops<6);
        assert(rtr_ev->getAdpRouted());
        assert(rtr_ev->val_route_pos > 0);
        assert(port >= params.p);
        assert(td_ev->dest.mid_router>=0 && td_ev->dest.mid_router<params.a);

        if(router_id == td_ev->dest.mid_router){
            assert(current_vc == start_vc + 1 || current_vc == start_vc + 2);

            next_vc = start_vc + 3;
            next_port = choose_port_for_group(td_ev, td_ev->dest.group);
        }
        else{
            if(current_vc == start_vc + 1){
                assert(port > params.p + params.a - 2); // global port
                next_vc = start_vc + 2;
                next_port = port_for_router(td_ev->dest.mid_router);

                // for debug purpose
                assert(rtr_ev->val_route_pos ==1 || rtr_ev->val_route_pos ==2);
                int minpathport = choose_port_for_group(td_ev, td_ev->dest.group);
                if(minpathport != next_port){
                    rtr_ev->val_route_pos += 2;
                }
            }
            else{
                assert(current_vc == start_vc + 3);
                assert(port >= params.p && port < params.p + params.a - 1); //local port

                next_vc = start_vc + 3;
                next_port = choose_port_for_group(td_ev, td_ev->dest.group);
            }
        }
    }

    assert(next_vc >=0);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

void topo_dragonfly::q_adaptive(int port, int vc, internal_router_event* ev){

    int next_vc = -1;
    int next_port = -1;
    int min_path_port = -1;
    bool first_mid_group_rtr=false;

    const int vn = ev->getVN();
    RouteAlgo qalg =  vns[vn].algorithm;
    assert( qalg == Q1);

    const int start_vc = vns[vn].start_vc;
    const int end_vc = start_vc + vns[vn].num_vcs;

    if(qalg == Q1 && src_mid_group_q){
        assert(vns[vn].num_vcs == 5);
    }

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);
    RtrEvent* rtr_ev = td_ev->getEncapsulatedEvent();
    int link_hops = rtr_ev->getNumHops();

    int current_vc = td_ev->getVC();
    assert(current_vc == vc);

    //init estimate to -1. Force qrouting to set propoer new estimate value in later steps
    td_ev->new_estimate = -1;

    //find minimal routing path port
    if ( td_ev->dest.group != group_id ) {
        min_path_port = choose_port_for_group(td_ev, td_ev->dest.group);
    }
    else if ( td_ev->dest.router != router_id ) {
        min_path_port = port_for_router(td_ev->dest.router);
    }
    else {
        min_path_port = td_ev->dest.host;
    }
    assert(min_path_port>=0);

    // Increment VC per hop
    assert(link_hops <= vns[vn].num_vcs);
    assert(current_vc == start_vc || current_vc == start_vc + link_hops -1 );
    // when forwarding to compute node, no vc increment 
    if(min_path_port >= params.p){
        next_vc = start_vc + link_hops;
    }
    else{
        // in dest group, either reached it or message src group == dest group
        next_vc = (link_hops-1 >= 0) ?  (start_vc + link_hops-1) : (start_vc + 0);
    }


    // msg in dest_group, min route to dest
    if(td_ev->dest.group == group_id) {
        if(qtable_row_type == "g" || qtable_row_type == "destG_srcN" ){
            td_ev->new_estimate = 0;
        } 
        else {
            assert(qtable_row_type == "r" || qtable_row_type == "n" );
            uint32_t dest_rtr = td_ev->dest.group * params.a + td_ev->dest.router;
            assert(dest_rtr == td_ev->dest.router_global);
            
            if(td_ev->dest.router_global == router_id_global) {
                td_ev->new_estimate = 0;
            }
            else {
                int row_index = qtable_get_row_idx(td_ev);
                auto best_pair = qtable_get_min_port(port, min_path_port, row_index, td_ev, params.p, params.k);
                int best_port = best_pair.first; 
                int64_t best_est = best_pair.second; 
                td_ev->new_estimate = best_est;
            }
        }
        
        uint64_t time2save = getCurrentSimTimeNano();
        if( save_qtable && 
            // td_ev->getTraceID() >= rtr_ev->getSpecialIndex()
            (time2save >= save_qtable_time) ){
            save_qtable_toFile();
            save_qtable = false;
        }

        next_port = min_path_port;
        
        assert((next_vc >= start_vc) && (next_vc < end_vc));
        assert(next_port >= 0);
        assert(td_ev->new_estimate >= 0);

        td_ev->setVC(next_vc);
        td_ev->setNextPort(next_port);   
        return;
    }

    //In Src or Mid group, Not in Dest group
    int row_index = qtable_get_row_idx(td_ev);
    auto best_pair =  qtable_get_min_port(port, min_path_port, row_index, td_ev, params.p, params.k);
    int best_port = best_pair.first;
    int64_t best_est = best_pair.second;
    // qrouting: update qvalue using best Est no matter which action is actually taken
    td_ev->new_estimate = best_est;  

    //need to set next_port explicitly
    assert(next_port == -1);

    //Set non-min routed counter, val_route_pos, 
    //only meaninful for q-adp routing: 
    //  0 -- minimally forwarded, 
    //  1 -- non_min decision made at src rtr, 
    //  2 -- non_min decision made at 1st rtr in intermediate group
    if(group_id == td_ev->src_group){
        if(src_group_q){
            assert(!rtr_ev->getAdpRouted());
        }    
    }
    else if (group_id != td_ev->src_group && group_id != td_ev->dest.group){
        if(!rtr_ev->getAdpRouted()){
            //1st time reach intermediate group
            assert(is_port_global(port));
            assert(rtr_ev->val_route_pos == 0);

            first_mid_group_rtr=true;
            rtr_ev->setAdpRouted();
            rtr_ev->val_route_pos = 1;
            //record midgroup info
            assert(rtr_ev->midgroup == -1);
            rtr_ev->midgroup = group_id;
        }
    }
    else{
        output.fatal(CALL_INFO, -1, "ERROR DFrtr %d: Msg not in Src, intermediate, Dest group\n", router_id_global);
    }

    //  1 -- forward minimally, 
    //  0 -- read q-table,
    int go_min = 0;
    if(src_group_q){
        if(link_hops >= max_hops || group_id != td_ev->src_group){
            go_min = 1;
        }
    }
    else if (src_mid_group_q){   
        if( qalg == Q1 ){
            if(group_id == td_ev->src_group && link_hops>0){
                go_min = 1;
            }

            if(group_id != td_ev->src_group && group_id != td_ev->dest.group && link_hops > 1){
                assert(!first_mid_group_rtr);
                go_min = 1;
            }
        }
        else{
            assert( 0 );
        }
    }
    else{
        if (link_hops >= max_hops ) {
            go_min = 1;
        }
    }

    if(go_min){
        next_port = min_path_port;
        assert((next_vc >= start_vc) && (next_vc < end_vc));
        assert(next_port >= 0);
        assert(td_ev->new_estimate >= 0);

        td_ev->setVC(next_vc);
        td_ev->setNextPort(next_port);   
        return;
    }

    // let rl choose port. We can arrive here only when following cases are satisfied
    if(src_group_q){
        assert(!rtr_ev->getAdpRouted());
        assert(group_id == td_ev->src_group);
        assert(link_hops < max_hops);
    }
    else if (src_mid_group_q){
        if(group_id == td_ev->src_group){
            assert(port<params.p);
            assert(link_hops == 0);
        }
        else{
            assert(port >= params.p + params.a -1);
            assert(link_hops < 3);
            assert(first_mid_group_rtr);
        }
    }
    else{
        assert(link_hops < max_hops);
    }

    double q_thld = -1.0;
    if(first_mid_group_rtr && src_mid_group_q){
        assert(next_port==-1); // next port is not decided yet
        assert(min_path_port >= params.p); //it must not be a host port
        if(min_path_port < params.p + params.a -1){
            uint32_t tmp_router = params.a+10;
            do {
               tmp_router = rng_q->generateNextUInt32() % params.a;
            }
            while (tmp_router == router_id);
            assert(tmp_router>=0 && tmp_router<params.a);

            best_port = port_for_router(tmp_router);
            assert(best_port>=params.p && best_port<params.a + params.p -1);
            best_est = qtable_get_est(best_port, td_ev, row_index); 

            q_thld = q_threshold2;
        }
        else{
            // has direct global link to dest group, min routing to dest
            next_port = min_path_port;            
            assert(next_vc >= start_vc && next_vc < end_vc);
            assert(next_port >= 0);
            assert(td_ev->new_estimate >= 0);

            td_ev->setVC(next_vc);
            td_ev->setNextPort(next_port);   
            return;
        }
    }
    else{
        q_thld = q_threshold1;
    }

    assert(q_thld>=0);
    if(best_port != min_path_port) {
        int64_t min_path_est = qtable_get_est(min_path_port, td_ev, row_index); 
        double est_diff = double(min_path_est - best_est) / double (min_path_est);

        // q-adp routing change best est as a random port est
        if(est_diff < 0){
            assert(first_mid_group_rtr && src_mid_group_q );
        }

        if(est_diff < q_thld){
            next_port = min_path_port;
        }
        else{
            if(first_mid_group_rtr){
                assert(rtr_ev->val_route_pos == 1);
                rtr_ev->val_route_pos = 2;
            }
            next_port = best_port;
        }
    }
    else{
        next_port = best_port;
    }

    assert(next_port>=params.p && next_port < params.k);

    //exploration
    if(qtable_bcast == NOBCAST){
        double prob = rng->nextUniform();
        if(prob < epsilon) {
            int random_port = ( rng->generateNextUInt32() % (params.k - params.p) ) + params.p;
            next_port = random_port;
            assert(next_port>=params.p && next_port < params.k);
        }
    } 

    // in qlearning: action policy != value policy
    assert(next_vc >= start_vc && next_vc < end_vc);
    assert(next_port >= 0);
    assert(td_ev->new_estimate >= 0);
 
    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);   
    return; 
}

void
topo_dragonfly::setOutput2NbrCreditArray(int const* array, int vcs)
{
    output_2nbr_credits = array;
    num_vcs = vcs;
}

void
topo_dragonfly::setOutputUsedCreditArray(int const* array, int vcs)
{
    output_used_credits = array;
    num_vcs = vcs;
}

// ATTENTION, destgroup number is [0 -- g-1], need to exclude the self group id, thus dest_group -=1 if dest_group >= group_id
// port index [0 -- number_local + global -1], thus return value needs to add number of host as the correct port
void
topo_dragonfly::setQtable()
{   
    if(!useQrouting) return;
    assert(link_latency_global > 0);
    assert(link_latency_local > 0);

    // uint32_t row = params.g - 1;
    uint32_t row = qtable_rows;
    uint32_t col = params.k - params.p;

    assert( (params.a-1 + params.h) == col );

    qtable = new int64_t[row * col];
    // qtable.resize(row * col);

    if(pathToQtableFile.empty()){
        if(router_id_global == 0) printf("\n\nRtr 0 setting qtable\n");

        if(qtable_row_type == "g"){     
            for(int g=0; g< row; g++){
                uint32_t dest_g = g;
                if (dest_g >= group_id) dest_g++;
                uint32_t port_to_destg = port_for_group(dest_g, 0);
                for(int p=0; p<col; p++){
                    if( p >= params.a-1 && p + params.p == port_to_destg)
                        qtable[g*col+p] = link_latency_global;
                    else
                        qtable[g*col+p] = link_latency_global + link_latency_local;
                }
            }
        }
        else if (qtable_row_type == "r"){
            for(int r=0; r< row; r++){
                
                int dest_r = r;
                if(dest_r >= router_id_global) dest_r++;
                int dest_group = dest_r/params.a;

                for(int p=0; p<col; p++){
                    if(dest_group == group_id){
                        if(p >= params.a-1){ // global port
                            qtable[r*col+p] = 2*link_latency_global;
                        } else {
                            qtable[r*col+p] = 2*link_latency_local;
                        }

                    }
                    else{
                        qtable[r*col+p] = link_latency_global + 2*link_latency_local;
                    }
                }
            }
        }
        else if (qtable_row_type == "n") {
            for(int n=0; n<row; n++){
                for(int p=0; p<col; p++){
                    qtable[n*col+p] = link_latency_global;
                }
            }
        }

        else if (qtable_row_type == "destG_srcN") {

            for(int g=0; g< params.g -1; g++){
                uint32_t dest_g = g;
                if (dest_g >= group_id) dest_g++;
                uint32_t port_to_destg = port_for_group(dest_g, 0);
                for(int n=0; n<params.p; n++){ 

                    int row_idx = (g * params.p) + n;
                    for(int p=0; p<col; p++){
                        assert(row_idx*col+p<(row*col));
                        if( p >= params.a-1 && p + params.p == port_to_destg)
                            qtable[row_idx*col+p] = link_latency_global;
                        else
                            qtable[row_idx*col+p] = link_latency_global + link_latency_local;
                    }
                }
            }
        }

        else {
            output.fatal(CALL_INFO,-1,"DFtopology: unknown qtable type %s\n", qtable_row_type.c_str());
        }

    } else {
    // load an existing qtable            
        //io load qtable
        std::ifstream file(pathToQtableFile + qtablefile);
        int count = -1; 
        std::string line;
        while (std::getline(file, line))
        {   
            if(count == -1){
                if (router_id_global == 0)
                    printf("Rtr %d loading qtable from %s w/ following info:\n%s\n", router_id_global, (pathToQtableFile + qtablefile).c_str(), line.c_str() );
            }
            else{
                qtable[count] = (int64_t) std::stoi(line);
            }
            count++;
        }
        
        assert(count == (row * col));
    }

    if(router_id_global == 0){
        dumpQtable();
    }
}

void topo_dragonfly::dumpQtable(){
    // uint32_t row = params.g - 1;
    uint32_t row = qtable_rows;
    uint32_t col = params.k - params.p;
    printf("\nRtr %d, Qtable:\n", router_id_global);

    if(qtable_row_type == "g"){
        for ( uint32_t i = 0; i < row ; i++ ) {
            if(i>=group_id) 
                printf("To group %d: ", i+1);
            else 
                printf("To group %d: ", i);
            
            for( uint32_t j = 0; j<col; j++){
                printf("(port %d: %ldns) ", j+params.p, qtable[i*col + j]);
            }     
            printf("\n");
        }

    }
    else if(qtable_row_type == "r"){
        for ( uint32_t i = 0; i < row ; i++ ) {
            if(i>=router_id_global) 
                printf("To Rtr %d: ", i+1);
            else 
                printf("To Rtr %d: ", i);
            
            for( uint32_t j = 0; j<col; j++){
                printf("(port %d: %ldns) ", j+params.p, qtable[i*col + j]);
            }     
            printf("\n");
        }
    }
    else if(qtable_row_type == "n"){
        for ( uint32_t i = 0; i < row ; i++ ) {
            printf("To Node %d: ", i);

            for( uint32_t j = 0; j<col; j++){
                printf("(port %d: %ldns) ", j+params.p, qtable[i*col + j]);
            }     
            printf("\n");
        }
    }
    else if(qtable_row_type == "destG_srcN"){

        for(int g=0; g< params.g -1; g++){
            uint32_t dest_g = g;
            if (dest_g >= group_id) dest_g++;

            printf("To Group %d\n", dest_g);

            for(int n=0; n<params.p; n++){ 

                printf("From node %d: ", n);
                int row_idx = (g * params.p) + n;
                for(int p=0; p<col; p++){
                    printf("(port %d: %ldns) ", p+params.p, qtable[row_idx*col+p]);
                }
                printf("\n");
            }
        }

    }

    else {
        output.fatal(CALL_INFO,-1,"DF dumpQtable unknown qtable type %s\n", qtable_row_type.c_str());
    }

    printf("\n\n\n");
}

int topo_dragonfly::qtable_get_row_idx(topo_dragonfly_event *td_ev){
    dgnflyAddr dest = td_ev->dest;
    int target_id =-1;
    if(qtable_row_type == "g"){
        assert(dest.group != group_id);
        target_id = dest.group;
        if (dest.group > group_id) target_id--;
    } 
    else if(qtable_row_type == "r")
    {
        
        assert(dest.router_global != router_id_global);
        target_id = dest.router_global;
        if (dest.router_global > router_id_global) target_id--;
    }
    else if(qtable_row_type == "n"){
        target_id = td_ev->getDest();
    }
    else if(qtable_row_type == "destG_srcN"){
        assert(dest.group != group_id);
        int src_node = td_ev->getSrc();
        src_node %= params.p;

        int destg = dest.group;
        if (dest.group > group_id) destg--;
        target_id = destg *  params.p + src_node;
    }
    else{
        assert(0);
    }

    assert(target_id>=0 && target_id < qtable_rows);
    return target_id;
}

int64_t topo_dragonfly::qtable_get_est(int target_port, topo_dragonfly_event* td_ev, int row_idx){

    assert(target_port>=params.p && target_port < params.k);

    const int vn = td_ev->getVN();
    const RouteAlgo ralg = vns[vn].algorithm;
    const int vc_start = vns[vn].start_vc;
    const int vc_end = vc_start + vns[vn].num_vcs;
    
    uint32_t col = params.k - params.p;
    int port_number_intable = target_port - params.p;

    int64_t min_est = qtable[row_idx * col + port_number_intable];

    return min_est;
}

std::pair<int,int64_t> topo_dragonfly::qtable_get_min_port(int port, int min_path_port, int row_idx, topo_dragonfly_event* td_ev, uint32_t port_start, uint32_t port_end)
{   
    assert(port_start >= params.p);
    assert(port_end <= params.k);
    assert(port_end > port_start);
    
    // so far consider all ports
    assert(port_start == params.p);
    assert(port_end == params.k);

    const int vn = td_ev->getVN();
    const RouteAlgo ralg = vns[vn].algorithm;
    assert(ralg == Q1);

    const int vc_start = vns[vn].start_vc;
    const int vc_end = vc_start + vns[vn].num_vcs;

    uint32_t col = params.k - params.p;
    int64_t min_est = INT64_MAX;
    int64_t min_est_w_buffer = INT64_MAX;
    int best_port = -1; 
    std::vector<int> index;
    
    for(int check_port = port_start - params.p ; check_port <  port_end - params.p; check_port++){
        int64_t est = qtable[row_idx * col + check_port];

        if(ralg == Q1){
            if(est < min_est) index.clear(); 
            if(est <= min_est) {
                index.push_back(check_port + params.p);
                min_est = est;
            }
        }  
        else{
            assert(0);
        }      
    }
    
    int port_min = 0;
    assert( port_min < index.size());
    
    assert(index[port_min]>=params.p && index[port_min]<params.k);
    assert( min_est <= min_est_w_buffer );

    std::pair<int, int64_t> ret_pair = std::make_pair(index[port_min], min_est);    
    return ret_pair;
}

//node id to group
uint32_t topo_dragonfly::idToGroup(int id){
    uint32_t hosts_per_group = params.p * params.a;
    uint32_t group = id / hosts_per_group;
    return group;
}
//node id to router
uint32_t topo_dragonfly::idToRouter(int id){
    uint32_t hosts_per_group = params.p * params.a;
    uint32_t group = id / hosts_per_group;
    uint32_t router = (id % hosts_per_group) / params.p; //relative to group

    router =  params.a * group + router;
    return router;
}

qtable_event* topo_dragonfly::create_qtable_event(internal_router_event* ev, int port_number, bool bcast){
    
    if(!useQrouting) return NULL;

    RouteAlgo ralg = vns[ev->getVN()].algorithm;
    if(ralg!=Q1){
        return NULL;
    }

    assert(port_number>=params.p && port_number<params.k); //is not a host port, assured by portcontrol
    uint32_t dest_g = idToGroup(ev->getDest());
    uint32_t dest_rtr = idToRouter(ev->getDest());

    //remember to del once receive
    qtable_event* qtable_ev = NULL;

    if(qtable_row_type == "g") {
        if(dest_g == group_id){
            if(port_number < params.p + params.a -1 && (ralg==Q1)  ) {
                // q1 routing local link no need to generate qevent
                return NULL; 
            } 
            //q1 global link , q_event only with queueing time
            qtable_ev = new qtable_event(dest_g, ev->queueing_time, 0, ev->getVN());
        }
        else{
            assert(ev->new_estimate != 0);
            qtable_ev = new qtable_event(dest_g, ev->queueing_time, ev->new_estimate, ev->getVN());
        } 
    } 
    else if (qtable_row_type == "r")
    {   
        if(dest_rtr == router_id_global){
            qtable_ev = new qtable_event(dest_rtr, ev->queueing_time, 0, ev->getVN());
        }
        else{
            assert(ev->new_estimate != 0);
            qtable_ev = new qtable_event(dest_rtr, ev->queueing_time, ev->new_estimate, ev->getVN());
        }
    }
    else if (qtable_row_type == "n")
    {   
        if(dest_rtr == router_id_global){
            qtable_ev = new qtable_event(ev->getDest(), ev->queueing_time, 0, ev->getVN());
        }
        else{
            assert(ev->new_estimate != 0);
            qtable_ev = new qtable_event(ev->getDest(), ev->queueing_time, ev->new_estimate, ev->getVN());
        }
    }
    else if (qtable_row_type == "destG_srcN"){

        int src_nid = ev->getSrc();
        if(dest_g == group_id){
            if(port_number < params.p + params.a -1 && (ralg==Q1)) {
                // q1 routing local link no need to generate qevent
                return NULL; 
            } 
            //q1 global link, q_event only with queueing time
            qtable_ev = new qtable_event(dest_g, ev->queueing_time, 0, ev->getVN(), src_nid);
        }
        else{
            assert(ev->new_estimate != 0);
            qtable_ev = new qtable_event(dest_g, ev->queueing_time, ev->new_estimate, ev->getVN(), src_nid);
        } 
    }
    else{
        assert(0);
    }

    assert(qtable_ev);
    return qtable_ev;
}

void topo_dragonfly::updateQtable(qtable_event* qe, int port){
    assert(qe->queueing_time>=0);
    assert(qe->new_estimate>=0);
    // uint32_t destg = qe->dest_group;
    uint32_t target_row = qe->target_row_in_table;
    RouteAlgo ralg = vns[qe->appvn].algorithm;

    const int target_port = port - params.p;
    uint32_t col = params.k - params.p;

    if(!qe->qBcast){
        // This is a normal qtable_event 
        if(qtable_row_type == "g") {
            uint32_t destg = target_row;
            if(destg == group_id){
                if(ralg == Q1) {
                    //no need to update qtable for my own group
                    output.fatal(CALL_INFO,-1,"q1 no same group qtable update should be guaranteed by portcontrol and create_qtable_event()");
                    // return;
                } else {
                    assert(0);
                }
            } else {
                if(destg > group_id){
                    destg--;
                } 
                uint32_t target_index = destg * col + target_port;
                int64_t newq=-10;
                int64_t deltaq=0;

                if(ralg == Q1) {
                    deltaq = qe->new_estimate + qe->queueing_time - qtable[target_index];
                } else {
                    assert(0);
                }

                //Hysteretic q-learning
                if(deltaq < 0) {
                    // This is a positive case, that Est. get smaller. 
                    // Update with lr1
                    newq = qtable[target_index] + learning_rate * deltaq;
                } 
                else {
                    // This is a negative case, that penalty arrives. 
                    // Update with lr2
                    newq = qtable[target_index] + learning_rate2 * deltaq;
                }
                assert(newq>=0); 
                qtable[target_index] = newq;
            }
        } 
        else if(qtable_row_type == "r")
        {
            uint32_t dest_rtr = target_row;
            assert(dest_rtr != router_id_global);

            if(dest_rtr > router_id_global){
                dest_rtr--;
            } 

            uint32_t target_index = dest_rtr * col + target_port;
            int64_t newq=-10;
            int64_t deltaq=0;

            if(ralg == Q1) {
                // newq = qtable[target_index] + learning_rate * ( qe->new_estimate + qe->queueing_time - qtable[target_index]);
                deltaq = qe->new_estimate + qe->queueing_time - qtable[target_index];

            } else {
                assert(0);
            }  

            //Hysteretic q-learning
            if(deltaq < 0) {
                // This is a positive case, that Est. get smaller. 
                // Update with lr1
                newq = qtable[target_index] + learning_rate * deltaq;
            } 
            else {
                // This is a negative case, that penalty arrives. 
                // Update with lr2
                newq = qtable[target_index] + learning_rate2 * deltaq;
            }
            assert(newq>=0); 
            qtable[target_index] = newq;
        }

        else if(qtable_row_type == "n"){
            assert(qtable_row_type == "n");
            
            uint32_t target_index = target_row * col + target_port;
            
            int64_t newq=-10;
            int64_t deltaq=0;
            uint32_t dest_rtr = target_row / params.p;

            if(dest_rtr == router_id_global){
                printf("DFrtr %d, updating for node %u, through port %d \n", router_id_global, target_row, port );
            }
            assert(dest_rtr != router_id_global);

            if(ralg == Q1) {
                deltaq = qe->new_estimate + qe->queueing_time - qtable[target_index];

            } else {
                assert(0);
            }  

            if(deltaq < 0) {
                // This is a positive case, that Est. get smaller. 
                // Update with lr1
                newq = qtable[target_index] + learning_rate * deltaq;
            } 
            else {
                // This is a negative case, that penalty arrives. 
                // Update with lr2
                newq = qtable[target_index] + learning_rate2 * deltaq;
            }
            assert(newq>=0); 
            qtable[target_index] = newq;
        }

        else if(qtable_row_type == "destG_srcN"){
            assert(qe->src_nid>=0);
            assert(ralg == Q1);

            uint32_t destg = target_row;

            if(destg == group_id){
                //no need to update qtable for my own group
                output.fatal(CALL_INFO,-1,"q1 no same group qtable update should be guaranteed by portcontrol and create_qtable_event()");

            } else {

                if(destg > group_id){
                    destg--;
                } 

                int src_node_on_rtr =  qe->src_nid % params.p;

                uint32_t target_index = (destg * params.p + src_node_on_rtr) * col + target_port;
                int64_t newq=-10;
                int64_t deltaq=0;
                deltaq = qe->new_estimate + qe->queueing_time - qtable[target_index];
                
                //Hysteretic q-learning
                if(deltaq < 0) {
                    // This is a positive case, that Est. get smaller. 
                    // Update with lr1
                    newq = qtable[target_index] + learning_rate * deltaq;

                } 
                else {
                    newq = qtable[target_index] + learning_rate2 * deltaq;
                }
                assert(newq>=0); 
                qtable[target_index] = newq;
            }
        }
        else{
            assert(0);
        }
    } else { 
        printf("BCASTING IS UNDER DEVELOP\n");
    }
}

void topo_dragonfly::save_qtable_toFile(){
    // uint32_t row = params.g - 1;
    uint32_t row = qtable_rows;
    uint32_t col = params.k - params.p;

    char *dir_mk = const_cast<char*>(qtableFileDir.c_str());

    mkdir(dir_mk ,0755);

    uint64_t time2save = getCurrentSimTimeNano();

    std::string table_s = "Table saved at time ";
    table_s += std::to_string(time2save/1000);
    table_s += "us \n";

    for ( uint32_t i = 0; i < row * col; i++ ) {
        // out2file.output("%lu\n", qtable[i]);
        table_s += std::to_string(qtable[i]);
        table_s += "\n";
    }
    char* table_sc = const_cast<char*>(table_s.c_str());
    
    // out2file.output(table_sc);

    std::ofstream myfile;
    myfile.open(qtablefile);
    myfile << table_sc;
    myfile.close();

    if(router_id_global == 0){
        printf("Rtr 0 qtable saved to file %s, time %ld (us)\n", qtablefile.c_str(), getCurrentSimTimeNano()/1000);
        dumpQtable();
    }
}

void topo_dragonfly::set_parent(Router* router){
    parent = router;
}

int topo_dragonfly::get_port_remote_group_id(int outport){
    //is must be a global port
    assert(outport>=params.p+ params.a-1);
    assert(outport<params.k);
    //only work with absolute so far
    assert(global_route_mode == ABSOLUTE);

    int remote_gid = router_id * params.h +  (outport - (params.p+ params.a-1) ) ;

    if(remote_gid >= group_id)
        remote_gid++;

    assert(remote_gid<params.g);

    return remote_gid;
}

bool topo_dragonfly::perid_funct_handler(Cycle_t cycle){

    int target_rtr = 2;

    assert(router_id_global == target_rtr );


    printf("DFrtr %d @ %luns, xbar crt, output queue, output crt, used crt:\n", router_id_global, getCurrentSimTimeNano());

    for(int prt = 0; prt < params.k; prt++){

        printf("Port %d: ", prt);

        for(int vc = 0; vc<num_vcs; vc++){
            printf("vc_%d[%d %d %d %d] ", vc,
                output_credits[prt * num_vcs + vc],
                output_queue_lengths[prt * num_vcs + vc],
                output_2nbr_credits[prt * num_vcs + vc],
                output_used_credits[prt * num_vcs + vc] );
        }

        printf("\n");
    }

    return false;
}

void topo_dragonfly::set_output_queue_length(int64_t obs){
        out_queue_len = obs;
}