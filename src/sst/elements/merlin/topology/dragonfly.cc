// Copyright 2009-2020 NTESS. Under the terms
// of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.
// 
// Copyright (c) 2009-2020, NTESS
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
#include <sst/core/sharedRegion.h>
#include "sst/core/rng/xorshift.h"

#include "dragonfly.h"

#include <stdlib.h>
#include <sstream>

using namespace SST::Merlin;



void
RouteToGroup::init(SharedRegion* sr, size_t g, size_t r)
{
    region = sr;
    data = sr->getPtr<const RouterPortPair*>();
    groups = g;
    routes = r;
    
}

const RouterPortPair&
RouteToGroup::getRouterPortPair(int group, int route_number)
{
    // data = static_cast<RouterPortPair*>(region->getRawPtr());
    return data[group*routes + route_number];
}

void
RouteToGroup::setRouterPortPair(int group, int route_number, const RouterPortPair& pair) {
    // output.output("%d, %d, %d, %d\n",group,route_number,pair.router,pair.port);
    region->modifyArray(group*routes+route_number,pair);
}


topo_dragonfly::topo_dragonfly(ComponentId_t cid, Params &p, int num_ports, int rtr_id, int num_vns) :
    Topology(cid),
    num_vns(num_vns),
    router_id_global(rtr_id)
{

    qtablefile = qtableFileDir+"rtr"+ std::to_string(router_id_global) +"_qtable"; 
    out2file.init("", 0, 0, Output::FILE, qtablefile);

    bool found;
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

    load_qtable = p.find<bool>("load_qtable", false, found);
    max_hops = p.find<int>("max_hops", 1, found);
    assert(max_hops >= 0);

    src_group_q = p.find<bool>("src_group_q", false, found);
    src_mid_group_q = p.find<bool>("src_mid_group_q", true, found);

    pathToQtableFile = p.find<std::string>("pathToQtableFile","",found);
    if(load_qtable && pathToQtableFile.empty()){
        output.fatal(CALL_INFO, -1, "Need to specify the qtable file for loading\n");

    }

    params.p = p.find<uint32_t>("hosts_per_router");
    params.a = p.find<uint32_t>("routers_per_group");
    params.k = num_ports;
    params.h = p.find<uint32_t>("intergroup_per_router");
    params.g = p.find<uint32_t>("num_groups");
    params.n = p.find<uint32_t>("intergroup_links");

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
    
    adaptive_threshold = p.find<double>("adaptive_threshold",2.0, found);

    // Get the global link map
    std::vector<int64_t> global_link_map;
    p.find_array<int64_t>("global_link_map", global_link_map);
    
    // Get a shared region
    SharedRegion* sr = Simulation::getSharedRegionManager()->getGlobalSharedRegion("group_to_global_port",
                                                                                  ((params.g-1) * params.n) * sizeof(RouterPortPair),
                                                                                   new SharedRegionMerger());
    // Set up the RouteToGroup object
    group_to_global_port.init(sr, params.g, params.n);

    // Fill in the shared region using the RouteToGroupObject (if
    // vector for param dragonfly:global_link_map is empty, then
    // nothing will be intialized.
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
        group_to_global_port.setRouterPortPair(group, route_num, rpp);
    }

    
    // Publish the shared region to make sure everyone has the data.
    sr->publish();


    bool show_adprouting = false;
    useQrouting = 0;

    // Setup the routing algorithms
    int curr_vc = 0;
    for ( int i = 0; i < num_vns; ++i ) {
        vns[i].start_vc = curr_vc;
        if ( !vn_route_algos[i].compare("VALn") ) {
            if ( params.g <= 2 ) {
                /* 2 or less groups... no point in valiant */
                vns[i].algorithm = MINIMAL;
                vns[i].num_vcs = 2;
            } else {
                //VALn routing
                vns[i].algorithm = VALIANT;
                vns[i].num_vcs = 4;
            }
        }
        // else if ( !vn_route_algos[i].compare("adaptive-local") ) {
        else if ( !vn_route_algos[i].compare("ugal-g") ) {    
            vns[i].algorithm = UGAL_3VC;
            vns[i].num_vcs = 3;
            show_adprouting = true;
        }
        else if ( !vn_route_algos[i].compare("ugal-n") ) {    
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
        else if ( !vn_route_algos[i].compare("q-adaptive") ) {
            vns[i].algorithm = Q1;
            if(src_group_q){
                assert(!src_mid_group_q);
                vns[i].num_vcs = max_hops + 3;
            }
            else if (src_mid_group_q) {
                assert(!src_group_q);
                assert(max_hops == 1);
                vns[i].num_vcs = 5;
            }
            else{
                vns[i].num_vcs = max_hops + 3;
            }
            useQrouting = 1;
        }
        else if ( !vn_route_algos[i].compare("q-adaptive2") ) {
            output.output("WARNING: q-adaptive2 is experimental, use with caution\n");
            vns[i].algorithm = Q2;
            vns[i].num_vcs = max_hops + 3;
            useQrouting = 2;
        }
        else {
            fatal(CALL_INFO_LONG,1,"ERROR: Unknown routing algorithm specified: %s\n",vn_route_algos[i].c_str());
        }
        curr_vc += vns[i].num_vcs;
    }
    
    group_id = rtr_id / params.a;
    router_id = rtr_id % params.a;

    rng = new RNG::XORShiftRNG(rtr_id+1);

    output.verbose(CALL_INFO, 1, 1, "%u:%u:  ID: %u   Params:  p = %u  a = %u  k = %u  h = %u  g = %u\n",
            group_id, router_id, rtr_id, params.p, params.a, params.k, params.h, params.g);


    q_threshold1 = p.find<double>("q_threshold1", 0.0, found);
    assert(q_threshold1>=0);
    q_threshold2 = p.find<double>("q_threshold2", 0.0, found);
    assert(q_threshold2>=0);

    qtable_row_type = p.find<std::string>("qtable_row_type", "destG_srcN");
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

    std::string qtable_bcast_s = p.find<std::string>("qtable_bcast", "nobcast", found);
    qbcastThsld = p.find<float>("qbcastThsld", 10, found);
    UnitAlgebra qbcastPerid_ua = p.find<UnitAlgebra>("qbcastPerid","1000ns");
    if ( !qbcastPerid_ua.hasUnits("s") ) {
        output.fatal(CALL_INFO,-1,"bcast period must specified in seconds");
    }
    qbcastPerid = (qbcastPerid_ua / UnitAlgebra("1ns")).getRoundedValue();
    bcast_itr = 0;

    if (useQrouting > 0){
        setQtable();
        for(int i = params.a-1; i<params.k - params.p; i++){
            //global port index - num_host_port
            global_port_idx.push_back(i);
        }
    }
    if (useQrouting == 2){ //used by q2
        set_t2nbrTable();
    }

    if (qtable_bcast_s == "nobcast"){
        qtable_bcast = NOBCAST;
    } 
    else if (qtable_bcast_s == "perid") {

        //Perid Bcasting designed to work with q1 rting, g or r type qtable for now
        assert( qtable_row_type == "r" || qtable_row_type == "g" );
        assert(useQrouting == 1);
        qtable_bcast = PERID; 
        set_qBcastTable(); //Both bcast thld and bcast perid need this
        set_tFromNbrTable();

        qPeridClock_handler = new Clock::Handler<topo_dragonfly>(this,&topo_dragonfly::perid_qBcast_handler);
        qPerid_tc = registerClock( qbcastPerid_ua, qPeridClock_handler, false);

    }
    else if (qtable_bcast_s == "thld") { 
        //Thld Bcasting designed to work with q2 rting, g or r type qtable for now

        assert(qtable_row_type == "r" || qtable_row_type == "g");
        assert(useQrouting == 2);
        qtable_bcast = THLD;
        set_qBcastTable();
    }
    else{
        output.fatal(CALL_INFO, -1, "Unknown Qtable Bcast: %s\n", qtable_bcast_s.c_str());
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

    if(rtr_id == 0 ){
        output.output("\n---------");
        output.output("Dragonfly Topology: %d nodes", params.p * params.a * params.g);
        output.output("---------\n");
        output.output("host/router %d\trouter/group %d\tnum_group %d\n", 
            params.p, params.a, params.g);

        output.output("num_ports %d\tglobal_link %d\tglobal_link/router %d\n", 
            params.k, params.n, params.h);
        
        output.output("link latency:local %ld\tglobal %ld\n", link_latency_local, link_latency_global);

        output.output("total VC num: %d\n", curr_vc);
        // output.output("periodic func %d\n", prid_func);

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
            output.output("saveqtable %d\tsaveqtable_time(us) %d\n", save_qtable, save_qtable_time/1000);
            // output.output("qtable Bcast %d\n", qtable_bcast);
            // output.output("global ports are ( - num_hosts): ");
            // for(auto x : global_port_idx){
            //     output.output("%d ",x);
            // }
            // output.output("\n");
            // if(qtable_bcast == PERID){
            //     output.output("Bcast period %d\t Bcast clock %s\n", qbcastPerid, qbcastPerid_ua.toStringBestSI().c_str());
            // }
            // if(qtable_bcast == THLD ){
            //     output.output("Bcast threshold %.4f\n", qbcastThsld);
            // }
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
    if (useQrouting ==2){
        delete[] t2nbrTable;
    }
}

void topo_dragonfly::route_nonadaptive(int port, int vc, internal_router_event* ev)
{
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
        td_ev->setVC(vc+1);
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

void topo_dragonfly::route_adaptive_local(int port, int vc, internal_router_event* ev)
{
    int vn = ev->getVN();
    if ( vns[vn].algorithm != ADAPTIVE_LOCAL ) return;

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
    }
    else { // Use direct route
        td_ev->dest.mid_group = td_ev->dest.group;
        td_ev->setNextPort(direct_route_port);
        td_ev->global_slice = direct_slice;
    }

}

int topo_dragonfly::choose_port_for_group(topo_dragonfly_event* td_ev, int dest_g){
    assert(dest_g>=0);
    int direct_slice1 = td_ev->global_slice_shadow;
    int direct_slice2 = (td_ev->global_slice_shadow + 1) % params.n;
    int direct_route_port1 = port_for_group(dest_g, direct_slice1, 0 );
    int direct_route_port2 = port_for_group(dest_g, direct_slice2, 1 );

    int direct_crt1 = 0;
    int direct_crt2 = 0;
    int port_to_g = -1;

    for(int i=0; i<num_vcs; i++){
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

//VALn routing
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

    assert(next_vc >=0);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

std::pair<int,int> topo_dragonfly::adp_select_port(topo_dragonfly_event *td_ev, int port, int vc_start, int vc_end, std::vector<int> vc_min, std::vector<int> vc_nonmin, std::string debuginfo)
{
    assert(vc_start >=0);
    assert(vc_end >= vc_start);
    assert(vc_min.size()>0);
    assert(vc_nonmin.size()>0);

    //as only 1 VN is used so far
    assert(vc_start == 0);
    assert(vc_end == num_vcs);

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
    const int ugla_vcs = vns[vn].num_vcs;
    assert(ugla_vcs == 3);

    assert(ugla_vcs == num_vcs); // only using 1 vn so far

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

            std::vector<int> vc_min {1};
            std::vector<int> vc_nonmin {0};

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, num_vcs, vc_min, vc_nonmin);

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

    assert(next_vc >=0);
    assert(next_port>=0);

    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

void topo_dragonfly::route_ugal_4vc(int port, int vc, internal_router_event* ev)
{   
    const int vn = ev->getVN();
    assert( vns[vn].algorithm == UGAL_4VC );
    const int start_vc = vns[vn].start_vc;

    const int ugla_vcs = vns[vn].num_vcs;
    assert(ugla_vcs == 4);

    assert(ugla_vcs == num_vcs); // only using 1 vn so far

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

            std::vector<int> vc_min {2};
            std::vector<int> vc_nonmin {0};

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, num_vcs, vc_min, vc_nonmin);

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

    assert(next_vc >=0);
    assert(next_port>=0);
    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);
}

// Progressive Adaptive Routing
void topo_dragonfly::route_PAR(int port, int vc, internal_router_event* ev)
{   
    const int vn = ev->getVN();
    assert( vns[vn].algorithm == PAR );
    const int start_vc = vns[vn].start_vc;
    const int par_vcs = vns[vn].num_vcs;
    
    RtrEvent* rtr_ev =  ev->getEncapsulatedEvent();
    int msgid = rtr_ev->getTraceID();
    int link_hops = rtr_ev->getNumHops();
    
    assert(link_hops < 8);
    assert(par_vcs == 5);
    assert(par_vcs == num_vcs); // only using 1 vn so far

    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);

    const int current_vc = td_ev->getVC();
    assert(current_vc == vc);

    int next_port = -1;
    int next_vc = -1;

    std::vector<int> vc_min {0};
    std::vector<int> vc_nonmin {1};

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

            std::pair<int,int> adp_result = adp_select_port(td_ev, port, start_vc, num_vcs, vc_min, vc_nonmin, "1st");

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
                    td_ev, port, start_vc, num_vcs, vc_min, vc_nonmin, "2nd");

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

// q-adaptive routing
void topo_dragonfly::q_adaptive(int port, int vc, internal_router_event* ev){

    int next_vc = -1;
    int next_port = -1;
    int min_path_port = -1;

    int vn = ev->getVN();
    int start_vc = vns[vn].start_vc;
    assert( vns[vn].algorithm == Q1 || vns[vn].algorithm == Q2);
    
    topo_dragonfly_event *td_ev = static_cast<topo_dragonfly_event*>(ev);
    RtrEvent* rtr_ev = td_ev->getEncapsulatedEvent();
    int link_hops = rtr_ev->getNumHops();

    assert(link_hops <= vns[vn].num_vcs);

    int current_vc = td_ev->getVC();
    assert(current_vc == vc);
    assert(current_vc == 0 || current_vc == link_hops -1 );

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

    // Increment vc at every hop
    // when forwarding to compute node, no vc increment 
    if(min_path_port >= params.p){
        next_vc = link_hops;
    }
    else{
        next_vc = (link_hops-1 >= 0) ? link_hops-1 : 0;
    }
    assert(next_vc < vns[vn].num_vcs);

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
                int best_port = qtable_get_min_port(port, row_index, vns[vn].algorithm, params.p, params.k);
                int64_t best_est = qtable_get_est(best_port, row_index, vns[vn].algorithm); 
                td_ev->new_estimate = best_est;
            }
        }
        
        uint64_t time2save = getCurrentSimTimeNano();
        if( save_qtable && 
            (time2save >= save_qtable_time) ){
            save_qtable_toFile();
            if(router_id_global == 0 && useQrouting>1){
                dump_t2nbrTable();
            }
            save_qtable = false;
        }

        next_port = min_path_port;

        assert(next_vc >= 0);
        assert(next_port >= 0);
        assert(td_ev->new_estimate >= 0);

        td_ev->setVC(next_vc);
        td_ev->setNextPort(next_port);   
        return;
    }

    //In Src or Mid group, Not in Dest group
    int row_index = qtable_get_row_idx(td_ev);
    int best_port = qtable_get_min_port(port, row_index, vns[vn].algorithm, params.p, params.k);
    int64_t best_est = qtable_get_est(best_port, row_index, vns[vn].algorithm); 

    // qrouting: update qvalue using best Est no matter which action is actually taken
    td_ev->new_estimate = best_est;  

    // tempororay set to the best: it may be changed later  
    next_port = best_port; 

    //Set non-min routed counter
    if(group_id == td_ev->src_group){
        if(src_group_q){
            assert(!rtr_ev->getAdpRouted());
        }    
    }
    else if (group_id != td_ev->src_group && group_id != td_ev->dest.group){
        if(!rtr_ev->getAdpRouted()){
            rtr_ev->setAdpRouted();
            assert(rtr_ev->val_route_pos == 0);
            rtr_ev->val_route_pos = 1;
        }
    }
    else{
        output.fatal(CALL_INFO, -1, "ERROR DFrtr %d: Msg not in Src, intermediate, Dest group\n", router_id_global);
    }

    // limit q routing by number of hops
    int go_min = 0;
    if(src_group_q){
        if(link_hops >= max_hops || group_id != td_ev->src_group){
            go_min = 1;
        }
    }
    else if (src_mid_group_q){        
        if(group_id == td_ev->src_group && link_hops>0){
            go_min = 1;
        }
        if(group_id != td_ev->src_group && group_id != td_ev->dest.group && link_hops > 1){
            go_min = 1;
        }
    }
    else{
        if (link_hops >= max_hops ) {
            go_min = 1;
        }
    }

    if(go_min){
        next_port = min_path_port;
        assert(next_vc >= 0);
        assert(next_port >= 0);
        assert(td_ev->new_estimate >= 0);

        td_ev->setVC(next_vc);
        td_ev->setNextPort(next_port);   
        return;
    }

    bool first_mid_router = false;
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
            assert(link_hops == 1);

            first_mid_router = true;
        }
    }
    else{
        assert(link_hops < max_hops);
    }

    double q_thld = -1.0;

    if(first_mid_router){
        assert(best_port == next_port);
        assert(next_port >= params.p);

        if(min_path_port < params.p + params.a -1){

            uint32_t tmp_router = params.a+10;
            do {
               tmp_router = rng->generateNextUInt32() % params.a;
            }
            while (tmp_router == router_id);
            assert(tmp_router>=0 && tmp_router<params.a);

            best_port = port_for_router(tmp_router);
            
            next_port = best_port;
            assert(next_port>=params.p && next_port<params.a + params.p -1);

            best_est = qtable_get_est(best_port, row_index, vns[vn].algorithm); 

            q_thld = q_threshold2;
        }
        else{
            // has direct global link to dest group, min routing to dest
            next_port = min_path_port;
            assert(next_vc >= 0);
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
        int64_t min_path_est = qtable_get_est(min_path_port, row_index, vns[vn].algorithm); 
        double est_diff = double(min_path_est - best_est) / double (min_path_est);

        if(est_diff < 0){
            assert(first_mid_router);
        }

        if(est_diff < q_thld){
            next_port = min_path_port;
        }
        else{
            if(first_mid_router){
                assert(rtr_ev->val_route_pos == 1);
                rtr_ev->val_route_pos = 2;
            }
        }
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

    assert(next_vc >= 0);
    assert(next_port >= 0);
    assert(td_ev->new_estimate >= 0);
 
    td_ev->setVC(next_vc);
    td_ev->setNextPort(next_port);   
     
    if(qtable_bcast == THLD){
        do_thld_qBcast(td_ev);
    }  
    return; 
}

void topo_dragonfly::route_packet(int port, int vc, internal_router_event* ev) {

    // route_nonadaptive(port,vc,ev);
    // route_adaptive_local(port,vc,ev);

    int vn = ev->getVN();
    if( vns[vn].algorithm == UGAL_3VC ) {
        route_ugal_3vc(port,vc,ev);
    }
    else if (vns[vn].algorithm == UGAL_4VC){
        route_ugal_4vc(port,vc,ev);
    }
    else if ( vns[vn].algorithm == Q1 || vns[vn].algorithm == Q2){
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
    case Q2:
    case VALIANT:
    case ADAPTIVE_LOCAL:
    case UGAL_3VC:
    case UGAL_4VC:
    case PAR:
        //random router/group are init in routing functions
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

    if ( td_ev->getTraceType() != SST::Interfaces::SimpleNetwork::Request::NONE ) {
        output.output("TRACE(%d): process_input():"
                      " mid_group_shadow = %u\n",
                      td_ev->getTraceID(),
                      td_ev->dest.mid_group_shadow);
    }
    return td_ev;
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
            next_port = port_for_group(td_ev->dest.group, td_ev->global_slice);
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
    if ( (uint32_t)port < params.p ) return R2N;
    else return R2R;
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


uint32_t topo_dragonfly::router_to_group(uint32_t group)
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

/* returns local router port if group can't be reached from this router */
uint32_t topo_dragonfly::port_for_group(uint32_t group, uint32_t slice, int id)
{   
    assert(group >=0 && group < params.g);
    // Look up global port to use
    switch ( global_route_mode ) {
    case ABSOLUTE:
        if ( group >= group_id ) group--;
        break;
    case RELATIVE:
        if ( group > group_id ) {
            group = group - group_id - 1;
        }
        else {
            group = params.g - group_id + group - 1;
        }
        break;
    default:
        break;
    }

    const RouterPortPair& pair = group_to_global_port.getRouterPortPair(group,slice);

    if ( pair.router == router_id ) {
        return pair.port;
    } else {
        return port_for_router(pair.router);
    }

}


uint32_t topo_dragonfly::port_for_router(uint32_t router)
{   
    assert(router >= 0 && router < params.a);
    uint32_t tgt = params.p + router;
    if ( router > router_id ) tgt--;
    return tgt;
}

//table index [0 -- num_ports-host_port] (local ports ... global ports ...)
void 
topo_dragonfly::set_t2nbrTable()
{   
    int table_size = params.k - params.p;
    assert(table_size == (params.a -1 + params.h ));

    t2nbrTable = new int64_t[table_size];
    for(int p = 0; p<table_size; p++){
        if (p < params.a-1){ // local ports
            t2nbrTable[p] = link_latency_local;
        } 
        else { //global ports
            t2nbrTable[p] = link_latency_global;
        }   
    }
    if(router_id_global == 0){
        dump_t2nbrTable();
    }

}

void 
topo_dragonfly::dump_t2nbrTable(){
    printf("\nRtr %d, t2nbrTable:\n", router_id_global);
    int table_size = params.k - params.p;
    for(int p = 0; p<table_size; p++){
        printf("%ld ", t2nbrTable[p]); 
    }
    printf("\n");
}


// destgroup number is [0 -- g-1], without itself
// port index [0 -- number_local + global -1]
void
topo_dragonfly::setQtable()
{   

    assert(link_latency_global > 0);
    assert(link_latency_local > 0);

    uint32_t row = qtable_rows;
    uint32_t col = params.k - params.p;

    assert( (params.a-1 + params.h) == col );

    qtable = new int64_t[row * col];

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
        if(router_id_global == 0) printf("Rtr 0 loading qtable from %s\n", (pathToQtableFile + qtablefile).c_str() );
        //io load qtable
        std::ifstream file(pathToQtableFile + qtablefile);
        int count = 0; 
        int64_t x;
        while (file >> x) 
            qtable[count++] = x;
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
            // printf("\n");
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

// return est delivery time from curr router to destination group
int64_t topo_dragonfly::qtable_get_est(int port, int row_idx, RouteAlgo ralg){

    assert(port>=params.p && port < params.k);
    uint32_t col = params.k - params.p;
    int port_number_intable = port - params.p;

    int64_t min_est = qtable[row_idx * col + port_number_intable];

    if(ralg == Q2){
        min_est += t2nbrTable[port_number_intable];
    }
    return min_est;
}

int topo_dragonfly::qtable_get_min_port(int port, int row_idx, RouteAlgo ralg, uint32_t port_start, uint32_t port_end)
{   
    assert(port_start >= params.p);
    assert(port_end <= params.k);
    assert(port_end > port_start);

    uint32_t col = params.k - params.p;
    int64_t min_est = INT64_MAX;
    std::vector<int> index; 
    
    for(int i = port_start - params.p ; i <  port_end - params.p; i++){
        int64_t est = qtable[row_idx * col + i];
        if(ralg == Q2){
            est += t2nbrTable[i];
        }

        if(est < min_est) index.clear(); 
        if(est <= min_est) {
            index.push_back(i);
            min_est = est;
        }
    }
    assert(index.size() > 0);

    int port_min = 0;
    port_min = index[port_min] + (int)params.p;
    assert(port_min>=params.p && port_min < params.k);
    return port_min;
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
    //is not a host port, assured by portcontrol
    assert(port_number>=params.p && port_number<params.k); 
    uint32_t dest_g = idToGroup(ev->getDest());
    uint32_t dest_rtr = idToRouter(ev->getDest());

    //remember to del once receive
    qtable_event* qtable_ev = NULL;

    if(qtable_row_type == "g") {
        if(dest_g == group_id){
            if(port_number < params.p + params.a -1 && useQrouting==1) {
                return NULL; 
            } 
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
            if(port_number < params.p + params.a -1 && useQrouting==1) {
                return NULL; 
            } 
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
    uint32_t target_row = qe->target_row_in_table;
    RouteAlgo ralg = vns[qe->appvn].algorithm;

    const int target_port = port - params.p;
    uint32_t col = params.k - params.p;

    if(!qe->qBcast){
        if(qtable_row_type == "g") {
            uint32_t destg = target_row;
            if(destg == group_id){
                if(ralg == Q1) {
                    output.fatal(CALL_INFO,-1,"q1 non same group qtable update should be guaranteed by portcontrol and create_qtable_event()");
                    // return;
                } else {
                    assert(ralg == Q2);
                    t2nbrTable[target_port] = qe->queueing_time; 
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
                    assert(ralg == Q2);
                    t2nbrTable[target_port] = qe->queueing_time; 
                    deltaq = qe->new_estimate - qtable[target_index];
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
                deltaq = qe->new_estimate + qe->queueing_time - qtable[target_index];

            } else {
                assert(ralg == Q2);
                t2nbrTable[target_port] = qe->queueing_time; 
                deltaq = qe->new_estimate - qtable[target_index];
            }  

            if(deltaq < 0) {
                newq = qtable[target_index] + learning_rate * deltaq;
            } 
            else {
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
                assert(ralg == Q2);
                output.fatal(CALL_INFO,-1,"q2 does not work with node type qtable for now");
            }  

            if(deltaq < 0) {
                newq = qtable[target_index] + learning_rate * deltaq;
            } 
            else {
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
            
                if(deltaq < 0) {
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
        assert( qe->queueing_time == 0);
        uint32_t current_row = -1;

        if(qtable_row_type == "g"){
            current_row = group_id;
        }
        else{
            assert(qtable_row_type == "r");
            current_row = router_id_global;
        }
        
        if(target_row == current_row){
            printf("DFrtr %d, in group %d, new est %ld\n", router_id_global, group_id, qe->new_estimate);
            output.fatal(CALL_INFO,-1,"q2 with bcast: no same group/router qtable update should be guaranteed by portcontrol and do_thld_qBcast()");
        }
        
        if(target_row > current_row){
            target_row--;
        } 
        uint32_t target_index = target_row * col + target_port;
        int64_t newq=-10;
        int64_t deltaq=qe->new_estimate - qtable[target_index];

        if(deltaq < 0) {
            newq = qtable[target_index] + learning_rate * deltaq;
        } 
        else {
            newq = qtable[target_index] + learning_rate2 * deltaq;
        }
        assert(newq>=0); 
        qtable[target_index] = newq;
    }
}
    
//used by periodic bcasting designed to work with q1 adaptive routing
void topo_dragonfly::updateWholeQtable(qBcast_event* qe, int port){
    assert(qe);
    int64_t q_time = qe->queueing_time;
    std::vector<int64_t> new_est = qe->qBcastTable; // table of size 1*group or 1* router/system
    int target_port = port - params.p;
    uint32_t col = params.k - params.p;

    int myid = -1;
    if(qtable_row_type == "g"){
        myid = group_id;
    } 
    else{
        assert(qtable_row_type == "r");
        myid = router_id_global;
    }

    for(int r=0; r<qtable_rows+1; r++){
        if(r == myid){
            continue;
        }

        int r_table = r;
        if(r>myid){
            r_table--;
        }
        uint32_t target_index = r_table * col + target_port; 
        int64_t newq = -10;
        int64_t deltaq=0;

        if(new_est[r] == -100){
            deltaq = 0 + q_time - qtable[target_index];
        } else {
            deltaq = new_est[r] + q_time - qtable[target_index];
        }

        if(deltaq < 0) {
            newq = qtable[target_index] + learning_rate * deltaq;
        } 
        else {
            newq = qtable[target_index] + learning_rate2 * deltaq;
        }
        assert(newq>=0); 
        qtable[target_index] = newq;
    }
}

void topo_dragonfly::save_qtable_toFile(){
    uint32_t row = qtable_rows;
    uint32_t col = params.k - params.p;

    char *dir_mk = const_cast<char*>(qtableFileDir.c_str());

    mkdir(dir_mk ,0755);

    uint64_t time2save = getCurrentSimTimeNano();

    std::string table_s = "Table saved at time ";
    table_s += std::to_string(time2save/1000);
    table_s += "us \n";

    for ( uint32_t i = 0; i < row * col; i++ ) {
        table_s += std::to_string(qtable[i]);
        table_s += "\n";
    }
    char* table_sc = const_cast<char*>(table_s.c_str());

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

void topo_dragonfly::set_qBcastTable(){
    uint32_t col = qtable_rows + 1; 
    qBcastTable.resize(1 * col);

    if(qtable_row_type == "g"){
        for(int g=0; g< col; g++){
            if(g == group_id){
                qBcastTable[g] = -100;
            }else{
                qBcastTable[g] = (link_latency_global + link_latency_local);
            }
        }
    }
    else{
        assert(qtable_row_type == "r");
        for(int r=0; r< col; r++){
            if(r == router_id_global){
                qBcastTable[r] = -100;
            }else{
                int dest_g = r/params.a;
                if(dest_g == group_id){
                    qBcastTable[r] = link_latency_local;
                }
                else{
                    qBcastTable[r] = link_latency_local * 2 + link_latency_global;
                }
            }
        }
    }

    if(router_id_global == 0){
        printf("\n\nRtr 0 qBcastTable set as:\n");
        for(int g=0; g< col; g++){
            printf("%ld ",  qBcastTable[g]);
        }
        printf("\n\n");
    } 
}

void topo_dragonfly::do_thld_qBcast(topo_dragonfly_event *td_ev){
    uint32_t dest_group = td_ev->dest.group;
    uint32_t dest_router = td_ev->dest.router_global;
    int64_t new_est = td_ev->new_estimate;
    int vn = td_ev->getVN(); 

    int target_row = -1;
    if(qtable_row_type == "g"){
        assert(dest_group != group_id);
        target_row = dest_group;
    } else {
        assert(qtable_row_type == "r");
        assert(dest_router != router_id_global);
        target_row = dest_router;
    }

    int64_t old_est = qBcastTable[target_row];
    assert(old_est>0);

    float diff = (float) new_est - (float) old_est;
    diff /= old_est;

    if(abs(diff) > qbcastThsld){
        qBcastTable[target_row] = new_est;
        std::vector<int> bcast_ports;
        for(int p = params.p; p<params.k; p++){
            PortInterface* curr_port = parent->get_rtr_port(p);
            int remote_rtr_id = curr_port->get_remote_rtr_id();
            int rmt_group = remote_rtr_id/params.a;
            
            if(qtable_row_type == "g" && rmt_group != target_row) {
                bcast_ports.push_back(p);
            }

            if(qtable_row_type == "r" && remote_rtr_id != target_row) {
                bcast_ports.push_back(p);
            }

        }
        assert(!bcast_ports.empty());
        parent->bcast_qvalue_thld(target_row, new_est, bcast_ports, vn);
    }
}

//table index [0 -- num_ports-host_port] (loccal ports ... global ports ...)
void 
topo_dragonfly::set_tFromNbrTable()
{   

    int table_size = params.k - params.p;
    assert(table_size == (params.a -1 + params.h ));

    tFromNbrTable.resize(table_size);

    for(int p = 0; p<table_size; p++){
        if (p < params.a-1){ // local ports
            tFromNbrTable[p] = link_latency_local;
        } 
        else { //global ports
            tFromNbrTable[p] = link_latency_global;
        }   
    }
    if(router_id_global == 0){
        dump_tFromNbrTable();
    }
}

void 
topo_dragonfly::dump_tFromNbrTable(){
    printf("\nRtr %d, tFromNbrTable:\n", router_id_global);
    int table_size = params.k - params.p;
    for(int p = 0; p<table_size; p++){
        printf("%ld ", tFromNbrTable[p]); 
    }
    printf("\n");
}

//this function is only useful for bcasting PERID
void topo_dragonfly::update_tFromNbrTable(int port, internal_router_event* ev){

    if(qtable_bcast != PERID) return;

    assert(port>=params.p && port<params.k); // it should not be a host port
    tFromNbrTable[port - params.p] = ev->queueing_time;
}

bool topo_dragonfly::perid_qBcast_handler(Cycle_t cycle){

    assert(useQrouting == 1 && qtable_bcast == PERID);
    std::vector<int> bcast_ports;

    // construct bcast table; fromNbr table is accurate already
    int myid = -1;
    if(qtable_row_type == "g"){
        myid = group_id;
    } 
    else{
        assert(qtable_row_type == "r");
        myid = router_id_global;
    }

    for(int r=0; r<qtable_rows+1; r++){
        if(r == myid){
            assert(qBcastTable[r] == -100);
        } else {
            int r_intable = r;
            if(r > myid){
                r_intable--;
            }
            int p = qtable_get_min_port(0, r_intable, Q1, params.p, params.k); //0 -- just a placeholder
            qBcastTable[r] = qtable_get_est(p, r_intable, Q1);
        }
    }

    for(int p = params.p; p<params.k; p++){
        bcast_ports.push_back(p);
    }
    assert(!bcast_ports.empty());

    parent->bcast_qvalue_perid(bcast_ports, params.p, qBcastTable, tFromNbrTable);

    return false;
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