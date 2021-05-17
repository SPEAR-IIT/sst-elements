// -*- mode: c++ -*-

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


#ifndef COMPONENTS_MERLIN_TOPOLOGY_DRAGONFLY_H
#define COMPONENTS_MERLIN_TOPOLOGY_DRAGONFLY_H

#include <sst/core/event.h>
#include <sst/core/link.h>
#include <sst/core/params.h>
#include <sst/core/rng/sstrng.h>

#include "sst/elements/merlin/router.h"

#include <fstream>
#include <sys/stat.h>

namespace SST {
class SharedRegion;

namespace Merlin {

class topo_dragonfly_event;

struct RouterPortPair {
    uint16_t router;
    uint16_t port;

    RouterPortPair(int router, int port) :
        router(router),
        port(port)
        {}

    RouterPortPair() {}
};

class RouteToGroup {
private:
    const RouterPortPair* data;
    SharedRegion* region;
    size_t groups;
    size_t routes;


public:
    RouteToGroup() {}

    void init(SharedRegion* sr, size_t g, size_t r);

    const RouterPortPair& getRouterPortPair(int group, int route_number);

    void setRouterPortPair(int group, int route_number, const RouterPortPair& pair);
};


class topo_dragonfly: public Topology {

public:

    SST_ELI_REGISTER_SUBCOMPONENT_DERIVED(
        topo_dragonfly,
        "merlin",
        "dragonfly",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "Dragonfly topology object.  Implements a dragonfly with a single all to all pattern within the group.",
        SST::Merlin::Topology)

    SST_ELI_DOCUMENT_PARAMS(
        {"dragonfly:hosts_per_router",      "Number of hosts connected to each router."},
        {"dragonfly:routers_per_group",     "Number of links used to connect to routers in same group."},
        {"dragonfly:intergroup_per_router", "Number of links per router connected to other groups."},
        {"dragonfly:intergroup_links",      "Number of links between each pair of groups."},
        {"dragonfly:num_groups",            "Number of groups in network."},
        {"dragonfly:algorithm",             "Routing algorithm to use [minmal (default) | valiant].", "minimal"},
        {"dragonfly:adaptive_threshold",    "Threshold to use when make adaptive routing decisions.", "2.0"},
        {"dragonfly:global_link_map",       "Array specifying connectivity of global links in each dragonfly group."},
        {"dragonfly:global_route_mode",     "Mode for intepreting global link map [absolute (default) | relative].","absolute"},

        {"hosts_per_router",      "Number of hosts connected to each router."},
        {"routers_per_group",     "Number of links used to connect to routers in same group."},
        {"intergroup_per_router", "Number of links per router connected to other groups."},
        {"intergroup_links",      "Number of links between each pair of groups."},
        {"num_groups",            "Number of groups in network."},
        {"algorithm",             "Routing algorithm to use [minmal (default) | valiant].", "minimal"},
        {"adaptive_threshold",    "Threshold to use when make adaptive routing decisions.", "2.0"},
        {"global_link_map",       "Array specifying connectivity of global links in each dragonfly group."},
        {"global_route_mode",     "Mode for intepreting global link map [absolute (default) | relative].","absolute"},

        {"dragonfly:link_lat_global",    "link latency", "100ns"},
        {"link_lat_global",    "link latency", "100ns"},

        {"dragonfly:link_lat_local",    "link latency", "100ns"},
        {"link_lat_local",    "link latency", "100ns"},

        {"dragonfly:learning_rate",    "learning rate for Qrouting", "0.5"},
        {"learning_rate",     "learning rate for Qrouting", "0.5"},

        {"dragonfly:pathToQtableFile",    "path to qtable file for loading", ""},
        {"pathToQtableFile",     "path to qtable file for loading", ""},

        {"dragonfly:save_qtable",    "save ?", "false"},
        {"save_qtable",     "save?", "false"},

        {"dragonfly:max_hops",    "maximum hops allowed for q", "6"},
        {"max_hops",     "maximum hops allowed for q", "6"},

        {"dragonfly:epsilon",    "probability for qrouting exploration instead of best", "0.1"},
        {"epsilon",    "probability for qrouting exploration instead of best", "0.1"},

        {"load_qtable",     "yes?", "false"},

        // broadcasting q-table to neighbors
        {"qtable_bcast",     "nobcast, perid, thld", "nobcast"}, 
        {"qbcastThsld",    "bcast qtable by Threshold", "0.1"}, //0.0 -- 1.0
        {"qbcastPerid",    "1000ns", "1000ns"},

        // qtable num of rows: groups or all routers
        {"qtable_row_type",     "g or r or n or destG_srcN", "g"}, 
        {"learning_rate2",     "learning rate beta for hysteretic q-learning for Qrouting", "learning_rate"},
        {"src_group_q",     "only use q-routing in src group", "false"},
        {"src_mid_group_q",     "Only the first router in the source group or in the first mid group do q-routing. src_group_q and src_mid_q can not be true at the same time. max_hops became irrelevant", "false"},
        {"save_qtable_time",    "when to save qtable (us)", "1000us"},
        {"prid_func",     "yes?", "false"},
        {"q_threshold1",    "Threshold for q-routing. If (min_path_est - min_est) / min_path_est < threshold, min_path will be selected", "0.0"},
        {"q_threshold2",    "if src_mid_group_q is enabled, this threshold is used for intermediate group, in group qrouting. if src_mid_group_q is not enabled, this params has no effect", "0.0"},
    )

    /****************************************************************
    q-adaptive  == q-adaptive-routing-1
    
    q-adaptive2 == q-adaptive-routing-2 // under development

    q2 adopts the idea of near-end delivery time (from current router to neighbor router) and far-end
    delivery time (from neighbor to dest GROUP). 
    Estimation used for routing is the sum of the two.
    in this case, near-end table is updated by replacemnt with new value
    far-end table is updated with a learning rate without immediate reward
    ************************************************************************/


    /* Assumed connectivity of each router:
     * ports [0, p-1]:      Hosts
     * ports [p, p+a-2]:    Intra-group
     * ports [p+a-1, k-1]:  Inter-group
     */

    struct dgnflyParams {
        uint32_t p;  /* # of hosts / router */
        uint32_t a;  /* # of routers / group */
        uint32_t k;  /* Router Radix */
        uint32_t h;  /* # of ports / router to connect to other groups */
        uint32_t g;  /* # of Groups */
        uint32_t n;  /* # of links between groups in a pair */
    };

    enum RouteAlgo {
        MINIMAL = 0,
        VALIANT,
        ADAPTIVE_LOCAL,
        Q1,
        Q2,
        PAR,
        UGAL_3VC,
        UGAL_4VC,
    };

    enum QBcastType {
        NOBCAST = 0,   // this may trigger exploration is set
        PERID,      // periodic bcasting => only works with q1
        THLD,       // threshold bcasting => only works with q2
    };

    RouteToGroup group_to_global_port;

    struct dgnflyParams params;
    double adaptive_threshold;
    uint32_t group_id;
    uint32_t router_id; 
    
    //global rtr id
    uint32_t router_id_global; 

    RNG::SSTRandom* rng;

    int const* output_credits;       
    int const* output_queue_lengths;

    //credit of output port track remote router input port free slot number
    int const* output_2nbr_credits; 
    int const* output_used_credits;

    int num_vcs;
    int num_vns;

    enum global_route_mode_t { ABSOLUTE, RELATIVE };
    global_route_mode_t global_route_mode;

    int64_t link_latency_global;
    int64_t link_latency_local;
    float learning_rate;    //alpha 
    float learning_rate2;   //beta
    float epsilon;
    bool save_qtable;
    int save_qtable_time;

    std::string qtablefile;
    std::string pathToQtableFile;
    std::string qtableFileDir="qtables/";
    int max_hops;
    Output out2file; //for io q-table
    int64_t* qtable;
    int64_t* t2nbrTable;
    std::string qtable_row_type;
    int qtable_rows;

    int t2nbrTable_size;
    bool load_qtable;
    QBcastType qtable_bcast;
    float qbcastThsld;

    Router* parent;
    std::vector<int64_t> qBcastTable; // this table stores the min est time from current router to dest group. Used for qvalue broadcasting

    int qbcastPerid;
    std::vector<int64_t> tFromNbrTable; 

    int bcast_itr; // used by perid bcast, record current bcasat iteration

    std::string outfile;

    bool src_group_q; 
    bool src_mid_group_q; 

    std::vector<int> global_port_idx; 

    double q_threshold1;
    double q_threshold2;

    bool prid_func;

public:
    struct dgnflyAddr {
        uint32_t group;
        uint32_t mid_group;
        uint32_t mid_group_shadow;
        uint32_t router;
        uint32_t host;

        uint32_t router_global;
        uint32_t mid_router; //relative to group
    };

    topo_dragonfly(ComponentId_t cid, Params& p, int num_ports, int rtr_id, int num_vns);
    ~topo_dragonfly();

    virtual void route_packet(int port, int vc, internal_router_event* ev);

    virtual internal_router_event* process_input(RtrEvent* ev);

    virtual PortState getPortState(int port) const;
    virtual std::string getPortLogicalGroup(int port) const;

    virtual void routeInitData(int port, internal_router_event* ev, std::vector<int> &outPorts);
    virtual internal_router_event* process_InitData_input(RtrEvent* ev);

    virtual void getVCsPerVN(std::vector<int>& vcs_per_vn) {
        for ( int i = 0; i < num_vns; ++i ) {
            vcs_per_vn[i] = vns[i].num_vcs;
        }
    }

    virtual int computeNumVCs(int vns) { return vns * 3; }
    virtual int getEndpointID(int port);
    virtual void setOutputBufferCreditArray(int const* array, int vcs);
    virtual void setOutputQueueLengthsArray(int const* array, int vcs);

    virtual void setOutput2NbrCreditArray(int const* array, int vcs);
    virtual void setOutputUsedCreditArray(int const* array, int vcs);
    virtual void setQtable();
    virtual qtable_event* create_qtable_event(internal_router_event* ev, int port_number, bool bcast);
    virtual void updateQtable(qtable_event* qe, int port);
    virtual void updateWholeQtable(qBcast_event* qe, int port);
    virtual void set_t2nbrTable();
    virtual void set_tFromNbrTable();
    virtual void set_parent(Router* router);
    virtual void set_qBcastTable(); 
    virtual void update_tFromNbrTable(int port, internal_router_event* ev);
    
private:
    void idToLocation(int id, dgnflyAddr *location);
    uint32_t router_to_group(uint32_t group);
    uint32_t port_for_router(uint32_t router);
    uint32_t port_for_group(uint32_t group, uint32_t global_slice, int id = -1);

    struct vn_info {
        int start_vc;
        int num_vcs;
        RouteAlgo algorithm;
    };

    vn_info* vns;

    void route_nonadaptive(int port, int vc, internal_router_event* ev);
    void route_adaptive_local(int port, int vc, internal_router_event* ev);

    int choose_port_for_group(topo_dragonfly_event* td_ev, int dest_g);

    std::pair<int,int> adp_select_port(topo_dragonfly_event *td_ev, int port, int vc_start, int vc_end, std::vector<int> vc_min, std::vector<int> vc_nonmin, std::string debuginfo="debugMsg");

    void route_ugal_3vc(int port, int vc, internal_router_event* ev);
    void route_ugal_4vc(int port, int vc, internal_router_event* ev);
    void route_minimal(int port, int vc, internal_router_event* ev);
    void route_valiant(int port, int vc, internal_router_event* ev);
    void route_PAR(int port, int vc, internal_router_event* ev);
    void q_adaptive(int port, int vc, internal_router_event* ev);
    uint32_t idToGroup(int id); 
    uint32_t idToRouter(int id); 

    int qtable_get_min_port(int port, int row_idx, RouteAlgo ralg, uint32_t port_start, uint32_t port_end);
    int64_t qtable_get_est(int port, int row_idx, RouteAlgo ralg);
    void dumpQtable();
    void save_qtable_toFile();

    int useQrouting; // 0(not q routing), 1, 2
    void dump_t2nbrTable();
    void dump_tFromNbrTable();
    int qtable_get_row_idx(topo_dragonfly_event *td_ev);
    void do_thld_qBcast(topo_dragonfly_event *td_ev);
    void update_qBcastTable(topo_dragonfly_event *td_ev);

    Clock::Handler<topo_dragonfly>* qPeridClock_handler; //qtable periodic update clock
    TimeConverter* qPerid_tc;
    UnitAlgebra bcast_clock;

    bool perid_qBcast_handler(Cycle_t cycle);
    void route_static(int port, int vc, internal_router_event* ev);
    int get_port_remote_group_id(int outport);

    Clock::Handler<topo_dragonfly>* peridClock_handler; //qtable periodic update clock
    TimeConverter* perid_tc;

    bool perid_funct_handler(Cycle_t cycle);
};


class topo_dragonfly_event : public internal_router_event {

public:
    uint32_t src_group;
    topo_dragonfly::dgnflyAddr dest;
    uint16_t global_slice;
    uint16_t global_slice_shadow;

    topo_dragonfly_event() { }
    topo_dragonfly_event(const topo_dragonfly::dgnflyAddr &dest) :
        dest(dest), global_slice(0)
        {}
    ~topo_dragonfly_event() { }

    virtual internal_router_event *clone(void) override
    {
        return new topo_dragonfly_event(*this);
    }

    void serialize_order(SST::Core::Serialization::serializer &ser)  override {
        internal_router_event::serialize_order(ser);
        ser & src_group;
        ser & dest.group;
        ser & dest.mid_group;
        ser & dest.mid_group_shadow;
        ser & dest.router;
        ser & dest.host;
        ser & dest.router_global; 
        ser & dest.mid_router;
        ser & global_slice;
        ser & global_slice_shadow;
    }

private:
    ImplementSerializable(SST::Merlin::topo_dragonfly_event)

};

}
}

#endif // COMPONENTS_MERLIN_TOPOLOGY_DRAGONFLY_H
