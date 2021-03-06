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


#ifndef COMPONENTS_MERLIN_ROUTER_H
#define COMPONENTS_MERLIN_ROUTER_H

#include <sst/core/component.h>
#include <sst/core/event.h>
#include <sst/core/link.h>
#include <sst/core/subcomponent.h>
#include <sst/core/timeConverter.h>
#include <sst/core/unitAlgebra.h>
#include <sst/core/interfaces/simpleNetwork.h>

#include <queue>

namespace SST {
namespace Merlin {

#define VERIFY_DECLOCKING 0
    
const int INIT_BROADCAST_ADDR = -1;

class TopologyEvent;

class PortInterface;
    
class Router : public Component {
    private:
        bool requestNotifyOnEvent;

        Router() :
        	Component(),
        	requestNotifyOnEvent(false),
        	vcs_with_data(0)
        {}

    protected:
        inline void setRequestNotifyOnEvent(bool state)
        { requestNotifyOnEvent = state; }

        int vcs_with_data;
        
    public:

        Router(ComponentId_t id) :
            Component(id),
            requestNotifyOnEvent(false),
            vcs_with_data(0)
        {}

        virtual ~Router() {}
        
        inline bool getRequestNotifyOnEvent() { return requestNotifyOnEvent; }
       
        virtual void notifyEvent() {}

        inline void inc_vcs_with_data() { vcs_with_data++; }
        inline void dec_vcs_with_data() { vcs_with_data--; }
        inline int get_vcs_with_data() { return vcs_with_data; }

        virtual void bcast_qvalue_thld(uint32_t dest_group, int64_t new_est, std::vector<int> bcast_ports, int vn) = 0;
        virtual PortInterface* get_rtr_port(int id) = 0;
        virtual void bcast_qvalue_perid(std::vector<int> bcast_ports, uint32_t num_host, std::vector<int64_t> qBcastTable, std::vector<int64_t> tFromNbrTable) = 0;

    virtual int const* getOutputBufferCredits() = 0;
    virtual void sendTopologyEvent(int port, TopologyEvent* ev) = 0;
    virtual void recvTopologyEvent(int port, TopologyEvent* ev) = 0;
};

#define MERLIN_ENABLE_TRACE


class BaseRtrEvent : public Event {

    public:
        enum RtrEventType {CREDIT, PACKET, INTERNAL, TOPOLOGY, INITIALIZATION, QTABLE, QBCAST};

        inline RtrEventType getType() const { return type; }

        void serialize_order(SST::Core::Serialization::serializer &ser)  override {
            Event::serialize_order(ser);
            ser & type;
        }    
        
    protected:
        BaseRtrEvent(RtrEventType type) :
            Event(),
            type(type)
        {}

    private:
        BaseRtrEvent()  {} // For Serialization only
        RtrEventType type;

        ImplementSerializable(SST::Merlin::BaseRtrEvent);
};
    

class RtrEvent : public BaseRtrEvent {

    friend class internal_router_event;

public:
    
    RtrEvent() :
        BaseRtrEvent(BaseRtrEvent::PACKET),
        injectionTime(0)
    {}

    RtrEvent(SST::Interfaces::SimpleNetwork::Request* req, SST::Interfaces::SimpleNetwork::nid_t trusted_src, int route_vn) :
        BaseRtrEvent(BaseRtrEvent::PACKET),
        request(req),
        trusted_src(trusted_src),
        route_vn(route_vn),
        injectionTime(0)
    {   
        job_id = -1;
        val_route_pos = 0;
        midgroup = -1;
    }

    RtrEvent(SST::Interfaces::SimpleNetwork::Request* req, SST::Interfaces::SimpleNetwork::nid_t trusted_src, int route_vn, int job_id) :
        BaseRtrEvent(BaseRtrEvent::PACKET),
        request(req),
        trusted_src(trusted_src),
        route_vn(route_vn),
        injectionTime(0),
        job_id(job_id)
    {
        val_route_pos = 0;
        midgroup = -1;
    }

    
    ~RtrEvent()
    {
        if (request) delete request;
    }
    
    inline void setInjectionTime(SimTime_t time) {injectionTime = time;}
    // inline void setTraceID(int id) {traceID = id;}
    // inline void setTraceType(TraceType type) {trace = type;}
    virtual RtrEvent* clone(void)  override {
        RtrEvent *ret = new RtrEvent(*this);
        ret->request = this->request->clone();
        return ret;
    }

    inline SimTime_t getInjectionTime(void) const { return injectionTime; }
    inline SST::Interfaces::SimpleNetwork::Request::TraceType getTraceType() const {return request->getTraceType();}
    inline int getTraceID() const {return request->getTraceID();}
    
    inline void computeSizeInFlits(int flit_size ) {size_in_flits = (request->size_in_bits + flit_size - 1) / flit_size; }
    inline int getSizeInFlits() { return size_in_flits; }
    inline int getSizeInBits() { return request->size_in_bits; }

    inline SST::Interfaces::SimpleNetwork::nid_t getDest() const {return request->dest;}
    
    inline SST::Interfaces::SimpleNetwork::nid_t getTrustedSrc() { return trusted_src; }
    inline int getRouteVN() { return route_vn; }
    inline int getLogicalVN() { return request->vn; }
    SST::Interfaces::SimpleNetwork::Request* takeRequest() {
        auto ret = request;
        request = nullptr;
        return ret;
    }
    
    virtual void print(const std::string& header, Output &out) const  override {
        out.output("%s RtrEvent to be delivered at %" PRIu64 " with priority %d. src = %ld (logical: %ld), dest = %ld\n",
                   header.c_str(), getDeliveryTime(), getPriority(), trusted_src, request->src, request->dest);
        if ( request->inspectPayload() != NULL) request->inspectPayload()->print("  -> ", out);
    }

    inline int getNumHops() {return request->num_hops; }
    inline void incNumHops() {request->num_hops++; }
    inline int getSpecialIndex() {return request->special_index; }
    inline void setAdpRouted() { assert(!request->adp_routed); request->adp_routed = true; }
    inline bool getAdpRouted() {return request->adp_routed; }
    inline SST::Interfaces::SimpleNetwork::nid_t getLogicalSrc() { return request->src; }
    inline int getJobId() {return job_id; }

    void serialize_order(SST::Core::Serialization::serializer &ser)  override {
        BaseRtrEvent::serialize_order(ser);
        ser & request;
        ser & trusted_src;
        ser & route_vn;
        ser & size_in_flits;
        ser & injectionTime;
        ser & job_id;
        ser & route_path;
        //where non-min route decision is made:
        // 0 - min route; 1 src router, 2 src group 2nd router
        // 3 - mid group reroute from 1
        // 4 - mid group reroute from 2
        ser & val_route_pos; 
        ser & midgroup; 
    }
    
private:
    SST::Interfaces::SimpleNetwork::Request* request;

    SST::Interfaces::SimpleNetwork::nid_t trusted_src;
    int route_vn;
    SimTime_t injectionTime;
    int size_in_flits;

    int job_id;   
    ImplementSerializable(SST::Merlin::RtrEvent)

public:
    std::vector<int> route_path;
    //0 - min route, 1-src router non min route, 2- src_group 2nd router nonmin route, 3- mid group reroute
    int val_route_pos; 
    int midgroup;
};

class TopologyEvent : public BaseRtrEvent {
        // Allows Topology events to consume bandwidth.  If this is set to
        // zero, then no bandwidth is consumed.
        int size_in_flits;
        
    public:
        TopologyEvent(int size_in_flits) :
    	BaseRtrEvent(BaseRtrEvent::TOPOLOGY),
    	size_in_flits(size_in_flits)
        {}

        TopologyEvent() :
    	BaseRtrEvent(BaseRtrEvent::TOPOLOGY),
    	size_in_flits(0)
        {}

        inline void setSizeInFlits(int size) { size_in_flits = size; }
        inline int getSizeInFlits() { return size_in_flits; }

        virtual void print(const std::string& header, Output &out) const  override {
            out.output("%s TopologyEvent to be delivered at %" PRIu64 " with priority %d\n",
                    header.c_str(), getDeliveryTime(), getPriority());
        }

        void serialize_order(SST::Core::Serialization::serializer &ser)  override {
            BaseRtrEvent::serialize_order(ser);
            ser & size_in_flits;
        }
        
        ImplementSerializable(SST::Merlin::TopologyEvent);
};

class credit_event : public BaseRtrEvent {
    public:
        int vc;
        int credits;

        credit_event() :
    	BaseRtrEvent(BaseRtrEvent::CREDIT)
        {}

        credit_event(int vc, int credits) :
    	BaseRtrEvent(BaseRtrEvent::CREDIT),
    	vc(vc),
    	credits(credits)
        {}

        virtual void print(const std::string& header, Output &out) const  override {
            out.output("%s credit_event to be delivered at %" PRIu64 " with priority %d\n",
                    header.c_str(), getDeliveryTime(), getPriority());
        }

        void serialize_order(SST::Core::Serialization::serializer &ser)  override {
            BaseRtrEvent::serialize_order(ser);
            ser & vc;
            ser & credits;
        }
        
    private:

        ImplementSerializable(SST::Merlin::credit_event)
    
};

class RtrInitEvent : public BaseRtrEvent {
public:

    enum Commands { REQUEST_VNS, SET_VNS, REPORT_ID, REPORT_BW, REPORT_FLIT_SIZE, REPORT_PORT };

    // int num_vns;
    // int id;

    Commands command;
    int int_value;
    UnitAlgebra ua_value;

    RtrInitEvent() :
        BaseRtrEvent(BaseRtrEvent::INITIALIZATION)
    {}

    virtual void print(const std::string& header, Output &out) const  override {
        out.output("%s RtrInitEvent to be delivered at %" PRIu64 " with priority %d\n",
                header.c_str(), getDeliveryTime(), getPriority());
        out.output("%s     command: %d, int_value = %d, ua_value = %s\n",
                   header.c_str(), command, int_value, ua_value.toStringBestSI().c_str());
    }

    void serialize_order(SST::Core::Serialization::serializer &ser)  override {
        BaseRtrEvent::serialize_order(ser);
        ser & command;
        ser & int_value;
        ser & ua_value;
    }
    
    
private:
    ImplementSerializable(SST::Merlin::RtrInitEvent)
};

class internal_router_event : public BaseRtrEvent {
    int next_port;
    int next_vc;
    int vc;
    int credit_return_vc;
    RtrEvent* encap_ev;

public:
    uint64_t previous_router_arrive_time;
    int64_t queueing_time;
    int64_t new_estimate;

    internal_router_event() :
        BaseRtrEvent(BaseRtrEvent::INTERNAL)
    {
        encap_ev = NULL;
        previous_router_arrive_time = 0;
        queueing_time = -1;
        new_estimate = -1;
    }
    internal_router_event(RtrEvent* ev) :
        BaseRtrEvent(BaseRtrEvent::INTERNAL)
    {
        encap_ev = ev;
        previous_router_arrive_time = 0;
        queueing_time = -1;
        new_estimate = -1;
    }

    virtual ~internal_router_event() {
        if ( encap_ev != NULL ) delete encap_ev;
    }

    virtual internal_router_event* clone(void) override
    {
        return new internal_router_event(*this);
    };

    inline void setCreditReturnVC(int vc) {credit_return_vc = vc; return;}
    inline int getCreditReturnVC() {return credit_return_vc;}

    inline void setNextPort(int np) {next_port = np; return;}
    inline int getNextPort() {return next_port;}

    // inline void setNextVC(int vc) {next_vc = vc; return;}
    // inline int getNextVC() {return next_vc;}

    inline void setVC(int vc_in) {vc = vc_in; return;}
    inline int getVC() {return vc;}

    // inline void setVN(int vn) {encap_ev->setVN(vn); return;}
    inline int getVN() {return encap_ev->route_vn;}

    inline int getFlitCount() {return encap_ev->getSizeInFlits();}

    inline void setEncapsulatedEvent(RtrEvent* ev) {encap_ev = ev;}
    inline RtrEvent* getEncapsulatedEvent() {return encap_ev;}

    inline SST::Interfaces::SimpleNetwork::Request* inspectRequest() { return encap_ev->request; }

    inline int getDest() const {return encap_ev->request->dest;}
    inline int getSrc() const {return encap_ev->getTrustedSrc();}

    // request is priviate -> not working anymore
    inline int64_t getReqSrc() const {return encap_ev->request->src;}

    inline SST::Interfaces::SimpleNetwork::Request::TraceType getTraceType() {return encap_ev->getTraceType();}
    inline int getTraceID() {return encap_ev->getTraceID();}

    virtual void print(const std::string& header, Output &out) const  override {
        out.output("%s internal_router_event to be delivered at %" PRIu64 " with priority %d.  src = %d, dest = %d\n",
                   header.c_str(), getDeliveryTime(), getPriority(), getSrc(), getDest());
        if ( encap_ev != NULL ) encap_ev->print(header + std::string("-> "),out);
    }

    void serialize_order(SST::Core::Serialization::serializer &ser)  override {
        BaseRtrEvent::serialize_order(ser);
        ser & next_port;
        ser & next_vc;
        ser & vc;
        ser & credit_return_vc;
        ser & encap_ev;

        ser & previous_router_arrive_time;
        ser & queueing_time;
        ser & new_estimate;
    }
        
    private:
        ImplementSerializable(SST::Merlin::internal_router_event)
};


class qtable_event : public BaseRtrEvent {
    public:
        uint32_t target_row_in_table;
        int64_t queueing_time; // currrent router arrrive time - previous router arrive time
        int64_t new_estimate;
        int appvn;
        bool qBcast;
        int src_nid;

        qtable_event():
            BaseRtrEvent(BaseRtrEvent::QTABLE)
        {
            target_row_in_table = 0;
            queueing_time = -1;
            new_estimate = -1;
            appvn = 0;
            qBcast = false;
            src_nid = -1;
        }

        qtable_event(uint32_t row, int64_t q_time, int64_t new_est, int vn):
            BaseRtrEvent(BaseRtrEvent::QTABLE)
        {
            target_row_in_table = row;
            queueing_time = q_time;
            new_estimate = new_est;
            appvn = vn;
            qBcast = false;
            src_nid = -1;
        }

        qtable_event(uint32_t row, int64_t q_time, int64_t new_est, int vn, bool bcast):
            BaseRtrEvent(BaseRtrEvent::QTABLE)
        {
            target_row_in_table = row;
            queueing_time = q_time;
            new_estimate = new_est;
            appvn = vn;
            qBcast = bcast;
            src_nid = -1;
        }

        qtable_event(uint32_t row, int64_t q_time, int64_t new_est, int vn, int srcnode):
            BaseRtrEvent(BaseRtrEvent::QTABLE)
        {
            target_row_in_table = row;
            queueing_time = q_time;
            new_estimate = new_est;
            appvn = vn;
            qBcast = false;
            src_nid = srcnode;
        }

        void serialize_order(SST::Core::Serialization::serializer &ser)  override {
            BaseRtrEvent::serialize_order(ser);
            ser & target_row_in_table;
            ser & queueing_time;
            ser & new_estimate;
            ser & appvn;
            ser & qBcast;
            ser & src_nid;

        }
        
    private:
        ImplementSerializable(SST::Merlin::qtable_event)
};


class qBcast_event : public BaseRtrEvent {
    public:
        int64_t queueing_time; // currrent router arrrive time - previous router arrive time
        std::vector<int64_t> qBcastTable;

        qBcast_event():
        BaseRtrEvent(BaseRtrEvent::QBCAST)
        {
            queueing_time = -1;
        }

        qBcast_event(int64_t q_time, std::vector<int64_t> new_est):
            BaseRtrEvent(BaseRtrEvent::QBCAST)
        {
            queueing_time = q_time;
            qBcastTable = new_est;
        }

        void serialize_order(SST::Core::Serialization::serializer &ser)  override {
            BaseRtrEvent::serialize_order(ser);
            ser & queueing_time;
            ser & qBcastTable;
        }
        
    private:
        ImplementSerializable(SST::Merlin::qBcast_event)
};


class Topology : public SubComponent {
public:

    // Parameters are:  num_ports, id, num_vns
    SST_ELI_REGISTER_SUBCOMPONENT_API(SST::Merlin::Topology, int, int, int)
    
    enum PortState {R2R, R2N, UNCONNECTED};
    Topology(ComponentId_t cid) : SubComponent(cid), output(Simulation::getSimulation()->getSimulationOutput()) {}
    virtual ~Topology() {}

    virtual void route(int port, int vc, internal_router_event* ev) __attribute__ ((deprecated("route() is deprecated and will be removed in SST 11. Please use route_packet(), which is now called when a packet reaches the head of the input queue."))) { }

    virtual void reroute(int port, int vc, internal_router_event* ev) __attribute__ ((deprecated("reroute() is deprecated and will be removed in SST 11. Please use route_packet(), which is now called when a packet reaches the head of the input queue."))) {
        DISABLE_WARN_DEPRECATED_DECLARATION
        route(port,vc,ev);
        REENABLE_WARNING
    }

    virtual void route_packet(int port, int vc, internal_router_event* ev)  {
        DISABLE_WARN_DEPRECATED_DECLARATION 
        route(port,vc,ev);
        reroute(port,vc,ev);
        REENABLE_WARNING
    }
    virtual internal_router_event* process_input(RtrEvent* ev) = 0;
    
    // Returns whether the port is a router to router, router to nic, or unconnected
    virtual PortState getPortState(int port) const = 0;
    inline bool isHostPort(int port) const { return getPortState(port) == R2N; }
    virtual std::string getPortLogicalGroup(int port) const {return "";}
    
    // Methods used during init phase to route init messages
    virtual void routeInitData(int port, internal_router_event* ev, std::vector<int> &outPorts) = 0;
    virtual internal_router_event* process_InitData_input(RtrEvent* ev) = 0;

    // Method used for autodiscovery of VC/VN
    virtual int computeNumVCs(int vns) __attribute__ ((deprecated("computeNumVCs() is deprecated and will be removed in SST 11. Please use getVCsPerVN() instead."))) {return vns;}

    // Gets the number of VCs per VN for each VN.  Vector that is
    // passed in must have it's size set to num_vns before making this
    // call.
    virtual void getVCsPerVN(std::vector<int>& vns_per_vn) {
        DISABLE_WARN_DEPRECATED_DECLARATION
        int vcs = computeNumVCs(1);
        REENABLE_WARNING
        for ( int& val : vns_per_vn ) val = vcs;
    }
    // Method used to set endpoint ID
    virtual int getEndpointID(int port) {return -1;}
    
    // Sets the array that holds the credit values for all the output
    // buffers.  Format is:
    // For port=n, VC=x, location in array is n*num_vcs + x.

    // If topology does not need this information, then default
    // version will ignore it.  If topology needs the information, it
    // will need to overload function to store it.
    virtual void setOutputBufferCreditArray(int const* array, int vcs) {};
    virtual void setOutputQueueLengthsArray(int const* array, int vcs) {};

    // For adaptive routing
    virtual void setOutput2NbrCreditArray(int const* array, int vcs) {};
    virtual void setOutputUsedCreditArray(int const* array, int vcs) {};
    

    // When TopologyEvents arrive, they are sent directly to the
    // topology object for the router
    virtual void recvTopologyEvent(int port, TopologyEvent* ev) {};

    virtual void setQtable() {};
    virtual void set_t2nbrTable() {};
    // virtual std::string getTopoName() {return "";};
    // virtual uint32_t idToGroup(int id) {return 0;};
    // virtual uint64_t qtable_get_min(int dest_group, bool get_port=false) {return 0;};
    virtual qtable_event* create_qtable_event(internal_router_event* ev, int port_number, bool bcast) {return NULL;};
    virtual void updateQtable(qtable_event* qe, int port) {};

    virtual void updateWholeQtable(qBcast_event* qe, int port) {};


    // virtual void set_ports(PortInterface const** array) {};
    virtual void set_parent(Router* router) {};
    
    // virtual void check_perid_qBcast() {};

    virtual void update_tFromNbrTable(int port, internal_router_event* ev) {};

protected:
    Output &output;
};



// Class to manage link between NIC and router.  A single NIC can have
// more than one link_control (and thus link to router).
class PortInterface : public SubComponent{

public:

    // params are: parent router, router id, port number, topology object
    SST_ELI_REGISTER_SUBCOMPONENT_API(SST::Merlin::PortInterface, Router*, int, int, Topology*)

    typedef std::queue<internal_router_event*> port_queue_t;
    typedef std::queue<TopologyEvent*> topo_queue_t;

    virtual void sendTopologyEvent(TopologyEvent* ev) = 0;
    // Returns true if there is space in the output buffer and false
    // otherwise.
    virtual void send(internal_router_event* ev, int vc) = 0;
    // Returns true if there is space in the output buffer and false
    // otherwise.
    virtual bool spaceToSend(int vc, int flits) = 0;
    // Returns NULL if no event in input_buf[vc]. Otherwise, returns
    // the next event.
    virtual internal_router_event* recv(int vc) = 0;
    virtual internal_router_event** getVCHeads() = 0;
    
    // time_base is a frequency which represents the bandwidth of the link in flits/second.
    PortInterface(ComponentId_t cid) :
        SubComponent(cid)
        {}

    // virtual void initVCs(int vns, int* vcs_per_vn, internal_router_event** vc_heads, int* xbar_in_credits, int* output_queue_lengths) = 0;
    virtual void initVCs(int vns, int* vcs_per_vn, internal_router_event** vc_heads_in, int* xbar_in_credits_in, int* output_queue_lengths_in, int* output_credits_in, int* output_used_credits_in) = 0;

    virtual ~PortInterface() {}
    // void setup();
    // void finish();
    // void init(unsigned int phase);
    // void complete(unsigned int phase);
    

    virtual void sendInitData(Event *ev) = 0;
    virtual Event* recvInitData() = 0;
    virtual void sendUntimedData(Event *ev) = 0;
    virtual Event* recvUntimedData() = 0;
    
    virtual void dumpState(std::ostream& stream) {}
    virtual void printStatus(Output& out, int out_port_busy, int in_port_busy) {}
    
    // void setupVCs(int vcs, internal_router_event** vc_heads
	virtual bool decreaseLinkWidth() = 0;
	virtual bool increaseLinkWidth() = 0; 

    virtual void link_send_qevent(uint32_t dest_group, int64_t new_est, int vn) = 0 ; //used to let hr_router be able to let port send event
    virtual int get_remote_rtr_id() = 0; 
    virtual int get_port_number() = 0;

    virtual void link_send_bcastEvent(int64_t queueing_t, std::vector<int64_t> new_est) = 0; // used for bcast2


    class OutputArbitration : public SubComponent {
    public:

        SST_ELI_REGISTER_SUBCOMPONENT_API(SST::Merlin::PortInterface::OutputArbitration)
    
        OutputArbitration(ComponentId_t cid) :
            SubComponent(cid)
        {}
        virtual ~OutputArbitration() {}

        virtual void setVCs(int num_vns, int* vcs_per_vn) = 0;
        virtual int arbitrate(Cycle_t cycle, PortInterface::port_queue_t* out_q, int* port_out_credits, bool isHostPort, bool& have_packets) = 0;
        virtual void dumpState(std::ostream& stream) {};
    };

};


class XbarArbitration : public SubComponent {
public:

    SST_ELI_REGISTER_SUBCOMPONENT_API(SST::Merlin::XbarArbitration)
    
    XbarArbitration(ComponentId_t cid) :
        SubComponent(cid)
    {}
    virtual ~XbarArbitration() {}

#if VERIFY_DECLOCKING
    virtual void arbitrate(PortInterface** ports, int* port_busy, int* out_port_busy, int* progress_vc, bool clocking) = 0;
#else
    virtual void arbitrate(PortInterface** ports, int* port_busy, int* out_port_busy, int* progress_vc) = 0;
#endif
    virtual void setPorts(int num_ports, int num_vcs) = 0;
    virtual bool isOkayToPauseClock() { return true; }
    virtual void reportSkippedCycles(Cycle_t cycles) {};
    virtual void dumpState(std::ostream& stream) {};
	
};

}
}

#endif // COMPONENTS_MERLIN_ROUTER_H
