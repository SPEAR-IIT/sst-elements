// Copyright 2009-2014 Sandia Corporation. Under the terms
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
// Government retains certain rights in this software.
//
// Copyright (c) 2009-2014, Sandia Corporation
// All rights reserved.
//
// This file is part of the SST software package. For license
// information, see the LICENSE file in the top level directory of the
// distribution.

#ifndef _H_SST_MEM_H_REQ_GEN_CPU
#define _H_SST_MEM_H_REQ_GEN_CPU

#include <sst/core/component.h>
#include <sst/core/interfaces/simpleMem.h>

#include "mirandaGenerator.h"

using namespace SST::Interfaces;

namespace SST {
namespace Miranda {

class RequestGenCPU : public SST::Component {
public:

	RequestGenCPU(SST::ComponentId_t id, SST::Params& params);
	void finish();
	void init(unsigned int phase);

private:
	RequestGenCPU();  // for serialization only
	RequestGenCPU(const RequestGenCPU&); // do not implement
	void operator=(const RequestGenCPU&); // do not implement
	~RequestGenCPU();

	void handleEvent( SimpleMem::Request* ev );
	bool clockTick( SST::Cycle_t );
	void issueRequest(RequestGeneratorRequest* req);

    	Output* out;

	RequestGenerator* reqGen;
	std::map<SimpleMem::Request::id_t, SimTime_t> requestsInFlight;
    	SimpleMem* cache_link;

	std::queue<RequestGeneratorRequest*> pendingRequests;

	uint32_t maxRequestsPending;
	uint32_t requestsPending;
	uint32_t reqMaxPerCycle;
	uint64_t cacheLine;

	uint64_t readRequestsIssued;
	uint64_t splitReadRequestsIssued;
	uint64_t writeRequestsIssued;
	uint64_t splitWriteRequestsIssued;
	uint64_t cyclesWithoutIssue;
	uint64_t cyclesWithIssue;
	uint64_t bytesRead;
	uint64_t bytesWritten;

	uint64_t reqLatency;

	bool printStats;

};

}
}
#endif /* _RequestGenCPU_H */