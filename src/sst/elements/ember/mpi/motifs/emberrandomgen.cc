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

#include <sst/core/rng/xorshift.h>
#include "emberrandomgen.h"

using namespace SST::Ember;
using namespace SST::RNG;

EmberRandomTrafficGenerator::EmberRandomTrafficGenerator(SST::ComponentId_t id, Params& params) :
	EmberMessagePassingGenerator(id, params),
	m_startTime(0),
	m_stopTime(0),
	m_commStartTime(0),
	m_isend_time(0)
{

	msgSize = (uint32_t) params.find("arg.messagesize", 1);
	maxIterations = (uint32_t) params.find("arg.iterations", 1);
	m_wait2start = (uint64_t) params.find("arg.wait2start", 1);
	iteration = 0;
	
	//
	m_compute = (uint32_t) params.find("arg.compute", 0);
	waitunblock = (int) params.find("arg.waitunblock", 1);
	report_step = (int) params.find("arg.report_step", 1);
	subiterations = (int) params.find("arg.subiterations", 1);

	mkdir("./ember_stats/" ,0755);
    outfile = "./ember_stats/emberRandom_rank" + std::to_string(rank()) + "_" + std::to_string(size()) + ".stats";
    FILE *file = fopen(outfile.c_str(), "w");
    fprintf(file, "rank,ite,start(ns),stop(ns),comm(ns),comm.start\n");
    fclose(file);

	// init ostring
    ostring.str("");

	if(rank()==0){
		printf("=============================\n");
		printf("EmberRandom: iterations %u \n", maxIterations);
		printf("EmberRandom: subiterations %d \n", subiterations);
		printf("EmberRandom: msgsize (double) %u\n", msgSize);
		printf("EmberRandom: m_wait2start? %lu \n", m_wait2start);
		printf("EmberRandom: computetime %lu \n", m_compute);
		printf("EmberRandom: waitunblock? %d \n", waitunblock);
		printf("EmberRandom: progress speed %d \n", report_step);
	}
}

bool EmberRandomTrafficGenerator::generate( std::queue<EmberEvent*>& evQ ) {
	//progresss print
	uint64_t comm_time = m_stopTime - m_commStartTime;
	if(rank()==0 && iteration%report_step==0){
		printf("\tEmberRandom rank 0, ite %u, start %lu, comm.start %lu, isend time %lu, stop %lu, comm time %lu\n", iteration, m_startTime, m_commStartTime, m_isend_time, m_stopTime, comm_time);
	}
	//IO
	if(iteration>0){
		ostring <<rank()<<","<<iteration<<","<<m_startTime<<","<<m_stopTime<<","<<comm_time<<","<<m_commStartTime<<"\n";
		// IO every 10MB or at end of sim
		if(ostring.str().size() > 10000000 || iteration == maxIterations){
			std::ofstream iofile;
			iofile.open(outfile.c_str(), std::ios::app);
			if (iofile.is_open())
			{
				iofile << ostring.rdbuf();
				iofile.close();
				ostring.str("");
			}
			else assert(0);
		}
	}

	if(iteration == maxIterations) {
		//
		if(rank()==0) printf("\t\tEmberRandom rank 0 finished %u ites\n", maxIterations);
		return true;
	} 
	else {
		if(iteration==0){
			enQ_compute( evQ, m_wait2start );
		}

		// iteration starts
		enQ_getTime( evQ, &m_startTime );
		enQ_compute( evQ, m_compute );

		const uint32_t worldSize = size();
		// loop over empty vector is legit
		std::vector<uint32_t> send_to_list;
		std::vector<uint32_t> recv_from_list;

		for(int subite=0; subite<subiterations; subite++){
			XORShiftRNG* rng = new XORShiftRNG(size() + iteration);	
			uint32_t* rankArray = (uint32_t*) malloc(sizeof(uint32_t) * worldSize);

			// Set up initial passing (everyone passes to themselves)
			for(uint32_t i = 0; i < worldSize; i++) {
				rankArray[i] = i;
			}

			for(uint32_t i = 0; i < worldSize; i++) {
				uint32_t swapPartner = rng->generateNextUInt32() % worldSize;
				uint32_t oldValue = rankArray[i];
				rankArray[i] = rankArray[swapPartner];
				rankArray[swapPartner] = oldValue;
			}

			// void* allocBuffer = memAlloc( msgSize * 8 );
			if(rankArray[rank()] != rank()) {
				send_to_list.push_back(rankArray[rank()]);
				// MessageRequest*  req  = new MessageRequest();
                // msgRequests.push_back(req);
				// enQ_isend( evQ, allocBuffer, msgSize, DOUBLE, rankArray[rank()],
				// 	0, GroupWorld, req);
			}
			// enQ_getTime( evQ, &m_isend_time );

			// void* recvBuffer = memAlloc( msgSize * 8 );
			for(uint32_t i = 0; i < worldSize; i++) {
				if(rank() == rankArray[i] && (i != rank())) {
					recv_from_list.push_back(i);
					// MessageRequest*  req  = new MessageRequest();
                	// msgRequests.push_back(req);
					// enQ_irecv( evQ, recvBuffer, msgSize, DOUBLE, i, 0,
						// GroupWorld, req);
					break;
				}
			}
			delete rng;
		}
		assert(send_to_list.size() == recv_from_list.size());
		int neighbor_size = send_to_list.size();
		msgRequests.resize( neighbor_size * 2 );

		enQ_getTime( evQ, &m_commStartTime );

		for(int sendid=0; sendid<neighbor_size; sendid++){
			void* allocBuffer = memAlloc( msgSize * 8 );
			enQ_isend( evQ, allocBuffer, msgSize, DOUBLE, send_to_list[sendid],
					0, GroupWorld, &msgRequests[sendid]);
		}

		enQ_getTime( evQ, &m_isend_time );

		for(int recvid=0; recvid<neighbor_size; recvid++){
			void* recvBuffer = memAlloc( msgSize * 8 );
			enQ_irecv( evQ, recvBuffer, msgSize, DOUBLE, recv_from_list[recvid], 0,
				GroupWorld, &msgRequests[neighbor_size+recvid]);
		}

		if(waitunblock && neighbor_size>0 ){
			enQ_waitall( evQ, neighbor_size*2, &msgRequests[0], NULL );
			// for(uint32_t i = 0; i < msgRequests.size(); ++i) {
            // 	enQ_wait( evQ, msgRequests[i]);
        	// }
		}
		enQ_getTime( evQ, &m_stopTime );
		// msgRequests.clear();
		iteration++;
		return false;
	}
}

