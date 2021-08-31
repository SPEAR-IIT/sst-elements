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
#include "emberallreduce.h"

using namespace SST::Ember;

static void test(void* a, void* b, int* len, PayloadDataType* ) {

	printf("%s() len=%d\n",__func__,*len);
}

EmberAllreduceGenerator::EmberAllreduceGenerator(SST::ComponentId_t id,
                                            Params& params) :
	EmberMessagePassingGenerator(id, params, "Allreduce"),
    m_loopIndex(0),
    m_startTime(0),
    m_stopTime(0),
    comm_start_time(0)
{
	m_iterations = (uint32_t) params.find("arg.iterations", 1);
	m_compute    = (uint32_t) params.find("arg.compute", 0);
	m_count      = (uint32_t) params.find("arg.count", 1);

    m_wait2start = (uint64_t) params.find("arg.wait2start", 1);

	if ( params.find<bool>("arg.doUserFunc", false )  )   {
		m_op = op_create( test, 0 );
	} else {
		m_op = Hermes::MP::SUM;
	}

    dobarrier=params.find<bool>("arg.dobarrier", false );

    mkdir("./ember_stats/" ,0755);
    outfile = "./ember_stats/emberAllreduce_" + std::to_string(m_iterations) + "_" + std::to_string(m_count) + "_" + std::to_string(m_compute) + "_rank" + std::to_string(rank()) + "_" + std::to_string(size()) + ".stats";
    FILE *file = fopen(outfile.c_str(), "w");
    fprintf(file, "rank,ite,start(ns),stop(ns),comm(ns),comm_start\n");
    fclose(file);

    if(0 == rank()){
        output( "EmberAllreduce %d ranks\n", size());
        output( "EmberAllreduce %u its\n", m_iterations);
        output( "EmberAllreduce compute time: %lu \n", m_compute);
        output( "EmberAllreduce data count %u \n", m_count);
        output( "EmberAllreduce wait2start: %lu ns \n", m_wait2start);
        output( "EmberAllreduce dobarrier?: %d\n", (bool)dobarrier);
    }
}

bool EmberAllreduceGenerator::generate( std::queue<EmberEvent*>& evQ) {

    //time unit is nanosec
    if(rank()==0){
        // progress msg
        if( (m_iterations/5 !=0) && ( m_loopIndex %  (m_iterations/5) == 0)) {
            printf("\t\tEmberAllreduce rank%d, ite %d, start %lu, end %lu @ %lu\n", rank(), m_loopIndex, m_startTime, m_stopTime, getCurrentSimTimeNano() );
        }
    }

    if(m_loopIndex>0){
        uint64_t comm_time = m_stopTime - comm_start_time;

        std::ofstream iofile;
        iofile.open(outfile.c_str(), std::ios::app);
        if (iofile.is_open())
        {
            iofile <<rank()<<","<<m_loopIndex<<","<<m_startTime<<","<<m_stopTime<<","<<comm_time<<","<<comm_start_time<<"\n";

            iofile.close();
        }
        else assert(0);
    }

    if ( m_loopIndex == m_iterations ) {
        if ( 0 == rank() ) {
            double latency = (double)(m_stopTime-m_startTime)/(double)m_iterations;
            latency /= 1000000000.0;
            output( "%s: ranks %d, loop %d, %d double(s), latency %.3f us\n",
                    getMotifName().c_str(), size(), m_iterations, m_count, latency * 1000000.0  );
        }
        return true;
    }
    if ( 0 == m_loopIndex ) {

        enQ_compute( evQ, m_wait2start );

		memSetBacked();
		m_sendBuf = memAlloc(sizeofDataType(DOUBLE)*m_count);
		m_recvBuf = memAlloc(sizeofDataType(DOUBLE)*m_count);
    }

    enQ_getTime( evQ, &m_startTime );
    enQ_compute( evQ, m_compute );

    enQ_getTime( evQ, &comm_start_time );

    // printf("Emberallreduce, rand %d, ite %d\n", rank(), m_loopIndex);
    enQ_allreduce( evQ, m_sendBuf, m_recvBuf, m_count, DOUBLE, m_op, GroupWorld );
    
    if(dobarrier){
        // printf("rank %d ember allreduce called barrier\n", rank());
        enQ_barrier( evQ, GroupWorld );
    }

    enQ_getTime( evQ, &m_stopTime );

    ++m_loopIndex;

    return false;
}
