#include <sst_config.h>
#include "emberlulesh.h"

using namespace SST::Ember;

EmberLuleshGenerator::EmberLuleshGenerator(SST::ComponentId_t id,
                                    Params& params) :
	EmberMessagePassingGenerator(id, params, "Bcast"),
    m_loopIndex(0),
    m_startTime(0),
    m_stopTime(0),
    curr_phase(0)
{   
    // 0 means No
    dobcast = (params.find("arg.dobcast", 1) != 0);
    doreduce = (params.find("arg.doreduce", 1) != 0);
    do3dnn26 = (params.find("arg.do3dnn26", 1) != 0);
    do3dsweep = (params.find("arg.do3dsweep", 1) != 0);

    msgscale = (double) params.find<double>("arg.msgscale", 1);

	m_iterations = (uint32_t) params.find("arg.iterations", 1);
    m_compute    = (uint32_t) params.find("arg.compute", 0);
    m_wait2start = params.find<uint64_t>("arg.wait2start", 0);
	
    m_root    = 0;
    m_sendBuf = NULL;

    px = (int32_t) params.find<int32_t>("arg.pep", 0);
	py = (int32_t) params.find<int32_t>("arg.peq", 0);
	kba = (uint32_t) params.find<uint32_t>("arg.kba", 1);

    	// Check that we are using all the processors or else lock up will happen :(.
	if( (px * py) != (signed) size() ) {
		fatal(CALL_INFO, -1, "Error: Sweep 3D motif checked "
			"processor decomposition: %" PRIu32 "x%" PRIu32 " != MPI World %"
			PRIu32 "\n", px, py, size());
	}

    //3d26NN
    config3d26(params);

    // print 26 neighbors
    // int32_t nn26[] = {xface_down, xface_up, yface_down,	yface_up, zface_down,	zface_up, line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h,	line_i,	line_j,	line_k,	line_l, corner_a, corner_b,	corner_c, corner_d, corner_e, corner_f,	corner_g, corner_h,};
    // if(rank()==0){
    //     output("3dnn 26 neighbors\n");
    // }
    // output("%d,", rank());
    // for (int32_t x : nn26){
    //     output("%d,", x);
    // }
    // output("\n");

    // comm. quantity B
    int typeSize = sizeofDataType(DOUBLE);
    bcast_vol = msgscale * 51984 / (6*6*6 * typeSize);  
    reduce_vol = msgscale * 51992 / (6*6*6 * typeSize);  
    nn3d26_vol = msgscale * 290279024 / (3880);  // for a 6x6x6 cude, 3880 data exchange in total
    sweep3d_vol = msgscale * 299785872 / (402*kba); // number data exchange if mapped to a 12x18 processes array 

    motifname = "emberLulesh";
    mkdir("./ember_stats/" ,0755);
	outfile = "./ember_stats/" + motifname + "_rank" + std::to_string(rank()) + "_" + std::to_string(size()) + ".stats";
	FILE *file = fopen(outfile.c_str(), "w");
	fprintf(file, "rank,ite,start(ns),stop(ns),comm(ns),phase(bcast-reduce-nn3d-sweep3d)\n");
	fclose(file);

    if(0 == rank()) {
		output("\n");
        output(" Lulesh: process %d: %ux%ux%u\n", size(), peX, peY, peZ);
        output(" Lulesh: iterations %u\n", m_iterations);
        output(" Lulesh: msg scale %f\n", msgscale);
        output(" Lulesh wait2start: %lu ns\n", m_wait2start);	
        output(" Lulesh compute time: %lu ns\n", m_compute);	
        output(" Lulesh: dobcast %d | doreduce %d | do3dnn26 %d | do3dsweep %d \n", dobcast, doreduce, do3dnn26, do3dsweep);
        output(" Lulesh bcast %u x DOUBLE\n", bcast_vol);
        output(" Lulesh reduce %u x DOUBLE\n", reduce_vol);
        output(" Lulesh nn3d26 %u\n", nn3d26_vol);
        output(" Lulesh sweep3d %u\n", sweep3d_vol);
        output(" Lulesh sweep3d process array %ux%u\n", px, py);
        output(" Lulesh sweep3d slice count (kba) %u\n", kba);
		output("\n");
	}
}

bool EmberLuleshGenerator::generate( std::queue<EmberEvent*>& evQ) {

    uint64_t comm_time = m_stopTime - m_startTime - m_compute;

	if(rank()==0 ){
	    output(" %s rank%d, ite %d, start %lu, end %lu, comm.time %lu @ %lu, phase %d\n", motifname.c_str(), rank(), m_loopIndex, m_startTime, m_stopTime, comm_time, getCurrentSimTimeNano(), curr_phase);
	}

    if(rank()==size()-1 ){
		output(" \t%s rank%d, ite %d, start %lu, end %lu, comm.time %lu @ %lu, phase %d\n", motifname.c_str(), rank(), m_loopIndex, m_startTime, m_stopTime, comm_time, getCurrentSimTimeNano(), curr_phase);
	}

	//yao IO
	if(m_loopIndex>0 || curr_phase > BCAST){
        int the_phase = 0;
        int the_ite = 0;
        // get previous loop info
        if (curr_phase == 0){
            the_phase = DONE - 1;
            the_ite = m_loopIndex;
        }
        else{
            the_ite = m_loopIndex + 1;
            the_phase = curr_phase - 1;
        }

		std::ofstream iofile;
		iofile.open(outfile.c_str(), std::ios::app);
		if (iofile.is_open())
		{	
			iofile <<rank()<<","<<the_ite<<","<<m_startTime<<","<<m_stopTime<<","<<comm_time<<","<<the_phase<<"\n";
			iofile.close();
		}
		else assert(0);
	}

	//change place of exit assertion so that last iteration is record
	if ( m_loopIndex == m_iterations) {
		//yao
        assert(curr_phase == BCAST );
		if(rank()==0 || rank()==size()-1){
			output("%s rank %d finished %d iterations\n", motifname.c_str(), rank(), m_iterations);
		}
		return true;
	}

    //yao wait 2 start
    if ( 0 == m_loopIndex && curr_phase == BCAST ){
        enQ_compute( evQ, m_wait2start );
    }

    enQ_getTime( evQ, &m_startTime );
    
    switch(curr_phase) {
        case BCAST:
            if(dobcast){
                enQ_bcast( evQ, m_sendBuf, bcast_vol, DOUBLE, m_root, GroupWorld );
                // if(rank()==0 || rank()==size()-1){
                //     output("     LULESH rank %d, bcast done @ %lu\n", rank(), getCurrentSimTimeNano());
                // }
            }
            enQ_compute( evQ, m_compute );
            break;

        case REDUCE:
            if(doreduce){
                m_sendBuf = NULL;
                m_recvBuf = NULL;
                enQ_reduce( evQ, m_sendBuf, m_recvBuf, reduce_vol, DOUBLE,
                        MP::SUM, m_root, GroupWorld );
                // if(rank()==0 || rank()==size()-1){
                //     output("     LULESH rank %d, reduce done @ %lu\n", rank(), getCurrentSimTimeNano());
                // }
            }
            enQ_compute( evQ, m_compute );
            break;

        case NN3D26:
            if(do3dnn26){
                // 3dNN26
                issue3D26stencil(evQ);
                // if(rank()==0 || rank()==size()-1){
                //     output("     LULESH rank %d, nn3d26 done @ %lu\n", rank(), getCurrentSimTimeNano());
                // }
            }
            enQ_compute( evQ, m_compute );
            break;
        
        case SWEEP3D:
            if(do3dsweep){
                configsweep();
                issue3dsweep(evQ);
                // if(rank()==0 || rank()==size()-1){
                //     output("     LULESH rank %d, sweep3d done @ %lu\n", rank(), getCurrentSimTimeNano());
                // }
            }
            enQ_compute( evQ, m_compute );
            break;

        default:
            output(" ERROR!!! LULESH rank %d, phase %d\n", rank(), curr_phase);
            assert(0);
    }
    enQ_getTime( evQ, &m_stopTime);

    curr_phase++;
    if (curr_phase == DONE){
        curr_phase=0;
        m_loopIndex++;
    }

    assert(curr_phase<DONE);
    return false;
}


void EmberLuleshGenerator::config3d26(Params& params) {
    peX = (uint32_t) params.find("arg.pex", 0);
	peY = (uint32_t) params.find("arg.pey", 0);
	peZ = (uint32_t) params.find("arg.pez", 0);

    assert( peX * peY * peZ == (unsigned) size() );

    xface_down = -1;
    xface_up = -1;

    yface_down = -1;
    yface_up = -1;

    zface_down = -1;
    zface_up = -1;

	line_a = -1;
        line_b = -1;
        line_c = -1;
        line_d = -1;
        line_f = -1;
	line_e = -1;
        line_g = -1;
        line_h = -1;
        line_i = -1;
        line_j = -1;
        line_k = -1;
	line_l = -1;

        corner_a = -1;
        corner_b = -1;
        corner_c = -1;
        corner_d = -1;
        corner_e = -1;
        corner_f = -1;
        corner_g = -1;
	corner_h = -1;

    	int32_t my_Z = 0;
	int32_t my_Y = 0;
	int32_t my_X = 0;

	getPosition(rank(), peX, peY, peZ, &my_X, &my_Y, &my_Z);

	xface_down = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y, my_Z);
        xface_up = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y, my_Z);

        yface_down = convertPositionToRank(peX, peY, peZ, my_X, my_Y - 1, my_Z);
        yface_up = convertPositionToRank(peX, peY, peZ, my_X, my_Y + 1, my_Z);

        zface_down = convertPositionToRank(peX, peY, peZ, my_X, my_Y, my_Z - 1);
        zface_up = convertPositionToRank(peX, peY, peZ, my_X, my_Y, my_Z + 1);

	line_a = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y - 1, my_Z);
        line_b = convertPositionToRank(peX, peY, peZ, my_X, my_Y - 1, my_Z - 1);
        line_c = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y - 1, my_Z);
        line_d = convertPositionToRank(peX, peY, peZ, my_X, my_Y - 1, my_Z + 1);
        line_e = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y, my_Z + 1);
        line_f = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y, my_Z + 1);
        line_g = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y, my_Z - 1);
        line_h = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y, my_Z - 1);
        line_i = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y + 1, my_Z);
        line_j = convertPositionToRank(peX, peY, peZ, my_X, my_Y + 1, my_Z + 1);
        line_k = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y + 1, my_Z);
        line_l = convertPositionToRank(peX, peY, peZ, my_X, my_Y + 1, my_Z - 1);

        corner_a = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y - 1, my_Z + 1);
        corner_b = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y - 1, my_Z + 1);
        corner_c = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y - 1, my_Z - 1);
        corner_d = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y - 1, my_Z - 1);
        corner_e = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y + 1, my_Z + 1);
        corner_f = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y + 1, my_Z + 1);
        corner_g = convertPositionToRank(peX, peY, peZ, my_X - 1, my_Y + 1, my_Z - 1);
	corner_h = convertPositionToRank(peX, peY, peZ, my_X + 1, my_Y + 1, my_Z - 1);

	// Count up the requests which could be posted (this is a high water mark)
	char requestLength = 0;

	requestLength =  (xface_down > -1) ? 1 : 0;
	requestLength += (xface_up > -1) ? 1 : 0;
	requestLength += (yface_down > -1) ? 1 : 0;
	requestLength += (yface_up > -1) ? 1 : 0;
	requestLength += (zface_up > -1) ? 1 : 0;
	requestLength += (zface_down > -1) ? 1 : 0;

	requestLength += (line_a > -1) ? 1 : 0;
	requestLength += (line_b > -1) ? 1 : 0;
	requestLength += (line_c > -1) ? 1 : 0;
	requestLength += (line_d > -1) ? 1 : 0;
	requestLength += (line_e > -1) ? 1 : 0;
	requestLength += (line_f > -1) ? 1 : 0;
	requestLength += (line_g > -1) ? 1 : 0;
	requestLength += (line_h > -1) ? 1 : 0;
	requestLength += (line_i > -1) ? 1 : 0;
	requestLength += (line_j > -1) ? 1 : 0;
	requestLength += (line_k > -1) ? 1 : 0;
	requestLength += (line_l > -1) ? 1 : 0;

	requestLength += (corner_a > -1) ? 1 : 0;
	requestLength += (corner_b > -1) ? 1 : 0;
	requestLength += (corner_c > -1) ? 1 : 0;
	requestLength += (corner_d > -1) ? 1 : 0;
	requestLength += (corner_e > -1) ? 1 : 0;
	requestLength += (corner_f > -1) ? 1 : 0;
	requestLength += (corner_g > -1) ? 1 : 0;
	requestLength += (corner_h > -1) ? 1 : 0;

    requests.resize( requestLength * 2 );
}

void EmberLuleshGenerator::configsweep()
{
	int32_t myX = 0;
	int32_t myY = 0;

	// Get our position in a 2D processor array
	getPosition(rank(), px, py, &myX, &myY);

	x_up   = (myX != (px - 1)) ? rank() + 1 : -1;
	x_down = (myX != 0) ? rank() - 1 : -1;

	y_up   = (myY != (py - 1)) ? rank() + px : -1;
	y_down = (myY != 0) ? rank() - px : -1;
}

void EmberLuleshGenerator::issue3dsweep(std::queue<EmberEvent*>& evQ){

    // Sweep from (0, 0) outwards towards (Px, Py)

    for(uint32_t i = 0; i < kba; i+= 1) {
        if(x_down >= 0) {
            enQ_recv(evQ, x_down, sweep3d_vol, 1000, GroupWorld );
        }

        if(y_down >= 0) {
            enQ_recv(evQ, y_down, sweep3d_vol, 1000, GroupWorld );
        }

        // enQ_compute( evQ, nsCompute );

        if(x_up >= 0) {
            enQ_send( evQ, x_up, sweep3d_vol, 1000, GroupWorld );
        }

        if(y_up >= 0) {
            enQ_send( evQ, y_up, sweep3d_vol, 1000, GroupWorld );
        }
    }
}

void EmberLuleshGenerator::issue3D26stencil(std::queue<EmberEvent*>& evQ)
{
    int nextRequest = 0;

    if(xface_down > -1) {
        enQ_irecv( evQ, xface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(xface_up > -1) {
        enQ_irecv( evQ, xface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(yface_down > -1) {
        enQ_irecv( evQ, yface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(yface_up > -1) {
        enQ_irecv( evQ, yface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(zface_down > -1) {
        enQ_irecv( evQ, zface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(zface_up > -1) {
        enQ_irecv( evQ, zface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_a > -1) {
        enQ_irecv( evQ, line_a, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_b > -1) {
        enQ_irecv( evQ, line_b, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_c > -1) {
        enQ_irecv( evQ, line_c, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_d > -1) {
        enQ_irecv( evQ, line_d, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_e > -1) {
        enQ_irecv( evQ, line_e, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_f > -1) {
        enQ_irecv( evQ, line_f, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_g > -1) {
        enQ_irecv( evQ, line_g, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_h > -1) {
        enQ_irecv( evQ, line_h, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_i > -1) {
        enQ_irecv( evQ, line_i, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_j > -1) {
        enQ_irecv( evQ, line_j, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_k > -1) {
        enQ_irecv( evQ, line_k, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_l > -1) {
        enQ_irecv( evQ, line_l, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_a > -1) {
        enQ_irecv( evQ, corner_a, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_b > -1) {
        enQ_irecv( evQ, corner_b, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_c > -1) {
        enQ_irecv( evQ, corner_c, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_d > -1) {
        enQ_irecv( evQ, corner_d, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_e > -1) {
        enQ_irecv( evQ, corner_e, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_f > -1) {
        enQ_irecv( evQ, corner_f, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_g > -1) {
        enQ_irecv( evQ, corner_g, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_h > -1) {
        enQ_irecv( evQ, corner_h, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }


    // Enqueue the sends
    if(xface_down > -1) {
        enQ_isend( evQ, xface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(xface_up > -1) {
        enQ_isend( evQ, xface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(yface_down > -1) {
        enQ_isend( evQ, yface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(yface_up > -1) {
        enQ_isend( evQ, yface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(zface_down > -1) {
        enQ_isend( evQ, zface_down, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(zface_up > -1) {
        enQ_isend( evQ, zface_up, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_a > -1) {
        enQ_isend( evQ, line_a, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_b > -1) {
        enQ_isend( evQ, line_b, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_c > -1) {
        enQ_isend( evQ, line_c, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_d > -1) {
        enQ_isend( evQ, line_d, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_e > -1) {
        enQ_isend( evQ, line_e, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_f > -1) {
        enQ_isend( evQ, line_f, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_g > -1) {
        enQ_isend( evQ, line_g, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_h > -1) {
        enQ_isend( evQ, line_h, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_i > -1) {
        enQ_isend( evQ, line_i, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_j > -1) {
        enQ_isend( evQ, line_j, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_k > -1) {
        enQ_isend( evQ, line_k, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(line_l > -1) {
        enQ_isend( evQ, line_l, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_a > -1) {
        enQ_isend( evQ, corner_a, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_b > -1) {
        enQ_isend( evQ, corner_b, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_c > -1) {
        enQ_isend( evQ, corner_c, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_d > -1) {
        enQ_isend( evQ, corner_d, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_e > -1) {
        enQ_isend( evQ, corner_e, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_f > -1) {
        enQ_isend( evQ, corner_f, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_g > -1) {
        enQ_isend( evQ, corner_g, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    if(corner_h > -1) {
        enQ_isend( evQ, corner_h, nn3d26_vol, 0, GroupWorld, &requests[nextRequest]);
        nextRequest++;
    }

    // Enqueue a wait all for all the communications we have set up
    enQ_waitall( evQ, nextRequest, &requests[0], NULL );

}