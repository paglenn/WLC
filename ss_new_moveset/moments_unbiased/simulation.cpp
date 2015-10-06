#include "engine.h"
#include<vector>
#include<ctime>
using std::cout ;
using std::endl  ;

int main(int argc, char* argv[]) {

	//cout << "Hi there " << endl ; 
	//
	/*
	if (argc < 4 ) { 
		cerr << "usage: ./wlc RP TP RPTP " << endl ; 
		return 0 ; 
	}
	*/
	int start_time = time(0); 
	double rp_in, tp_in, rptp_in ; 
	rp_in = tp_in = rptp_in = 0.0 ; 
	//double rp_in = atof( argv[1]) ; 
	//double tp_in = atof(argv[2]) ; 
	//double rptp_in = atof(argv[3]) ; 


	init(); 
	//adjustTP(tp_in) ; 
	//adjustRP(rp_in) ; 
	//cout << numSweeps << endl ;  
	//return 0 ; 
	
	char progress[75] ; 
	//alignRPTP(rptp_in) ;

	double pct_acc = 0; 
	sprintf(progress, "RP = %g \t TP = %g \t RPTP = %g ", getRP(),getTP(),getRPTP() ) ; 
	fputs(progress,progressFile) ;


	int blockLength; 
	for(int j = 0; j < numSweeps; j++) { 

		// decide resolution 
		if(j == 0) {
			blockLength = equilibrationTime;  
		}	else { 
			blockLength = sampleRate; 
		}
		
		// carry out a block of steps before recording data 
		for(int i = 0 ; i < blockLength; i++) {
			int acc = mc_step(); 
			pct_acc += acc / (float) sampleRate ; 
		}
		
		
		WriteEventData(j);	
		
		
		// Write progress report 
		if(j%progressRate == 0 ) { 
			
			sprintf(progress,"step %g/%g\n", (float) j , (float) numSweeps);
			fputs(progress,progressFile);
		}
		

	} // end simulation 

	checkNorms(); 

	// write statistics and close up 
	writeHistograms();
	write_metadata(); 
	GoldstoneModes() ;
	
	sprintf(progress, "Z = %g \t RP = %g \t TP = %g \t RPTP = %g ", getZ() ,getRP(),getTP(),getRPTP() ) ; 
	fputs(progress,progressFile) ;

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	pct_acc = pct_acc / (double) numSweeps; 
	sprintf(summary,"%.2f%% steps accepted\nrunning time: %.2f minutes\n",100*(pct_acc),tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



