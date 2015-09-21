#include "engine.h"
#include<vector>
#include<ctime>
using std::cout ;
using std::endl  ;

int main(int argc, char* argv[]) {

	//cout << "Hi there " << endl ; 
	if (argc < 4 ) { 
		cerr << "usage: ./wlc RP TP RPTP " << endl ; 
		return 0 ; 
	}

	int start_time = time(0); 
	double rp_in = atof( argv[1]) ; 
	double tp_in = atof(argv[2]) ; 
	double rptp_in = atof(argv[3]) ; 


	init(); 
	adjustTP(tp_in) ; 
	adjustRP(rp_in) ; 
	
	char progress[75] ; 
	alignRPTP(rptp_in) ;

	double pct_acc = 0; 
	sprintf(progress, "RP = %g \t TP = %g \t RPTP = %g ", getRP(),getTP(),getRPTP() ) ; 
	fputs(progress,progressFile) ;


	for(int wi = 0; wi < numWindows; wi ++ ) {

		for(int wpass = 0; wpass < numPasses; wpass ++) {


			if( numWindows > 1 ) {
				int max_reached = adjustZ(wi); 
				set_zero() ; 
				if ( max_reached) {
					wi = numWindows; 
					wpass = numPasses; 
					break ; 
				}
			}

			for(int j = 0; j < numSteps; j++) { 


				int acc = umbrella_mc_step(wi); 
				pct_acc += acc / float(numFrames) ; 
				
				
				if(j > equilibrationTime and j%sampleRate == 0  ) WriteEventData(j,wi);	
				
				
				// Write progress report 
				if(j%progressRate == 0 ) { 
					
					sprintf(progress,"window %d/%d\tpass %g/%g\tstep %g/%g\n",wi+1,numWindows,(float) wpass+1,(float) numPasses,(float) j , (float) numSteps);
					fputs(progress,progressFile);
				}
			checkNorms(); 

			} // end simulation 
		
		} // end pass 

	}	// end window 


	// write statistics and close up 
	writeHistograms();
	write_metadata(); 
	
	sprintf(progress, "Z = %g \t RP = %g \t TP = %g \t RPTP = %g ", getZ() ,getRP(),getTP(),getRPTP() ) ; 
	fputs(progress,progressFile) ;

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",100*(pct_acc),tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



