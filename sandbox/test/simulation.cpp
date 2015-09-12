#include "engine.h"
#include<vector>
#include<ctime>
using std::cout ;
using std::endl  ;

int main(int argc, char* argv[]) {

	cout << "Hi there " << endl ; 
	int start_time = time(0); 
	//char* input_rp = argv[1]; 
	//char* input_tp = argv[2] ;
	//char* input_rptp = argv[3] ;
	double rp_in = atof( argv[1]) ; 
	double tp_in = atof(argv[2]) ; 
	double rptp_in = atof(argv[3]) ; 

	cout << rp_in << " " << tp_in << " " << rptp_in  << endl ; 

	for(int i = 0 ; i < 10000 ; i++ ) cout << gaussian(0,1) << endl ;   
	return 0 ; 

	init(); 
	adjustTP(tp_in) ; 
	adjustRP(rp_in) ; 
	
	char progress[50] ; 
	printf("RP = %g \t TP = %g \t RPTP = %g \n", getRP(),getTP(),getRPTP() ) ; 
	alignRPTP(rptp_in) ;
	printf("RP = %g \t TP = %g \t RPTP = %g \n", getRP(),getTP(),getRPTP() ) ; 

	double pct_acc = 0; 
	sprintf(progress, "RP = %g \t TP = %g \t RPTP = %g ", getRP(),getTP(),getRPTP() ) ; 
	fputs(progress,progressFile) ;
	cout << "meow" << endl ; 
	

	set_zero() ; 

	for(int wi = 0; wi < numWindows; wi ++ ) {
		cout << "sim starting " << endl ; 

		for(int wpass = 0; wpass < numPasses; wpass ++) {
			cout << "sim starting " << endl ; 

			if( numWindows > 1) adjustZ(wi); 

			for(int j = 0; j < numSteps; j++) { 
				cout << "sim starting " << endl ; 
				return 0 ; 

				int acc = umbrella_mc_step(wi); 
				pct_acc += acc / float(numFrames) ; 
				
				
				if(j > equilibrationTime and j%sampleRate == 0  ) WriteEventData(j,wi);	
				
				
				// Write progress report 
				if(j%progressRate == 0 ) { 
					
					//char progress[50]; 
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



