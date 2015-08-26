#include "engine.h"
#include<vector>
#include<ctime>
using std::cout ;

int main() {

	int start_time = time(0); 
	init(); 
	double pct_acc = 0; 
	
	for(int wi = 0; wi < numWindows; wi ++ ) {

		for(int wpass = 0; wpass < numPasses; wpass ++) {
			adjustZ(wi); 

			for(int j = 0; j < numSteps; j++) { 

				int acc = umbrella_mc_step(wi); 
				pct_acc += acc / float(numFrames) ; 
				//if ( pct_acc < 49) gaussVar *= 0.9; 
				//else if (pct_acc > 51) gaussVar *= 1.1;  
				
				
				if(j > equilibrationTime and j%sampleRate == 0  ) WriteEventData(j,wi);	
				
				
				// Write progress report 
				if(j%progressRate == 0 ) { 
					
					char progress[50]; 
					//sprintf(progress,"step %g/%g\tpass %g/%g\twindow%d/%d\n",(float) j,(float) numSteps,(float) wpass,(float) numPasses,wi,numWindows);
					sprintf(progress,"window %d/%d\tpass %g/%g\tstep %g/%g\n",wi,numWindows,(float) wpass,(float) numPasses,(float) j , (float) numSteps);
					fputs(progress,progressFile);
				}
			checkNorms(); 

			} // end simulation 
		
		} // end pass 

	}	// end window 


	// write statistics and close up 
	writeHistograms();
	write_metadata(); 

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",100*(pct_acc),tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



