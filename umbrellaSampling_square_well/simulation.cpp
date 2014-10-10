#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); 
	init(); 
	double pct_acc = 0; 
	
	for(int wi = 0; wi < numWindows; wi ++ ) {
		reset(); 

		for(int wpass = 0; wpass < numPasses; wpass ++) {
			adjustZ(wi); 

			for(int j = 0; j < numSteps; j++) { 

				int acc = umbrella_mc_step(); 
				pct_acc += acc / float(numFrames) ; 
				
				zTimeSeries->Fill(j, getZ()) ; 
				
				
				if(j%sampleRate == 0 & j > equilibrationTime ) WriteEventData(j,wi);	
				
				
				// Write progress report 
				if(j%progressRate == 0 ) { 
					
					char progress[50]; 
					sprintf(progress,"step %g/%g\tpass %g/%g\twindow%d/%d\n",(float) j,(float) numSteps,(float) wpass,(float) numPasses,wi,numWindows);
					fputs(progress,progressFile);
				}

			} // end simulation 
		
		} // end pass 

	}	// end window 


	// write statistics and close up 
	writeHistograms();

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",100*(pct_acc),tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



