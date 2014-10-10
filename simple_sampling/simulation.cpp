#include "engine.h"
#include<vector>
#include<ctime>
#include<climits>

int main() {

	int start_time = time(0); 
	init(); 
	double pct_acc = 0; 

	//if ( numSteps < 0 ) std::cerr<< " numerical overflow " << std::endl; 
	for(unsigned long long int j = 0; j < numSteps; j++) { 

		int acc = umbrella_mc_step(); 
		pct_acc += 100*acc / float(numSteps) ; 
		
		zTimeSeries->Fill(j, getZ()) ; 
		
		
		if(j%sampleRate == 0 & j > equilibrationTime ) WriteEventData(j);	
		
		
		// Write progress report 
		char progress[50]; 
		sprintf(progress,"step %g/%g\n",(float) j,(float) numSteps);
		if(j%progressRate == 0 ) fputs(progress,progressFile);

	} // end simulation 
	if ( pct_acc < 0 ) std::cerr<< " numerical underflow " << std::endl; 

	// write statistics and close up 
	writeHistograms();

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",pct_acc,tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



