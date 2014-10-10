#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); 
	init(); 
	readInParameters();
	double pct_acc = 0; 
	
	for(int j = 0; j < numSteps; j++) { 
		int acc = umbrella_mc_step(); 
		pct_acc += acc / float(numSteps) ; 
		//std::cout<<acc<<std::endl;
		zTimeSeries->Fill(j, getZ()) ; 
		if(j%sampleRate == 0 & j > equilibrationTime ) WriteEventData(j);	
		// Write progress report 
		char progress[50]; 
		sprintf(progress,"step %g/%g\n",(float) j,(float) numSteps);
		if(j%progressRate == 0 ) fputs(progress,progressFile);
	} // end simulation 

	// write statistics and close up 
	writeHistograms();
	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // in minutes  
	//sprintf(summary,"%.1f%% steps accepted\n",100*pct_acc/float(numFrames));
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",100*(pct_acc),tdiff);
	fputs(summary,progressFile); 
	return cleanup(); 
}



