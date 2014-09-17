#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); 
	init(); 
	int num_acc = 0; 
	
	for(int j = 0; j < numSteps; j++) { 
		int acc = umbrella_mc_step(); 
		num_acc += acc; 
		//std::cout<<acc<<std::endl;
		WriteEventData();	
		// Write progress report 
		char progress[50]; 
		sprintf(progress,"step %d/%d\n",j,numSteps);
		fputs(progress,progressFile);
	} // end simulation 

	// write statistics and close up 
	writeHistograms();
	char ratio[50];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // in minutes  
	//sprintf(ratio,"%.1f%% steps accepted\n",100*num_acc/float(numFrames));
	sprintf(ratio,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",100*(num_acc/float(numSteps)),tdiff);
	fputs(ratio,progressFile); 
	return cleanup(); 
}



