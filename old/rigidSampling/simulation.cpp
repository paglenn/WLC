#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); // track running time 
	init(); 
	double num_acc = 0.;
	
	for(int wi = 0; wi < numWindows; wi++) {
		reset();

		for(int wpass = 0; wpass < numPasses; wpass++) {
			adjustZ(wi,wpass); 
			writeZ(wi,0); 

			for(int step = 1; step <= numSteps; step++) { 
				int acc = umbrella_mc_step(wi); 
				num_acc += 100.*acc/numFrames;  

				if(step%sampleRate==0) writeZ(wi,step); 
				//if(step%sampleRate==0) writeZFile(wi); 

			} // end run 
			char progress[50]; 
			sprintf(progress,"window %d/%d\tpass %d/%d\tstep %d\n",wi+1,numWindows,wpass+1,numPasses,numSteps);
			fputs(progress,progressFile); 

		} //end window 
	
	} // end sim 
	char ratio[50];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // in minutes  
	//sprintf(ratio,"%.1f%% steps accepted\n",100*num_acc/float(numFrames));
	sprintf(ratio,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",num_acc,tdiff);
	fputs(ratio,progressFile); 
	
	//write_metadata(); 
	
	cleanup(); 
	return 0; 
}



