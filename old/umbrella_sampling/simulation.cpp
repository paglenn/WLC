#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); // track running time 
	init(); 
	double pct_acc = 0; // track num accepted moves 
	
	for(int wi = 0; wi < numWindows; wi++) {
		reset();

		for(int wpass = 0; wpass < numPasses; wpass++) {
			adjustZ(wi,wpass); 
			writeZ(wi,0); 
			writeZFile(wi);

			for(int step = 1; step <= numSteps; step++) { 
				int acc = umbrella_mc_step(wi); 
				pct_acc += acc*100./numFrames;  
				//writeCosines();
				if(step > equilTime and step%corrTime==0) writeZ(wi,step); 
				if(step > equilTime and step%mbarFreq ==0) writeZFile(wi); 
				//writeZFile(wi);

				if(step%progressFreq == 0) { 
					char progress[50]; 
					sprintf(progress,"window %d/%d\tpass %d/%d\tstep %g\n",wi+1,numWindows,wpass+1,numPasses,(double) step);
					fputs(progress,progressFile); 
				}

			} // end run 

		} //end window 
	
	} // end sim 
	char ratio[50];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // in minutes  
	//sprintf(ratio,"%.1f%% steps accepted\n",100*pct_acc/float(numFrames));
	sprintf(ratio,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",pct_acc,tdiff);
	fputs(ratio,progressFile); 
	
	write_metadata(); 
	
	cleanup(); 
	return 0; 
}



