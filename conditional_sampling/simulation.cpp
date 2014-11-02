#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); 
	double pct_acc = 0; // track percentage of accepted moves 
	unsigned long long int j ; // number of trials gets very high... 
	int acc; 
	char progress[50]; 

	init(); 

	for( j = 0; j < numSteps; j++) { 

		acc = mc_step(); // acc == 0 if move rejected, 1 otherwise 
		pct_acc += 100*acc / float(numSteps) ;
		
		//if (pct_acc < 49) gaussVar /= 10.; 
		//if( pct_acc > 51) gaussVar *= 10.; 
		//if (j% progressRate == 0 ) std::cout<< pct_acc <<" gv" <<  gaussVar << std::endl;
		//zTimeSeries->Fill(j, getZ()) ; 
		
		
		if(j%sampleRate == 0 & j > equilibrationTime ) WriteEventData(j);	
		
		
		// Write progress report 
		sprintf(progress,"step %g/%g\n",(float) j,(float) numSteps);
		if(j%progressRate == 0 ) fputs(progress,progressFile);

	} // end simulation 

	// write statistics and close up 
	writeHistograms();

	char summary[70];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // simulation time in minutes  
	sprintf(summary,"%.1f%% steps accepted\nrunning time: %.2f minutes\n",pct_acc,tdiff);
	fputs(summary,progressFile); 

	return cleanup(); 
}



