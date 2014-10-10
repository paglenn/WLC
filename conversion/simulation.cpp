#include "engine.h"
#include<vector>
#include<ctime>

int main() {

	int start_time = time(0); // track running time 
	init(); 
	for(int i = 0 ; i < numWindows; i++) { 
		zvals[i] = new TH1D("","",10,0.,1.); 
		zvals[i]->SetBit(TH1::kCanRebin); 
	}
	long double pct_acc = 0.; // track num accepted moves 
	

	for(int win = 0; win < numWindows; win++) {
		reset();

		for(int wpass = 0; wpass < numPasses; wpass++) {
			//std::cout<<"adjusting"<<std::endl;
			adjustZ(win); 
			//std::cout<<"adjusted"<<std::endl;

			for(int step = 1; step <= numSteps; step++) { 
				int acc = umbrella_mc_step(win); 
				pct_acc += acc*100./numFrames;  
				//writeCosines();
				if(step > equilTime and step%corrTime==0) writeZ(win,step); 
				//if(step > equilTime and step%mbarFreq ==0) writeZFile(win); 
				if(step > equilTime and step%mbarFreq ==0) zvals[win] -> Fill(getZ());  
				if(step%checkSumFreq == 0) checkNorms(); // make sure |t| == 1 for all t[j]    

				if(step%progressFreq == 0) { 
					char progress[50]; 
					sprintf(progress,"window %d/%d\tpass %d/%d\tstep %g\n",win+1,numWindows,wpass+1,numPasses,(double) step);
					fputs(progress,progressFile); 
				}

			} // end run 

		} //end window 
	
	} // end sim 

	// write simulation summary
	char ratio[50];
	int end_time = time(0); 
	double tdiff = (end_time - start_time)/float(60); // in minutes  
	//sprintf(ratio,"%.1f%% steps accepted\n",100*pct_acc/float(numFrames));
	sprintf(ratio,"%.1Lf%% steps accepted\nrunning time: %.2f minutes\n",pct_acc,tdiff);
	fputs(ratio,progressFile); 
	write_metadata(); 
	
	return cleanup(); 
}



