#include "engine.h"
#include<vector>
#include<ctime>
using std::cout ;

int main() {

	int num_conditional_bins = 10; 
	int start_time = time(0); 
	init(num_conditional_bins); 
	double pct_acc = 0; 
	int numWindows = nw; 
		
	// modulate RP
	for (int r = 0; r < num_conditional_bins; r++ ) { 
		
		target = (r+0.5) / num_conditional_bins; 
		adjustRP( target) ; 
		
		// modulate TP 
		for(int t = 0; t <  num_conditional_bins; t++ ) { 

			adjustTP((t+0.5)/ num_conditional_bins) ; 	

			// Modulate RPTP 
			for( int rt = 0; rt < num_conditional_bins; rt++) { 	

				align(( rt+0.5) / num_conditional_bins) ; 
				// note: we'll have to make sure to preserve this orientation when adjusting RP. 
	

				for(int wi = 0; wi < numWindows; wi ++ ) {
					reset(); 

					adjustZ(wi); 

						for(int j = 0; j < numSweeps; j++) { 

							for(int k = 0; k < N; k++ ) { 

								int acc = umbrella_mc_step(wi,target ); 
								pct_acc += acc / float(N) ; 
								
							}
							
						if(j > equilibrationTime and j%sampleRate == 0  ) WriteEventData(j,wi, r , t, rt);	
							
						// Write progress report 
						if(j%progressRate == 0 ) { 
							char progress[50]; 
							sprintf(progress,"bins %d %d %d;  window %d/%dsweep %g/%g \n",r, t, rt ,wi, numWindows, (float) j , (float) numSweeps);
							fputs(progress,progressFile);
						}

					} // end simulation 

				}	// end window 

			} // end rptp 

		} // end tp 

	} // end rp  
	pct_acc /= (numSweeps * numWindows * nb * nb * nb ); 


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



