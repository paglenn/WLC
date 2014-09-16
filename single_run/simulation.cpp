#include "engine.h"
#include<vector>

int main() {

	init(); 
	int num_acc = 0; 
	
	for(int j = 0; j < numSteps; j++) { 
		int acc = umbrella_mc_step(numWindows-1); 
		if(acc) num_acc += 1; 
		write_cosines();
		write_RP();
		write_TP();
		write_RPTP();
	}
	std::cout<<100*num_acc/float(numSteps)<< "%"<<std::endl;

	return 0; 
}



