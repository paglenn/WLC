#include<iostream>
#include<cmath>
#include<vector>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<random>
#include<ctime>
//using namespace std; 

int j, step, w, wpass; 
int seed = time(0); 
std::default_random_engine rng (seed);
const int numSteps = 50000; 
double L = 0.1; 
const int N = 100; 
double delta = L/double(N); 
double gauss_var = pow(5*delta,2); 
double tx[N],ty[N],tz[N]; 

int numPasses = 5000; 
const int numWindows = 10; 
int numBins = 100; 
double binWidth = 1./numBins; 
double binOverlap = 3*binWidth; 
double q = 1./delta; 
double k[numWindows];
int numFrames = numWindows * numPasses * numSteps;
double winMin[numWindows],winMax[numWindows],zmin[numWindows];
FILE * zFile;
FILE * windowFiles[numWindows];
FILE * progressFile; 
//FILE * cosFile; 
FILE * metaFile;


double wmax,wmin,width, wmean, z ;
double target, tol; 
int random_index ; 
std::vector<double> t0(3,0.),t_old(3,0.); 
double kdE; 

double rnd() {

	double randnum = rand()/float(RAND_MAX); 
	return 2*randnum - 1 ; 
}

void init() {

	t0.at(2) = 1.; 
	zFile = fopen("uwham.dat","w");
	srand(seed); 

	progressFile = fopen("progress.out","w"); 
	//cosFile = fopen("cos.dat","w"); 

	for(int j = 0; j < N; j++) {

		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}

	for(int j = 0.; j < numWindows; j++ ) {

		winMin[j] = double(j)/numWindows; 
		winMax[j] = double(j+1)/numWindows + binOverlap;
		if(j +1==numWindows) winMax[j] = 1.0;  
		zmin[j] = 0.5*(winMin[j] + winMax[j]);
		if( j < 5) { k[j] = delta/10.; }
		else { k[j] = delta;  }
	}
	
	for(int j = 0; j < numWindows; j++) { 
		std::string fileName = "window_" + std::to_string(j);
		windowFiles[j] = fopen(fileName.c_str(),"w"); 
	}

}

void reset() {
	
	for(int j = 0; j < N; j++) {

		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}
}

void normalize(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}

void perturb(int index) { 

	t_old[0] = tx[index]; 
	t_old[1] = ty[index]; 
	t_old[2] = tz[index]; 

	std::normal_distribution<double> gaus(0.0,gauss_var);

	tx[index] += gaus(rng);
	gaus.reset();
	ty[index] += gaus(rng);
	gaus.reset();
	tz[index] += gaus(rng);
	gaus.reset();

	normalize(index);
}

double dot(double x1, double y1, double z1, double x2, double y2, double z2) {

	double sum = x1*x2 + y1*y2 + z1*z2; 
	return sum ; 
}

double getZ() {
	double z = 0.;
	for(int j = 0; j < N; j++) z += delta*tz[j] ;
	return z/L; 
}

void check() {

	std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	std::cout<<"norm: "<< (tx[j]*tx[j] + ty[j]*ty[j] + tz[j]*tz[j]) << std::endl;

}

double deltaE(int index, int w_index) {
	double t0x,t0y,t0z; 
	t0x = t_old[0]; t0y = t_old[1]; t0z = t_old[2]; 
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= dot(tx[index+1],ty[index+1],tz[index+1],t0x,t0y,t0z) ;
		E_new -= dot(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= dot(tx[index-1],ty[index-1],tz[index-1],t0x,t0y,t0z) ;
		E_new -= dot(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}
	
	double dE = (E_new - E_old) ; 

	z = getZ(); 
	dE += 0.5*k[w_index]*pow(z - zmin[w_index],2.) ; 

	return q*dE;

}

double revert(int index) {

	tx[index] = t_old[0]; 
	ty[index] = t_old[1]; 
	tz[index] = t_old[2]; 

	return 0; 

}

void adjustZ(int w_index) {

	wmax = winMax[w_index] ; 
	wmin = winMin[w_index] ;
	width = wmax - wmin ;
	wmean = 0.5*(wmax + wmin) ; 
	z = 0.; 
	target = wmean + width/2. * rnd(); // (-1,1)
	//std::cout<<(target-wmean)*2./width<<std::endl;
	tol = width/8.;

	while(fabs(z - target) > tol) {

		random_index = 1+rand()%(N-1) ; 
		perturb(random_index); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target && dz > 0 ) dz = revert(random_index); 
		if(z < target && dz < 0 ) dz = revert(random_index); 

		if(dz!=0) z = z_new; 
	}
}


bool umbrella_mc_step(int w_index) {

	//std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	kdE = deltaE(random_index, w_index); 
	//std::cout<< kdE << std::endl;

	if(kdE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-kdE) ) {

		revert(random_index); 
		return false;
	}

	return true; 

}
/*
void writeCosines() {

	for(int j =0; j< N-1; j++) {
		double cos = dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]);  
		char cosString[10]; 
		sprintf(cosString,"%.5f\n",cos);
		fputs(cosString,cosFile); 
	}
}
*/
void writeZ(int window,int timeStep) {

	double z = getZ(); 
	char zString[20];
	sprintf(zString,"%d\t%f\n",timeStep,z);
	fputs(zString,windowFiles[window]);
}

void write_metadata() {
	metaFile = fopen("metadata.dat","w") ;
	fputs("#window_file\tz_min\tk\n",metaFile);
	for(int j = 0; j < numWindows; j++) {

		char data[100];
		sprintf(data,"/Users/paulglen/github/actin/conversion/window_%d\t%f\t%f\n",j,zmin[j],k[j]);
		fputs(data,metaFile);
	}
	fclose(metaFile); 
	
}

void cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	//fclose(cosFile);
	for(int j = 0; j < numWindows; j++) fclose(windowFiles[j]);
}

