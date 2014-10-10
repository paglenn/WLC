// modules for actin monte carlo sim
#include<iostream>
#include<cmath>
#include<vector>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<random>
#include<ctime>

int j, step, w, wpass; 
const int numSteps = 1e7; 
int corrTime = 100.; //100; 
int equilTime = 2e4; 
int progressFreq = 1e4; 
int mbarFreq = 1e5; // *actual correlation time* 
double L = 0.1; 
const int N = 100; 
double delta = L/double(N); 
double ds = delta;
double tx[N],ty[N],tz[N]; 
int seed = time(0); 
double gauss_var = 0.05;// 0.02*L;
std::default_random_engine rng (seed);

int numPasses = 1; 
int numBins = 1000; 
const int numWindows = 250; 
double binWidth = 1/float(numBins); 
double binOverlap =0; //1*binWidth;  
double q = 1./ds; 
//double k0 = L;
double K = 3.5e2; 
double k[numWindows];
long int numFrames = numWindows * numPasses * numSteps;
double winMin[numWindows],winMax[numWindows],zmin[numWindows];
FILE * zFile; 
FILE * windowFiles[numWindows];
FILE * progressFile; 
FILE * parameterFile;
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

double getZ() {
	double z = 0.;
	for(int j = 0; j < N; j++) z += tz[j] ;
    //if(z < 0) { std::cerr<<"z is negative!"<<std::endl; exit(1);}
	return z/float(N); 
}

void init() {

	t0.at(2) = 1.; 
	t_old.at(2) = 1.; 

	zFile = fopen("uwham.dat","w");
    int seed = time(0);
	srand(seed); 

	progressFile = fopen("progress.out","w"); 
    parameterFile = fopen("parameters.dat","w");
	//cosFile = fopen("cos.dat","w"); 

	for(int j = 0; j < N; j++) {

		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}

    char params[100]; 
    sprintf(params,"%d windows\n%d passes/window\n%d steps/pass\n",numWindows,numPasses,numSteps);
    fputs(params,parameterFile);
	for(int j = 0.; j < numWindows; j++ ) {

		winMin[j] = double(j)/numWindows; 
		winMax[j] = double(j+1)/numWindows;
		//if(j +1==numWindows) winMax[j] = 1.0; // note: this may be important for wham  
		zmin[j] = 0.5*(winMin[j] + winMax[j]) ;
       	k[j] = K; 
        //if(j +1 == numWindows) k[j] = 0.5*k[j]; 
        char winParams[100]; 
        sprintf(winParams,"window %d: winMin = %f\twinMax = %f\tzmin = %f\tk = %f\n",j,winMin[j],winMax[j],zmin[j],k[j]);
        fputs(winParams,parameterFile);
	}
    fclose(parameterFile);
	
	for(int j = 0; j < numWindows; j++) { 
		std::string fileName = "window_" + std::to_string(j);
		windowFiles[j] = fopen(fileName.c_str(),"w"); 
	}

}


void normalize(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}

void perturb(int index, double var = gauss_var) { 

	t_old[0] = tx[index]; 
	t_old[1] = ty[index]; 
	t_old[2] = tz[index]; 

	std::normal_distribution<double> gaus(0.0,var);

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


void check() {

	std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	std::cout<<"norm: "<< (tx[j]*tx[j] + ty[j]*ty[j] + tz[j]*tz[j]) << std::endl;

}

double deltaE(int index, int w_index, double z_prev) {
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

	E_old += 0.5*k[w_index]* pow(z_prev - zmin[w_index],2.); 
	E_new += 0.5*k[w_index]*pow(getZ() - zmin[w_index],2.);
	double dE = E_new - E_old ; 



	//std::cout<<dE/delta<<std::endl;
	//dE += 0.5*k[w_index]*pow(getZ() - zmin[w_index],2.) ; 
    //std::cout<<w_index<<std::endl;
	//std::cout<<dE/delta<<std::endl;
	//exit(1); 
	return dE;

}

double revert(int index) {

	tx[index] = t_old[0]; 
	ty[index] = t_old[1]; 
	tz[index] = t_old[2]; 

	return 0; 

}

void adjustZ(int w_index, int pass) {

	wmax = winMax[w_index] ; 
	wmin = winMin[w_index] ;
	width = wmax - wmin ;
	wmean = 0.5*(wmax + wmin) ; 
	z = getZ() ;
	target = wmean; //+ rnd() * (width/2. + binOverlap); 
	tol = width/10.;
	while(target>1) target = wmean + rnd()*(width/2.) ;
	//std::cout<<(target-wmean)*2./width<<std::endl;

	while(fabs(z - target) > tol) {
       // std::cout<<z<<std::endl;

		random_index = 1+rand()%(N-1) ; 
        if(random_index<1) { std::cerr<<"random number gen is off"<<std::endl; exit(1); }
		perturb(random_index,1); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target & dz > 0 ) dz = revert(random_index); 
		if(z < target & dz < 0 ) dz = revert(random_index); 

	    if(dz!=0) z = z_new; 
	}
}


bool umbrella_mc_step(int w_index) {

	double z_old = getZ();

	//std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 
    if(random_index<1) { std::cerr<<"random number gen is off"<<std::endl; exit(1); }

	kdE = q*deltaE(random_index, w_index, z_old ); 
	//std::cout<< kdE << std::endl;

	if(kdE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-kdE) ) {

		revert(random_index); 
		return false;
	}
    //std::cout<<"t[0][z]: "<<tz[0] << std::endl;

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
float getE() { 
	float E = 0. ; 
	for (int j = 0; j+1<N; j++) {
		E -= dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]); 
	}
	return E; 
}

void writeZFile(int wi) {

    double z = getZ() ; 
	//std::cout<<z<<std::endl;
    char zString[20]; 
	sprintf(zString,"%d ",wi);
    fputs(zString,zFile);
	float E = getE();
    for(int n = 0; n < numWindows; n++) {
		double V = q*E + q*0.5*k[n]*pow(z-zmin[n],2);
		//std::cout<< V << std::endl; exit(1);
        char u_kln[10]; 
        sprintf(u_kln,"%f ",V);
        fputs(u_kln,zFile); 
    }
    char newline[3];
    sprintf(newline,"\n");
    fputs(newline,zFile); 
}


void writeZ(int window,int timeStep) {

	double z = getZ(); 
	char zString[20];
	sprintf(zString,"%d\t%.4f\n",timeStep,z);
	fputs(zString,windowFiles[window]);
}

void write_metadata() {
	metaFile = fopen("metadata.dat","w") ;
	fputs("#window_file\tz_min\tk\n",metaFile);
	for(int j = 0; j < numWindows; j++) {

		char data[100];
		//sprintf(data,"/home/pglenn/test_mk3/window_%d\t%f\t%f\n",j,zmin[j],k[j]);
		sprintf(data,"/Users/paulglen/github/actin/umbrella_sampling/window_%d\t%f\t%f\n",j,zmin[j],k[j]);
		fputs(data,metaFile);
	}
	fclose(metaFile); 
	
}

void reset() {
    int seed = time(0); 
    std::default_random_engine rng (seed);
    srand(seed);
	
	for(int j = 0; j < N; j++) {
		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}
}

void cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	//fclose(cosFile);
	for(int j = 0; j < numWindows; j++) fclose(windowFiles[j]);
}


