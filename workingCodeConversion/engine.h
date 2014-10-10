// Umbrella sampling code rewritten from simple sampling , 
// since that was known to be working. Hopefully, if there is an error, 
// WE WILL FIND IT! 
#include<iostream>
#include<cmath>
#include<vector>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<random>
#include<TH1D.h>
#include<TFile.h>
//#include<TCanvas.h>
//#include<TF1.h>
//#include<TPad.h>
//#include<TProfile.h>
//using namespace std; 

int j, step, w, wpass; 

int seed = time(0); 
std::default_random_engine rng (seed);

// simulation parameters 
const int numSteps = 7e6; // 7.5e5; 
int progressFreq = numSteps/10; 
int checkSumFreq = numSteps/5;
int corrTime = 2e4; 
int equilTime = 4e4; 
int mbarFreq = corrTime; 
int numPasses = 1; 
const int numWindows = 100; 
long int numFrames = numWindows * numPasses * numSteps;
double equilibrationTime = 0. ; 
double k[numWindows]; 
double zmin[numWindows]; 


// system parameters 
//double L = 0.1; 
double L = 3./5; 
const int N = 100; 
double ds = L/float(N); // segment length 
double beta = 1./ds; 
double stepVariance = 0.15;  // for sphere sampling 
double adjustingVariance = 0.05; 
//double stepVariance = 0.5;
//double K = 3.25e2/(N*N); 
double K = 2.2e2/(N*N); 
//double K = 10./(N*N); 
double tx[N],ty[N],tz[N]; 

// Data files 
FILE * windowFiles[numWindows];
FILE * mbarFile_u_kln;
FILE * mbarFile_u_kn; 
FILE * progressFile; 
FILE * metaFile; 
//FILE * cosFile; 
//FILE * tpFile; 
//FILE * rpFile;
//FILE * rptpFile; 
TFile DataFile("data.root","recreate");

/* // ROOT histograms
TH1D* zHist = new TH1D(); 
TH1D* rpHist = new TH1D(); 
TH1D* tpHist = new TH1D(); 
TH1D* rptpHist = new TH1D(); 
TH1D* cosHist = new TH1D(); 
TProfile* corr = new TProfile(); 
*/
TH1D * zvals[numWindows]; 
TH1D * normHist = new TH1D(); 
const double PI = acos(-1.);

double wmax,wmin,width, wmean, z ;
double target, tol; 
int random_index ; 
std::vector<double> t0(3,0.),t_old(3,0.); 

/*
void initHists() {
	cosHist->SetBit(TH1::kCanRebin);
	rpHist->SetBit(TH1::kCanRebin);
	tpHist->SetBit(TH1::kCanRebin);
	rptpHist->SetBit(TH1::kCanRebin);
	zHist->SetBit(TH1::kCanRebin);
    corr->SetBit(TH1::kCanRebin); 
	cosHist->SetBins(100,0.9,1.0);
	rpHist->SetBins(100,-0.01,0.01);
	tpHist->SetBins(100,-0.01,0.01);
	rptpHist->SetBins(100,-0.01,0.01);
	zHist->SetBins(100,0.8,1.0) ;
    corr->SetBins(100,0.95,1.0); 
}
*/

// all initialization is here 
void init() {

	//initHists(); 
	normHist->SetBit(TH1::kCanRebin); 
	normHist->SetBins(100,0.95,1.05); 
	t_old[2] = 1. ; 
	t0.at(2) = 1.; 
	mbarFile_u_kln = fopen("ukln.dat","w");
	mbarFile_u_kn = fopen("ukn.dat","w"); 
	srand(seed); 
	//rpFile = fopen("rpvals.dat","w");
	//tpFile = fopen("tpvals.dat","w");
	//rptpFile = fopen("rptpvals.dat","w");
	//cosFile = fopen("cos.dat","w"); 

	progressFile = fopen("progress.out","w"); 

	 char params[100];
	 sprintf(params,"%d windows\n%d passes/window\n%d steps/pass\n",numWindows,numPasses,numSteps);
	 fputs(params,progressFile); 

	for(int j = 0; j < N; j++) {

		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}

	
	for(int j = 0.; j < numWindows; j++ ) {
		zmin[j] = N*(j + 0.5)/float(numWindows);
		k[j] = K; 
		if(j +1 == numWindows ) k[j] = K*5.; 
		std::string fileName = "window_" + std::to_string(j); 
		windowFiles[j] = fopen(fileName.c_str(), "w"); 
	}
	

}
 
void normalizeVector(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}

// for rotation of a segment 
void perturb(int index, double var = stepVariance ) { 

	t_old[0] = tx[index]; 
	t_old[1] = ty[index]; 
	t_old[2] = tz[index]; 

	std::normal_distribution<double> gaus(0.0,var);

	tx[index] += gaus(rng);
	gaus.reset(); // needed for independent rv's 
	ty[index] += gaus(rng);
	gaus.reset();
	tz[index] += gaus(rng);
	gaus.reset();

	normalizeVector(index);
}

// define dot product 
double dot(double x1, double y1, double z1, double x2, double y2, double z2) {

	// correct for any numerical imprecision
	//double chk1 = sqrt(x1*x1+y1*y1+z1*z1); 
	//double chk2 = sqrt(x2*x2+y2*y2+z2*z2); 

	return x1*x2 + y1*y2 + z1*z2 ; 
}

double getZ() {
	// height of polymer along z-axis  
	double z = 0.;
	for(int j = 0; j < N; j++) z += tz[j] ;
	return z; 
}

double getRP() {
	// projection of eed vector into the plane 
	double rpx = 0.,rpy = 0.;
	for(int j = 0; j < N; j++) {
		rpx += tx[j]; 
		rpy += ty[j];
	}
	return sqrt(rpx*rpx+rpy*rpy);
}

double getRPTP() {
	// dot product of perpendicular component of final segment and 
	// projection of eed vector into the plane 
	double rx = 0.,ry = 0.; 
	double tpx=0.,tpy=0.;
	tpx = tx[N-1]; tpy = ty[N-1]; 
	for(int i = 1; i < N; i++) {
		rx += tx[i];
		ry += ty[i];
	}

	return (rx*tpx + ry*tpy); 
}

	
double getTP() { 
	//perpendicular component of final segment 
	double tp = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]); 
	return tp; 
}

// bias potential 
double V(double z, int wj) {
	return 0.5*k[wj]*pow(z - zmin[wj],2.);
}


double deltaE(int index, int w_index, double z_prev ) {
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= dot(tx[index+1],ty[index+1],tz[index+1],t_old[0],t_old[1],t_old[2]) ;
		E_new -= dot(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= dot(tx[index-1],ty[index-1],tz[index-1],t_old[0],t_old[1],t_old[2]) ;
		E_new -= dot(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}
	//std::cout<<E_new - E_old << std::endl; 
	
	E_old += V(z_prev,w_index); 
	E_new += V(getZ(),w_index); 
	//std::cout<<E_new - E_old << std::endl; 
	//exit(1); 

	return E_new - E_old ;

}

// reversion for acceptance criterion failure 
double revert(int index) {

	tx[index] = (float) t_old[0]; 
	ty[index] = (float) t_old[1]; 
	tz[index] = (float) t_old[2]; 

	return 0; 

}

// adjust down to z window 
void adjustZ(int w_index) {

	double z = getZ(); 
	double tol = 1.;

	while(fabs(z - zmin[w_index]) > tol) {

		random_index = 1+rand()%(N-1) ; 
		perturb(random_index,adjustingVariance); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target & dz > 0 ) dz = revert(random_index); 
		if(z < target & dz < 0 ) dz = revert(random_index); 

		if(dz!=0) z = z_new; 
	}
}


/*
void write_cosines() {
	for(int j =0; j< N-1; j++) {
		double tcos = dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]);  
		cosHist->Fill(tcos); 
	}
}

void write_TP() {

	double tp = getTP(); 
	tpHist->Fill(tp);
}

void write_RP() {

	double rp = getRP(); 
	rpHist->Fill(rp);
}

void write_RPTP() {

	rptpHist->Fill( getRPTP() );
}

void write_Z() {
	zHist->Fill(getZ()); 
}

void writeCorrelator() {
    int rmax = 10 ; 
	for(int r = 0; r <= rmax; r++) {

		for(int i = 0; i < N-r; i++) {
            
			double c_r = dot(tx[i],ty[i],tz[i],tx[i+r],ty[i+r],tz[i+r]);
            corr->Fill(r,c_r);
        }
    }
}
*/

// after each window , reset all the segment orientations before adjusting down again (decorrelation) 
void reset() {
	srand(time(0)); 
	for(int j = 1; j < N; j++) {
		tx[j] = 0.; 
		ty[j] = 0.; 
		tz[j] = 1.;
	}
}

// Metropolis acceptance routine 
bool umbrella_mc_step(int w_index) {

	double z_old = getZ(); 
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index,z_old); 

	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-beta*dE) ) {
		return (bool) revert(random_index); 
	}
	

	return true; 

}

// get total energy for MBAR (unbiased) 
double getE() { 
	float E = 0. ; 
	for (int j = 0; j+1<N; j++) {
		E -= dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]); 
	}
	return E; 
}

double getBiasedE(int wi) { 

	return getE() + V(getZ(),wi) ; 
}


// write input for MBAR/uwham 
void writeZFile(int wi) {

    double z = getZ() ; 
    char zString[20]; 
	sprintf(zString,"%d ",wi);
    fputs(zString,mbarFile_u_kln);

	double E = getE();

    // compute u_kln 
	for(int w = 0; w < numWindows; w++) {
        char u_kln[10]; 
        sprintf(u_kln,"%g ",beta*(E + V(z,w)) );
        fputs(u_kln,mbarFile_u_kln); 
    }

    char newline[3];
    sprintf(newline,"\n");
    fputs(newline,mbarFile_u_kln); 

	char u_kn[10]; 
	sprintf(u_kn,"%g\n",beta*V(z,wi)); 
	fputs(u_kn,mbarFile_u_kn); 
}


void writeZ(int currWin,int timeStep) {
	double z = getZ(); 
	char zString[20];
	sprintf(zString,"%d\t%f\t%f\n",timeStep,z,getBiasedE(currWin));
	fputs(zString,windowFiles[currWin]);
}

void write_metadata() {
	metaFile = fopen("metadata.dat","w") ;
	//fputs("#window_file\tz_min\tk\n",metaFile);
	for(int j = 0; j < numWindows; j++) {
		char data[100];
		//sprintf(data,"/home/pglenn/test_mk3/window_%d\t%f\t%f\n",j,zmin[j],k[j]);
		sprintf(data,"/Users/paulglen/github/actin/workingCodeConversion/window_%d\t%f\t%f\t%d\t%f\t\n",j,zmin[j],k[j],corrTime,ds);
		// window harmonic_min harmonic_const. correl_time temperature 
		fputs(data,metaFile);
	}
	fclose(metaFile); 
	
}

void checkNorms() { 
	
	for ( int jj = 0; jj < N; jj++ ) {
		normHist->Fill(tx[j]*tx[j] + ty[j]*ty[j] + tz[j]*tz[j]); 
	}
}

/*
void WriteEventData() {

	write_cosines();
    writeCorrelator();
	write_RP();
	write_TP();
	write_RPTP();
	write_Z(); 
}

void writeHistograms() {
	// normalize all histograms  
	rpHist -> Scale(1./rpHist->Integral("width")) ; 
	tpHist -> Scale(1./tpHist->Integral("width")) ; 
	rptpHist -> Scale(1./rptpHist->Integral("width") ) ; 
	cosHist -> Scale(1./cosHist->Integral("width")) ; 
	zHist->Scale(1./zHist->Integral("width")) ;
	
	rpHist->Write("rp");
	rptpHist->Write("rptp");
	tpHist->Write("tp");
	cosHist->Write("dcos");
	zHist->Write("z"); 
    corr->Write("corr");
}
*/
int cleanup() {
	fclose(mbarFile_u_kln); 
	fclose(mbarFile_u_kn); 
	fclose(progressFile);
	//fclose(cosFile);
	//DataFile.Write(); 
	normHist->Write("norms");

	for(int j = 0; j < numWindows; j++ ) { 
		std::string name = "zvals_" + std::to_string(j); 
		zvals[j]->Write(name.c_str()) ;  
		fclose(windowFiles[j]); 
	}
	DataFile.Close();
	return 0 ; 
}


/*
 * checksum to make sure the code was working 
void check() {

	std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	std::cout<<"norm: "<< (tx[j]*tx[j] + ty[j]*ty[j] + tz[j]*tz[j]) << std::endl;

}
*/
