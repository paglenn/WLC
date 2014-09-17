#include<iostream>
#include<cmath>
#include<vector>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<random>
#include<ctime>
#include<TH1D.h>
#include<TFile.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TPad.h>
//using namespace std; 

int j, step, w, wpass; 
int seed = time(0); 
std::default_random_engine rng (seed);
const int numSteps =1000000000; 
double L = 0.1; 
const int N = 100; 
double delta = L/double(N); 
double gauss_var = pow(1.5*L,2); 
double tx[N],ty[N],tz[N]; 

int numPasses = 1; 
const int numWindows = 1; 
int numBins = 100; 
double binWidth = 1./numBins; 
double binOverlap = 5*binWidth; 
int numFrames = numWindows * numPasses * numSteps;
double winMin[numWindows],winMax[numWindows],zmin[numWindows];
double k = 1./delta; 
FILE * zFile;
FILE * progressFile; 
FILE * cosFile; 
FILE * tpFile; 
FILE * rpFile;
FILE * rptpFile; 
TFile DataFile("data.root","recreate");
TH1D* rpHist = new TH1D(); 
TH1D* tpHist = new TH1D(); 
TH1D* rptpHist = new TH1D(); 
TH1D* cosHist = new TH1D(); 
const double PI = acos(-1.);


double wmax,wmin,width, wmean, z ;
double target, tol; 
int random_index ; 
std::vector<double> t0(3,0.),t_old(3,0.); 

void init() {

	t0.at(2) = 1.; 
	zFile = fopen("uwham.dat","w");
	//rpFile = fopen("rpvals.dat","w");
	//tpFile = fopen("tpvals.dat","w");
	//rptpFile = fopen("rptpvals.dat","w");
	//cosFile = fopen("cos.dat","w"); 
	rpHist->SetTitle("rp");
	tpHist->SetTitle("tp");
	rptpHist->SetTitle("<R_{#perp}, T_{#perp}>");
	cosHist->SetTitle("cosines");
	cosHist->SetBit(TH1::kCanRebin);
	rpHist->SetBit(TH1::kCanRebin);
	tpHist->SetBit(TH1::kCanRebin);
	rptpHist->SetBit(TH1::kCanRebin);
	cosHist->SetBins(100,0.9,1.0);
	rpHist->SetBins(100,-0.01,0.01);
	tpHist->SetBins(100,-0.01,0.01);
	rptpHist->SetBins(100,-0.01,0.01);
	srand(seed); 

	progressFile = fopen("progress.out","w"); 

	for(int j = 0; j < N; j++) {

		tx[j] = 0.; 
		ty[j] = 0.;
		tz[j] = 1.;
	}

	for(int j = 0.; j < numWindows; j++ ) {

		winMin[j] = double(j)/numWindows; 
		winMax[j] = double(j+1)/numWindows + binOverlap;
		zmin[j] = 0.5*(winMin[j] + winMax[j]);
		if(j +1==numWindows) winMax[j] = 1.0; 
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
	double chk1 = sqrt(x1*x1+y1*y1+z1*z1); 
	double chk2 = sqrt(x2*x2+y2*y2+z2*z2); 
	return sum/(chk1*chk2) ; 
}

double getZ() {
	double z = 0.;
	for(int j = 0; j < N; j++) z += delta*tz[j] ;
	return z/L; 
}

double getRP() {
	double rpx = 0.,rpy = 0.;
	for(int j = 0; j < N; j++) {
		rpx += tx[j]; 
		rpy += ty[j];
	}
	double rp = delta*sqrt(rpx*rpx+rpy*rpy);
	return rp/L;
}

double getRPTP() {
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
	double tp = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]); 
	return tp; 
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

	//z = getZ(); 
	//dE += 0.5*k*pow(z - zmin[w_index],2.) ; 

	return dE;

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
	target = wmean + width/2. * (2*rand()/double(RAND_MAX) - 1.); // (-1,1)
	tol = width/4.;

	while(fabs(z - target) > tol) {

		random_index = 1+rand()%(N-1) ; 
		perturb(random_index); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target & dz > 0 ) dz = revert(random_index); 
		if(z < target & dz < 0 ) dz = revert(random_index); 

		if(dz!=0) z = z_new; 
	}
}



void write_cosines() {

	for(int j =0; j< N-1; j++) {
		double tcos = dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]);  
		/*
		char cosString[50]; 
		sprintf(cosString,"%f\n",cos);
		fputs(cosString,cosFile); 
		*/
		cosHist->Fill(tcos); 
	}
}

void write_TP() {

	double tp = getTP(); 
	/*
	char tpString[50]; 
	sprintf(tpString,"%f\n",tp);
	fputs(tpString,tpFile); 
	*/
	tpHist->Fill(tp);
}

void write_RP() {

	double rp = getRP(); 
	/*
	char rpString[50]; 
	sprintf(rpString,"%f\n",rp);
	fputs(rpString,rpFile); 
	*/
	rpHist->Fill(rp);
}

void write_RPTP() {

	double rptp = getRPTP(); 
	/*
	char rptpString[50]; 
	sprintf(rptpString,"%f\n",rptp);
	fputs(rptpString,rptpFile); 
	*/
	rptpHist->Fill(rptp);
}

bool umbrella_mc_step(int w_index = numWindows-1) {

	//std::cout<<tx[1]<<ty[1]<<tz[1]<<std::endl;
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index); 
	//std::cout<< kdE << std::endl;

	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE/delta) ) {

		revert(random_index); 
		return false;
	}
	

	return true; 

}

void WriteEventData() {

	write_cosines();
	write_RP();
	write_TP();
	write_RPTP();
}

void writeHistograms() {
	// normalize all 
	rpHist -> Scale(1./rpHist->Integral("width")) ; 
	tpHist -> Scale(1./tpHist->Integral("width")) ; 
	rptpHist -> Scale(1./rptpHist->Integral("width") ) ; 
	cosHist -> Scale(1./cosHist->Integral("width")) ; 
	
	rpHist->Write("rp");
	rptpHist->Write("rptp");
	tpHist->Write("tp");
	cosHist->Write("dcos");
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	//fclose(cosFile);
	//DataFile.Write(); 
	DataFile.Close();
	return 0 ; 
}





