#include<iostream>
#include<cmath>
#include<vector>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<random>
#include<TH1D.h>
#include<TFile.h>
#include<TProfile.h>
using std::cout; 
using std::endl; 
using std::vector; 
//using namespace std; 

// system parameters 
double L, delta, q; 
int N; 
vector<double> tx, ty, tz, t0, t_old ; 

// simulation parameters 
unsigned long long int numSteps , numSweeps, numFrames; 
int numPasses, numWindows; 
unsigned long int sampleRate, progressRate, equilibrationTime; 
int numBins ; 
double binWidth, binOverlap; 
double gaussVar; 
vector<double> winMin, winMax, zmin; 
double wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
TFile* DataFile;
TH1D* zHist;
TH1D* rpxHist; 
TH1D* rpyHist; 
TH1D* tpxHist; 
TH1D* tpyHist; 
TH1D* rptpHist; 
TH1D* cosHist ; 
TH1D* normsHist; 
TProfile* corr; 
TProfile* zTimeSeries; 

int seed = time(0); 
std::default_random_engine rng (seed);

void readInParameters() {

	inFile = fopen("parameters.txt","r"); 
	if (inFile == NULL) { 
		std::cerr << "parameter file not found!" << std:: endl; 
		exit(1); 
	}
	
	while(!feof(inFile)) {
		char str[50]; 
		fscanf(inFile,"%s",str);
		if ( *str == '#')  fgets(str,100,inFile) ; 
		else if ( *str == 'L') fscanf(inFile,"%lf",&L);
		else if ( *str == 'N') fscanf(inFile,"%d", &N); 
		else if (!strcmp(str,"numSweeps")) fscanf(inFile,"%lld",&numSweeps); 
		else if (!strcmp(str,"sampleRate")) fscanf(inFile,"%ld",&sampleRate); 
		else if (!strcmp(str,"progressRate")) fscanf(inFile,"%ld",&progressRate); 
		else if (!strcmp(str,"equilibrationTime")) fscanf(inFile,"%ld",&equilibrationTime); 
		else if (!strcmp(str,"gaussVar")) fscanf(inFile,"%lf",&gaussVar); 
		else if (!strcmp(str,"numWindows")) fscanf(inFile,"%d",&numWindows); 
		else if (!strcmp(str,"numPasses")) fscanf(inFile,"%d",&numPasses); 
		else if (!strcmp(str,"numBins")) fscanf(inFile,"%d",&numBins); 
		//else std::cout<<"whoami" << str<<std::endl;
	}

}

void createHistograms() {
	DataFile = new TFile("data.root","recreate"); 
	zHist = new TH1D(); 
	rpxHist = new TH1D(); 
	rpyHist = new TH1D(); 
	tpxHist = new TH1D(); 
	tpyHist = new TH1D(); 
	rptpHist = new TH1D(); 
	cosHist = new TH1D(); 
	normsHist = new TH1D(); 

	corr = new TProfile(); 
	zTimeSeries = new TProfile() ; 

	cosHist->SetBit(TH1::kCanRebin);
	rpxHist->SetBit(TH1::kCanRebin);
	rpyHist->SetBit(TH1::kCanRebin);
	tpxHist->SetBit(TH1::kCanRebin);
	tpyHist->SetBit(TH1::kCanRebin);
	rptpHist->SetBit(TH1::kCanRebin);
	zHist->SetBit(TH1::kCanRebin);
    corr->SetBit(TH1::kCanRebin); 
	zTimeSeries->SetBit(TH1::kCanRebin); 
	cosHist->SetBins(100,0.9,1.0);
	rpxHist->SetBins(100,-0.01,0.01);
	rpyHist->SetBins(100,-0.01,0.01);
	tpyHist->SetBins(100,-0.01,0.01);
	tpxHist->SetBins(100,-0.01,0.01);
	rptpHist->SetBins(100,-0.01,0.01);
	zHist->SetBins(numBins,0.7,1.0) ;
    corr->SetBins(100,0.95,1.0); 
}

void init() {

	readInParameters();
	createHistograms(); 

	delta = L/double(N); 
	q = 1./delta; 
	numSteps = numSweeps * N; 
	gaussVar *= delta; 
	numFrames = numWindows * numPasses * numSteps ; 
	binWidth = 1./numBins; 
	binOverlap = 5 * binWidth; 
	srand(seed); 

	zFile = fopen("fe_z.dat","w");

	progressFile = fopen("progress.out","w"); 

	t0.assign(3,0.); 
	t0[2] = 1.; 
	t_old = t0; 
	tx.assign(N,0.); 
	ty.assign(N,0.); 
	tz.assign(N,1.);

	
	for(int j = 0.; j < numWindows; j++ ) {
		wmax = double(j+1)/numWindows + binOverlap;
		wmin = double(j)/numWindows;
		winMin.push_back( wmin ) ; 
		winMax.push_back( wmax ) ;
		double min_loc = 0.5*(wmax + wmin);
		zmin.push_back(min_loc); 
		if(j +1==numWindows) winMax[j] = 1.0; 
	}
	
}

void normalize(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}

vector<double> perturb(int index) { 

	t_old[0] = tx[index]; 
	t_old[1] = ty[index]; 
	t_old[2] = tz[index]; 

	std::normal_distribution<double> gaus(0.0,gaussVar);
	
	vector<double> v; 
	for(int i =0; i < 3; i++) {
		v.push_back(gaus(rng)); 
		gaus.reset(); 
	}
	tx[index] += v[0];
	ty[index] += v[1];
	tz[index] += v[2];

	//normalize(index);
	return v ; 

}

void compensate(vector<double> v, int i_exc) {

	for ( int i = 1 ; i < N-1; i++ ) {
		if ( i != i_exc) { 

			tx[i]  -= v[0] / (N-3) ; 
			ty[i] -= v[1] / (N-3 ) ; 
			tz[i] -= v[2] / (N-3 ) ; 
		} 
	}
}



double dot(double x1, double y1, double z1, double x2, double y2, double z2) {

	double sum = x1*x2 + y1*y2 + z1*z2; 

	return sum; 
}

double t_dot(int index_1, int index_2 = -1 ) {
	
	double sum; 
	if (index_2 == -1) {
		sum = tx[index_1]*t_old[0] + ty[index_1]*t_old[1] + tz[index_1] * t_old[2] ;  
	} else {
		sum = tx[index_1]*tx[index_2] + ty[index_1]*ty[index_2] + tz[index_1] * tz[index_2] ; 
	}
	return sum; 
}

double getZ() {
	// height of polymer along z-axis  
	double zsum  = 0.;
	for(int j = 0; j < N; j++) zsum += tz[j] ;
	return zsum/float(N); 
}

double getRP(int elem) {
	// projection of eed vector into the plane 
	double rpx = 0.,rpy = 0.;
	for(int j = 0; j < N; j++) {
		rpx += tx[j]; 
		rpy += ty[j];
	}
	double ans = (elem == 0) ? rpx : rpy; 
	return ans/float(N); 
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

	return (rx*tpx + ry*tpy)/float(N); 
}

	
double getTP(int elem) { 
	//perpendicular component of final segment 
    double ans = (elem == 0) ? tx[N-1] : ty[N-1] ;
	return ans; 
}


double getE() { 
	double E = 0; 
	for( int i = 0; i +1 < N; i++ ) {

		E -= t_dot(i, i+1) ; 
	}

	return q*E; 
}

double revert(int index) {

	tx[index] = t_old[0]; 
	ty[index] = t_old[1]; 
	tz[index] = t_old[2]; 

	return 0; 

}

void adjustZ(int w_index) {

	double wmax = winMax[w_index] ; 
	double wmin = winMin[w_index] ;
	double width = wmax - wmin ;
	double wmean = 0.5*(wmax + wmin) ; 
	double z = 0.; 
	double target = wmean + width/2. * (2*rand()/double(RAND_MAX) - 1.); // (-1,1)
	double tol = width/4.;

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
		//double tcos = dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]);  
		double tcos = t_dot(j,j+1); 
		cosHist->Fill(tcos); 
	}
}

void computeNorms() { 

	for (int ii = 1; ii < N; ii++ ) {
		normsHist->Fill(tx[ii]*tx[ii] + ty[ii] * ty[ii] + tz[ii]*tz[ii]) ;   
	}
}

void writeTP() {

	double tpx = getTP(0); 
	double tpy = getTP(1); 
	tpxHist->Fill(tpx);
	tpyHist->Fill(tpy);
}

void writeRP() {

	double rpx = getRP(0); 
    double rpy = getRP(1); 
    
	rpxHist->Fill(rpx);
    rpyHist->Fill(rpy); 
}

void writeRPTP() {

	double rptp = getRPTP(); 
	rptpHist->Fill(rptp);
}

void writeCorrelator() {
    int rmax = N/10 ; 
	double c_r ; 
	for(int r = 0; r <= rmax; r++) {
		for(int i = 0; i < N-r; i++) {
			//c_r = dot(tx[i],ty[i],tz[i],tx[i+r],ty[i+r],tz[i+r]);
			c_r = t_dot(i,i+r); 
            corr->Fill(r,c_r);
        }
    }
}

bool mc_step(int w_index = 0) {

	random_index = 1 + rand() % (N-2); // random int (1,N-2)  

	double E_old = getE() ; 
	double rp_old_x = getRP(0); 
	double rp_old_y = getRP(1); 
	vector<double> dv = perturb(random_index); 

	compensate( dv, random_index) ; 
	
	for(int i = 0; i +1 < N; i++ ) normalize(i) ; 
	
	double diff_tol = 1e-4; 
	if(fabs ( getRP(0) - rp_old_x) > diff_tol | fabs( getRP(1) - rp_old_y) > diff_tol ) {

		cout <<" for x "; 
		cout << "old: " << rp_old_x; 
		cout << "new: " << getRP(0) ;  
		cout << " for y: " ; 
		cout << "old : " << rp_old_y; 
		cout << "new: " << getRP(1) ; 
		cout << endl; 
		exit(1); 
	}
	
	double dE = getE() - E_old;  

	// Metropolis criterion 
	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE) ) 
	{
		revert(random_index); 
		return false;
	}

	return true; 

}

void writeZ(int timesteps ) { 
	zHist->Fill(getZ()); 
	zTimeSeries->Fill(timesteps,getZ()); 
}

void WriteEventData(int timesteps) {

	write_cosines();
    writeCorrelator();
	writeRP();
	writeTP();
	writeRPTP();
	writeZ(timesteps); 
}

void writeZHist() {
	int numBins = zHist->GetNbinsX(); 
	for ( int i = 0; i < numBins; i++ ) {

		char histVals[50]; 
		if (zHist->GetBinContent(i) != 0) {

			double binCenter = zHist->GetBinCenter(i); 
			double logP = - log( zHist->GetBinContent(i)) ; 
			sprintf(histVals,"%f %f\n",binCenter,logP) ; 
			fputs(histVals,zFile) ; 
		}

	}

}

void writeHistograms() {
	// normalize all histograms  
	computeNorms(); 

	rpxHist -> Scale(1./rpxHist->Integral("width")) ; 
	rpyHist -> Scale(1./rpyHist->Integral("width")) ; 
	tpxHist -> Scale(1./tpxHist->Integral("width")) ; 
	tpyHist -> Scale(1./tpyHist->Integral("width")) ; 
	rptpHist -> Scale(1./rptpHist->Integral("width") ) ; 
	cosHist -> Scale(1./cosHist->Integral("width")) ; 
	zHist->Scale(1./zHist->Integral("width")) ; 

	
	rpxHist->Write("rpx");
	rpyHist->Write("rpy");
	rptpHist->Write("rptp");
	tpxHist->Write("tpx");
	tpyHist->Write("tpy");
	cosHist->Write("dcos");
	zHist->Write("z"); 
    corr->Write("corr");
	zTimeSeries->Write("z_ts"); 
	normsHist->Write("t_norm"); 

	writeZHist(); // do AFTER z hist has been scaled down to a probability dist 

}

int cleanup() {
	// housekeeping 
	fclose(zFile); 
	fclose(progressFile);
	DataFile->Close();
	return 0 ; 
}

