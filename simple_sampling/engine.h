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
//using namespace std; 

// system parameters 
double L, delta; 
int N; 
std::vector<double> tx, ty, tz, t0, t_old ; 

// simulation parameters 
unsigned long long int numSteps , numSweeps, numFrames; 
int numPasses, numWindows; 
unsigned long int sampleRate, progressRate, equilibrationTime; 
int numBins ; 
double binWidth, binOverlap; 
double gaussVar; 
std::vector<double> winMin, winMax, zmin; 
double wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
TFile* DataFile;
TH1D* zHist;
TH1D* rpHist; 
TH1D* tpHist; 
TH1D* rptpHist; 
TH1D* cosHist ; 
TH1D* normsHist; 
TProfile* corr; 
TProfile* zTimeSeries; 

const double PI = acos(-1.);
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
	rpHist = new TH1D(); 
	tpHist = new TH1D(); 
	rptpHist = new TH1D(); 
	cosHist = new TH1D(); 
	normsHist = new TH1D(); 

	corr = new TProfile(); 
	zTimeSeries = new TProfile() ; 

	cosHist->SetBit(TH1::kCanRebin);
	rpHist->SetBit(TH1::kCanRebin);
	tpHist->SetBit(TH1::kCanRebin);
	rptpHist->SetBit(TH1::kCanRebin);
	zHist->SetBit(TH1::kCanRebin);
    corr->SetBit(TH1::kCanRebin); 
	zTimeSeries->SetBit(TH1::kCanRebin); 
	cosHist->SetBins(100,0.9,1.0);
	rpHist->SetBins(100,-0.01,0.01);
	tpHist->SetBins(100,-0.01,0.01);
	rptpHist->SetBins(100,-0.01,0.01);
	zHist->SetBins(numBins,0.7,1.0) ;
    corr->SetBins(100,0.95,1.0); 
}

void init() {

	readInParameters();
	createHistograms(); 

	delta = L/double(N); 
	numSteps = numSweeps * N; 
	//std::cout<<numSweeps << std::endl; 
	//std::cout << N << std::endl;
	//exit(1); 
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

void perturb(int index) { 

	t_old[0] = tx[index]; 
	t_old[1] = ty[index]; 
	t_old[2] = tz[index]; 

	std::normal_distribution<double> gaus(0.0,gaussVar);

	tx[index] += gaus(rng);
	gaus.reset(); // needed for independent rv's 
	ty[index] += gaus(rng);
	gaus.reset();
	tz[index] += gaus(rng);
	gaus.reset();

	normalize(index);
}

double dot(double x1, double y1, double z1, double x2, double y2, double z2) {

	double sum = x1*x2 + y1*y2 + z1*z2; 

	return sum; 
}

double getZ() {
	// height of polymer along z-axis  
	double zsum  = 0.;
	for(int j = 0; j < N; j++) zsum += tz[j] ;
	return zsum/float(N); 
}

double getRP() {
	// projection of eed vector into the plane 
	double rpx = 0.,rpy = 0.;
	for(int j = 0; j < N; j++) {
		rpx += tx[j]; 
		rpy += ty[j];
	}
	double rp = sqrt(rpx*rpx+rpy*rpy);
	return rp/float(N);
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

	
double getTP() { 
	//perpendicular component of final segment 
	double tp = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]); 
	return tp; 
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
	
	//z = getZ(); 
	//dE += 0.5*k*pow(z - zmin[w_index],2.) ; 

	return E_new - E_old ;

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
		double tcos = dot(tx[j],ty[j],tz[j],tx[j+1],ty[j+1],tz[j+1]);  
		cosHist->Fill(tcos); 
	}
}

void computeNorms() { 

	for (int ii = 1; ii < N; ii++ ) {
		normsHist->Fill(tx[ii]*tx[ii] + ty[ii] * ty[ii] + tz[ii]*tz[ii]) ;   
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

	double rptp = getRPTP(); 
	rptpHist->Fill(rptp);
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

bool umbrella_mc_step(int w_index = numWindows-1) {

	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index); 
	//std::cout<<dE<<std::endl;

	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE/delta) ) {

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
	write_RP();
	write_TP();
	write_RPTP();
	writeZ(timesteps); 
}

void writeZHist() {
	int numBins = zHist->GetNbinsX(); 
	for ( int i = 0; i < numBins; i++ ) {

		char histVals[50]; 
		if (zHist->GetBinContent(i) != 0) {
			sprintf(histVals,"%f %f\n",zHist->GetBinCenter(i),-log(zHist->GetBinContent(i))) ; 
			fputs(histVals,zFile) ; 
		}
	}
}

void writeHistograms() {
	// normalize all histograms  
	computeNorms(); 

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
	zTimeSeries->Write("z_ts"); 
	normsHist->Write("t_norm"); 

	writeZHist(); 
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	DataFile->Close();
	return 0 ; 
}

