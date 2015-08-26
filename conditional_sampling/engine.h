#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<cstdlib>
#include<random>
#include<TH1D.h>
#include<TFile.h>
#include<ctime>
using std::vector; 
using std::endl; 
using std::cout;
using std::cerr; 
using std::normal_distribution;
using std::default_random_engine;  

// system parameters 
double L, delta, q; 
int N; 
vector<double> tx, ty, tz, t0x,t0y,t0z, tox,toy,toz ; 

// simulation parameters 
int numSweeps, numSteps; 
int numPasses, numWindows, numFrames; 
int sampleRate, progressRate, equilibrationTime; 
int numBins ; 
double binWidth, binOverlap; 
double gaussVar, bias ; 
vector<double> winMin, winMax, zmin,K ; 
double zstart, wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;

// hard conditions on certain variables
double RP, TP, RPTP ; 

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
FILE * metaFile; 
vector<FILE*> windowFiles; 
vector<FILE*> histograms; 
vector<TH1D*> zHists; 
TFile* DataFile;
TH1D* zHist;
TH1D* rpHist; 
TH1D* tpHist; 
TH1D* rptpHist; 
TH1D* cosHist ; 
TH1D* normsHist; 

const double PI = acos(-1.);
int iseed ; 
default_random_engine generator; 

void readInParameters() {

	inFile = fopen("parameters.txt","r"); 
	if (inFile == NULL) { 
		cerr << "parameter file not found!" << endl; 
		exit(1); 
	}
	
	while(!feof(inFile)) {
		char str[50]; 
		fscanf(inFile,"%s",str);
		if ( *str == '#')  fgets(str,100,inFile) ; 
		else if ( *str == 'L') fscanf(inFile,"%lf",&L);
		else if ( *str == 'N') fscanf(inFile,"%d", &N); 
		else if (!strcmp(str,"numSweeps")) fscanf(inFile,"%d",&numSweeps); 
		else if (!strcmp(str,"sampleRate")) fscanf(inFile,"%d",&sampleRate); 
		else if (!strcmp(str,"progressRate")) fscanf(inFile,"%d",&progressRate); 
		else if (!strcmp(str,"equilibrationTime")) fscanf(inFile,"%d",&equilibrationTime); 
		else if (!strcmp(str,"iseed"))	fscanf(inFile, "%d", &iseed ) ;  
		else if (!strcmp(str,"gaussVar")) fscanf(inFile,"%lf",&gaussVar); 
		else if (!strcmp(str,"numWindows")) fscanf(inFile,"%d",&numWindows); 
		else if (!strcmp(str,"numPasses")) fscanf(inFile,"%d",&numPasses); 
		else if (!strcmp(str,"zstart")) fscanf(inFile,"%lf",&zstart); 
		else if (!strcmp(str,"bias")) fscanf(inFile,"%lf",&bias); 
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


	cosHist->SetBit(TH1::kCanRebin);
	rpHist->SetBit(TH1::kCanRebin);
	tpHist->SetBit(TH1::kCanRebin);
	rptpHist->SetBit(TH1::kCanRebin);
	zHist->SetBit(TH1::kCanRebin);
	cosHist->SetBins(100,0.9,1.0);
	rpHist->SetBins(100,-0.01,0.01);
	tpHist->SetBins(100,-0.01,0.01);
	rptpHist->SetBins(100,-0.01,0.01);
	zHist->SetBins(numBins,0.8,1.0) ;

	//std::cout<<"Initialized" << std::endl;
	for(int i = 0.; i < numWindows; i++) {
		char hname[10]; 
		sprintf(hname,"z_%d",i);
		zHists.push_back( new TH1D(hname,hname,numBins,0.,N) ); 
		zHists[i]->SetBit(TH1::kCanRebin); 
	}
}

void normalize(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}

double getZ() {
	// height of polymer along z-axis  
	double zsum  = 0.;
	for(int j = 0; j < N; j++) zsum += tz[j] ;
	return zsum; 
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


void init() {

	readInParameters();
//	std::cout<<"Initialized" << std::endl;

	delta = L/double(N); 
	q  = 1./delta; 
	numSteps = numSweeps * N; 
	gaussVar *= delta; 
	numFrames = numWindows * numPasses * numSteps ; 
	binWidth = (1. - zstart)/numBins; 
	binOverlap = 1 * binWidth; 
	srand(iseed); 
	generator.seed(iseed); 

	zFile = fopen("uwham.dat","w");

	progressFile = fopen("progress.out","w"); 
	if ( progressRate > numFrames/1000 ) progressRate = numFrames / 1000 ; 

	tx.assign(N,0.); 
	ty.assign(N,0.); 
	tz.assign(N,1.);
	t0x.assign(N,0.); 
	t0y.assign(N,0.); 
	t0z.assign(N,1.);
	tox.assign(N,0.); 
	toy.assign(N,0.); 
	toz.assign(N,1.);
	// initialize to a random configuration 
	for (int index = 0 ; index < N ; index++ ) { 

		normal_distribution<double> gaus(0.0,1);

		tx[index] = gaus(generator);
		gaus.reset(); // needed for independent rv's 
		ty[index] = gaus(generator);
		gaus.reset();
		tz[index] = gaus(generator);
		gaus.reset();

		normalize(index);
		t0x[index] = tx[index]; 
		tox[index] = tx[index];
		t0y[index] = ty[index]; 
		toy[index] = ty[index]; 
		t0z[index] = tz[index];
		toz[index] = tz[index];
	}

	RP = getRP() ; 
	RPTP = getRPTP() ; 
	TP = getTP() ; 

	
	for(int j = 0.; j < numWindows; j++ ) {
		double wmax =N 	* ( zstart + (1-zstart) * double(j+1)/numWindows) ;
		double wmin =N	* ( zstart + (1-zstart) * double(j)/numWindows) ;
		//binOverlap = (wmax - wmin) / numBins; 
		//std::cout<<wmin<<"\t"<<wmax<<std::endl; 
		winMin.push_back( wmin ) ; 
		winMax.push_back( wmax ) ;
		double min_loc = 0.5*(wmax + wmin);
		zmin.push_back(min_loc); 
		//if(j +1==numWindows) winMax[j] = N; 

		K.push_back(bias); 

		char fileName[15]; 
		sprintf(fileName,"window_%d",j); 
		windowFiles.push_back(fopen(fileName,"w")); 
		
		sprintf(fileName,"hist_%d",j); 
		histograms.push_back(fopen(fileName,"w")); 
	}

	createHistograms(); 
	
}

double V_bias(double z, int wj ) { 
	return 0.5*K[wj]*pow(z - zmin[wj], 2.);
}
/*
void normalize(int index) {

	double norm = sqrt(tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] );
	tx[index] /= norm; 
	ty[index] /= norm; 
	tz[index] /= norm; 

}
*/
void perturb(int index) { 

	tox[index] = tx[index]; 
	toy[index] = ty[index]; 
	toz[index] = tz[index]; 

	normal_distribution<double> gaus(0.0,gaussVar);

	tx[index] += gaus(generator);
	gaus.reset(); // needed for independent rv's 
	ty[index] += gaus(generator);
	gaus.reset();
	tz[index] += gaus(generator);
	gaus.reset();

	normalize(index);
}

double dot(double x1, double y1, double z1, double x2, double y2, double z2) {

	double sum = x1*x2 + y1*y2 + z1*z2; 

	return sum; 
}


double deltaE(int index, int w_index, double z_prev) {
	
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= dot(tx[index+1],ty[index+1],tz[index+1],tox[index],toy[index],toz[index]) ;
		E_new -= dot(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= dot(tx[index-1],ty[index-1],tz[index-1],tox[index],toy[index],toz[index]) ;
		E_new -= dot(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}

	double z_new = getZ(); 
	//double dVz = ( z_new < winMin[w_index] || z_new > winMax[w_index] ) ? bigNum : 0.;
	double dVz  = V_bias(z_new,w_index) - V_bias(z_prev,w_index) ; 
	//std::cout<<z_new<<z_prev<<std::endl; exit(1);
	//std::cout<<"dE: " << ( E_new - E_old) /delta << std::endl;
	//std::cout<<"V : "<<  dVz << std::endl; 
	//std::cout<< bias << std:: endl ; exit(1);


	return q*E_new- q*E_old + dVz;

}

double revert(int index) {

	tx[index] = tox[index]; 
	ty[index] = toy[index]; 
	tz[index] = toz[index]; 

	return 0; 

}

void adjustZ(int w_index) {

	double wmean = zmin[w_index] ; 
	double z = 0.; 
	double target = wmean ; // (-1,1)
	double tol = delta;

	while(fabs(z - target) > tol) {

		random_index = 1+rand()%(N-1) ; 
		perturb(random_index); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target && dz > 0 ) dz = revert(random_index); 
		if(z < target && dz < 0 ) dz = revert(random_index); 

		if(dz!=0) z = z_new; 
	}

	RP = getRP(); 
	TP = getTP(); 
	RPTP = getRPTP(); 
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

bool diff(double a, double b) { 
	if (fabs(a-b)  > 1e-5) return true; 
	else return false; 
}

// monte carlo move and acceptance/rejection criteria 
bool umbrella_mc_step(int w_index = numWindows-1) {

	double z_prev = getZ() ; 
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index, z_prev ); 

	// check if things have been changed 
	if (diff(RP,getRP()) || diff(TP,getTP()) || diff(RPTP,getRPTP()) ) {
			revert(random_index) ; 
			return false; 
	
	}

	// Metropolis part 
	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE) ) {

		revert(random_index); 
		return false;
	}
	

	return true; 

}


double getE()  { 
	double E = 0. ; 
	for ( int i = 0; i+1 < N; i++) {
		E -= dot(tx[i],ty[i],tz[i],tx[i+1],ty[i+1],tz[i+1]); 
	}
	return q*E ; 
}

double getBiasedE(int wi ) { 

	return getE() + V_bias(getZ(),wi) ; 
}

void writeZ(int step, int wi ) { 
	//cout<<"hi"<<endl;
	double z = getZ(); 
	//double E = getE(); 
	//cout<<"hi"<<endl;
	zHists[wi]->Fill(z); 

	//char zval[30]; 
	//sprintf(zval,"%d\t%f\t%f\n",step,z,E); 
	//fputs(zval, windowFiles[wi]); 
	//cout<<"hi"<<endl;
}

void writeZFile(int step, int wi) {
	double z = getZ(); 
	double E = getE(); 
	char zString[40];
	sprintf(zString,"%d\t%f\t%f\n",step,z,E); 
	fputs(zString,windowFiles[wi]);

}

void WriteEventData(int step, int wi ) {

	//write_cosines();
	//write_RP();
	//write_TP();
	//write_RPTP();
	writeZ(step,wi); 
	writeZFile(step,wi); 
}

void writeZHist() {
	int numBins = zHist->GetNbinsX(); 
	for ( int i = 0; i < numBins; i++ ) {

		char histVals[50]; 
		if (zHist->GetBinContent(i) != 0) {
			float binCenter = zHist->GetBinCenter(i); 
			float logP = log ( zHist->GetBinContent(i) ) ; 
			sprintf(histVals,"%f %f\n",binCenter,-logP) ; 
			fputs(histVals,zFile) ; 
		}
	}
}
/*
void reset() { 

	// must re-seed after a few billion steps 
	iseed = time(0); 
	srand(iseed);
	generator.seed(iseed);
	tx.assign(N,0); 
	ty.assign(N,0); 
	tz.assign(N,1); 

	if(getZ() < 1.0 ) {
		cout<< "can't normalize? "<< endl; 
		exit(1); 
	}

}
*/
void writeHistograms() {

	for(int j = 0; j < numWindows; j++ ) {

		for(int i = 1; i <= zHists[j]->GetNbinsX(); i++) {

			double binCenter = zHists[j]->GetBinCenter(i);
			double binContent = zHists[j]->GetBinContent(i); 

			char histVals[25]; 
			if(binContent != 0 ) {
				sprintf(histVals,"%f\t%f\n",binCenter,binContent); 
				fputs(histVals,histograms[j]) ;

			}

		}
		zHists[j]->Write(); 
	}
}

void write_metadata() {

	metaFile = fopen("metadata.dat","w") ; 

	for (int j = 0 ; j < numWindows; j++ ) {
		char data[100]; 
		sprintf(data, "/Users/paulglen/github/WLC/umbrellaSampling_harmonic/window_%d\t%f\t%f\t%d\t%f\t\n",j,zmin[j],K[j],0,delta) ; 
		fputs(data,metaFile) ; 
	}
	fclose(metaFile) ;

}

void checkNorms() { 
	// enforce normality of the t_i to within 1e-4 %  
	for(int i = 0; i < N; i++) {
		
		if(fabs(dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) - 1.0) > 1e-6) {

			cerr << " vector " << i << "failed to be normal " << std::endl; 
			cout << "instead has magnitude " << dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) << endl;
			exit(1); 
		}
	}
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	DataFile->Close();
	return 0 ; 
}

