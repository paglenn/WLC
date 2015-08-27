#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<cstdlib>
#include<random>
#include<TH1D.h>
//#include<TFile.h>
#include<ctime>
#include"VectorMethods.h"
using std::vector; 
using std::endl; 
using std::cout;
using std::cerr; 
using std::normal_distribution;
using std::default_random_engine;  

// system parameters 
double L, delta, q; 
int N; 
vector<double> tx, ty, tz, t0, t_old ; 

// simulation parameters 
int numSweeps, numSteps; 
int nw; 
int sampleRate, progressRate, equilibrationTime; 
int numBins ; 
double binWidth, binOverlap; 
double gaussVar, bias ; 
vector<double> winMin, winMax, zmin,K ; 
double zstart, wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;
int nb; 

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
FILE * metaFile; 
vector<FILE*> windowFiles; 
vector<FILE*> histograms; 
vector<TH1D*> zHists; 
//TFile* DataFile;
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
		char str[100]; 
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
		else if (!strcmp(str,"numWindows")) fscanf(inFile,"%d",&nw); 
		else if (!strcmp(str,"zstart")) fscanf(inFile,"%lf",&zstart); 
		else if (!strcmp(str,"bias")) fscanf(inFile,"%lf",&bias); 
		else if (!strcmp(str,"numBins")) fscanf(inFile,"%d",&numBins); 
		//else std::cout<<"whoami" << str<<std::endl;
	}

}

void createHistograms() {
	//DataFile = new TFile("data.root","recreate"); 
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
	for(int _wi = 0.; _wi < nw; _wi++) {
		for (int _rp = 0; _rp < nb; _rp++) {  
			for( int _tp = 0 ; _tp < nb ; _tp++ ) { 
				for (int _rptp = 0 ; _rptp < nb ; _rptp++ ) {
					//cout << "inner loop" << endl; 
					char hname[50]; 
					sprintf(hname,"z_%d_%d_%d_%d",_wi,_rp, _tp, _rptp );
	//				cout << zHists.size() << endl;
	//				cout << _wi*nb*nb*nb + _rp*nb*nb + _tp *nb + _rptp << endl; 
	//				cout << "printed" << endl; 
					zHists.push_back( new TH1D(hname,hname,numBins,0.,N) ); 
					zHists[_wi*nb*nb*nb + _rp*nb*nb + _tp * nb + _rptp ]->SetBit(TH1::kCanRebin); 
				}
			}
		}
	}
	//exit(1);
}

void init(int num_conditional_bins ) {

	readInParameters();
//	std::cout<<"Initialized" << std::endl;
	//cout << "parameters read" << endl; 

	delta = L/double(N); 
	q  = 1./delta; 
	numSteps = numSweeps * N; 
	gaussVar *= delta; 
	binWidth = (1. - zstart)/numBins; 
	binOverlap = 1 * binWidth; 
	srand(iseed); 
	generator.seed(iseed); 
	nb = num_conditional_bins ; 

	zFile = fopen("uwham.dat","w");

	progressFile = fopen("progress.out","w"); 

	t0.assign(3,0.); 
	t0[2] = 1.; 
	t_old = t0;// shallow copy -- will need to change, given class structure  
	tx.assign(N,0.); 
	ty.assign(N,0.); 
	tz.assign(N,1.);

	for(int j = 0.; j < nw; j++ ) {
		double wmax =N 	* ( zstart + (1-zstart) * double(j+1)/nw) ;
		double wmin =N	* ( zstart + (1-zstart) * double(j)/nw) ;
		//binOverlap = (wmax - wmin) / numBins; 
		//std::cout<<wmin<<"\t"<<wmax<<std::endl; 
		winMin.push_back( wmin ) ; 
		winMax.push_back( wmax ) ;
		double min_loc = 0.5*(wmax + wmin);
		zmin.push_back(min_loc); 
		//if(j +1==nw) winMax[j] = N; 

		K.push_back(bias); 
		for (int _rp = 0; _rp  < nb ; _rp++) { 
			for (int _tp = 0; _tp  < nb ; _tp++) { 
				for (int _rptp = 0; _rptp  < nb ; _rptp++) { 
					char fileName[30]; 
					sprintf(fileName,"window_%d_%d_%d_%d",j, _rp, _tp, _rptp); 
					windowFiles.push_back(fopen(fileName,"w")); 
					
					sprintf(fileName,"hist_%d_%d_%d_%d",j, _rp, _tp, _rptp); 
					histograms.push_back(fopen(fileName,"w")); 
				}
			}
		}
	}

	createHistograms(); 
	
}

double V_bias(double z, int wj ) { 
	return 0.5*K[wj]*pow(z - zmin[wj], 2.);
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

double getZ() {
	// height of polymer along z-axis  
	double zsum  = 0.;
	for(int j = 0; j < N; j++) zsum += tz[j] ;
	return zsum; 
}

double getRP(int elem) {
	// projection of eed vector into the plane 
	double rpx = 0.,rpy = 0.;
	
	for(int j = 0; j < N; j++) {
		//rpx += tx[j]; 
		//rpy += ty[j];
		rpx += tx[j]; 
		rpy += ty[j]; 
	}
	double ans = (elem == 0) ? rpx : rpy; 
	return ans/ N; 
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

double getTP(int elem) { 
	//perpendicular component of final segment 
    double ans = (elem == 0) ? tx[N-1] : ty[N-1] ;
	return ans; 
}

	
double getTP() { 
	//perpendicular component of final segment 
	double tp = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]); 
	return tp; 
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

double revert(int index) {

	tx[index] = t_old[0]; 
	ty[index] = t_old[1]; 
	tz[index] = t_old[2]; 

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

void rotate_T(double angle, int dim = 3) {

	vector<double> v(3,0.); 
	v[0] = tx[N-1] ; 
	v[1] = ty[N-1] ; 
	v[2] = tz[N-1] ; 

	if (dim == 3) v = rotate_phi(v,angle) ;
	else v = rotate_2D(v,angle) ; 

	tx[N-1]  = v[0] ; 
	ty[N-1]  = v[1] ; 
	tz[N-1]  = v[2] ; 

}

void adjustRP(double target) {
	
	double rp = -1.; 
	double tol = 1./N;

	while(fabs(rp - target) > tol) {

		random_index = 1+rand()%(N-2) ; 
		perturb(random_index); 
		double rp_new = getRP(); 
		double dRP = rp_new - rp; 

		if(rp > target & dRP > 0 ) dRP = revert(random_index); 
		if(rp < target & dRP < 0 ) dRP = revert(random_index); 

		if(dRP!=0) rp = rp_new; 
	}
}

void adjustTP(double target) {

	double t_z =  sqrt(1. - target*target)  ; 
	double tp = target/sqrt(2); 
	double targetV[3] = {tp,tp,t_z } ; 
	tx[N-1] = targetV[0] ; 
	ty[N-1] = targetV[1]; 
	tz[N-1] = targetV[2]; 
}

void align(double target) {

	double rpx = getRP(0) ; double rpy = getRP(1) ; 
	double tpx = getTP(0) ; double tpy = getTP(1) ; 
	double norm = sqrt( (rpx*rpx + rpy * rpy) * ( tpx*tpx + tpy*tpy) ) ;
	double theta_curr = acos( getRPTP() / norm )  ;
	double target_theta = acos(target) ;
	double dtheta = target_theta - theta_curr ; 
	rotate_T(dtheta, 2 ) ; 
}


bool umbrella_mc_step(int _win, double targetRP) {
	// the _item implies the index of a bin 

	random_index = 1 + rand() % (N-2); // random int (1,N-2) (conditioning on last t )   
	double E_old = getBiasedE(_win) ;  

	perturb(random_index); 
	adjustRP(targetRP) ; // this way adjustment of RP is part of moveset 

	double dE = getBiasedE(_win) - E_old; 


	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE) ) {

		revert(random_index); 
		return false;
	}
	
	return true; 

}



void writeZ(int step, int _wi, int _rp, int _tp, int _rptp ) { 
	double z = getZ(); 
	zHists[_wi*nb*nb*nb + _rp*nb*nb + _tp*nb + _rptp ]->Fill(z); 

}

void writeZFile(int sweep, int _wi, int _rp, int _tp, int _rptp) {
	double z = getZ(); 
	double E = getE(); 
	char zString[40];
	sprintf(zString,"%d\t%.6f\t%.6f\n",sweep,z,E); 
	fputs(zString,windowFiles[_wi*nb*nb*nb + nb *nb* _rp + nb * _tp + _rptp ]);

}

void WriteEventData(int step, int _wi, int _rp, int _tp, int _rptp ) {

	writeZ(step,_wi, _rp, _tp, _rptp); 
	writeZFile(step,_wi, _rp, _tp, _rptp ); 
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

void writeHistograms() {

	for(int j = 0; j < nw; j++ ) {

		for(int i = 1; i <= zHists[j]->GetNbinsX(); i++) {

			double binCenter = zHists[j]->GetBinCenter(i);
			double binContent = zHists[j]->GetBinContent(i); 

			char histVals[25]; 
			if(binContent != 0 ) {
				sprintf(histVals,"%f\t%f\n",binCenter,binContent); 
				fputs(histVals,histograms[j]) ;

			}

		}
		//zHists[j]->Write(); 
	}
}

void write_metadata() {

	metaFile = fopen("metadata.dat","w") ; 

	for (int j = 0 ; j < nw; j++ ) {
		char data[100]; 
		sprintf(data, "/Users/paulglen/github/WLC/umbrellaSampling_harmonic/window_%d\t%f\t%f\t%d\t%f\t\n",j,zmin[j],K[j],0,delta) ; 
		fputs(data,metaFile) ; 
	}
	fclose(metaFile) ;

}

void checkNorms() { 
	// enforce normality of the t_i to within 1e-4 %  
	for(int i = 0; i < N; i++) {
		
		if(abs(dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) - 1.0) > 1e-6) {

			cerr << " vector " << i << "failed to be normal " << std::endl; 
			cout << "instead has magnitude " << dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) << endl;
			exit(1); 
		}
	}
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	//DataFile->Close();
	return 0 ; 
}
