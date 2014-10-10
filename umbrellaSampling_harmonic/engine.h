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

// system parameters 
double L, delta; 
int N; 
std::vector<double> tx, ty, tz, t0, t_old ; 

// simulation parameters 
int numSweeps, numSteps; 
int numPasses, numWindows, numFrames; 
int sampleRate, progressRate, equilibrationTime; 
int numBins ; 
double binWidth, binOverlap; 
double gaussVar, bias ; 
std::vector<double> winMin, winMax, zmin,K ; 
double zstart, wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
FILE * metaFile; 
std::vector<FILE*> windowFiles; 
std::vector<FILE*> histograms; 
std::vector<TH1D*> zHists; 
TFile* DataFile;
TH1D* zHist;
TH1D* rpHist; 
TH1D* tpHist; 
TH1D* rptpHist; 
TH1D* cosHist ; 
TH1D* normsHist; 
TProfile* corr; 

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
		else if (!strcmp(str,"numSweeps")) fscanf(inFile,"%d",&numSweeps); 
		else if (!strcmp(str,"sampleRate")) fscanf(inFile,"%d",&sampleRate); 
		else if (!strcmp(str,"progressRate")) fscanf(inFile,"%d",&progressRate); 
		else if (!strcmp(str,"equilibrationTime")) fscanf(inFile,"%d",&equilibrationTime); 
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

	corr = new TProfile(); 

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
	zHist->SetBins(numBins,0.8,1.0) ;
    corr->SetBins(100,0.95,1.0); 

	//std::cout<<"Initialized" << std::endl;
	for(int i = 0.; i < numWindows; i++) {
		char hname[10]; 
		sprintf(hname,"z_%d",i);
		zHists.push_back( new TH1D(hname,hname,numBins,0.,N) ); 
		zHists[i]->SetBit(TH1::kCanRebin); 
	}
}

void init() {

	readInParameters();
//	std::cout<<"Initialized" << std::endl;

	delta = L/double(N); 
	numSteps = numSweeps * N; 
	gaussVar *= delta; 
	numFrames = numWindows * numPasses * numSteps ; 
	binWidth = (1. - zstart)/numBins; 
	binOverlap = 1 * binWidth; 
	srand(seed); 

	zFile = fopen("uwham.dat","w");

	progressFile = fopen("progress.out","w"); 
	if ( progressRate > numFrames/1000 ) progressRate = numFrames / 1000 ; 

	t0.assign(3,0.); 
	t0[2] = 1.; 
	t_old = t0;// shallow copy -- will need to change, given class structure  
	tx.assign(N,0.); 
	ty.assign(N,0.); 
	tz.assign(N,1.);

	
	for(int j = 0.; j < numWindows; j++ ) {
		double wmax =N 	* ( zstart + (1-zstart) * double(j+1)/numWindows) ;
		double wmin =N	* ( zstart + (1-zstart) * double(j)/numWindows) ;
		//binOverlap = (wmax - wmin) / numBins; 
		//std::cout<<wmin<<"\t"<<wmax<<std::endl; 
		winMin.push_back( wmin ) ; 
		winMax.push_back( wmax ) ;
		double min_loc = 0.5*(wmax + wmin);
		zmin.push_back(min_loc); 
		if(j +1==numWindows) winMax[j] = N; 

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


double deltaE(int index, int w_index, double z_prev) {
	
	double tx_old,ty_old,tz_old; 
	tx_old = t_old[0]; ty_old = t_old[1]; tz_old = t_old[2]; 
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= dot(tx[index+1],ty[index+1],tz[index+1],tx_old,ty_old,tz_old) ;
		E_new -= dot(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= dot(tx[index-1],ty[index-1],tz[index-1],tx_old,ty_old,tz_old) ;
		E_new -= dot(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}

	double z_new = getZ(); 
	//double dVz = ( z_new < winMin[w_index] || z_new > winMax[w_index] ) ? bigNum : 0.;
	double dVz  = V_bias(z_new,w_index) - V_bias(z_prev,w_index) ; 
	//std::cout<<z_new<<z_prev<<std::endl; exit(1);
	//std::cout<<"dE: " << ( E_new - E_old) /delta << std::endl;
	//std::cout<<"V : "<<  dVz << std::endl; 
	//std::cout<< bias << std:: endl ; exit(1);


	return (E_new - E_old)/delta + dVz;

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

	double z_prev = getZ() ; 
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index, z_prev ); 

	if(dE <= 0.) return true; 
	else if (rand()/double(RAND_MAX) > exp(-dE) ) {

		revert(random_index); 
		return false;
	}
	

	return true; 

}

void writeZ(int step, int wi ) { 
	double z = getZ(); 
	zHists[wi]->Fill(z); 

	//char zval[10]; 
	//sprintf(zval,"%d\t%f\n",step,z); 
	//fputs(zval, windowFiles[wi]); 
}

double getE()  { 
	double E = 0. ; 
	for ( int i = 0; i+1 < N; i++) {
		E -= dot(tx[i],ty[i],tz[i],tx[i+1],ty[i+1],tz[i+1]); 
	}
	return E/delta ; 
}

double getBiasedE(int wi ) { 

	return getE() + V_bias(getZ(),wi) ; 
}

void writeZFile(int wi) {
	double z = getZ(); 
	double E = getE(); 
	char zString[20];
	sprintf(zString,"%f\t%f\n",z,E);
	fputs(zString,windowFiles[wi]);

}

void WriteEventData(int step, int wi ) {

	//write_cosines();
    //writeCorrelator();
	//write_RP();
	//write_TP();
	//write_RPTP();
	writeZ(step,wi); 
	writeZFile(wi); 
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

void reset() { 
	seed = time(0); 
	std::default_random_engine rng (seed); 
	srand(seed); 
	tx.assign(N,0); 
	ty.assign(N,0); 
	tz.assign(N,1); 

	if(getZ() < 1.0 ) {
		std::cout<< "can't normalize? "<< std::endl; exit(1); 
	}

}

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

			//std::cout<<binCenter<<binContent<<std::endl;
		}
		zHists[j]->Write(); 
	}
}

void write_metadata() {

	metaFile = fopen("metadata.dat","w") ; 

	for (int j = 0 ; j < numWindows; j++ ) {
		char data[100]; 
		sprintf(data, "/Users/paulglen/github/actin/umbrellaSampling_harmonic/window_%d\t%f\t%f\t%d\t%f\t\n",j,zmin[j],K[j],0,delta) ; 
		fputs(data,metaFile) ; 
	}
	fclose(metaFile) ;

}

void checkNorms() { 
	for(int i = 0; i < N; i++) {
		
		if(abs(dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) - 1.0) > 1e-6) {

			std::cerr << " vector " << i << "failed to be normal " << std::endl; 
			std::cout << "instead has magnitude " << dot(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]);
			std::cout << std::endl; 
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

