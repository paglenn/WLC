#include<iostream>
#include<cmath>
#include<math.h>
#include<vector>
#include<string>
#include<cstdlib>
#include<ctime>
#include"VectorMethods.h"
#include<stdio.h>
#include<cstring>
#include "MersenneTwister.h"
#include<fstream>
#include<iomanip>
using std::vector; 
using std::endl; 
using std::cout; 
using std::cerr; 
using std::ofstream;
using std::setw;

// system parameters 
double L,_L,ds, _ds, _N ; 
int N; 
vector<double> tx, ty, tz, t0x,t0y,t0z, tox,toy,toz ; 
double tp[2], rp[2] ; 

// simulation parameters 
int numSweeps, numSteps; 
int numPasses, numWindows, numFrames; 
int sampleRate, progressRate, equilibrationTime; 
int nbins ; 
double binWidth, binOverlap; 
double gaussVar, bias ; 
vector<double> winMin, winMax, zmin,K ; 
double zstart, wmax,wmin,width, wmean, z ;
double target, tol;
int random_index ;
int totalCounts; 

const int nMoments = 50 ; 
double zMoments[nMoments]; 


// hard conditions on certain variables
bool flag ; 
double RP, TP, RPTP ; 

// Data files 
FILE * inFile; 
FILE * zFile;
FILE * progressFile; 
FILE * metaFile; 
vector<FILE*> windowFiles; 
vector<FILE*> histFiles; 
vector< vector<double> > zHists; 
ofstream momentsFile ; 
//ofstream zFiles[20];

const double PI = acos(-1.);
const double sqrt2 = sqrt(2) ; 
const double eps = 1e-6 ; // for singular things or comparison 
const long int MAX = pow(2,32) -1 ; 
int iseed ; 
//std::mt19937 generator; 
MTRand generator; 


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
		else if (!strcmp(str,"numWindows")) fscanf(inFile,"%d",&numWindows); 
		else if (!strcmp(str,"numPasses")) fscanf(inFile,"%d",&numPasses); 
		else if (!strcmp(str,"zstart")) fscanf(inFile,"%lf",&zstart); 
		else if (!strcmp(str,"bias")) fscanf(inFile,"%lf",&bias); 
		else if (!strcmp(str,"nbins")) fscanf(inFile,"%d",&nbins); 
		//else std:://cout<<"whoami" << str<<std::endl;
	}

	// MF correction / rescaling
	_L = 1./L ; 
	L = 1./(_L + 2) ; 
}

void createHistograms() {

	totalCounts = 0 ; 
	for(int i = 0.; i < numWindows; i++) {
		zHists.push_back( vector<double>(nbins, 0.0) ) ; 
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
	if (isnan(zsum) ) {
		cerr << "z is nan!" << endl ; 
		exit(1) ; 
	}

	return zsum * _N; 
}

double getRP() {
	// projection of eed vector into the plane 
	double ans = sqrt( rp[0]*rp[0] + rp[1]*rp[1]) ; 
	return ans/ (float) N ; 
}

double getTP() { 
	//perpendicular component of final segment 
	double ans = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]);
	return ans ; 
}

double getRPTP() {
	// innerProduct product of perpendicular component of final segment and 
	// projection of eed vector into the plane 

	return (rp[0]*tp[0] + rp[1]*tp[1])/ N ; 
}


double ranf() { 

	double u = 10;
	while(u >= 1 || u < eps ) u = generator.rand() ;
	return u ; 
}

/*
double gaussian (double mean, double var) { 
	
	double ans, u , v; 
	u = ranf() ; 
	v = ranf() ; 

	if (u < eps ) u += eps ; 
	ans = mean + sqrt( - 2*var* log(u) ) * cos(2*PI*v) ; 
	if (isnan(ans)) { 
		cerr <<" Something in gaussian calculation screwed up! " << endl ; 
		exit(1) ; 
	}


	// Box-Muller polar form 
	double x1, x2, w = 1.0; 
	double ans ; 
	while(w >= 1.0) {

		x1 = 2.0 * ranf() - 1.0;
		x2 = 2.0 * ranf() - 1.0;
		 w = x1 * x1 + x2 * x2;
	}
	w = sqrt( (-2.0 * log( w ) ) / w ) ; 
	ans = x1 * w ; 

	while (isnan(ans)) {
		// just in case 
		ans = gaussian(mean, var) ; 
	}

	// note that x2 * w provides another gaussian random number... 
	return ans ; 
	
}
*/

void init() {

	readInParameters();
//	std:://cout<<"Initialized" << std::endl;
//
	_N = 1./N ; 
	ds = L * _N ; 
	_ds  = 1./ds; 
	numSteps = numSweeps * N; 
	////cout << numSteps << endl ;
	gaussVar *= ds; 
	numFrames = numWindows * numPasses * numSteps ; 
	binWidth = (1. - zstart)/nbins; 
	binOverlap = 1 * binWidth; 
	generator.seed(); 

	zFile = fopen("uwham.dat","w");
	momentsFile.open("moments.dat"); 

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
	tp[0] = tp[1] = 0.0 ; 
	rp[0] = rp[1] = 0.0 ; 

	// initialize to a random configuration 
	// don't change first ti from (0,0,1) 
	// limit to ti > 0 at first 
	z = -1;
	while ( z  < 0 ) { 	
	for (int index = 1 ; index < N ; index++ ) { 

		

		tx[index] = generator.randNorm(0.0,1.0);
		ty[index] = generator.randNorm(0.0,1.0);
		while( tz[index] <= 0. ) { 
			tz[index] = generator.randNorm(0.0,1.0) ;
		}

		normalize(index);
		t0x[index] = tx[index]; 
		tox[index] = tx[index];
		t0y[index] = ty[index]; 
		toy[index] = ty[index]; 
		t0z[index] = tz[index];
		toz[index] = tz[index];

		rp[0] += tx[index] ;
		rp[1] += ty[index] ; 

	}
		z = getZ(); 
	}

	tp[0] = tx[N-1] ;
	tp[1] = ty[N-1] ; 

	RP = getRP() ; 
	RPTP = getRPTP() ; 
	TP = getTP() ; 
	////cout << "Starting Z: " << z << endl; 
	////cout << "RP " << RP << endl ; 
	////cout << "TP " << TP << endl; 
	////cout << "RPTP " << RPTP << endl ; 

	
	for(int j = 0.; j < numWindows; j++ ) {
		double wmax = ( zstart + (1-zstart) * double(j+1)/numWindows) ;
		double wmin = ( zstart + (1-zstart) * double(j)/numWindows) ;
		winMin.push_back( wmin ) ; 
		winMax.push_back( wmax ) ;
		double min_loc = 0.5*(wmax + wmin);
		zmin.push_back(min_loc); 

		K.push_back(bias); 

		char fileName[15]; 
		sprintf(fileName,"window_%d",j); 
		windowFiles.push_back(fopen(fileName,"w")); 
		
		sprintf(fileName,"hist_%d",j); 
		histFiles.push_back(fopen(fileName,"w")); 
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

	tx[index] += generator.randNorm(0.0,gaussVar);
	 // needed for independent rv's 
	ty[index] += generator.randNorm(0.0,gaussVar);
	
	tz[index] += generator.randNorm(0.0,gaussVar);
	


	normalize(index);
	rp[0] += (tx[index] - tox[index]) ; 
	rp[1] += (ty[index] - toy[index]) ; 
	
	if (index == N-1 ) { 
		tp[0] = tx[index] ; 
		tp[1] = ty[index] ;
	}
}

double innerProduct(double x1, double y1, double z1, double x2, double y2, double z2) {
	double ans = x1*x2 + y1*y2 + z1*z2; 
	return ans; 
}


double deltaE(int index) {
	
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= innerProduct(tx[index+1],ty[index+1],tz[index+1],tox[index],toy[index],toz[index]) ;
		E_new -= innerProduct(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= innerProduct(tx[index-1],ty[index-1],tz[index-1],tox[index],toy[index],toz[index]) ;
		E_new -= innerProduct(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}

	return _ds*( E_new - E_old) ;

}

double revert(int index) {

	rp[0] += (tox[index] - tx[index]) ; 
	rp[1] += (toy[index] - ty[index]) ; 

	tx[index] = tox[index]; 
	ty[index] = toy[index]; 
	tz[index] = toz[index]; 

	if (index == N-1 ) { 
		tp[0] = tx[index] ; 
		tp[1] = ty[index] ;
	}

	return 0 ; 

}

bool diff(double a, double b) { 
	if (fabs(a-b)  > eps) return true; 
	else return false; 
}

bool conditional_check() {  
	bool ans = false  ; 
	ans = ans || diff(RP,getRP()) ; 
	ans = ans || diff(TP,getTP()) ; 
	ans = ans || diff(RPTP, getRPTP()) ; 
	return false; 
	return ans ; 

}

int adjustZ(int w_index) {

	double wmean = zmin[w_index] ; 
	double z = 0.; 
	double target = wmean ; // (-1,1)
	double tol = 0.25 / numWindows;

	// this is bound to be horribly slow 
	// PROBLEM
	// if this happens then the value of Z is too large to be reached for this
	// RP bin 
	if ( getRP()*getRP() + wmean*wmean > 1 ) return 1 ;  


	while(fabs(z - target) > tol || conditional_check() ) {

		random_index = 1+generator.randInt(N-1) ; 
		perturb(random_index); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target && dz > 0 ) tz[random_index] *= -1 ;  
		if(z < target && dz < 0 ) tz[random_index] *= -1 ; 

		if(dz!=0) z = z_new; 

	}


	return 0 ; 

}

void rotate_T(double angle, int dim = 3) {
	////cout << "Angle " << angle << endl ; 

	vector<double> v(3,0.); 
	v[0] = tx[N-1] ; 
	v[1] = ty[N-1] ; 
	v[2] = tz[N-1] ; 
	//double norm2 = tx[N-1]*tx[N-1]  + ty[N-1] *ty[N-1] ; 
	////cout << "norm before rotation: " << norm2 << endl ; 
//	//cout <<"v: "<<  v[0] << " " << v[1] << endl;

	//if (dim == 3) v = rotate_phi(v,angle) ;
	v = rotate_2D(v,angle) ; 

	//update RP 
	rp[0] += (v[0] - tx[N-1]) ; 
	rp[1] += (v[1] - ty[N-1]) ; 

	tx[N-1]  = v[0] ; 
	ty[N-1]  = v[1] ; 
	tz[N-1]  = v[2] ; 
	//norm2 = tx[N-1]*tx[N-1]  + ty[N-1] *ty[N-1] ; 
	////cout << "norm after rotation: " << norm2 << endl ; 

	tp[0] = tx[N-1] ; 
	tp[1] = ty[N-1] ; 

}

void GoldstoneModes() { 
	// output gradient of angle as a function of position
	ofstream fout; 
	fout.open("modes.dat") ; 
	double dtx, dty, dtz, dt2 ; 
	for(int i = 0 ; i < N-1 ; i++ ) { 

		dtx = tx[i+1] - tx[i];
		dty = ty[i+1] - ty[i];
		dtz = tz[i+1] - tz[i];
		dt2 = dtx*dtx + dty*dty + dtz*dtz ; 
		fout << i ; 
		fout << setw(12) << dt2; 
		fout << endl ; 
	}
}
		


void alignRPTP(double target) {
	// target is mapped to an angle between 0 and PI 
	//double tpx = tp[0] ; double tpy = tp[1] ; 
	double norm = sqrt( (rp[0]*rp[0] + rp[1] * rp[1]) * ( tp[0]*tp[0] + tp[1]*tp[1]) )/N  ;
	double theta_0 ; 
	if (norm < eps)  return; 
	else { 
		double TP_norm = getTP() ; 
		if (rp[0] != 0.0 ) {
			theta_0 = atan(rp[1]/rp[0] ) ; 
			if (rp[0] < 0) {
				theta_0 += PI; 
			}
		} else { 
			if (rp[1] > 0 ) { theta_0 = PI/2.; } 
			else { theta_0 = 3*PI/2. ;  } 
		}
		double theta_t = PI * target ;
		tp[0] = TP_norm * cos(theta_0 + theta_t) ; 
		tp[1] = sqrt(TP_norm*TP_norm - tp[0]*tp[0])  ;
		//tp[1] = TP_norm * sin(theta_0 + theta_t)  ; 	
		tx[N-1] = tp[0] ; 
		ty[N-1] = tp[1] ;
		RP = getRP() ; 
		RPTP = getRPTP() ; 
		TP = getTP() ; 
	}
}

void adjustRP(double target) {
	
	double tol = 10*eps;
	double myRP = getRP() ; 

	while(fabs(myRP- target) > tol) {

		random_index = 1+generator.randInt(N-2) ; 
		perturb(random_index); 
		double rp_new = getRP(); 
		double dRP = rp_new - myRP; 

		if(( myRP > target)  & (dRP > 0)  ) dRP = revert(random_index); 
		if(( myRP < target) & ( dRP < 0) ) dRP = revert(random_index); 

		if(dRP!=0) myRP = rp_new; 
	}

}

void adjustTP(double target) {

	double tp_ip = target/sqrt(2); 
	double t_z =  sqrt(1. - target*target)  ; 
	tp[0] = tp_ip  ; 
	tp[1] = tp_ip  ; 

	rp[0] += (tp[0] - tx[N-1] ) ; 
	rp[1] += (tp[1] - ty[N-1] ) ; 

	tx[N-1] = tp[0] ; 
	ty[N-1] = tp[1] ; 
	tz[N-1] = t_z ; 
}


void set_zero() {

	for(int index = 0 ; index < N ; index++) { 
		
		t0x[index] = tx[index]; 
		tox[index] = tx[index];
		t0y[index] = ty[index]; 
		toy[index] = ty[index]; 
		t0z[index] = tz[index];
		toz[index] = tz[index];
	
	}

	
}

// monte carlo move and acceptance/rejection criteria 
bool mc_step() {

	double z_prev = getZ() ; 
	random_index = 1 + generator.randInt(N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index); 

	// keeping RP, TP, RPTP constant here 
	// possible PROBLEM but maybe not. sampling still happens. 
	if (conditional_check()  ) {
			revert(random_index) ; 
			return false; 
	
	}

	// Metropolis part 
	if(dE <= 0.) return true; 
	else {
		double urn = generator.rand() ; 
		if ( urn > exp(-dE) ) {
			revert(random_index); 
			return false;
		}
	}
	

	return true; 

}


double getE()  { 
	double E = 0. ; 
	for ( int i = 0; i+1 < N; i++) {
		E -= innerProduct(tx[i],ty[i],tz[i],tx[i+1],ty[i+1],tz[i+1]); 
	}
	return _ds*E ; 
}

double getBiasedE(int wi ) { 

	return getE() + V_bias(getZ(),wi) ; 
}

void writeZ(int step) { 
	int wi = 0 ; 
	double z_curr = getZ() ; 
	double bin_index = floor( z_curr * nbins) ; 

	zHists[wi][bin_index] += 1 ; 
	totalCounts += 1 ; 

	for(int i = 0 ; i< nMoments; i++ ) {
		zMoments[i] += pow( z_curr, i+1) ; 
	}
}

void writeZMoments() { 
	int wi = 0 ; 
	for(int i = 0 ; i < nMoments; i++) { 
		double TC = totalCounts ; 
		zMoments[i] /=  TC ; 
		if (TC > 100 ) {
			momentsFile << i+1 ; 
			momentsFile << setw(15) << zMoments[i] ; 
			momentsFile << endl ; 
		}
	}
}
			

void writeZFile(int step) {
	int wi = 0 ; 
	char zString[40];
	sprintf(zString,"%d\t%f\t%f\n",step, getZ() , getE()); 
	//fputs(zString,windowFiles[wi]);

}

void WriteEventData(int step) {

	writeZ(step); 
	writeZFile(step); 
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
		//cout<< "can't normalize? "<< endl; 
		exit(1); 
	}

}
*/
void writeHistograms() {

	writeZMoments() ;

	for(int j = 0; j < numWindows; j++ ) {
		double TC = totalCounts ; 

		for(int i = 0; i < nbins; i++) {
			double binLowEdge = i/(double) nbins ; 
			double binCenter = binLowEdge + 1/(2.*nbins) ; 
			double binContent = zHists[j][i]; 
			double P = binContent / TC ; 

			char histVals[25]; 
			if(binContent >100 ) {
				sprintf(histVals,"%f\t%f\n",binCenter, P ); 
				fputs(histVals,histFiles[j]) ;

			}

		}

	}
}

void write_metadata() {

	metaFile = fopen("metadata.dat","w") ; 

	for (int j = 0 ; j < numWindows; j++ ) {
		char data[100]; 
		sprintf(data, "/Users/paulglen/github/WLC/umbrellaSampling_harmonic/window_%d\t%f\t%f\t%d\t%f\t\n",j,zmin[j],K[j],0,ds) ; 
		fputs(data,metaFile) ; 
	}
	fclose(metaFile) ;

}

void checkNorms() { 
	// enforce normality of the t_i to within 1e-4 %  
	for(int i = 0; i < N; i++) {
		double norm = innerProduct(tx[i],ty[i],tz[i],tx[i],ty[i],tz[i]) ; 	
		if(fabs(norm - 1.0) > 1e-6) {
			cerr << " vector " << i << "failed to be normal " << std::endl; 
			//cout << "instead has magnitude " << norm  << endl;
			exit(1); 
		}
	}
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	momentsFile.close(); 
	return 0 ; 
}

