#include<iostream>
#include<cassert>
#include<cmath>
#include<math.h>
#include<vector>
#include<string>
#include<cstdlib>
#include<ctime>
#include"VectorMethods.h"
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
double L,_L,ds, _ds, _N , _Nm2; 
int N; 
vector<long double> tx, ty, tz, t0x,t0y,t0z, tox,toy,toz ; 
long double tx_new , ty_new, tz_new ;  
double xi ,yi,zi ; 
double tp[2], rp[2] ; 
double urn; 
long double norm2; 
long double snorm2, _snorm2;
double zmax ; 
double z_mp ;

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

const int nMoments = 100; 
double zMoments[nMoments]; 

double thetaMax ; 
double r2; 
double sx, sy; 
double ans; 


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
vector<double> cosHist; 
int rmax = 15 ;
vector<double> corr(rmax,0.0); // for average correlator  
vector<int> TC_corr(rmax,0); // for total counts  
ofstream fout ; 
ofstream coorFile ; 

//ofstream zFiles[20];

const double PI = acos(-1.);
const double sqrt2 = sqrt(2) ; 
const double eps = 1e-4; // for singular things or comparison 
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
		else if (!strcmp(str,"thetaMax")) fscanf(inFile,"%lf",&thetaMax); 
		else if (!strcmp(str,"numWindows")) fscanf(inFile,"%d",&numWindows); 
		else if (!strcmp(str,"numPasses")) fscanf(inFile,"%d",&numPasses); 
		else if (!strcmp(str,"zstart")) fscanf(inFile,"%lf",&zstart); 
		else if (!strcmp(str,"bias")) fscanf(inFile,"%lf",&bias); 
		else if (!strcmp(str,"nbins")) fscanf(inFile,"%d",&nbins); 
		//else std:://cout<<"whoami" << str<<std::endl;
	}

	// MF correction / rescaling
	cout << L << endl ; 
	_L = 1./L ; 
	L = 1./(_L + 2) ; 
	cout << "Rescaled: " << L << endl ; 
	fclose(inFile) ;
}

bool diff(double a, double b) { 
	// percent difference 
	//if (fabs(a-b)  > eps) return true; 
	if (a == 0 || b ==0 ) { 
		if (fabs(a-b) > eps) return true ; 
	} else { 
		double abs_diff = fabs(a -b ) ;	
		if (abs_diff > 1e-4 ) return true ; 
	}
	return false; 
}

void checkNorms() { 
	// enforce normality of the t_i to within 1e-4 %  
	for(int i = 1; i +1 < N; i++) {
		norm2 = tx.at(i)*tx.at(i) + ty.at(i)*ty.at(i) + tz.at(i)*tz.at(i) ; 
		if(diff(norm2,1.0)) {
			cerr << " vector " << i << "failed to be normal " << std::endl; 
			cerr << "instead has magnitude " << norm2  << endl;
			exit(1); 
		}
	}
}

double innerProduct(double x1, double y1, double z1, double x2, double y2, double z2) {
	ans = x1*x2 + y1*y2 + z1*z2; 
	return ans; 
}


void createHistograms() {

	totalCounts = 0 ; 
	for(int i = 0.; i < numWindows; i++) {
		zHists.push_back( vector<double>(nbins, 0.0) ) ; 
	}
	cosHist = vector<double>(nbins,0.0) ; 
}

void normalize(int index) {

	norm2  = tx[index]*tx[index] + ty[index]*ty[index] + tz[index]*tz[index] ;  
	long double _norm = 1./sqrt( (long double) norm2) ; 
	tx[index] *= _norm; 
	ty[index] *= _norm; 
	tz[index] *= _norm; 

}

double getZ() {
	// height of polymer along z-axis  
	double zsum  = 0.;
	for(int j = 0; j < N; j++) zsum += tz.at(j) ;
	assert(!isnan(zsum)) ;
	return zsum * _N; 
}
/*
void getCorrelator() {
	for(int r = 0; r < rmax; r++) {

		for(int i = 0; i < N-r; i++) {

			corr.at(r) += innerProduct(tx[i],ty[i],tz[i],tx[i+r],ty[i+r],tz[i+r]);
			TC_corr.at(r) += 1 ; 
		}
	}
}
*/

double getRP() {
	// projection of eed vector into the plane 
	//double ans = sqrt( rp[0]*rp[0] + rp[1]*rp[1]) ; 
	//return ans/ (float) N ; 

	double rpf = 0.0;  
	sx = 0.0 ; sy = 0.0 ;  
	for(int i = 1 ; i < N ; i++) { 
		sx += tx.at(i) ; 
		sy += ty.at(i) ; 
	}

	rpf = sqrt(sx*sx + sy*sy ) ;
	return rpf * _N ; 
}

double getTP() { 
	//perpendicular component of final segment 
	ans = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]);
	return ans ; 
}

double getRPTP() {
	// innerProduct product of perpendicular component of final segment and 
	// projection of eed vector into the plane 
	sx = 0.0 ; sy = 0.0 ;  
	for(int i = 1 ; i < N ; i++) { 
		sx += tx.at(i) ; 
		sy += ty.at(i) ; 
	}

	return (sx*tp[0] + sy*tp[1])* _N ; 
}


double ranf() { 

	double u = 10;
	while(u >= 1 || u < eps ) u = generator.rand() ;
	return u ; 
}


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
	zmax = 0.0 ; 

	//coorFile.open("coor.xyz") ; 
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

	//tx[index] += generator.randNorm(0.0,gaussVar);
	// needed for independent rv's 
	//ty[index] += generator.randNorm(0.0,gaussVar);

	//tz[index] += generator.randNorm(0.0,gaussVar);

	//compute the t_new and store separately 
	tx_new =  tx[index] + generator.randNorm(0.0,gaussVar);
	ty_new =  ty[index] + generator.randNorm(0.0,gaussVar);
	tz_new =  tz[index] + generator.randNorm(0.0,gaussVar);

	norm2 = tx_new * tx_new + ty_new *ty_new + tz_new * tz_new ; 
	_snorm2 = 1./sqrt(norm2) ; 
	tx.at(index) = tx_new * _snorm2;
	ty.at(index) = ty_new * _snorm2;
	tz.at(index) = tz_new * _snorm2; 


	//normalize(index);
	rp[0] += (tx[index] - tox[index]) ; 
	rp[1] += (ty[index] - toy[index]) ; 

	if (index == N-1 ) { 
		tp[0] = tx[index] ; 
		tp[1] = ty[index] ;
	}
}


double getE()  { 
	double E = 0. ; 
	for ( int i = 0; i+1 < N; i++) {
		E -= innerProduct(tx[i],ty[i],tz[i],tx[i+1],ty[i+1],tz[i+1]); 
	}
	return _ds * E ; 
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

void revert() { 
	for(int i = 1 ; i+1 < N ; i++) {
		revert(i); 
	}
}



bool conditional_check() {  
	bool ans = false  ; 
	ans = ans || diff(RP,getRP()) ; 
	ans = ans || diff(TP,getTP()) ; 
	ans = ans || diff(RPTP, getRPTP()) ; 
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
	//if (dim == 3) v = rotate_phi(v,angle) ;
	v = rotate_2D(v,angle) ; 

	//update RP 
	rp[0] += (v[0] - tx[N-1]) ; 
	rp[1] += (v[1] - ty[N-1]) ; 

	tx[N-1]  = v[0] ; 
	ty[N-1]  = v[1] ; 
	tz[N-1]  = v[2] ; 
	//norm2 = tx[N-1]*tx[N-1]  + ty[N-1] *ty[N-1] ; 

	tp[0] = tx[N-1] ; 
	tp[1] = ty[N-1] ; 

}

void GoldstoneModes() { 
	// output gradient of angle as a function of position
	//cout << "opening modes file: " << endl ;
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
	//cout << "closing modes file: " << endl ;
	fout.close(); 
	fout.clear();
	//cout << "closed." << endl ; 
}



void alignRPTP(double target) {
	// target is mapped to an angle between 0 and PI 
	//double tpx = tp[0] ; double tpy = tp[1] ; 
	double norm = RP * TP   ; 
	double theta_0 ; 
	cout << "RP before: " << RP << " : " << getRP() << endl ;
	cout << "TP: " << TP << ":" << getTP() << endl ;
	if (norm < eps)  return; 
	else { 
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
		tp[0] = TP * cos(theta_0 + theta_t) ; 
		tp[1] = TP * sin(theta_0 + theta_t)  ; 	
		rp[0] += tp[0] - TP / sqrt2 ; 
		rp[1] += tp[1] - TP / sqrt2 ; 
		tx[N-1] = tp[0] ; 
		ty[N-1] = tp[1] ;
		// this changes RP a little bit 
		RP = getRP() ;


	}
}

void adjustRP(double target) {
	// @param: target RP/N val (in [0,1] )
	_Nm2 = 1./(N-2) ; 

	double txi = N * _Nm2 *target /  sqrt2 - tp[0]*_Nm2 ;

	rp[0] = 0.0 ; 
	rp[1] = 0.0 ; 
	for(int i = 1 ; i +1 < N ; i++) {

		tx.at(i) = txi ; 
		ty.at(i) = txi ;
		tz.at(i) = sqrt(1 - 2*txi*txi) ; 

		rp[0] += txi ; 
		rp[1] += txi ; 
	}

	rp[0] += tp[0] ; 
	rp[1] += tp[1] ; 


	RP = getRP() ; 
	zmax = sqrt( 1 - RP*RP) ;  
	//cout << "target: " << target << endl ; 
	checkNorms(); 

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

	TP = getTP() ;
	//cout << TP << " (actual):(tgt) " << target << endl ;
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

void do_crankshaft() { 

	int i = 1 + generator.randInt(N-6) ; 
	int j = i + 4 - generator.randInt( 2) ; 
	double u, v , w, uvw; 
	//cout << i << ": " << j << endl ;
	u = tx.at(j) - tx.at(i) ; 
	v = ty.at(j) - ty.at(i) ; 
	w = tz.at(j) - tz.at(i) ; 
	uvw = sqrt(u*u + v*v + w*w) ; 
	u /= uvw; 
	v /= uvw; 
	w /= uvw;


	double phi = 2*(generator.rand() - 0.5) *thetaMax * PI ; // thetaMax*[-pi,pi] 

	double c,s ; 
	c = cos(phi) ; 
	s = sin(phi);

	vector<double> tnew; 
	//cout << i << " :" << j << endl ; 
	for(int k = i+1 ; k < j ; k++) {
		//cout << "chosen: " << k << endl; 
		//cout <<"old: " <<  norm2 << endl ;
		tnew = fast_rotate(tx[k],ty[k], tz[k], u, v, w , c,s) ; 
		tx[k] = tnew[0] ;
		ty[k] = tnew[1] ; 
		tz[k] = tnew[2] ; 

		rp[0] += (tx[k] - tox[k]) ; 
		rp[1] += (ty[k] - toy[k]) ; 

		if (k == N-1 ) { 
			tp[0] = tx[k] ; 
			tp[1] = ty[k] ;
		}
	}
}

void do_pivot() { 
	int i = generator.randInt(N-2) ; 
	int j = N; 
	double u, v, w; 
	u = generator.randNorm(0.0,1.0); 
	v = generator.randNorm(0.0,1.0); 
	w = generator.randNorm(0.0,1.0); 
	double phi = 2*(generator.rand() - 0.5) * thetaMax * PI ; // [-pi,pi] 
	double c,s ; 

	double uvw = sqrt(u*u + v*v + w*w) ; 
	u /= uvw; 
	v /= uvw; 
	w /= uvw;


	c = cos(phi) ; 
	s = sin(phi);

	vector<double> tnew; 
	//cout << i << " :" << j << endl ; 
	for(int k = i+1 ; k < j ; k++) {
		tnew = fast_rotate(tx[k],ty[k], tz[k], u, v, w , c,s) ; 
		tx[k] = tnew[0] ;
		ty[k] = tnew[1] ; 
		tz[k] = tnew[2] ;

		rp[0] += (tx[k] - tox[k]) ; 
		rp[1] += (ty[k] - toy[k]) ; 

		if (k == N-1 ) { 
			tp[0] = tx[k] ; 
			tp[1] = ty[k] ;
		}
	}
}

void savestate() {

	for(int index = 0 ; index < N; index++) {
		tox[index] = tx[index]; 
		toy[index] = ty[index]; 
		toz[index] = tz[index]; 
	}
}

int sgn(double val) { 
	return (val > 0 ) - (val < 0) ; 
}

int move() { 
	
	savestate() ; 
	int ip  = 1 + generator.randExc() * (N-2) ; 
	perturb(ip) ; 

	int jp =  1 + generator.randExc() * (N-2) ; 
	urn = generator.rand() - 0.5 ;  // generate sign of proposed tz[j]
	tx.at(jp) = tx.at(jp) - (tx.at(ip) - tox.at(ip)); 
	ty.at(jp) = ty.at(jp) - (ty.at(ip) - toy.at(ip)); 
	r2 = tx.at(jp) * tx.at(jp) +  ty.at(jp)*ty.at(jp) ;
	if ( r2 > 1 )  return 1  ; // rejection flag 

	tz.at(jp) = sgn(urn) * sqrt( (long double ) (1. - r2) ) ; 
	rp[0] += (tx.at(jp) - tox.at(jp)) ; 
	rp[1] += (ty.at(jp) - toy.at(jp)) ; 

	//if (ip == N-1 || jp == N-1) cout << "something went wrong" << endl; 
	assert(ip != N-1 && jp != N-1 ) ; 
	return 0 ; 

}


// monte carlo move acceptance/rejection criteria 
int mc_step() {

	//double z_prev = getZ() ; 
	double z_curr = getZ() ; 
	if (z_curr > zmax) zmax = z_curr ; 
	double E0 = getE() ; 
	//random_index = 1 + generator.randInt(N-1); // random int (1,N-1)  

	//perturb(random_index); 
	//cout << "before move: " << endl ; 
	int notOK = move() ;
	if(notOK) {
		revert() ; 
		return 0 ; 
	}

	double dE = getE() - E0 ; 

	if(false /*conditional_check() */) { 
		cout << "new rp: " << getRP() ; 
		cout << "target : " << RP  ; 
		cout << " new tp: " << getTP() ; 
		cout << endl ;
		revert() ; 
		return 0 ; 
	}

	// Metropolis part 
	if(dE <= 0.) return 1; 
	else {
		urn = generator.rand() ; 
		if ( urn > exp(-dE) ) {
			revert(); 
			return 0;
		}
	}

	return 1; 

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
		zMoments[i] += pow( N*z_curr, i+1) ; 
	}
}

void write_cosines() {
	for(int j =0; j< N-1; j++) {
		double tcos = innerProduct(tx.at(j),ty.at(j),tz.at(j),tx[j+1],ty[j+1],tz.at(j+1));  
		int bin_index = floor((0.5 + 0.5*tcos )* nbins) ; 
		if (bin_index > 0) {
			cosHist[bin_index] += 1; 
		}
	}
}

void writeZMoments() { 
	int wi = 0 ; 
	fout.open("moments.dat"); 
	for(int i = 0 ; i < nMoments; i++) { 
		double TC = totalCounts ; 
		zMoments[i] /=  TC ; 
		if (TC > 100 ) {
			fout << i+1 ; 
			fout << setw(15) << zMoments[i] ; 
			fout << endl ; 
		}
	}
	fout.close(); 
	fout.clear();
}


void writeZFile(int step) {
	int wi = 0 ; 
	char zString[40];
	sprintf(zString,"%d\t%f\t%f\n",step, getZ() , getE()); 
	fputs(zString,windowFiles[wi]);

}

void writeCoordinates( ) {

	xi = 0.0 ; 
	yi = 0.0 ; 
	zi = 0.0 ; 
	coorFile << N << endl ; 
	coorFile << endl ; 
	for(int i = 0 ; i < N ; i++) {

		xi += tx.at(i) ; 
		yi += ty.at(i) ; 
		zi += tz.at(i) ; 
		coorFile << 1  ;
		coorFile << " " ; 
		coorFile << xi ; 
		coorFile << " " ; 
		coorFile << yi ; 
		coorFile << " " ; 
		coorFile << zi ; 
		coorFile << endl ; 
	}
}


void WriteEventData(int step) {

	writeZ(step); 
	//writeZFile(step); 
	//getCorrelator() ; 
	//writeCoordinates() ; 
}

void writeHistograms(int wi = 0 ) {

	writeZMoments() ;
	double P, binContent;
	//fout.open("cosines.dat"); 

	double TC = totalCounts ; 

	double maxP = 0.0  ;
	for(int i = 0; i < nbins; i++) {
		double binLowEdge = i/(double) nbins ; 
		double binCenter = binLowEdge + 1/(2.*nbins) ; 
		binContent = zHists[wi][i]; 
		P = binContent / TC ; 
		if (P > maxP ) {
			maxP = P ; 
			z_mp = binCenter ; 
		}

		char histVals[25]; 
		if(binContent >100 ) {
			sprintf(histVals,"%f\t%f\n",binCenter, - log(P) ); 
			fputs(histVals,histFiles[0]) ;

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

void writeLogFile() { 

	fout.open("log.dat"); 
	fout << "RP " << setw(15) << RP << endl ; 
	fout << "TP " << setw(15) << TP << endl ; 
	fout << "RPTP " << setw(15) << RPTP << endl ; 
	fout << "zmax " << setw(15) << zmax << endl ; 
	fout << "Zmp " << setw(15) << z_mp << endl; 
	fout.close() ; 
	fout.clear() ; 
}


int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	//coorFile.close() ; 
	cout << "all files closed " << endl ; 
	return 0 ; 
}

