#include<iostream>
#include<cmath>
#include<math.h>
#include<vector>
#include<string>
#include<cstdlib>
//#include<random>
#include<ctime>
#include"VectorMethods.h"
#include<stdio.h>
#include<cstring>
#include "MersenneTwister.h"
using std::vector; 
using std::endl; 
using std::cout;
using std::cerr; 
using std::normal_distribution ; 

// system parameters 
double L,_L,ds, _ds; 
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
vector<int> totalCounts; 

const double PI = acos(-1.);
const double sqrt2 = sqrt(2) ; 
const int MAX = pow(2,32) -1 ; 
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
		else if (!strcmp(str,"nbins")) fscanf(inFile,"%d",&nbins); 
		//else std::cout<<"whoami" << str<<std::endl;
	}

	// MF correction / rescaling
	_L = 1./L ; 
	L = 1./(_L + 2) ; 
}

void createHistograms() {

	totalCounts = vector<int>(numWindows,0) ; 
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
	return zsum/(float) N; 
}

double getRP() {
	// projection of eed vector into the plane 
	double ans = sqrt( rp[0]*rp[0] + rp[1]*rp[1]) ; 
	return ans/ (float) N ; 
}

double getTP() { 
	//perpendicular component of final segment 
    //double ans = (elem == 0) ? tx[N-1] : ty[N-1] ;
	double ans = sqrt(tx[N-1]*tx[N-1] + ty[N-1]*ty[N-1]);
	return ans ; 
}

double getRPTP() {
	// innerProduct product of perpendicular component of final segment and 
	// projection of eed vector into the plane 
	/*
	double rx = 0.,ry = 0.; 
	for(int i = 0; i < N; i++) {
		rx += tx[i];
		ry += ty[i];
	}
	*/
	//cout << tp[0] << ":" << tx[N-1] << endl ; 
	//cout << tp[1] << ":" << ty[N-1] << endl ; 

	return (rp[0]*tp[0] + rp[1]*tp[1])/ N ; 
}


double gaussian (double mean, double var) { 
	
	double ans, u , v; 
	u = generator.randInt()%MAX / (double) MAX ;  
	v = generator.randInt()%MAX / (double) MAX ;
	ans = mean + sqrt( - 2* var * log(u) ) * cos(2*PI*v) ; 
	return ans ; 
}

void init() {

	readInParameters();
//	std::cout<<"Initialized" << std::endl;

	ds = L/ (double) N ; 
	_ds  = 1./ds; 
	numSteps = numSweeps * N; 
	//cout << numSteps << endl ;
	gaussVar *= ds; 
	numFrames = numWindows * numPasses * numSteps ; 
	binWidth = (1. - zstart)/nbins; 
	binOverlap = 1 * binWidth; 
	iseed = time(0) ; 
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
	tp[0] = tp[1] = 0.0 ; 
	rp[0] = rp[1] = 0.0 ; 

	// initialize to a random configuration 
	// don't change first ti from (0,0,1) 
	// limit to ti > 0 at first 
	z = -1;
	while ( z  < 0 ) { 	
	for (int index = 1 ; index < N ; index++ ) { 

		

		tx[index] = gaussian(0.0,1.0);
		ty[index] = gaussian(0.0,1.0);
		while( tz[index] <= 0. ) { 
			tz[index] = gaussian(0.0,1.0) ;
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
	cout << "Starting Z: " << z << endl; 
	cout << "RP " << RP << endl ; 
	cout << "TP " << TP << endl; 
	cout << "RPTP " << RPTP << endl ; 

	
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

	tx[index] += gaussian(0.0,gaussVar);
	 // needed for independent rv's 
	ty[index] += gaussian(0.0,gaussVar);
	
	tz[index] += gaussian(0.0,gaussVar);
	


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


double deltaE(int index, int w_index, double z_prev) {
	
	double E_old = 0, E_new = 0  ; 
	
	if( index+1 != N) {
		E_old -= innerProduct(tx[index+1],ty[index+1],tz[index+1],tox[index],toy[index],toz[index]) ;
		E_new -= innerProduct(tx[index+1],ty[index+1],tz[index+1],tx[index],ty[index],tz[index]) ;
	}
	
	if( index != 0) {
		E_old -= innerProduct(tx[index-1],ty[index-1],tz[index-1],tox[index],toy[index],toz[index]) ;
		E_new -= innerProduct(tx[index-1],ty[index-1],tz[index-1],tx[index],ty[index],tz[index]) ;
	}

	double z_new = getZ(); 
	//double dVz = ( z_new < winMin[w_index] || z_new > winMax[w_index] ) ? bigNum : 0.;
	double dVz  = V_bias(z_new,w_index) - V_bias(z_prev,w_index) ; 


	return _ds*( E_new - E_old) + dVz;

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
	if (fabs(a-b)  > 1e-6) return true; 
	else return false; 
}

void adjustZ(int w_index) {

	double wmean = zmin[w_index] ; 
	double z = 0.; 
	double target = wmean ; // (-1,1)
	double tol = ds;

	// this is bound to be horribly slow 
	// PROBLEM
	while(fabs(z - target) > tol || diff(RP,getRP()) || diff(TP,getTP()) || diff(RPTP, getRPTP() ) ) {

		random_index = 1+rand()%(N-1) ; 
		cout << random_index << endl ; 
		perturb(random_index); 
		double z_new = getZ(); 
		double dz = z_new - z; 

		if(z > target && dz > 0 ) dz = revert(random_index); 
		if(z < target && dz < 0 ) dz = revert(random_index); 

		if(dz!=0) z = z_new; 
	}

}

void rotate_T(double angle, int dim = 3) {
	//cout << "Angle " << angle << endl ; 

	vector<double> v(3,0.); 
	v[0] = tx[N-1] ; 
	v[1] = ty[N-1] ; 
	v[2] = tz[N-1] ; 
	double norm2 = tx[N-1]*tx[N-1]  + ty[N-1] *ty[N-1] ; 
	//cout << "norm before rotation: " << norm2 << endl ; 
//	cout <<"v: "<<  v[0] << " " << v[1] << endl;

	//if (dim == 3) v = rotate_phi(v,angle) ;
	v = rotate_2D(v,angle) ; 

	//update RP 
//	cout <<"Before: "<<  rp[0] << " " << rp[1] << endl; 
	rp[0] += (v[0] - tx[N-1]) ; 
	rp[1] += (v[1] - ty[N-1]) ; 
//	cout <<"v: "<<  v[0] << " " << v[1] << endl;
//	cout <<"After: "<<  rp[0] << " " << rp[1] << endl; 

	tx[N-1]  = v[0] ; 
	ty[N-1]  = v[1] ; 
	tz[N-1]  = v[2] ; 
	norm2 = tx[N-1]*tx[N-1]  + ty[N-1] *ty[N-1] ; 
	//cout << "norm after rotation: " << norm2 << endl ; 

	tp[0] = tx[N-1] ; 
	tp[1] = ty[N-1] ; 

}

void alignRPTP(double target) {
	if (rp[0] == 0 ) return ; 
	double K = getTP() ; 
	double theta_0 = atan(rp[1]/rp[0]) ; 
	cout << "init theta "<< theta_0 << endl ; 
	double theta_t = acos(target / (getRP()*getTP() )) ; 
	cout << " rp*tp " << getRP() * getTP() << endl; 
	cout << "target " << target << endl ; 

	cout << "target_theta: "<< theta_t << endl ; 
	tp[0] = K  * cos(theta_0+ theta_t) ; 
	tp[1] = K * sin(theta_0+theta_t)  ; 
	//double tpx = tp[0] ; double tpy = tp[1] ; 
	double norm = sqrt( (rp[0]*rp[0] + rp[1] * rp[1]) * ( tp[0]*tp[0] + tp[1]*tp[1]) ) ;
	cout << " norm before: " << norm << endl ; 
	if (norm < 1e-6)  return; 
	else { 
		double theta_curr = acos( (rp[0]*tp[0] + rp[1]*tp[1]) / norm )  ;
		double theta_target  = acos(target / norm ) ; 
		cout << "target_theta: "<< theta_target << endl ; 
		double dtheta = theta_target - theta_curr ; 
		if (dtheta > 0 ) {
			//rotate_T( dtheta, 2 ) ; 
		} else { 
			//rotate_T( -dtheta, 2 ) ;
		}

		// check to see if the right vector was rotated 
		//norm = sqrt( (rp[0]*rp[0] + rp[1] * rp[1]) * ( tp[0]*tp[0] + tp[1]*tp[1]) ) ;
		double theta_new = acos( (rp[0]*tp[0] + rp[1]*tp[1]) / norm )  ;
		cout << "new theta (1)" << theta_new  << endl ;  
		/*
		if (fabs(theta_new - theta_target) > fabs(dtheta)) {
			rotate_T(-2*dtheta, 2) ; 
		}
		*/
		//cout << "new theta (2)" << acos( (rp[0]*tp[0] + rp[1]*tp[1]) / norm)  << endl ;  
		cout << "rptp " <<  (rp[0]*tp[0] + rp[1]*tp[1]) << endl ; 
		//cout << tx[N-1] << " " << tp[0] << endl ; 
		//cout << ty[N-1] << " " << tp[1] << endl ; 
		//cout << " norm after: " << norm << endl ; 
	}
}

void adjustRP(double target) {
	
	double tol = 1./(2*N);
	double myRP = getRP() ; 

	while(fabs(myRP- target) > tol) {

		random_index = 1+rand()%(N-2) ; 
		perturb(random_index); 
		double rp_new = getRP(); 
		double dRP = rp_new - myRP; 

		if(myRP > target & dRP > 0 ) dRP = revert(random_index); 
		if(myRP < target & dRP < 0 ) dRP = revert(random_index); 

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
	cout << "t:" << tx[N-1] << " " << ty[N-1] << " " << tz[N-1] << endl ;
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
bool umbrella_mc_step(int w_index = numWindows-1) {

	double z_prev = getZ() ; 
	random_index = 1 + rand() % (N-1); // random int (1,N-1)  

	perturb(random_index); 

	double dE = deltaE(random_index, w_index, z_prev ); 

	// keeping RP, TP, RPTP constant here 
	// possible PROBLEM but maybe not. sampling still happens. 
	/*
	if (diff(RP,getRP()) || diff(TP,getTP()) || diff(RPTP,getRPTP()) ) {
			revert(random_index) ; 
			return false; 
	
	}
	*/
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
		E -= innerProduct(tx[i],ty[i],tz[i],tx[i+1],ty[i+1],tz[i+1]); 
	}
	return _ds*E ; 
}

double getBiasedE(int wi ) { 

	return getE() + V_bias(getZ(),wi) ; 
}

void writeZ(int step, int wi ) { 
	double bin_index = floor( getZ() * 1000) ; 

	zHists[wi][bin_index] += 1 ; 
	totalCounts[wi] += 1 ; 
}

void writeZFile(int step, int wi) {
	char zString[40];
	sprintf(zString,"%d\t%f\t%f\n",step, getZ() , getE()); 
	fputs(zString,windowFiles[wi]);

}

void WriteEventData(int step, int wi ) {

	writeZ(step,wi); 
	writeZFile(step,wi); 
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
		double TC = totalCounts[j] ; 

		for(int i = 0; i < nbins; i++) {
			double binLowEdge = i/(double) nbins ; 
			double binCenter = binLowEdge + 1/(2.*nbins) ; 
			double binContent = zHists[j][i]; 
			double P = binContent / TC ; 

			char histVals[25]; 
			if(binContent != 0 ) {
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
			cout << "instead has magnitude " << norm  << endl;
			exit(1); 
		}
	}
}

int cleanup() {
	fclose(zFile); 
	fclose(progressFile);
	return 0 ; 
}

