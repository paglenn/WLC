#include<vector> 
#include<cassert>
#include<iostream>
#include<cmath>
using std::vector; 
using std::ostream; 
using std::cout; 

// first allow printing of vectors to std::ostream 
ostream& operator<<(ostream& os , const vector<double>& v) {
	for(int i = 0; i < v.size() ; i++ ) os << v[i] << "  " ; 
	os << std::endl; 
	return os; 
}

// check compatibility of two vectors for addition/subtraction 
void check_compatibility(int s1, int s2) {
	if (s1 != s2) { 
		std::cout << "Vector sizes do not match!" << std::endl; 
		exit(1);  
	} 
}

// element-wise addition of vectors . 
vector<double> operator+(vector<double> v1, vector<double> v2)  { 
	int length = v1.size(); 

	check_compatibility(length, v2.size()) ; 
	vector<double> s(v1.size(),0.) ; 
	for(int i = 0; i < v1.size() ; i++ ) s[i] = v1[i] + v2[i] ; 
	return s; 
}

// element-wise subtraction of vectors. 
vector<double> operator-(vector<double> v1, vector<double> v2)  { 
	int length = v1.size(); 

	check_compatibility(length, v2.size()) ; 

	vector<double> s(v1.size(),0.) ; 
	for(int i = 0; i < v1.size() ; i++ ) s[i] = v1[i] - v2[i] ; 
	return s; 
}

// multiplication of a vector by a scalar 
vector<double> operator*(vector<double> v1, double scalar) { 
	vector<double> P(v1.size() , 0.) ; 
	for( int i = 0; i < v1.size() ; i++ ) P[i] = v1[i] * scalar; 	
	return P; 
}

// division of a vector by a scalar 
vector<double> operator/(vector<double> v1, double divisor) { 
	vector<double> Q(v1.size() , 0.) ; 
	if (divisor != 0 ) {
		for( int i = 0; i < v1.size() ; i++ ) Q[i] = v1[i] / divisor; 	
	} else {
		std::cout << "Warning : division by zero " << std::endl; 
		exit(1); 
	}
	return Q; 
}


// scalar product of two vectors 
double dot ( vector<double> v1, vector<double> v2 ) {

	int length = v1.size(); 
	check_compatibility(length, v2.size()) ; 

	double sum = 0.0; 
	for (int i = 0;  i < v1.size() ; i++ ) sum += v1[i]*v2[i]; 

	return sum; 
}

// length of a vector (the Frobenius norm)  
double norm(vector<double> v1 ) {  
	return sqrt ( dot(v1,v1) ) ; 
}

// rotation in the xy-plane 
vector<double> rotate_2D(vector<double> v, double phi ) {
	// 
	vector<double>  vprime(3, 0.) ; 
	double vx, vy ; 
	vx = v[0] ; vy = v[1] ; 
	vprime[0] = vx * cos(phi) + vy * sin(phi) ; 
	vprime[1] = -vx * sin(phi) + vy * cos(phi) ; 
	vprime[2] = v[2] ; 
	return vprime ; 

}

// rotation about the x-axis , in the phi direction 
vector<double> rotate_phi(vector<double> v, double phi) { 
	vector<double> v_new(3,0.) ; 
	double vx, vy, vz ; 
	vx = v[0]; vy = v[1]; vz = v[2] ;  
	v_new[0] = vx; 
	v_new[1] = vy*cos(phi) - vz*sin(phi);
	v_new[2] = vy* sin(phi) + vz * cos(phi) ; 
	return v_new; 
}

// rotation about an arb. line passing through the origin 
vector<double> rotate_phi_true(vector<double> vec, vector<double> axis, double phi) 
{ 
	double uvw = dot(axis, axis) ; 
	double suvw = sqrt(uvw) ; 
	double u = axis[0]; 
	double v = axis[1]; 
	double w = axis[2]; 
	double x = vec[0], y = vec[1] , z = vec[2] ; 
	double c = cos(phi), s = sin(phi) ; 

	double ux = u*x, uy = u*y, uz = u*z ; 
	double vx = v*x, vy = v*y, vz = v*z ; 
	double wx = w*x, wy = w*y, wz = w*z ; 

	vector<double> ans(3,0.0); 
	ans[0] = u*(ux + vy + wz)*(1-c) + uvw*x*c + suvw * (-wy + vz ) * s  ; 
	ans[1] = v*(ux + vy + wz)*(1-c) + uvw*y*c + suvw * (wx - uz) * s 	; 
	ans[2] = w*(ux + vy + wz)*(1-c) + uvw*z*c + suvw * (-vx + uy) * s	;

	return ans; 

}
	
	

	
