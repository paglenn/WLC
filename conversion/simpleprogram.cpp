#include<iostream> 
using namespace std ; 


double * my; 

void assign() {
	int n = 2; 
	my = new double[2]; 
	my[0] = 1; 
}

int main() {
	assign(); 
	cout<<my[0] << endl; 
}
