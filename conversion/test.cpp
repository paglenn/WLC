#include<random>
#include<iostream>
using namespace std; 

int main() {

	srand(1); 
	for(int i = 0; i < 100; i++ ) {
		cout<<2*(rand()/float(RAND_MAX))-1<<endl;
	}

}
