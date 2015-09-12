#include <stdio.h>
#include <iostream>
#include "MersenneTwister.h"
#include <gsl/gsl_rng.h>

int main(int argc, char * argv[])
{
    if (argc != 2)
    {
        printf("Usage: ./a.out (integer, # to display)\n");
        exit(1);
    }
    
    // WARNING: Code will fail if input is not an integer!
    //   Type checking is an advanced feature and will not be
    //   done here.
    int input = atoi(argv[1]);
    printf("Hello World, %d\n", atoi(argv[1]));
    std::cout << "Hello World, " << argv[1] << std::endl;


    // Populate an array of size "len" with the index of that
    //   element of the array.
    int len = 10;
    int * data = new int [ len ];
    for (int i = 0 ; i < len ; i++ )
    {
        data[i] = i;
    }

    for (int i = 0 ; i < len ; i++ )
    {
        std::cout << data[i] << " ";
    }
    std::cout << std::endl;

    // Initialize the Mersenne Twister object with seed 90210
    MTRand rng(90210);
    for (int i = 0 ; i < len ; i++ )
    {
        data[i] = rng.randInt(input);
    }

    const gsl_rng_type * T;

    GSL_RNG_SEED=90210;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    delete data;
}
