//run hello world in parallel

#include <omp.h>
#include <iostream>

int main(int argc, char *argv[])
{
    //get number of threads
    int nthreads = omp_get_max_threads();
    //print number of threads
    std::cout << "Number of threads: " << nthreads << std::endl;
    //run hello world in parallel
    #pragma omp parallel
    {
        //generate hello world string in variable so that threads don't clash while writing to std
        std::string hello = "Hello, World! from thread " + std::to_string(omp_get_thread_num()) + "\n";
        //print hello world string
        std::cout << hello;
        
    }

return 0;
}