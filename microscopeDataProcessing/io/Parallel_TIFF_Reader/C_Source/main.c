#include "parallelReadTiff.h"
#include <time.h>

int main(){
    clock_t begin = clock();
    void* test = readTiffParallelWrapper("C:\\Users\\Matt\\Desktop\\testTiff\\test.tif");
    if(!test) return 1;
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f\n",time_spent);
    free(test);
    return 0;
}
