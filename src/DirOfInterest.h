#ifndef DIROFINTEREST_H_
#define DIROFINTEREST_H_


//CHECK OTHER CODE FOR u and other specified values
//Minimum value would be around 512, ideal is 1024 or 2048
int n = 5;
int k = 10; //Specify k number of points in the xy plane

double* getDirs(int k,int n);

double* getEquator(const int k); //Returned array is dynamically allocated, must free up
// all sub-lists as well

double* getRotationMat(double z[3],double u[3]); //Also dynamically allocated



#endif /* DIROFINTEREST_H_ */
