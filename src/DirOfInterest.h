#ifndef DIROFINTEREST_H_
#define DIROFINTEREST_H_


//CHECK OTHER CODE FOR u and other specified values
//Minimum value would be around 512, ideal is 1024 or 2048
//int n = 5;
//int k = 10; //Specify k number of points in the xy plane

//double* getDirs(int k,int n); Unused?

void getEquator(const int k, double ** C); //Initialize array and pass in C as argument
// Remember to free up all space taken

void getRotationMat(double* z,double* u, double** out); //Also dynamically allocated



#endif /* DIROFINTEREST_H_ */
