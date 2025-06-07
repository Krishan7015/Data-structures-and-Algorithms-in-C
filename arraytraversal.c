#include<stdio.h>
int main () {
    int a[]={3,4,5,6,7}; //initialize array
    int items=sizeof(a)/sizeof(a[0]); //calculate array size
    printf("Array elements are: \n");
    for(int i=0;i<items;i++) {
        printf("a[%d]=%d\n",i,a[i]); //print each element
    }
return 0;
}