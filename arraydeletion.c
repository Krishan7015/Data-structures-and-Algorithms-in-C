#include<stdio.h>
void main() {
    int a[]={2,5,6,8,9};
    int n=5,b;
    printf("The elements of array are:\n");
    for(int i=0;i<n;i++) {
        printf("a[%d]=%d\n",i,a[i]);
    }
    printf("Enter number of Elements to delete:\n");
    scanf("%d",&b);
    if(b>n){
        printf("ERROR! You cannot delete that much elements.");
        return;
    }
    for (int i=0;i<n-b;i++){
        a[i]=a[i+b];
        }
        n-=b;
     printf("The array elements after deletion :\n");
     for(int i=0;i<n;i++) {
        printf("a[%d]=%d\n",i,a[i]);
     }
}