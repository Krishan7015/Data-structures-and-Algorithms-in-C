#include<stdio.h>
void binary_search(int a[],int small, int big, int item){//small= lowest index of current search range
    int median;                                          //big= highest index of current search range
    median=(small+big)/2;                                 //item= element to be searched
    if(small<=big) {
        if(a[median]==item){ 
            printf("Element found at %d\n", median+1); //median+1 because array indexing is from 0
        }
            else if(item < a[median]){
            binary_search(a, small, median-1,item);
            }
        else if( item > a[median]){
        binary_search(a,median+1,big,item);
        }
    }
        else {
        printf("Unable to search!");
    }
}
int main(){
    int n=5,small=0,big=n-1,item;
    int a[]={2,3,7,9,12};
    printf("Enter element to search:");
    scanf("%d",&item);
    binary_search(a,small,big,item);
    return 0;
    }       
    
