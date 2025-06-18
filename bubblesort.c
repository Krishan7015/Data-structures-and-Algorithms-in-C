#include<stdio.h>
void bubble_sort(int a[],int size){
    for(int i=0;i<size;i++){
        int swaps=0;
        for(int j=0;j<size-i-1;j++){
            if(a[j]>a[j+1]){  //if lower indexed number is greater than higher indexed number                         
                int temp;    
                temp=a[j];     //swapping
                a[j]=a[j+1];
                a[j+1]=temp;
                swaps=1;
            }
        }
        if(swaps==0) //if no swaps array is already sorted
        break;
    }
}
int main(){
    int n=5;
    int array[n];  //n = size of array
    printf("Enter elements of array:\n");
    for(int i=0;i<n;i++){
        scanf("%d",&array[i]);
    }
    printf("Array before sorting:\n");
    for(int i=0;i<5;i++){
        printf("array[%d]=%d\n",i,array[i]);
    }
bubble_sort(array,n);  //function calling
printf("Array after sorting:\n");
for(int i=0;i<5;i++){
    printf("array[%d]=%d\n",i,array[i]);
}
return 0;
}
