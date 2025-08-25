#include<stdio.h>
void insertion_sort(int arr[], int size) {
    int temp, j;
    for(int i=1; i<size; i++) {
    temp = arr[i];
    j = i;
    while(j>0 && arr[j-1]>temp){
        arr[j] = arr[j-1];
        j--;
    }
    arr[j] = temp;   
  }
}
int main(){
    int n;
    n=7;
    int array[7] = {44,56,67,24,52,58,73};
    printf("Array before sorting:");
    for(int i = 0; i<n; i++){
        printf("array[%d] = %d\n",i,array[i]);
    }
    insertion_sort(array,n);
    printf("\nArray after insertion:");
    for(int i = 0; i<n; i++){
        printf("array[%d] = %d\n",i,array[i]);
    }
}
