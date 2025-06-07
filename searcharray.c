#include<stdio.h>
int main() {
int Array[]={3,5,7,8,2,1}; //initalize array
int a,found=0,items=6;
printf("Elements in array are: \n");
for(int i=0;i<items;i++){
printf("Array[%d]=%d\n",i,Array[i]); //print each array element
}
printf("Enter element to search: ");
scanf("%d",&a);
for(int i=0;i<items;i++) {
  if(Array[i]==a) {
        printf("Found element %d at position %d",a,i+1);
        found=1;
        break;
    }
}
    if (!found) {
    printf("Item not found\n");
    }
    return 0;
}
    
       
       
    
