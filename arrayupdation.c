#include<stdio.h>
int main() {
int position,element,i=0;
int a[]={2,4,6,8,10}; //array initialize 
printf("Array before update:\n");
for(i=0;i<5;i++){ //print original array
printf("a[%d]=%d\n",i,a[i]);
}
printf("\nInsert position where you want to update element: \n");
scanf("%d",&position);
if(position>5||position<1) {
    printf("Error! Invalid Position."); //stops program if invalid input
    return 1;
}
printf("Insert the new element: \n");
scanf("%d",&element);
a[position-1]=element; //position-1 because indexing starts from 0
printf("Array after update:\n");
for(i=0;i<5;i++){
printf("a[%d]=%d\n",i,a[i]);
   }
   return 0;
}