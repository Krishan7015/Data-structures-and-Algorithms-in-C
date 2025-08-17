#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#define Max 5
int queue[Max];
int front = 0;
int rear = -1;
int items = 0;
bool isFull() {
    return items == Max;
}
bool isEmpty() {
    return items == 0;
}
int peek(){
    return queue[front];
}
void insert(int data){
    if(!isFull()) {
        rear = (rear + 1)% Max;
        queue[rear] = data;
        items++;
    }
}
int removedata(){
    int data = queue[front++];
    if(front == Max) {
        front = 0;
    }
    items--;
    return data;
}
int main() {
    insert(12);
    insert(32);
    insert(14);
    insert(42);
    insert(90);
    printf("\nQueue after insertion:\n");
    for( int i = 0; i<items; i++) {
        printf("queue[%d] = %d\n",i,queue[i]);
    }
    int num = removedata();
    printf("\nElement removed: %d\n",num);
    if(!isEmpty()){
    printf("\nElement at front: %d\n",peek());
    }
        printf("\nUpdated Queue:\n");
    while(!isEmpty()) {
        int n = removedata();
        printf("%d\n",n);
    }
    if(!isEmpty()){
    printf("\nElement at front: %d\n",peek());
    }
    return 0;
    }
    