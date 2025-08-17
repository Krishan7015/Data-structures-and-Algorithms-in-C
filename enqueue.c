#include<stdio.h>
#include<stdlib.h>     
#include<string.h>   // dequeue from front and enqueue at rear
#include<stdbool.h>  // inc front pointer to remove and inc rear pointer to add
#define Max 5
int queue[Max];
int front = 0;
int rear = -1;
int items = 0;
bool isFull() {
    if( items == Max ){
        return true;
    }
    else {
        return false;
    }
}
bool isEmpty() {
    if (items == 0){
        return true;
    }
    else {
        return false;
    }
}
int removedata() {
    int data = queue[front++];
    if( front == Max ){ //forms circular queue
        front = 0;
    }
    items--;
    return data;
}
void insert(int data) {
    if(!isFull()){
        rear = (rear + 1)% Max;
        queue[rear] = data;
        items++;
    }
}
int main() {
    insert(12);
    insert(32);
    insert(14);
    insert(42);
    insert(15);
printf("Queue: \n");
while (!isEmpty()) {
int n = removedata();
printf("%d \n",n);
}
}
