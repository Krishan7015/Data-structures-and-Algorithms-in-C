#include<stdio.h>
int Max=5;
int stack[5];
int top=-1;
int isfull(){
    if (top==Max-1)
        return 1;
    else
        return 0;
}
int peek(){
    return stack[top];
}
int push(int data){
    if (!isfull()){
       top++;
       stack[top]=data;
       return 1;
    }
    else {
        printf("Could not insert data, stack is full\n");
        return 0;
    }
}
int main(){
    push(20);
    push(88);
    push(55);
    push(99);
    push(24);
    push(13); //Cannot be pushed as Stack is full
    printf("Top element is: %d \n",peek());
printf("Stack elements after pushing:\n");
    for( int i=0; i<=top; i++){
        printf("stack[%d]=%d\n",i,stack[i]);
    }
}