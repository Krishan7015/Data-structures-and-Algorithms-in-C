#include<stdio.h>
int maxsize=8;
int top=-1; 
int stack[8];
int isfull() {
    if(top==maxsize - 1) {
        return 1;
    } else {
        return 0;
    }
}
int push(int elements) {
     if(!isfull()){
        top=top+1; //increase pointer by 1
        stack[top]=elements;
        return 1;
    } else {
        printf("Stack is full, can not push elements\n");
        return 0;
    }    
}
int main() {
    push(3);
    push(5);
    push(8);
    push(13);
    push(21);
    printf("Elements in stack are:\n");
    for(int i=0;i<=top;i++) {
        printf("stack[%d]=%d\n",i,stack[i]);
    }
    return 0;
}


        
        
    
 
