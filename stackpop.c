#include<stdio.h>
int maxsize=8;
int top=-1; 
int stack[8];
int isempty(){ //check if stack is empty
    if (top==-1)
    return 1;
    else
    return 0;   
}
int isfull(){ //check if stack is full
    if (top==maxsize-1)
    return 1;
    else
    return 0;
}
int pop(){ //pop data from stack
    int data;
    if(!isempty()){ 
        data= stack[top];
        top=top-1;
        return data;
    }
    else{
        printf("Could not retreive data, Stack is empty:\n");
        return -1;
    }
}
int push(int data){ //push data into stack
    if(!isfull()){ 
        top=top+1;
        stack[top]=data;
        return 1;
    }
    else{
        printf("Could not insert data, stack is full:\n");
        return 0;
    }
}
int main(){
    push(10);
    push(20);
    push(30);
    push(40);
    push(50);
    printf("Stack Elements:\n"); //print stack elements
    for(int i = 0; i<=top; i++){
        printf("stack[%d]=%d\n",i,stack[i]);
    }
printf("Popped Elements:\n"); //print popped elements
while(!isempty()){
    int data = pop();
    printf("%d ",data);
}
return 0;
}


