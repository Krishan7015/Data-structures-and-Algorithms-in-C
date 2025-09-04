public class stackpush {
    final static int maxsize = 10;
    public static int stack[] = new int[maxsize];
    public static int top = -1;
    public static boolean isempty(){
        return top == -1;
    }
    public static boolean isfull(){
        return top == maxsize - 1;
    }
    public static int push(int data){
        if(!isfull()) {
            top++;
            stack[top] = data;
        }
        else {
            System.out.println("Stack is already Full.");
        }
        return data;
    }
    public static void main(String args []) {
        push(34);
        System.out.println("Stack Elements:");
        for(int i = 0; i <= top; i++) {
            System.out.println("Stack["+i+"] = "+stack[i]);
        }
    }
}
