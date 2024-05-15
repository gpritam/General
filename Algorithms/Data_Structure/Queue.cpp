#include <bits/stdc++.h>

struct node
{
    int data;
	
    node *next;
};

const int MaxSize = 10;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Print (node *front)
{
    while (front != nullptr)
    {
        std::cout << front->data << std::endl;
		
        front = front->next; 
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int Size (node *front)
{
	int count = 0;
	
	while (front != nullptr)
	{
		count++;
		
		front = front->next;
	}
	
	return count;
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
bool isEmpty (node* front)
{
	return (Size(front) == 0 ? true : false);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
bool isFull (node* front)
{
    return (Size(front) == MaxSize ? true : false);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Enqueue ( node **front, 
			   node **rear, 
			   int data )
{
    if (isFull(*front))
	{
        std::cout << "Queue is full!" << std::endl;
		
		exit(1);
    }
    else
	{
		node *newNode = new node;
		
        newNode->data = data;
        newNode->next = nullptr;
		
        if (*front == nullptr)
            *front = *rear = newNode;
        else
		{
            (*rear)->next = newNode;
			
            *rear = newNode;
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Dequeue (node **front)
{
    if (isEmpty(*front))
	{
        std::cout << "Queue is empty!" << std::endl;
		
		exit(1);
    }
    else
	{
		node *temp = *front;
		
        *front = (*front)->next;
        
        delete temp;
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int main()
{
	node *front = nullptr, *rear = nullptr;
	
    Enqueue(&front, &rear, 34);
    Enqueue(&front, &rear, 4);
    Enqueue(&front, &rear, 7);
    Enqueue(&front, &rear, 17);
	
	Print(front);
	
    Dequeue(&front);
    
	std::cout << std::endl;
	Print(front);
	
    return 0;
}
