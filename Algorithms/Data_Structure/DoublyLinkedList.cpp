#include <bits/stdc++.h>

struct Node
{
	double data;
	
	Node *NextNode, *PreviousNode;
};

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate (Node **A)
{
	Node *Temporary;
	
	while (*A != nullptr)
	{
		Temporary = *A;
		
		*A = (*A)->NextNode;
		
		delete Temporary;
	}
}

//________________________________________________________________________________________________
// Insert a node at k-th position
//________________________________________________________________________________________________
void Insert ( Node **A, 
			  const double data, 
			  const int k )
{
	Node *Temporary, *Previous = nullptr, *Current = nullptr;
	
	if (*A == nullptr)
	{
		Temporary = (Node *)malloc(sizeof(Node));
		Temporary->data = (k == 0 ? data : 0.0);
		Temporary->PreviousNode = nullptr;
		Temporary->NextNode = nullptr;
		
		*A = Temporary;
	}
	
	if (k > 0)
	{
		Previous = *A;
		Current = (*A)->NextNode;
		
		for (int i = 1; i < k; ++i)
		{
			if (Current == nullptr)
			{
				Temporary = (Node *)malloc(sizeof(Node));
				Temporary->data = 0.0;
				Temporary->PreviousNode = Previous;
				Temporary->NextNode = nullptr;
				
				Current = Temporary;
				Previous->NextNode = Current;
			}
			
			Previous = Current;
			Current = Current->NextNode;
		}
		
		Temporary = (Node *)malloc(sizeof(Node));
		Temporary->data = data;
		Temporary->PreviousNode = Previous;
		Temporary->NextNode = (Current == nullptr ? nullptr : Current);
		
		Previous->NextNode = Temporary;
		
		if (Current != nullptr)
		{
			Current = Current->NextNode;
			
			if (Current != nullptr)
				Current->PreviousNode = Temporary;
		}
		
		Current = Temporary;
	}
}

//________________________________________________________________________________________________
// Delete the k-th node
//________________________________________________________________________________________________
void Delete ( Node **A, 
			  const int k )
{
	if (*A != nullptr)
	{
		Node *Temporary = *A, *Previous, *Next;
		
		for (int i{}; i < k; ++i)
		{
			if (Temporary == nullptr)
			{
				std::cout << "Element is not present!" << std::endl;
				
				return;
			}
			
			Temporary = Temporary->NextNode;
		}
		
		Previous = Temporary->PreviousNode;
		Next = Temporary->NextNode;
		
		if (Previous != nullptr)
			Previous->NextNode = Next;
		
		if (Next != nullptr)
			Next->PreviousNode = Previous;
		
		if ((Previous == nullptr) && (Next == nullptr))
		{
			delete Temporary;
			
			*A = nullptr;
		}
		else if ((Previous == nullptr) && (Next != nullptr))
		{
			*A = (*A)->NextNode;
			
			delete Temporary;
		}
		else
			delete Temporary;
	}
}

//________________________________________________________________________________________________
// Modify data of the k-th node
//________________________________________________________________________________________________
void ModifyData ( Node **A, 
				  const double data, 
				  const int k )
{
	if (*A != nullptr)
	{
		Node *Temporary = *A;
		
		for (int i{}; i < k; ++i)
		{
			if (Temporary == nullptr)
			{
				std::cout << "Element is not present!" << std::endl;
				
				return;
			}
			
			Temporary = Temporary->NextNode;
		}
		
		Temporary->data = data;
	}
}

//________________________________________________________________________________________________
// Swap the data of i and j th element of the list A
//________________________________________________________________________________________________
void Swap ( Node **A, 
			const int i, 
			const int j )
{
	if (A != nullptr)
	{
		Node *First = *A;
		
		for (int k{}; k < i; ++k)
		{
			if (First == nullptr)
			{
				std::cout << "Element is not present!" << std::endl;
				
				return;
			}
			
			First = First->NextNode;
		}
		
		Node *Second = *A;
		
		for (int k{}; k < j; ++k)
		{
			if (Second == nullptr)
			{
				std::cout << "Element is not present!" << std::endl;
				
				return;
			}
			
			Second = Second->NextNode;
		}
		
		double data = First->data;
		First->data = Second->data;
		Second->data = data;
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void ReverseList (Node **A)
{
	if (*A != nullptr)
	{
		Node *Temp = *A;
		
		while (Temp != nullptr)
		{
			*A = Temp;
			Temp = (*A)->NextNode
			(*A)->NextNode = (*A)->PreviousNode;
			(*A)->PreviousNode = Temp;
		}
	}
}

//________________________________________________________________________________________________
// Print the data of a particular node
//________________________________________________________________________________________________
void PrintElement ( Node *A, 
					const int i )
{
	if (A != nullptr)
	{
		for (int k{}; k < i; ++k)
		{
			if (A == nullptr)
			{
				std::cout << "Element is not present!" << std::endl;
				
				return;
			}
			
			A = A->NextNode;
		}
		
		std::cout << A->data << std::endl;
	}
}

//________________________________________________________________________________________________
// Print all the data of the list (from head to tail)
//________________________________________________________________________________________________
void PrintList ( Node *A )
{
	while (A != nullptr)
	{
		std::cout << A->data << std::endl;
		
		A = A->NextNode;
	}
}

//________________________________________________________________________________________________
// Print all the data of the list (from tail to head)
//________________________________________________________________________________________________
void PrintReverseList ( Node *A )
{
	if (A != nullptr)
	{
		while (A->NextNode != nullptr)
			A = A->NextNode;
		
		while (A != nullptr)
		{
			std::cout << A->data << std::endl;
			
			A = A->PreviousNode;
		}
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int main()
{
	Node *A = nullptr;
	int choice;
	
	while (true)
	{
		std::cout << "Enter your choice : " << std::endl 
			  << "1 : Insert node" << std::endl
			  << "2 : Delete node" << std::endl
			  << "3 : Modify data of a node" << std::endl
			  << "4 : Swap data of two nodes" << std::endl
			  << "5 : Reverse data of this complete list" << std::endl
			  << "6 : Print data of a particular node" << std::endl
			  << "7 : Print the complete list" << std::endl
			  << "8 : Print the list from the tail to head node" << std::endl
			  << "9 : Quit" << std::endl;
		
		std::cin >> choice;
		
		switch (choice)
		{
			case 1:
				Insert(&A);
				break;
			
			case 2:
				Delete(&A);
				break;
			
			case 3:
				ModifyData(&A);
				break;
			
			case 4:
				Swap(&A);
				break;
			
			case 5:
				ReverseList(&A);
				break;
			
			case 6:
				PrintElement(A);
				break;
			
			case 7:
				PrintList(A);
				break;
			
			case 8:
				PrintReverseList(A);
				break;
			
			case 9:
				Deallocate(&A);
				exit(1);
			
			default:
				std::cout << "Not a valid input, try again!" << std::endl;
		}
	}
	
	return 0;
}
