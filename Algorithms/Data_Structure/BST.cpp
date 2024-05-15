#include <bits/stdc++.h>

struct node
{
	int data;
	
	node *leftChild, *rightChild;
};

node *root = nullptr;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
node* search ( int data )
{
	node *current = root, *parent;
	
	while (true)
	{
		parent = current;
		
		if (parent == nullptr)
			return parent;
		else
		{
			if (parent->data == data)
				return parent;
			else
				current = (data > parent->data ? parent->rightChild : parent->leftChild);
		}
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void insert (int data)
{
	node *tempNode = new node;
	
	tempNode->data = data;
	
	tempNode->leftChild = tempNode->rightChild = nullptr;
	
	if (root == nullptr)
	{
		root = tempNode;
		
		return;
	}
	
	node *current = root, *parent;
	
	while (true)
	{
		parent = current;
		
		if (data > parent->data)
		{
			current = parent->rightChild;
			
			if (current == nullptr)
			{
				parent->rightChild = tempNode;
				
				return;
			}
		}
		else
		{
			current = parent->leftChild;
			
			if (current == nullptr)
			{
				parent->leftChild = tempNode;
				
				return;
			}
		}
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void deleteNode (int data)
{
    // Check if the key is actually present in the BST.
	node *current = root, *parent = nullptr;
	
    while (current != nullptr)
	{
		if (current->data == data)
			break;
		
        parent = current;
		
        current =  (data < current->data ? current->leftChild : current->rightChild);
    }
	
    if (current == nullptr)
	{
        std::cout << "Value is not present in this BST" << std::endl;
		
        return;
    }
	
    if ( (current->leftChild == nullptr) || (current->rightChild == nullptr) )
	{
		// node to be deleted has atmost one child
		
        node* newCurrent = (current->leftChild == nullptr ? current->rightChild : current->leftChild);
		
        if (parent == nullptr)
            root = newCurrent;
		else
		{
		    if (current == parent->leftChild)
		        parent->leftChild = newCurrent;
		    else
		        parent->rightChild = newCurrent;
		}
		
        delete current;
		
		return;
    }
    else
	{
		// node to be deleted has two children 
		// [Find first the immediate inorder successor (left-most node of a BST whose root is current->rightChild) and successorParent]
		
		node *successorParent = current, *successor = current->rightChild;
		
		while (successor->leftChild != nullptr)
		{
			successorParent = successor;
			
			successor = successor->leftChild;
		}
		
		if (successorParent == current)
			successorParent->rightChild = successor->rightChild;
		else
			successorParent->leftChild = successor->rightChild;
		
		current->data = successor->data;
		
		delete successor;
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void inorder(node* Node)
{
	if (Node != nullptr)
	{
		inorder(Node->leftChild);
		
		std::cout << Node->data << " --> ";
		
		inorder(Node->rightChild);
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void preorder(node* Node)
{
	if (Node != nullptr)
	{
		std::cout << Node->data << " --> ";
		
		preorder(Node->leftChild);
		
		preorder(Node->rightChild);
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void postorder(node* Node)
{
	if (Node != nullptr)
	{
		postorder(Node->leftChild);
		
		postorder(Node->rightChild);
		
		std::cout << Node->data << " --> ";
	}
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int main()
{
    /* Let us create following BST
              50
           /     \
          30      70
         /  \    /  
       20   40  60    */
	
    insert(50);
    insert(30);
    insert(20);
    insert(40);
    insert(70);
    insert(60);
	
	inorder(root);
	std::cout << std::endl;
	preorder(root);
	std::cout << std::endl;
	postorder(root);
	std::cout << std::endl;
	
	node* k = search(30);
	
	if (k != nullptr)
		std::cout << std::endl << std::endl << "Element " << k->data << " found." << std::endl;
	else
		std::cout << std::endl << std::endl << "Element not found!" << std::endl;
	
	std::cout << std::endl << std::endl << "Delete a leaf node: 20" << std::endl;
    deleteNode(20);
    inorder(root);
	
	std::cout << std::endl << std::endl << "Delete Node with single child node: 70" << std::endl;
    deleteNode(70);
    inorder(root);
	
	std::cout << std::endl << std::endl << "Delete Node with both child nodes: 50" << std::endl;
    deleteNode(50);
    inorder(root);
	
	std::cout << std::endl << std::endl;
	
	return 0;
}
