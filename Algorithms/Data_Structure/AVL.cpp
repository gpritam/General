#include <bits/stdc++.h>

struct node
{
	int data;
	
	node *leftChild, *rightChild;
	
	int height;
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
int max (int a, int b)
{
	return (a > b ? a : b);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int height (node *a)
{
	return (a == nullptr ? 0 : a->height);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
struct Node *rightRotate(struct Node *y){
   struct Node *x = y->leftChild;
   struct Node *T2 = x->rightChild;
   x->rightChild = y;
   y->leftChild = T2;
   y->height = max(height(y->leftChild), height(y->rightChild)) + 1;
   x->height = max(height(x->leftChild), height(x->rightChild)) + 1;
   return x;
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
struct Node *leftRotate(struct Node *x){
   struct Node *y = x->rightChild;
   struct Node *T2 = y->leftChild;
   y->leftChild = x;
   x->rightChild = T2;
   x->height = max(height(x->leftChild), height(x->rightChild)) + 1;
   y->height = max(height(y->leftChild), height(y->rightChild)) + 1;
   return y;
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int getBalance(struct Node *N){
   if (N == NULL)
      return 0;
   return height(N->leftChild) - height(N->rightChild);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
struct Node *insertNode(struct Node *node, int data){
   if (node == NULL)
      return (newNode(data));
   if (data < node->data)
      node->leftChild = insertNode(node->leftChild, data);
   else if (data > node->data)
      node->rightChild = insertNode(node->rightChild, data);
   else
      return node;
   node->height = 1 + max(height(node->leftChild),
                     height(node->rightChild));
   int balance = getBalance(node);
   if (balance > 1 && data < node->leftChild->data)
      return rightRotate(node);
   if (balance < -1 && data > node->rightChild->data)
      return leftRotate(node);
   if (balance > 1 && data > node->leftChild->data) {
      node->leftChild = leftRotate(node->leftChild);
      return rightRotate(node);
   }
   if (balance < -1 && data < node->rightChild->data) {
      node->rightChild = rightRotate(node->rightChild);
      return leftRotate(node);
   }
   return node;
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
struct Node *deleteNode(struct Node *root, int data)
{
   if (root == NULL)
      return root;
   if (data < root->data)
      root->leftChild = deleteNode(root->leftChild, data);
   else if (data > root->data)
      root->rightChild = deleteNode(root->rightChild, data);
   else {
      if ((root->leftChild == NULL) || (root->rightChild == NULL)) {
         struct Node *temp = root->leftChild ? root->leftChild : root->rightChild;
         if (temp == NULL) {
            temp = root;
            root = NULL;
         } else
            *root = *temp;
         free(temp);
      } else {
         struct Node *temp = minValueNode(root->rightChild);
         root->data = temp->data;
         root->rightChild = deleteNode(root->rightChild, temp->data);
      }
   }
   if (root == NULL)
      return root;
   root->height = 1 + max(height(root->leftChild),
                     height(root->rightChild));
   int balance = getBalance(root);
   if (balance > 1 && getBalance(root->leftChild) >= 0)
      return rightRotate(root);
   if (balance > 1 && getBalance(root->leftChild) < 0) {
      root->leftChild = leftRotate(root->leftChild);
      return rightRotate(root);
   }
   if (balance < -1 && getBalance(root->rightChild) <= 0)
      return leftRotate(root);
   if (balance < -1 && getBalance(root->rightChild) > 0) {
      root->rightChild = rightRotate(root->rightChild);
      return leftRotate(root);
   }
   return root;
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
	insertNode(22);
	insertNode(14);
	insertNode(72);
	insertNode(44);
	insertNode(25);
	insertNode(63);
	insertNode(98);
	
	inorder(root);
	
	std::cout << std::endl << std::endl << "Delete Node: 25" << std::endl;
	deleteNode(root, 25);
	inorder(root);
	
	return 0;
}
