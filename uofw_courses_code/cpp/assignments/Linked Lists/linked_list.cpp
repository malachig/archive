#include <iostream.h>

#include "linked_list.h"

#define errMsg "<Error> LinkedList::"

/****
  Default constructor. Creates an empty list.
*****/
LinkedList::LinkedList()
{
  dm_size = 0;
  dm_list = NULL;
}


/****
  Copy constructor. Duplicates the nodes of given list, new nodes point the *same* data.
*****/
LinkedList::LinkedList(const LinkedList &l)
{
  dm_size = 0;
  dm_list = NULL;

  for ( int i=0; i<l.size(); i++ )
    add( (void*) l[i] );
}


/****
  Clears the list.
*****/
LinkedList::~LinkedList()
{
  clear_list();
}


/****
  Assignment operator. Duplicates the nodes of given list, new nodes point the *same* data.
*****/
const LinkedList& LinkedList::operator=(const LinkedList &l)
{
  if ( this == &l )
    return *this;		// LHS is the same as RHS

  clear_list();

  for ( int i=0; i<l.size(); i++ )
    add( (void*) l[i] );

  return *this;
}


/****
	add nodes to end of list, points to given data
*****/
void LinkedList::add( void *data )
{
  if ( data == NULL ) {
    cerr << errMsg << "append_entry() - NULL data."<<endl;
    return;
  }

  Node* n = new Node;
  n->data = data;
  n->link = NULL;
	
  if ( dm_size == 0 )
  {
    dm_list = n;
    dm_size = 1;
    return;
  }

  Node*	c = dm_list;	// current node

  // loop until link is null, then c is last node
  while ( c->link != NULL )
    c = c->link;

  c->link = n;
  dm_size++;
}

/*****
  Returns data at given position. Does bounds checking. Returns NULL if out of bounds
*****/
void* LinkedList::operator[]( int idx )
{
  if ( idx<0 ||  idx>=dm_size ) {
    cerr << errMsg << "operator[] " << idx << " out of bounds." << endl;
    return NULL;
  }

  Node* n = dm_list;

  // loop until reach index position
  for ( int i=0; i < idx; i++ )
    n = n->link;
		
  return n->data;
}

/*****
  Returns data at given position. Does bounds checking. Returns NULL if out of bounds
*****/
const void* LinkedList::operator[]( int idx ) const
{
  if ( idx<0 ||  idx>=dm_size ) {
    cerr << errMsg << "operator[] " << idx << " out of bounds." << endl;
    return NULL;
  }

  Node* n = dm_list;

  // loop until reach index position
  for ( int i=0; i < idx; i++ )
    n = n->link;

  return n->data;
}


/*****
  Removes node at given position. Returns data for that node.
  Does bounds checking. Returns NULL if out of bounds.
*****/
void* LinkedList::remove( int idx )
{
  if ( idx<0 ||  idx>=dm_size ) {
    cerr << errMsg << "remove_entry() " << idx << " out of bounds." << endl;
    return NULL;
  }

  void* data;
  Node* n = dm_list;  // node to remove

  if ( idx == 0 ) 
  {
    dm_list = dm_list->link;
    dm_size--;
    data = n->data;
    delete n;
    return data;
  }

  Node* p;    // previous node
  for ( int i=0; i < idx; i++ ) {
    p = n;
    n = n->link;
  }
		
  p->link = n->link;  // skip around node to remove
  dm_size--;

  data = n->data; // data at node
  delete n;       // delete node, no longer needed
	
  return data;
}


/*****
  Removes all the nodes from the list, does not delete the data.
*****/
void LinkedList::clear_list()
{
  Node* n;	// node to delete

  while( dm_list != NULL )
  {
    n = dm_list;
    dm_list = dm_list->link;

    delete n;
  }

  dm_list = NULL;
  dm_size = 0;
}