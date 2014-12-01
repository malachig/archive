#ifndef LINKED_LIST_HPP
#define LINKED_LIST_HPP


class LinkedList
{
public:
  LinkedList();
  LinkedList( const LinkedList &l );
  ~LinkedList();

  const LinkedList& operator=(const LinkedList &l);

  int		size() const { return dm_size; }	

  void*		operator[]( int idx );
  const void* operator[]( int idx ) const;

  void	add( void *data );
  void*	remove( int idx );
  void	clear_list();

private:

  struct Node {
    void* data;
    Node*	link;
  };
	
  int		dm_size;
  Node*	dm_list;
};


#endif