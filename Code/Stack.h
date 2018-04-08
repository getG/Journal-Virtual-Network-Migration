#ifndef STACK_H
#define STACK_H
#include <memory>

using namespace std;

template<typename T>
class Stack
{
private:
	T * m_elem;
	int m_size;
	int m_capacity;

public:
	// Constructor
	Stack(int size = 0, const T & t = T())
		:m_elem(NULL), m_size(0), m_capacity(0)
	{
		m_elem = (T*)malloc(sizeof(T)*size);
		uninitialized_fill(m_elem, m_elem + size, t);
		m_size = m_capacity = size;
	}

	// Destructor
	virtual ~Stack()
	{
		if (m_elem)
		{
			free(m_elem);
			m_elem = NULL;
		}
	}

	int size() const
	{
		return m_size;
	}

	int capacity() const
	{
		return m_capacity;
	}

	bool empty() const
	{
		return m_size == 0;
	}

	void clear()
	{
		m_size = 0;
	}

	const T & top() const
	{
		return m_elem[m_size - 1];
	}

	void pop()
	{
		m_size--;
	}

	T & operator[](int index)
	{
		return m_elem[index];
	}

	bool full() const
	{
		return m_size == m_capacity;
	}

	int next_capacity() const
	{
		if (m_capacity <= 0)
			return 1;
		else
			return 2 * m_capacity;
	}

	void reserve(int size)
	{
		if (size > m_capacity)
		{
			m_elem = (T*)realloc(m_elem, sizeof(T)*size);
			m_capacity = size;
		}
	}

	void push(const T & t)
	{
		if (full())
		{
			reserve(next_capacity());
		}
		new(m_elem + m_size)T(t);
		m_size++;
	}
};//end of class Stack
#endif
