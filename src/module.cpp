#include "module.hpp"

// Headers from libphysica
#include "Numerics.hpp"

// Headers from obscura
#include "Astronomy.hpp"

int fib(int n)
{ 
	if (n <= 1) 
		return n;
	else
		return fib(n-1) + fib(n-2); 
} 