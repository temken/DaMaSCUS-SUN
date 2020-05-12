#include "gtest/gtest.h"

#include "module.hpp"

TEST(TestFib, ReturnsCorrectFibonacciNumber)
{
	//ARRANGE
	//ACT
	//ASSERT
    EXPECT_EQ( fib(10),  55);
}
