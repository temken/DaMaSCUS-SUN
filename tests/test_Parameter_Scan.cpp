#include "Parameter_Scan.hpp"

#include "gtest/gtest.h"
#include <mpi.h>

// Headers from libphysica
#include "Natural_Units.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
{
	int result = 0;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}