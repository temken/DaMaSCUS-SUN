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

TEST(TestParameterScan, TestConfiguration)
{
	// ARRANGE
	// ACT & ASSERT
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	// ASSERT
	EXPECT_TRUE(cfg.compute_halo_constraints);
	EXPECT_EQ(cfg.sample_size, 50);
	EXPECT_DOUBLE_EQ(cfg.cross_section_min, 1.0e-35 * cm * cm);
	EXPECT_DOUBLE_EQ(cfg.cross_section_max, 1.0e-32 * cm * cm);
	EXPECT_EQ(cfg.cross_sections, 5);
	EXPECT_EQ(cfg.interpolation_points, 150);
	EXPECT_EQ(cfg.isoreflection_rings, 3);
}

TEST(TestParameterScan, TestConfigurationSummary)
{
	// ARRANGE
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	// ACT & ASSERT
	cfg.Print_Summary(0);
	// ASSERT
}
TEST(TestParameterScan, TestSTAScan)
{
	// ARRANGE
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	Solar_Model SSM;
	// ACT
	Parameter_Scan scan(cfg);
	scan.Perform_STA_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, 0);
	std::vector<std::vector<double>> limit_curve = scan.Limit_Curve();
	scan.Print_Grid();
	// ASSERT
	ASSERT_GT(limit_curve.size(), 0);
	// for(auto& row : scan.p_value_grid)
	// 	for(auto& entry : row)
	// 		ASSERT_GE(entry, 0.0);
}

TEST(TestParameterScan, TestFullScan)
{
	// ARRANGE
	Configuration cfg(PROJECT_DIR "tests/config_unittest.cfg", 1);
	Solar_Model SSM;
	// ACT
	Parameter_Scan scan(cfg);
	scan.Perform_Full_Scan(*cfg.DM, *cfg.DM_detector, SSM, *cfg.DM_distr, 1);
	std::vector<std::vector<double>> limit_curve = scan.Limit_Curve();
	scan.Print_Grid();
	// ASSERT
	ASSERT_GT(limit_curve.size(), 0);
	// for(auto& row : scan.p_value_grid)
	// 	for(auto& entry : row)
	// 		ASSERT_GE(entry, 0.0);
}