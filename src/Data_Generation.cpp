#include "Data_Generation.hpp"

#include <algorithm>
#include <chrono>
#include <mpi.h>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"

using namespace libphysica::natural_units;

Simulation_Data::Simulation_Data(unsigned int sample_size, double u_min, unsigned int iso_rings)
: min_sample_size(sample_size), minimum_speed(u_min), isoreflection_rings(iso_rings), number_of_trajectories(0), number_of_free_particles(0), number_of_reflected_particles(0), number_of_captured_particles(0), average_number_of_scatterings(0.0), computing_time(0.0), number_of_data_points(std::vector<unsigned long int>(iso_rings, 0)), data(iso_rings, std::vector<libphysica::DataPoint>())
{
}

void Simulation_Data::Configure(double initial_radius, unsigned int min_scattering, unsigned int max_scattering, unsigned long int max_free_steps)
{
	initial_and_final_radius	  = initial_radius;
	minimum_number_of_scatterings = min_scattering;
	maximum_number_of_scatterings = max_scattering;
	maximum_free_time_steps		  = max_free_steps;
}

void Simulation_Data::Generate_Data(obscura::DM_Particle& DM, Solar_Model& solar_model, int mpi_rank, unsigned int fixed_seed)
{
	auto time_start = std::chrono::system_clock::now();

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	unsigned int smallest_sample_size = 0;
	while(smallest_sample_size < min_sample_size)
	{
		Event IC = Initial_Conditions(SHM, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
		Trajectory_Result trajectory = simulator.Simulate(IC, DM);

		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		if(trajectory.Particle_Captured())
			number_of_captured_particles++;
		else
		{
			if(trajectory.Particle_Free())
				number_of_free_particles++;
			else if(trajectory.Particle_Reflected())
				number_of_reflected_particles++;
			Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
			double v_final = trajectory.final_event.Speed();
			if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > minimum_speed)
			{
				unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
				number_of_data_points[isoreflection_ring]++;
				data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
			}
			smallest_sample_size = *std::min_element(std::begin(number_of_data_points), std::end(number_of_data_points));
			libphysica::Print_Progress_Bar(1.0 * smallest_sample_size / min_sample_size, mpi_rank, 61);
		}
	}
	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
}

void Simulation_Data::Generate_Data_RMA(obscura::DM_Particle& DM, Solar_Model& solar_model, unsigned int fixed_seed)
{
	auto time_start = std::chrono::system_clock::now();

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	// 1. Generate the data in parallel

	// Open an MPI Remote Memory Access window for the data counters
	int mpi_rank, mpi_processes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Win mpi_window;
	// std::vector<unsigned long int> counters(isoreflection_rings, 0);
	unsigned long int data_counters[isoreflection_rings];
	MPI_Win_create(&data_counters, isoreflection_rings * sizeof(unsigned long int), sizeof(unsigned long int), MPI_INFO_NULL, MPI_COMM_WORLD, &mpi_window);
	int mpi_target_rank	  = 0;
	unsigned long int one = 1;

	unsigned int smallest_sample_size = 0;
	while(smallest_sample_size < min_sample_size)
	{
		Event IC = Initial_Conditions(SHM, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
		Trajectory_Result trajectory = simulator.Simulate(IC, DM);

		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		if(trajectory.Particle_Captured())
			number_of_captured_particles++;
		else
		{
			if(trajectory.Particle_Free())
				number_of_free_particles++;
			else if(trajectory.Particle_Reflected())
				number_of_reflected_particles++;
			Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
			double v_final = trajectory.final_event.Speed();
			if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > minimum_speed)
			{
				unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
				number_of_data_points[isoreflection_ring]++;
				MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, mpi_window);
				MPI_Accumulate(&one, 1, MPI_UNSIGNED_LONG, mpi_target_rank, isoreflection_ring, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_window);
				MPI_Win_unlock(0, mpi_window);
				data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
			}
			MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, mpi_window);
			MPI_Get(&data_counters, isoreflection_rings, MPI_UNSIGNED_LONG, mpi_target_rank, 0, isoreflection_rings, MPI_UNSIGNED_LONG, mpi_window);
			MPI_Win_unlock(0, mpi_window);

			smallest_sample_size = data_counters[0];
			for(unsigned int i = 1; i < isoreflection_rings; i++)
				if(data_counters[i] < smallest_sample_size)
					smallest_sample_size = data_counters[i];
			libphysica::Print_Progress_Bar(1.0 * smallest_sample_size / min_sample_size, mpi_rank, 61);
		}
	}
	MPI_Win_free(&mpi_window);

	// 2. Reduce the worker's data sets.
	average_number_of_scatterings *= number_of_trajectories;
	MPI_Allreduce(MPI_IN_PLACE, &number_of_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_free_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_reflected_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &average_number_of_scatterings, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings /= number_of_trajectories;

	MPI_Datatype mpi_datapoint;
	MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_datapoint);
	MPI_Type_commit(&mpi_datapoint);
	std::vector<std::vector<libphysica::DataPoint>> global_data;
	for(unsigned int i = 0; i < isoreflection_rings; i++)
	{
		// 1. How many data points did this worker gather?
		unsigned long int local_number_of_data_points = data[i].size();
		// 2. How many data points did all workers gather in total?
		unsigned long int total_number_of_data_points;
		MPI_Allreduce(&number_of_data_points[i], &total_number_of_data_points, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
		number_of_data_points[i] = total_number_of_data_points;
		// 3. Every worker needs to know how much every worker did.
		std::vector<unsigned long int> data_points_of_workers(mpi_processes);
		MPI_Allgather(&local_number_of_data_points, 1, MPI_UNSIGNED_LONG, data_points_of_workers.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		//4. Collect info on the data packages to be received.
		std::vector<int> receive_counter(mpi_processes);
		std::vector<int> receive_displacements(mpi_processes);
		for(int j = 0; j < mpi_processes; j++)
		{
			receive_counter[j]		 = data_points_of_workers[j];
			receive_displacements[j] = (j == 0) ? 0 : receive_displacements[j - 1] + data_points_of_workers[j - 1];
		}
		// 5.
		std::vector<libphysica::DataPoint> ring_data(total_number_of_data_points);
		MPI_Allgatherv(&data[i].front(), local_number_of_data_points, mpi_datapoint, &ring_data.front(), receive_counter.data(), receive_displacements.data(), mpi_datapoint, MPI_COMM_WORLD);
		global_data.push_back(ring_data);
	}
	data = global_data;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void Simulation_Data::Generate_Data_Ring(obscura::DM_Particle& DM, Solar_Model& solar_model, unsigned int fixed_seed)
{
	auto time_start = std::chrono::system_clock::now();

	//MPI ring communication
	int mpi_rank, mpi_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int mpi_source		= (mpi_rank == 0) ? mpi_processes - 1 : mpi_rank - 1;
	int mpi_destination = (mpi_rank == mpi_processes - 1) ? 0 : mpi_rank + 1;
	int mpi_tag			= 0;
	MPI_Status mpi_status;
	MPI_Request mpi_request;

	// Configure the simulator
	Trajectory_Simulator simulator(solar_model, maximum_free_time_steps, maximum_number_of_scatterings, initial_and_final_radius);
	if(fixed_seed != 0)
		simulator.Fix_PRNG_Seed(fixed_seed);

	obscura::Standard_Halo_Model SHM;
	SHM.Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));

	//Get the MPI ring communication started
	std::vector<unsigned long int> local_counter_new(isoreflection_rings, 0);
	if(mpi_rank == 0)
		MPI_Isend(&number_of_data_points.front(), isoreflection_rings, MPI_UNSIGNED_LONG, mpi_destination, mpi_tag, MPI_COMM_WORLD, &mpi_request);

	unsigned int smallest_sample_size = 0;
	while(smallest_sample_size < min_sample_size)
	{
		Event IC = Initial_Conditions(SHM, solar_model, simulator.PRNG);
		Hyperbolic_Kepler_Shift(IC, initial_and_final_radius);
		Trajectory_Result trajectory = simulator.Simulate(IC, DM);

		number_of_trajectories++;
		average_number_of_scatterings = 1.0 / number_of_trajectories * ((number_of_trajectories - 1) * average_number_of_scatterings + trajectory.number_of_scatterings);

		if(trajectory.Particle_Captured())
			number_of_captured_particles++;
		else
		{
			if(trajectory.Particle_Free())
				number_of_free_particles++;
			else if(trajectory.Particle_Reflected())
				number_of_reflected_particles++;
			Hyperbolic_Kepler_Shift(trajectory.final_event, 1.0 * AU);
			double v_final = trajectory.final_event.Speed();
			if(trajectory.number_of_scatterings >= minimum_number_of_scatterings && v_final > minimum_speed)
			{
				unsigned int isoreflection_ring = (isoreflection_rings == 1) ? 0 : trajectory.final_event.Isoreflection_Ring(obscura::Sun_Velocity(), isoreflection_rings);
				local_counter_new[isoreflection_ring]++;
				data[isoreflection_ring].push_back(libphysica::DataPoint(v_final));
			}
			//Check if MPI token arrived.
			int mpi_flag;
			MPI_Iprobe(mpi_source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_flag, &mpi_status);
			if(mpi_flag)
			{
				//Receive the tokens
				MPI_Recv(&number_of_data_points.front(), isoreflection_rings, MPI_UNSIGNED_LONG, mpi_source, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
				//Update the tokens
				unsigned long int nMin1 = *std::min_element(std::begin(number_of_data_points), std::end(number_of_data_points));
				for(int i = 0; i < isoreflection_rings; i++)
				{
					number_of_data_points[i] += local_counter_new[i];
					local_counter_new[i] = 0;
				}
				smallest_sample_size = *std::min_element(std::begin(number_of_data_points), std::end(number_of_data_points));

				//Check if we are done
				if(nMin1 < min_sample_size && smallest_sample_size >= min_sample_size)
					mpi_tag = mpi_source + 1;
				else if(nMin1 >= min_sample_size)
					mpi_tag = mpi_status.MPI_TAG;
				//Pass on the tokens, unless you are the very last process.
				if(mpi_tag != (mpi_rank + 1))
					MPI_Isend(&number_of_data_points.front(), isoreflection_rings, MPI_UNSIGNED_LONG, mpi_destination, mpi_tag, MPI_COMM_WORLD, &mpi_request);
				libphysica::Print_Progress_Bar(1.0 * smallest_sample_size / min_sample_size, mpi_rank, 61);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// 2. Reduce the worker's data sets.
	average_number_of_scatterings *= number_of_trajectories;
	MPI_Allreduce(MPI_IN_PLACE, &number_of_trajectories, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_free_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_reflected_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &number_of_captured_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &average_number_of_scatterings, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	average_number_of_scatterings /= number_of_trajectories;

	MPI_Datatype mpi_datapoint;
	MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_datapoint);
	MPI_Type_commit(&mpi_datapoint);
	std::vector<std::vector<libphysica::DataPoint>> global_data;
	for(unsigned int i = 0; i < isoreflection_rings; i++)
	{
		// 1. How many data points did this worker gather?
		unsigned long int local_number_of_data_points = data[i].size();
		// 2. How many data points did all workers gather in total?
		unsigned long int total_number_of_data_points;
		MPI_Allreduce(&number_of_data_points[i], &total_number_of_data_points, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
		number_of_data_points[i] = total_number_of_data_points;
		// 3. Every worker needs to know how much every worker did.
		std::vector<unsigned long int> data_points_of_workers(mpi_processes);
		MPI_Allgather(&local_number_of_data_points, 1, MPI_UNSIGNED_LONG, data_points_of_workers.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		//4. Collect info on the data packages to be received.
		std::vector<int> receive_counter(mpi_processes);
		std::vector<int> receive_displacements(mpi_processes);
		for(int j = 0; j < mpi_processes; j++)
		{
			receive_counter[j]		 = data_points_of_workers[j];
			receive_displacements[j] = (j == 0) ? 0 : receive_displacements[j - 1] + data_points_of_workers[j - 1];
		}
		// 5.
		std::vector<libphysica::DataPoint> ring_data(total_number_of_data_points);
		MPI_Allgatherv(&data[i].front(), local_number_of_data_points, mpi_datapoint, &ring_data.front(), receive_counter.data(), receive_displacements.data(), mpi_datapoint, MPI_COMM_WORLD);
		global_data.push_back(ring_data);
	}
	data = global_data;

	auto time_end  = std::chrono::system_clock::now();
	computing_time = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	MPI_Allreduce(MPI_IN_PLACE, &computing_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

double Simulation_Data::Free_Ratio()
{
	return 1.0 * number_of_free_particles / number_of_trajectories;
}
double Simulation_Data::Capture_Ratio()
{
	return 1.0 * number_of_captured_particles / number_of_trajectories;
}
double Simulation_Data::Reflection_Ratio()
{
	return 1.0 * number_of_reflected_particles / number_of_trajectories;
}

double Simulation_Data::Lowest_Speed(unsigned int iso_ring) const
{
	return (*std::min_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

double Simulation_Data::Highest_Speed(unsigned int iso_ring) const
{
	return (*std::max_element(data[iso_ring].begin(), data[iso_ring].end())).value;
}

void Simulation_Data::Print_Summary(unsigned int mpi_rank)
{
	if(mpi_rank == 0)
	{
		unsigned long int number_of_data_points_tot = std::accumulate(number_of_data_points.begin(), number_of_data_points.end(), 0);
		std::cout << SEPARATOR
				  << "Simulation data summary" << std::endl
				  << std::endl
				  << "Configuration:" << std::endl
				  << "Minimum DM speed [km/sec]:\t" << libphysica::Round(In_Units(minimum_speed, km / sec)) << std::endl
				  << "Minimum sample size:\t\t" << min_sample_size << std::endl
				  << "Isoreflection rings:\t\t" << isoreflection_rings << std::endl
				  << std::endl
				  << "Results:" << std::endl
				  << "Simulated trajectories:\t\t" << number_of_trajectories << std::endl
				  << "Generated data points (total):\t" << number_of_data_points_tot << std::endl
				  << "Average # of scatterings:\t" << libphysica::Round(average_number_of_scatterings) << std::endl
				  << "Free particles [%]:\t\t" << libphysica::Round(100.0 * Free_Ratio()) << std::endl
				  << "Reflected particles [%]:\t" << libphysica::Round(100.0 * Reflection_Ratio()) << std::endl
				  << "Captured particles [%]:\t\t" << libphysica::Round(100.0 * Capture_Ratio()) << std::endl;

		if(isoreflection_rings > 1)
		{
			std::cout << std::endl
					  << "Ring\tData points\tRelative [%]\t<u> [km/sec]\tu_max [km/sec]" << std::endl;
			for(unsigned int i = 0; i < isoreflection_rings; i++)
			{
				double rel_number_of_data_points = 100.0 * number_of_data_points[i] / number_of_data_points_tot;
				std::vector<double> u_average	 = libphysica::Weighted_Average(data[i]);
				double u_max					 = Highest_Speed(i);
				std::cout << i + 1 << "\t" << number_of_data_points[i] << "\t\t" << libphysica::Round(rel_number_of_data_points) << "\t\t" << libphysica::Round(In_Units(u_average[0], km / sec)) << " +- " << libphysica::Round(In_Units(u_average[1], km / sec)) << "\t" << libphysica::Round(In_Units(u_max, km / sec)) << std::endl;
			}
		}
		else
		{
			std::vector<double> u_average = libphysica::Weighted_Average(data[0]);
			double u_max				  = Highest_Speed();
			std::cout << "<u> [km/sec]:\t\t\t" << libphysica::Round(In_Units(u_average[0], km / sec)) << " +- " << libphysica::Round(In_Units(u_average[1], km / sec)) << std::endl
					  << "u_max [km/sec]:\t\t\t" << libphysica::Round(In_Units(u_max, km / sec)) << std::endl;
		}
		std::cout << std::endl
				  << "Trajectory rate [1/s]:\t\t" << libphysica::Round(1.0 * number_of_trajectories / computing_time) << std::endl
				  << "Data generation rate [1/s]:\t" << libphysica::Round(1.0 * number_of_data_points_tot / computing_time) << std::endl;

		if(computing_time > 60.0)
			std::cout
				<< "Simulation time:\t\t"
				<< "[" << floor(computing_time / 3600.0) << ":" << floor(fmod(computing_time / 60.0, 60.0)) << ":" << floor(fmod(computing_time, 60.0)) << "]." << std::endl;
		else
			std::cout << "Simulation time [s]:\t\t" << libphysica::Round(computing_time) << std::endl;

		std::cout << SEPARATOR << std::endl;
	}
}