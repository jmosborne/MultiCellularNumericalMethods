/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTNUMERICSCRYPT3D_HPP_
#define TESTNUMERICSCRYPT3D_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/lexical_cast.hpp>

//Misc
#include "CellBasedEventHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimpleWntUniformCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "CommandLineArguments.hpp"
#include "WntConcentrationModifier.hpp"
#include "PlaneBasedTimedCellKiller.hpp" 
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellBasedSimulationArchiver.hpp"

//Populations
#include "AbstractOffLatticeCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

//Forces etc
#include "GeneralisedLinearSpringForce.hpp"
#include "CryptBoundaryForce.hpp"

//Numerical Methods
#include "ForwardEulerNumericalMethod.hpp"
#include "BackwardEulerNumericalMethod.hpp"
#include "AdamsMoultonNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"
#include "RK3NumericalMethod.hpp"
#include "MidPointNumericalMethod.hpp"

#include "CellAgesWriter.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
* Tests different numerical method and cell population combinations. 
* Tracks the position of cell centres at regular intervals. 
*/

class TestNumericsCrypt3d : public AbstractCellBasedTestSuite
{

private:
 
	void ResetForNewReplicate(int seed)
	{    
        RandomNumberGenerator::Instance()->Reseed(seed);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellId::ResetMaxCellId();
        CellBasedEventHandler::Reset();
    };

	void SetBaseFilename(std::string& filename, std::string simulation_type, bool isStochasticCCM, std::string method, int rep)
	{
    	std::stringstream repAsString;
		repAsString << rep;

		filename = "NumericalConvergence/" 
		           + simulation_type 
				   + "/StochasticCCM_" + boost::lexical_cast<std::string>(isStochasticCCM) 
				   + "/" + method
				   + "/SteadyState"
				   + "/Run_" + repAsString.str();	
    };


    void SetFilename(std::string& filename, std::string simulation_type, bool isStochasticCCM, std::string method, int rep, int power, double amt)
	{
    	std::stringstream repAsString;
		repAsString << rep;
		std::stringstream powAsString;
		powAsString << power;
		std::stringstream amtAsString;
		amtAsString << amt;
		
		filename = "NumericalConvergence/" 
		           + simulation_type 
				   + "/StochasticCCM_" + boost::lexical_cast<std::string>(isStochasticCCM) 
				   + "/" + method
				   + "/AMT_" + amtAsString.str()
				   + "/Power_" + powAsString.str()
				   + "/Compression_1"
				   + "/Run_" + repAsString.str(); 
    };

    template<unsigned DIM>
    boost::shared_ptr<AbstractNumericalMethod<DIM> > MakeNumericalMethod(std::string method,bool is_adaptive)
	{
    	if(method=="FE")
		{
    		MAKE_PTR(ForwardEulerNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
    		return p_numerical_method;
		}
		else if(method=="MP")
		{
			MAKE_PTR(MidPointNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
    		return p_numerical_method;
		}
		else if(method=="RK3")
		{
			MAKE_PTR(RK3NumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
    		return p_numerical_method;
		}
		else if(method=="RK4")
		{
			MAKE_PTR(RK4NumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
    		return p_numerical_method;
		}
		else if(method=="BE_Tol_5")
		{
			MAKE_PTR(BackwardEulerNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
			p_numerical_method->SetSolverTolerance(1e-5);
    		return p_numerical_method;
		}
		
		else if(method=="BE_Tol_10")
		{
			MAKE_PTR(BackwardEulerNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
			p_numerical_method->SetSolverTolerance(1e-10);
    		return p_numerical_method;
		}
		else if(method=="AM2_Tol_5")
		{
			MAKE_PTR(AdamsMoultonNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
			p_numerical_method->SetSolverTolerance(1e-5);
    		return p_numerical_method;
		}
		else if(method=="AM2_Tol_10")
		{
			MAKE_PTR(AdamsMoultonNumericalMethod<DIM>, p_numerical_method);
			p_numerical_method->SetUseAdaptiveTimestep(is_adaptive);
			p_numerical_method->SetSolverTolerance(1e-10);
    		return p_numerical_method;
		}
		else
		{
			EXCEPTION("Unrecognised numerical method");
		}
		NEVER_REACHED;
		MAKE_PTR(ForwardEulerNumericalMethod<DIM>, defaultNm);
		return defaultNm;
    }

    template<unsigned DIM>
	void GenerateInitialCellsWithDistributedAges(std::vector<CellPtr>& rCells,
                                              const std::vector<unsigned> realCellIndices,
                                              double ccLength)
	{
    	assert(!realCellIndices.empty());
    	unsigned nCells = realCellIndices.size();

		double max_timestep = pow(2,-4);

    	rCells.clear();
    	rCells.reserve(nCells);
    	CellPropertyRegistry::Instance()->Clear();

    	for (unsigned i=0; i<nCells; i++)
    	{	
			double birthTime = 0.0;
		
			AbstractCellCycleModel* pCCM = new SimpleWntUniformCellCycleModel();
			
			dynamic_cast<SimpleWntUniformCellCycleModel*>(pCCM)->SetMinCellCycleDuration(ccLength-5);
			dynamic_cast<SimpleWntUniformCellCycleModel*>(pCCM)->SetMaxCellCycleDuration(ccLength+5); //i.e U(n-1,n+1)
			//if stochastic have random birth events but only at largest timesteps
			birthTime = - max_timestep*ceil(ccLength*RandomNumberGenerator::Instance()->ranf()/max_timestep);
			
    		pCCM->SetDimension(DIM);

    	    boost::shared_ptr<AbstractCellProperty> pState(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    	    CellPtr pCell(new Cell(pState, pCCM));

			pCell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

    	    pCell->SetBirthTime(birthTime);
    	    rCells.push_back(pCell);
    	}
	}


	template<unsigned DIM>
	void SetAbsoluteMovementThreshold(OffLatticeSimulation<DIM>* simulation, float AMT){
		AbstractCellPopulation<DIM,DIM>* pop = &(simulation->rGetCellPopulation());
        dynamic_cast<AbstractOffLatticeCellPopulation<DIM,DIM>*>(pop)->SetAbsoluteMovementThreshold((double)AMT);
	}


	template<unsigned DIM>
    void AddTracking(OffLatticeSimulation<DIM>* p_simulation, int interval)
	{
		// First clear the modifiers to make sure only one output step 
		p_simulation->GetSimulationModifiers()->clear();

    	MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<DIM>, trackingModifier,(interval,1));
    	p_simulation->AddSimulationModifier(trackingModifier);

		MAKE_PTR(WntConcentrationModifier<DIM>, p_wnt_modifier);
		p_wnt_modifier->SetCryptLength(2.0);
		p_simulation->AddSimulationModifier(p_wnt_modifier);
		
		boost::shared_ptr<CellDataItemWriter<DIM,DIM> > p_cell_data_item_writer(new CellDataItemWriter<DIM,DIM>("wnt_level"));
		p_simulation->rGetCellPopulation().AddCellWriter(p_cell_data_item_writer);

    }

	template<unsigned DIM>
    void CleanUpNodes(std::vector< Node<DIM>* > nodes)
    {
    	for(unsigned i=0; i<nodes.size(); i++)
		{
    		delete nodes[i];
    	}
    };


public:

	void TestWithPositionRecording() 
	{

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-step_range_lower"));
		unsigned step_range_lower = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-step_range_lower");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-step_range_upper"));
		unsigned step_range_upper = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-step_range_upper");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-min_amt"));
		unsigned min_amt = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-min_amt");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-seed_range_lower"));
		unsigned seed_range_lower = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed_range_lower");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-seed_range_upper"));
		unsigned seed_range_upper = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed_range_upper");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-method"));
		std::string method = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-method");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ccd"));
		double ccLength = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ccd");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-end_time"));
		double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-steady_state"));
		double steady_state_end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-steady_state");


		std::string simulation_type = "3d_crypt_steady";		
	    
		// int step_range_lower = 6; //FE and RK4 Fail for larger steps
		// int step_range_upper = 6; // 14;
		// int seed_range_lower = 0;
		// int seed_range_upper = 0;
		// std::string method = "FE";	
		// bool isAdaptive = true; 
		// double ccLength = 10.0;
		// double end_time = 50.0;
		// double steady_state_end_time = 100.0;
		
		bool isStochasticCCM = true; // Just used to keep folder structure sensible 
				
		// SETUP CONSTANTS

		double crypt_length = 4.0;
		double interactionCutoff = 5.0;
		// int nCells = 1;
        double dampingNormal = 1.0;

		unsigned num_amts = (unsigned)min_amt+1;
		// double AMTs[num_amts]; 
		// for (unsigned power = 0; power<num_amts; power++)
		// {
		// 	AMTs[power] = ((double)pow(10.0,-(double)power));
		// 	PRINT_3_VARIABLES(min_amt,power,AMTs[power]);
		// }
		bool isAdaptive = false;
		if (num_amts>1)
		{
			isAdaptive = true;	
		}
	
		for(unsigned random_seed = seed_range_lower; random_seed <= seed_range_upper; random_seed++)
		{
			//Run to steady state with RK4 and reasonable timestep 
			//double steady_state_end_time = 50;
			double default_power = 6;		
			double default_dt = 1.0/((double)pow(2,default_power));
			std::string default_method = "RK4";
			int default_update_interval = (int)pow(2,default_power-4); // Remesh etc every 2^-4 hours
			int default_output_interval = (int)pow(2,default_power); // Sample every hour
			double default_AMT = 1.0;
			
			ResetForNewReplicate(random_seed);

			std::string base_filename;
			SetBaseFilename(base_filename, simulation_type, isStochasticCCM, method, random_seed);
			std::cout << base_filename << std::endl;

			
			// Run to steady state with RK4
			boost::shared_ptr<AbstractNumericalMethod<3> > numericalMethod = MakeNumericalMethod<3>(default_method,isAdaptive);
			
			OffLatticeSimulation<3>* p_simulation;
			std::vector<Node<3>*> nodes;
			std::vector<unsigned> realCellIndices;
			std::vector<CellPtr> cells;
			AbstractOffLatticeCellPopulation<3,3>* p_cell_population;
			
			// Initialy start with one cell
			assert(nCells==1);
			nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
			realCellIndices.push_back(0);
						
			MAKE_PTR(NodesOnlyMesh<3>, pMesh);
			pMesh->ConstructNodesWithoutMesh(nodes, interactionCutoff);

			GenerateInitialCellsWithDistributedAges<3>(cells, realCellIndices, ccLength);       
			
			p_cell_population = new NodeBasedCellPopulation<3>(*pMesh, cells);
			p_cell_population->SetDampingConstantNormal(dampingNormal);
			p_cell_population->AddCellWriter<CellAgesWriter>();

			p_simulation = new OffLatticeSimulation<3>(*p_cell_population, false, true);
			p_simulation->SetNumericalMethod(numericalMethod);
			p_simulation->SetOutputDirectory(base_filename.c_str());
			p_simulation->SetDt(default_dt);
			p_simulation->SetSamplingTimestepMultiple(default_output_interval);
			p_simulation->SetUpdatingTimestepMultiple(default_update_interval);
			p_simulation->SetEndTime(steady_state_end_time);
			p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.

			// Force based boundary
			MAKE_PTR(CryptBoundaryForce<3>, p_force);
			p_simulation->AddForce(p_force);

			//Add cell killer
			MAKE_PTR_ARGS(PlaneBasedTimedCellKiller<3>, p_base_killer, (p_cell_population, crypt_length*unit_vector<double>(3,2), unit_vector<double>(3,2), default_update_interval, base_filename));
			p_simulation->AddCellKiller(p_base_killer);

			SetAbsoluteMovementThreshold<3>(p_simulation, default_AMT);
			
			MAKE_PTR(GeneralisedLinearSpringForce<3>, pForce);
			pForce->SetCutOffLength(interactionCutoff);
			pForce->SetMeinekeDivisionRestingSpringLength(1.0); 
			p_simulation->AddForce(pForce);

			AddTracking<3>(p_simulation, default_output_interval);
			
			Timer::Reset();
			p_simulation->Solve();
			PRINT_VARIABLE(Timer::GetElapsedTime());
			p_simulation->RemoveAllCellKillers();
			
			// Save simulation in steady state
			CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(p_simulation);

			// Avoid Memory Leaks
			CleanUpNodes<3>(nodes);
			delete p_simulation;
			delete p_cell_population;

			for(unsigned amt_index = 0; amt_index < num_amts; amt_index++)
			{
				//double AMT = AMTs[amt_index];
				double AMT = ((double)pow(10,-amt_index));
		
				for(unsigned power = step_range_lower; power <= step_range_upper; power+=2)
				{
					Timer::Reset();
					CellId::ResetMaxCellId();
        			CellBasedEventHandler::Reset();
					double dt = 1.0/((double)pow(2,power));
					assert(step_range_lower>=4);
					int update_interval = (int)pow(2,power-4); // Remesh every 2^-4 hours
					int output_interval = (int)pow(2,power); // Output every hour

					
					std::string filename;
					SetFilename(filename, simulation_type, isStochasticCCM, method, random_seed, power, AMT);
					std::cout << filename << std::endl;
				    OffLatticeSimulation<3>* p_loaded_simulation = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load(base_filename,steady_state_end_time);
					boost::shared_ptr<AbstractNumericalMethod<3> > numericalMethod = MakeNumericalMethod<3>(method,isAdaptive);
					
					p_loaded_simulation->SetNumericalMethod(numericalMethod);
					p_loaded_simulation->SetOutputDirectory(filename.c_str());
					p_loaded_simulation->SetDt(dt);
					p_loaded_simulation->SetSamplingTimestepMultiple(output_interval);
					p_loaded_simulation->SetUpdatingTimestepMultiple(update_interval);
					p_loaded_simulation->SetEndTime(steady_state_end_time+end_time);
						
					//Reset cell killer so at correct interval 
					p_loaded_simulation->RemoveAllCellKillers();
					MAKE_PTR_ARGS(PlaneBasedTimedCellKiller<3>, p_killer_2, (&(p_loaded_simulation->rGetCellPopulation()), crypt_length*unit_vector<double>(3,2), unit_vector<double>(3,2), update_interval, filename));
					
					p_loaded_simulation->AddCellKiller(p_killer_2);

					SetAbsoluteMovementThreshold<3>(p_loaded_simulation, AMT);
					AddTracking<3>(p_loaded_simulation,output_interval);

					Timer::Reset();

					try
					{
						p_loaded_simulation->Solve();
					}
					catch (Exception& e)
					{
						// If it throws then we report the error message and go to the next simulation
						WARNING("Simulation didnt run" << filename << ".");
						WARNING(e.GetMessage());
					}

					OutputFileHandler output_file_handler(filename, false);
					out_stream p_stream = output_file_handler.OpenOutputFile("timing.dat");
					*p_stream <<  Timer::GetElapsedTime();
					p_stream->close();

					PRINT_VARIABLE(Timer::GetElapsedTime());

					// Avoid Memory Leaks
					delete p_loaded_simulation;
				}

			}
		}
	};

};

#endif /*TESTNUMERICSCRYPT3D_HPP_*/