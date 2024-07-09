#ifndef TESTNUMERICS_HPP_
#define TESTNUMERICS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/lexical_cast.hpp>

//Misc
#include "CellBasedEventHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "CommandLineArguments.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "Debug.hpp"

//Populations
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

//Forces etc
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
//#include "StableRandomDirectionVertexBasedDivisionRule.hpp"

//Numerical Methods
#include "ForwardEulerNumericalMethod.hpp"
#include "BackwardEulerNumericalMethod.hpp"
#include "AdamsMoultonNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"
#include "RK3NumericalMethod.hpp"
#include "MidPointNumericalMethod.hpp"

#include "CellAgesWriter.hpp"
#include "Warnings.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
* Tests different numerical method and cell population combinations. 
* Tracks the position of cell centres at regular intervals. 
*/

class TestNumerics : public AbstractCellBasedTestSuite
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


    void SetFilename(std::string& filename, std::string simulation_type, bool isStochasticCCM, std::string method, int power, int rep, double amt, double compression)
	{
    	std::stringstream repAsString;
		repAsString << rep;
		std::stringstream powAsString;
		powAsString << power;
		std::stringstream amtAsString;
		amtAsString << amt;
	    std::stringstream compAsString;
		compAsString << compression;
		
		filename = "NumericalConvergence/" 
		           + simulation_type 
				   + "/StochasticCCM_" + boost::lexical_cast<std::string>(isStochasticCCM) 
				   + "/" + method
				   + "/AMT_" + amtAsString.str()
				   + "/Power_" + powAsString.str()
				   + "/Compression_" + compAsString.str()
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
                                              bool isStochasticCCM,
                                              double ccLength,
											  bool isWntDependent)
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

			AbstractCellCycleModel* pCCM;
		
			pCCM = new UniformCellCycleModel();
		
			if(isStochasticCCM)
			{
				assert(ccLength >5);
				dynamic_cast<UniformCellCycleModel*>(pCCM)->SetMinCellCycleDuration(ccLength-5);
				dynamic_cast<UniformCellCycleModel*>(pCCM)->SetMaxCellCycleDuration(ccLength+5); //i.e U(n-5,n+5)
				//if stochastic have random birth events but only at largest timesteps
				birthTime = - max_timestep*ceil(ccLength*RandomNumberGenerator::Instance()->ranf()/max_timestep);
			}
			else
			{
				assert(ccLength >=1);
				dynamic_cast<UniformCellCycleModel*>(pCCM)->SetMinCellCycleDuration(ccLength);
				dynamic_cast<UniformCellCycleModel*>(pCCM)->SetMaxCellCycleDuration(ccLength); // i.e fixed duration
				birthTime = -(double)i;
			}
			
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
    void AddTracking(OffLatticeSimulation<DIM>* simulation, int interval){
    	MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<DIM>, trackingModifier,(interval,1));
    	simulation->AddSimulationModifier(trackingModifier);
    }


    template<unsigned DIM>
    void AddForce(OffLatticeSimulation<DIM>* simulation, std::string simulation_type, double cutoff){
    	
    	if(simulation_type=="3d_node" || simulation_type=="2d_node" || simulation_type=="2d_mesh")
		{
    		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, pForce);
    		pForce->SetCutOffLength(cutoff);
			pForce->SetMeinekeDivisionRestingSpringLength(1.0); 
    		simulation->AddForce(pForce);
    	}
		else if(simulation_type=="1d_compression")
		{
    		MAKE_PTR(GeneralisedLinearSpringForce<DIM>, pForce);
    		pForce->SetCutOffLength(cutoff);
			pForce->SetMeinekeDivisionRestingSpringLength(1.0); //so can compare to exact solution
    		simulation->AddForce(pForce);
    	}
		else if(simulation_type=="2d_vertex")
		{
			MAKE_PTR(NagaiHondaForce<DIM>, pForce);
			// Change parameters so tissue forms rounded shape.
			pForce->SetNagaiHondaDeformationEnergyParameter(1.0); //1.0
        	pForce->SetNagaiHondaMembraneSurfaceEnergyParameter(0.1); // 0.1
        	pForce->SetNagaiHondaCellCellAdhesionEnergyParameter(0.1); //0.05
        	pForce->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.2); //0.1
			simulation->AddForce(pForce);
		}
		else
		{
			NEVER_REACHED;
		}
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

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-simulation_type"));
		// std::string simulation_type  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-simulation_type");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-step_range_lower"));
		// unsigned step_range_lower = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-step_range_lower");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-step_range_upper"));
		// unsigned step_range_upper = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-step_range_upper");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-min_amt"));
		// unsigned min_amt = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-min_amt");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-seed_range_lower"));
		// unsigned seed_range_lower = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed_range_lower");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-seed_range_upper"));
		// unsigned seed_range_upper = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed_range_upper");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-method"));
		// std::string method = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-method");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ccd"));
		// double ccLength = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ccd");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-stochastic"));
		// bool isStochasticCCM = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-stochastic");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-end_time"));
		// double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");

		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-compression"));
		// double compression = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-compression");


		std::string simulation_type = "1d_compression";		
	 	unsigned step_range_lower = 6; 
		unsigned step_range_upper = 20; // 14;
		unsigned min_amt = 0;
		unsigned seed_range_lower = 0;
		unsigned seed_range_upper = 0;
		std::string method = "RK4";	
		double ccLength = 12.0;
		bool isStochasticCCM = true;
		double end_time = 50.0;
		double compression = 1.0;
		
		// SETUP CONSTANTS

		double interactionCutoff = 20.0; //Large enough so no change in neighbours
		int nCellsPerSide = 2;
    	int nCells = nCellsPerSide*nCellsPerSide;
        double dampingNormal = 1;


		unsigned num_amts = (unsigned)min_amt+1;
		// double AMTs[num_amts]; 
		// for (unsigned power = 0; power<num_amts; power++)
		// {
		// 	AMTs[power] = ((double)pow(10,-power));
		// 	PRINT_2_VARIABLES(power,AMTs[power]);
		// }
		bool isAdaptive = false;
		if (num_amts>1)
		{
			isAdaptive = true;	
		}
		
		for(unsigned amt_index = 0; amt_index < num_amts; amt_index++)
		{
			double AMT = ((double)pow(10,-amt_index));
			//double AMT = AMTs[amt_index];

			for(unsigned power = step_range_lower; power <= step_range_upper; power+=2)
			{
				for(unsigned random_seed = seed_range_lower; random_seed <= seed_range_upper; random_seed++)
				{
					Timer::Reset();

					ResetForNewReplicate(random_seed);
					std::string filename;
					SetFilename(filename, simulation_type, isStochasticCCM, method, power, random_seed, AMT, compression);
					std::cout << filename << std::endl;

					double dt = 1.0/((double)pow(2,power));
					assert(step_range_lower>=6);
					int interval = (int)pow(2,power-4); // Remesh every 2^-4 hours
										
					// =================================================================================================== 
					// 1D Mesh Based compression SIMULATION ============================================================== 
					// =================================================================================================== 
					
					if(simulation_type=="1d_compression")
					{
						unsigned num_elements=nCellsPerSide-1; //20

						if(ccLength>=100) // i.e no proliferation
						{
							num_elements = 10;
						}

						boost::shared_ptr<AbstractNumericalMethod<1> > numericalMethod = MakeNumericalMethod<1>(method,isAdaptive);
					
						OffLatticeSimulation<1>* p_simulation;
						std::vector<unsigned> realCellIndices;
						std::vector<CellPtr> cells;
						AbstractCentreBasedCellPopulation<1,1>* p_cell_population;

						MutableMesh<1,1> mesh;
						mesh.ConstructRegularSlabMesh(1,num_elements);
						mesh.Translate(-(double)num_elements/2.0);
						mesh.Scale(compression);

						double num_cells = mesh.GetNumNodes();

						for(int n=0; n<num_cells; n++)
						{
							realCellIndices.push_back((unsigned)n);
						}

						GenerateInitialCellsWithDistributedAges<1>(cells, realCellIndices, isStochasticCCM, ccLength, (simulation_type=="1d_crypt")?(true):(false)); 
						
						p_cell_population = new MeshBasedCellPopulation<1>(mesh, cells);
						p_cell_population->SetDampingConstantNormal(dampingNormal);
						
						p_simulation = new OffLatticeSimulation<1>(*p_cell_population, false, true);
						p_simulation->SetNumericalMethod(numericalMethod);
						p_simulation->SetOutputDirectory(filename.c_str());
						p_simulation->SetDt(dt);
						p_simulation->SetSamplingTimestepMultiple((int)pow(2,power));
						p_simulation->SetEndTime(end_time);
						p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.
						p_simulation->SetUpdatingTimestepMultiple(interval);

						SetAbsoluteMovementThreshold<1>(p_simulation, AMT);
						AddForce<1>(p_simulation, simulation_type, interactionCutoff);
						assert(power>=6);
						if(ccLength>=100) // i.e no proliferation
						{
							AddTracking<1>(p_simulation, interval); // Output as much as possible
						}
						else
						{
							AddTracking<1>(p_simulation,(int)pow(2,power)); //Output every Hour
						}

						try
						{
							p_simulation->Solve();
						}
						catch (Exception& e)
						{
							// If it throws then we report the error message and go to the next simulation
							WARNING("Simulation didnt run" << filename << ".");
							WARNING(e.GetMessage());
						}

						// Avoid Memory Leaks
						delete p_simulation; 
						delete p_cell_population;
					}


					// =================================================================================================== 
					// 3D Node Simulation  ===============================================================================
					// =================================================================================================== 
					else if(simulation_type=="3d_node")
					{
						boost::shared_ptr<AbstractNumericalMethod<3> > numericalMethod = MakeNumericalMethod<3>(method,isAdaptive);
						
						OffLatticeSimulation<3>* p_simulation;
						std::vector<Node<3>*> nodes;
						std::vector<unsigned> realCellIndices;
						std::vector<CellPtr> cells;
						AbstractOffLatticeCellPopulation<3,3>* p_cell_population;
						
						assert(nCells==4); // If not need to change the below initialisation
						nodes.push_back(new Node<3>((unsigned)0, false, 0, 0, 0.0));
						nodes.push_back(new Node<3>((unsigned)1, false, 0, 1, 0.0));
						nodes.push_back(new Node<3>((unsigned)2, false, 1, 0, 0.0));
						nodes.push_back(new Node<3>((unsigned)3, false, 1, 1, 0.0));
						realCellIndices.push_back((unsigned)0);
						realCellIndices.push_back((unsigned)1);
						realCellIndices.push_back((unsigned)2);
						realCellIndices.push_back((unsigned)3);
						
						MAKE_PTR(NodesOnlyMesh<3>, pMesh);
						pMesh->ConstructNodesWithoutMesh(nodes, interactionCutoff);

						pMesh->Scale(compression,compression,compression);

						GenerateInitialCellsWithDistributedAges<3>(cells, realCellIndices, isStochasticCCM, ccLength, false);       
						
						p_cell_population = new NodeBasedCellPopulation<3>(*pMesh, cells);
						p_cell_population->SetDampingConstantNormal(dampingNormal);
						p_cell_population->AddCellWriter<CellAgesWriter>();
				
						p_simulation = new OffLatticeSimulation<3>(*p_cell_population, false, true);
						p_simulation->SetNumericalMethod(numericalMethod);
						p_simulation->SetOutputDirectory(filename.c_str());
						p_simulation->SetDt(dt);
						p_simulation->SetSamplingTimestepMultiple((int)pow(2,power));
						p_simulation->SetEndTime(end_time);
						p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.
						p_simulation->SetUpdatingTimestepMultiple(interval);

						SetAbsoluteMovementThreshold<3>(p_simulation, AMT);
						AddForce<3>(p_simulation, simulation_type, interactionCutoff);
						AddTracking<3>(p_simulation, (int)pow(2,power));

						try
						{
							p_simulation->Solve();
						}
						catch (Exception& e)
						{
							// If it throws then we report the error message and go to the next simulation
							WARNING("Simulation didnt run" << filename << ".");
							WARNING(e.GetMessage());
						}

						// Avoid Memory Leaks
						CleanUpNodes<3>(nodes);
						delete p_simulation;
						delete p_cell_population;
					}

					// =================================================================================================== 
					// 2D Node SIMULATION ================================================================================ 
					// =================================================================================================== 
					else if(simulation_type=="2d_node")
					{

						boost::shared_ptr<AbstractNumericalMethod<2> > numericalMethod = MakeNumericalMethod<2>(method,isAdaptive);
						
						OffLatticeSimulation<2>* p_simulation;
						std::vector<Node<2>*> nodes;
						std::vector<unsigned> realCellIndices;
						std::vector<CellPtr> cells;
						AbstractOffLatticeCellPopulation<2,2>* p_cell_population;

						for(int n=0; n<nCells; n++){
							nodes.push_back(new Node<2>((unsigned)n, false, n/nCellsPerSide, n%nCellsPerSide, 0.0));
							realCellIndices.push_back((unsigned)n);
						}
						
						MAKE_PTR(NodesOnlyMesh<2>, pMesh);
						pMesh->ConstructNodesWithoutMesh(nodes, interactionCutoff);

						pMesh->Scale(compression,compression);


						GenerateInitialCellsWithDistributedAges<2>(cells, realCellIndices, isStochasticCCM, ccLength, false);       
						
						p_cell_population = new NodeBasedCellPopulation<2>(*pMesh, cells);
						p_cell_population->SetDampingConstantNormal(dampingNormal);
						p_cell_population->AddCellWriter<CellAgesWriter>();

						p_simulation = new OffLatticeSimulation<2>(*p_cell_population, false, true);
						p_simulation->SetNumericalMethod(numericalMethod);
						p_simulation->SetOutputDirectory(filename.c_str());
						p_simulation->SetDt(dt);
						p_simulation->SetSamplingTimestepMultiple((int)pow(2,power));
						p_simulation->SetEndTime(end_time);
						p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.
						p_simulation->SetUpdatingTimestepMultiple(interval);
						
						SetAbsoluteMovementThreshold<2>(p_simulation, AMT);
						AddForce<2>(p_simulation, simulation_type, interactionCutoff);
						AddTracking<2>(p_simulation, (int)pow(2,power));

						try
						{
							p_simulation->Solve();
						}
						catch (Exception& e)
						{
							// If it throws then we report the error message and go to the next simulation
							WARNING("Simulation didnt run" << filename << ".");
							WARNING(e.GetMessage());
						}

						// Avoid Memory Leaks
						CleanUpNodes<2>(nodes);
						delete p_simulation;
						delete p_cell_population;
					}

					// =================================================================================================== 
					// 2D Mesh SIMULATION ================================================================================ 
					// =================================================================================================== 
					else if(simulation_type=="2d_mesh")
					{

						boost::shared_ptr<AbstractNumericalMethod<2> > numericalMethod = MakeNumericalMethod<2>(method,isAdaptive);
						
						OffLatticeSimulation<2>* p_simulation;
						std::vector<Node<2>*> nodes;
						std::vector<unsigned> realCellIndices;
						std::vector<CellPtr> cells;
						AbstractOffLatticeCellPopulation<2,2>* p_cell_population;

						HoneycombMeshGenerator generator(2, 2, 0);
						boost::shared_ptr<MutableMesh<2,2> > pMesh = generator.GetMesh();
						pMesh->Scale(compression,compression);
						realCellIndices = generator.GetCellLocationIndices();

						GenerateInitialCellsWithDistributedAges<2>(cells, realCellIndices, isStochasticCCM, ccLength, false);       
						
						p_cell_population = new MeshBasedCellPopulation<2>(*pMesh, cells);
					
						p_cell_population->SetDampingConstantNormal(dampingNormal);
						p_cell_population->AddCellWriter<CellAgesWriter>();
						
						static_cast<MeshBasedCellPopulation<2>*>(p_cell_population)->SetWriteVtkAsPoints(true);
						static_cast<MeshBasedCellPopulation<2>*>(p_cell_population)->SetBoundVoronoiTessellation(true);
						
        				p_cell_population->AddPopulationWriter<VoronoiDataWriter>();
						
						p_simulation = new OffLatticeSimulation<2>(*p_cell_population, false, true);
						p_simulation->SetNumericalMethod(numericalMethod);
						p_simulation->SetOutputDirectory(filename.c_str());
						p_simulation->SetDt(dt);
						p_simulation->SetSamplingTimestepMultiple((int)pow(2,power));
						p_simulation->SetEndTime(end_time);
						p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.
						p_simulation->SetUpdatingTimestepMultiple(interval);

						SetAbsoluteMovementThreshold<2>(p_simulation, AMT);
						AddForce<2>(p_simulation, simulation_type, interactionCutoff);
						AddTracking<2>(p_simulation, (int)pow(2,power));

						try
						{
							p_simulation->Solve();
						}
						catch (Exception& e)
						{
							// If it throws then we report the error message and go to the next simulation
							WARNING("Simulation didnt run" << filename << ".");
							WARNING(e.GetMessage());
						}

						// Avoid Memory Leaks
						CleanUpNodes<2>(nodes);
						delete p_simulation;
						delete p_cell_population;
					}

					// =================================================================================================== 
					// 2D Vertex SIMULATION ============================================================================== 
					// =================================================================================================== 
					else if(simulation_type=="2d_vertex")
					{
						boost::shared_ptr<AbstractNumericalMethod<2> > numericalMethod = MakeNumericalMethod<2>(method,isAdaptive);
						
						OffLatticeSimulation<2>* p_simulation;
						std::vector<unsigned> realCellIndices;
						std::vector<CellPtr> cells;
						AbstractOffLatticeCellPopulation<2,2>* p_cell_population;

						 // Create a simple 2D MutableVertexMesh
						HoneycombVertexMeshGenerator generator(nCellsPerSide, nCellsPerSide);
						boost::shared_ptr< MutableVertexMesh<2,2> > pMesh = generator.GetMesh();
						pMesh->Scale(compression,compression);

						double num_cells = pMesh->GetNumElements();
						for(int n=0; n<num_cells; n++)
						{
							realCellIndices.push_back((unsigned)n);
						}

						GenerateInitialCellsWithDistributedAges<2>(cells, realCellIndices, isStochasticCCM, ccLength, false);       
						
						p_cell_population = new VertexBasedCellPopulation<2>(*pMesh, cells);
					
						p_cell_population->SetDampingConstantNormal(dampingNormal);
						p_cell_population->AddCellWriter<CellAgesWriter>();
						// We dont use a growth model as tis would need to be updated to also work on substeps.
						p_cell_population->SetDataOnAllCells("target area",1.0); 
						
						// MAKE_PTR(StableRandomDirectionVertexBasedDivisionRule<2>, p_div_rule);
						// dynamic_cast<VertexBasedCellPopulation<2>*>(p_cell_population)->SetVertexBasedDivisionRule(p_div_rule);
						
						p_simulation = new OffLatticeSimulation<2>(*p_cell_population, false, true);
						p_simulation->SetNumericalMethod(numericalMethod);
						p_simulation->SetOutputDirectory(filename.c_str());
						p_simulation->SetDt(dt);
						p_simulation->SetSamplingTimestepMultiple((int)pow(2,power));
						p_simulation->SetEndTime(end_time);
						p_simulation->SetOutputDivisionLocations(true); // to check division events are at the same time.
						p_simulation->SetUpdatingTimestepMultiple(interval);

						SetAbsoluteMovementThreshold<2>(p_simulation,AMT);
						AddForce<2>(p_simulation, simulation_type, interactionCutoff);
						AddTracking<2>(p_simulation, (int)pow(2,power));

						try
						{
							p_simulation->Solve();
						}
						catch (Exception& e)
						{
							// If it throws then we report the error message and go to the next simulation
							WARNING("Simulation didn't run" << filename << ".");
							WARNING(e.GetMessage());
						}

						delete p_simulation;
						delete p_cell_population;
					}

					OutputFileHandler output_file_handler(filename, false);
					out_stream p_stream = output_file_handler.OpenOutputFile("timing.dat");
					*p_stream <<  Timer::GetElapsedTime();
					p_stream->close();

					PRINT_VARIABLE(Timer::GetElapsedTime());

				}	
			}
		}
	};

};

#endif /*TESTNUMERICS_HPP_*/