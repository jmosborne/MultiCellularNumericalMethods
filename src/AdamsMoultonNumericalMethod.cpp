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

#include "AdamsMoultonNumericalMethod.hpp"
#include "ReplicatableVector.hpp"
#include "PetscVecTools.hpp"
#include "ModifiedSimplePetscNonlinearSolver.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM> :: AdamsMoultonNumericalMethod():
  AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>(),
  mFeStepSize(0.01),
  mSolverTolerance(1e-8), 
  mSolverDt(DOUBLE_UNSET)
{	
    mpNonlinearSolver = new ModifiedSimplePetscNonlinearSolver();
    mpNonlinearSolver->SetTolerance(mSolverTolerance);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~AdamsMoultonNumericalMethod()
{
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{

    if(!this->mUseUpdateNodeLocation)
    {
        mSolverDt = dt;

        unsigned systemSize = this->mpCellPopulation->GetNumNodes()*SPACE_DIM;

        std::vector< c_vector<double,SPACE_DIM> > initialLocations = this->SaveCurrentLocations();
        std::vector< c_vector<double,SPACE_DIM> > initialF = this->ComputeForcesIncludingDamping();

        // Setup an initial guess consisting of the current node locations + one forward Euler step (note no bcs are applied here)
        Vec initialCondition = PetscTools::CreateAndSetVec(systemSize, 0.0);
        int index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> location = node_iter->rGetLocation();
            for(unsigned i=0; i<SPACE_DIM; i++)
            {
                PetscVecTools::SetElement(initialCondition, SPACE_DIM*index + i,  location[i] + mFeStepSize*initialF[index][i]); // Use a FE step as initial Guess
            }
        }
        
        // Call nonlinear solver
        Vec solnNextTimestep = mpNonlinearSolver->Solve( &ADAMSMOULTON_ComputeResidual<ELEMENT_DIM, SPACE_DIM>,
                                                      &SNESComputeJacobianDefault,  
                                                      initialCondition,   
                                                      UINT_MAX,          
                                                      this);              
                                                      
        // Unpack solution.
        ReplicatableVector solnNextTimestepRepl(solnNextTimestep);

        index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> oldLocation = initialLocations[index];
            c_vector<double, SPACE_DIM> solution;
            for(unsigned i=0; i<SPACE_DIM; i++)
            {
                solution[i] = solnNextTimestepRepl[SPACE_DIM * index + i];
            }

            c_vector<double, SPACE_DIM> displacement = solution - oldLocation;
            this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);

            c_vector<double, SPACE_DIM> newLocation = oldLocation + this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(oldLocation, oldLocation+displacement);
            this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);

            node_iter->ClearAppliedForce();
            double damping = this->mpCellPopulation->GetDampingConstant(node_iter->GetIndex());
            c_vector<double, SPACE_DIM> effectiveForce = (damping/dt)*displacement;
            node_iter->AddAppliedForceContribution(effectiveForce);
        }

        PetscTools::Destroy(initialCondition);
    }
    else
    {
      // If this type of cell population does not support the new numerical methods, delegate 
      // updating node positions to the population itself.
      NEVER_REACHED; 

    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>::SetFeStepSize(double stepSize)
{
    mFeStepSize = stepSize;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
const double& AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>::GetFeStepSize() const
{
    return mFeStepSize;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
const double& AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>::GetSolverTolerance() const
{
    return mSolverTolerance;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>::SetSolverTolerance(double solverTolerance)
{
    mpNonlinearSolver->SetTolerance(solverTolerance);
    mSolverTolerance = solverTolerance;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ADAMSMOULTONComputeResidual(const Vec currentGuess, Vec residualVector)
{    
    std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > current_node_locations = this->SaveCurrentNodeLocations();
    this->ImposeBoundaryConditions(current_node_locations); // Double check bcs are applied
    std::vector< c_vector<double, SPACE_DIM> > current_locations = this->SaveCurrentLocations();
    
    std::vector< c_vector<double, SPACE_DIM> > Fcurrent = this->ComputeForcesIncludingDamping();

    ReplicatableVector guessPositions(currentGuess);

    std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations = this->SaveCurrentNodeLocations();
    int index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {  
        c_vector<double, SPACE_DIM> guessLocation; 
        for(unsigned i=0; i<SPACE_DIM; i++){
            guessLocation[i] = guessPositions[SPACE_DIM * index + i];
        }
        node_iter->rGetModifiableLocation() = guessLocation;
    }

    // Get force at the guess locations
    this->ImposeBoundaryConditions(old_node_locations); // Apply any boundary conditions for half timestep 
    std::vector< c_vector<double, SPACE_DIM> > Fguess = this->ComputeForcesIncludingDamping();
    
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {     
        node_iter->rGetModifiableLocation() = current_locations[index];

        for(unsigned i=0; i<SPACE_DIM; i++){

            double residual_ith_cpt = guessPositions[SPACE_DIM * index + i] - current_locations[index][i] 
                                      - (mSolverDt/2.0) * ( Fguess[index][i] + Fcurrent[index][i] );

            PetscVecTools::SetElement(residualVector, SPACE_DIM * index + i, residual_ith_cpt);
        }
    } 
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<mFeStepSize>" << mFeStepSize << "</mFeStepSize>\n";
    *rParamsFile << "\t\t\t<SolverTolerance>" << mSolverTolerance << "</SolverTolerance>\n";

    // Call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
};

///////// Explicit instantiation
template class AdamsMoultonNumericalMethod<1,1>;
template class AdamsMoultonNumericalMethod<1,2>;
template class AdamsMoultonNumericalMethod<2,2>;
template class AdamsMoultonNumericalMethod<1,3>;
template class AdamsMoultonNumericalMethod<2,3>;
template class AdamsMoultonNumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AdamsMoultonNumericalMethod)
