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

#include "MidPointNumericalMethod.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
MidPointNumericalMethod<ELEMENT_DIM,SPACE_DIM> :: MidPointNumericalMethod():
  AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{	
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
MidPointNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~MidPointNumericalMethod()
{
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void MidPointNumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{
    if(!this->mUseUpdateNodeLocation)
    {
        // First apply Boundary conditions to fix any issues from proliferation
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations = this->SaveCurrentNodeLocations();
        this->ImposeBoundaryConditions(old_node_locations);
        old_node_locations = this->SaveCurrentNodeLocations();

        // Compute K1
        std::vector<c_vector<double, SPACE_DIM> > K1 = this->ComputeForcesIncludingDamping();

        // Update node positions by +(dt K1/2) and compute K2
        int index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * K1[index]/2.0;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);
        }
        this->ImposeBoundaryConditions(old_node_locations); // Apply any boundary conditions for half timestep 
        std::vector< c_vector<double, SPACE_DIM> > K2 = this->ComputeForcesIncludingDamping();

        // Final position update, according to the Midpoint numerical method
        index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> effectiveForce = K2[index];
            c_vector<double, SPACE_DIM> oldLocation = old_node_locations.find(&(*node_iter))->second; //Revert
            c_vector<double, SPACE_DIM> finalDisplacement =  dt * effectiveForce;

            this->DetectStepSizeExceptions(node_iter->GetIndex(), finalDisplacement, dt);

            c_vector<double, SPACE_DIM> newLocation = oldLocation + finalDisplacement;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), newLocation);

            // Ensure that each nodes holds an accurate applied force value,
            // incase it's accessed by some other class
            double damping = this->mpCellPopulation->GetDampingConstant(node_iter->GetIndex());
            node_iter->ClearAppliedForce();
            c_vector<double, SPACE_DIM> force = effectiveForce*damping;
            node_iter->AddAppliedForceContribution(force);
        }
    }
    else
    {
        // If this type of cell population does not support the new numerical methods, delegate
        // updating node positions to the population itself.
        NEVER_REACHED;
    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void MidPointNumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    // Nothing yet

    // Call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
};


///////// Explicit instantiation
template class MidPointNumericalMethod<1,1>;
template class MidPointNumericalMethod<1,2>;
template class MidPointNumericalMethod<2,2>;
template class MidPointNumericalMethod<1,3>;
template class MidPointNumericalMethod<2,3>;
template class MidPointNumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MidPointNumericalMethod)
