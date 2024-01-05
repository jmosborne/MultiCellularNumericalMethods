/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "CryptBoundaryForce.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
CryptBoundaryForce<DIM>::CryptBoundaryForce()
    : AbstractForce<DIM>(),
      mRepulsionStrength(100) // default to 100
{
}

template<unsigned DIM>
CryptBoundaryForce<DIM>::~CryptBoundaryForce()
{
}

template<unsigned DIM>
void CryptBoundaryForce<DIM>::SetRepulsionStrength(double newValue)
{
    assert(newValue > 0.0);
    mRepulsionStrength = newValue;
}


template<unsigned DIM>
void CryptBoundaryForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    double mRadiusOfSphere = 1.5;
    c_vector<double,DIM> mCentreOfSphere = mRadiusOfSphere*unit_vector<double>(DIM,DIM-1);

    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {   
        // Find the radial distance between this cell and the surface of the sphere
        const c_vector<double, DIM>& node_location = node_iter->rGetLocation();


        c_vector<double,DIM> centre_to_cell =  node_location - mCentreOfSphere;

        if (centre_to_cell(DIM-1)>0 )
        { 
          centre_to_cell(DIM-1)=0;
        }        

        double radius = norm_2(centre_to_cell);
        
        assert(radius != 0.0); //Can't project the centre to anywhere sensible

        double force_magnitude = mRepulsionStrength*(mRadiusOfSphere-radius);

        c_vector<double, DIM> force_contribution = force_magnitude * centre_to_cell/radius;

        node_iter->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void CryptBoundaryForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RepulsionStrength>" << mRepulsionStrength << "</RepulsionStrength> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CryptBoundaryForce<1>;
template class CryptBoundaryForce<2>;
template class CryptBoundaryForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptBoundaryForce)
