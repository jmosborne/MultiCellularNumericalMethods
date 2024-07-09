/*

Copyright (c) 2005-2020, University of Oxford.
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

#include "StableRandomDirectionVertexBasedDivisionRule.hpp"

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> StableRandomDirectionVertexBasedDivisionRule<SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> short_axis;
    unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
    short_axis = rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
    c_vector<double, SPACE_DIM> random_vector;
    if(SPACE_DIM==2)
    {
        random_vector(0)= short_axis(1);
        random_vector(1)= -short_axis(0);
    }
    else
    {
        NEVER_REACHED; // Not implemented for 1D or 3D
    }

    //Look for intersected edges and use the midpoint vector of theses as division axis
    MutableVertexMesh<SPACE_DIM,SPACE_DIM>* p_mesh  = static_cast<MutableVertexMesh<SPACE_DIM,SPACE_DIM>*>(&(rCellPopulation.rGetMesh()));

    VertexElement<SPACE_DIM,SPACE_DIM>* p_element = p_mesh->GetElement(rCellPopulation.GetLocationIndexUsingCell(pParentCell));
    
    // Get the centroid of the element
    c_vector<double, SPACE_DIM> centroid = p_mesh->GetCentroidOfElement(p_element->GetIndex());

    /*
     * Find which edges the axis of division crosses by finding any node
     * that lies on the opposite side of the axis of division to its next
     * neighbour.
     */
    unsigned num_nodes = p_element->GetNumNodes();
    std::vector<unsigned> intersecting_nodes;
    bool is_current_node_on_left = (inner_prod(p_mesh->GetVectorFromAtoB(p_element->GetNodeLocation(0), centroid), random_vector) >= 0);
    for (unsigned i=0; i<num_nodes; i++)
    {
        bool is_next_node_on_left = (inner_prod(p_mesh->GetVectorFromAtoB(p_element->GetNodeLocation((i+1)%num_nodes), centroid), random_vector) >= 0);
        if (is_current_node_on_left != is_next_node_on_left)
        {
            intersecting_nodes.push_back(i);
        }
        is_current_node_on_left = is_next_node_on_left;
    }

    // If the axis of division does not cross two edges then we cannot proceed
    if (intersecting_nodes.size() != 2)
    {
        EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    // Now make the division axis through the midpoints of these crosed edges
    

    random_vector = p_mesh->GetVectorFromAtoB(p_element->GetNodeLocation(intersecting_nodes[0]), p_element->GetNodeLocation((intersecting_nodes[1]+1)%num_nodes));
    random_vector = random_vector + p_mesh->GetVectorFromAtoB(p_element->GetNodeLocation((intersecting_nodes[0]+1)%num_nodes), p_element->GetNodeLocation(intersecting_nodes[1]));

    return random_vector;
}

// Explicit instantiation
template class StableRandomDirectionVertexBasedDivisionRule<1>;
template class StableRandomDirectionVertexBasedDivisionRule<2>;
template class StableRandomDirectionVertexBasedDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StableRandomDirectionVertexBasedDivisionRule)
