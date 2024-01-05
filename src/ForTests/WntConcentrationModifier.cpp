/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "WntConcentrationModifier.hpp"

template<unsigned DIM>
WntConcentrationModifier<DIM>::WntConcentrationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mCryptLength(DOUBLE_UNSET),
      mLengthSet(false)
{
}

template<unsigned DIM>
WntConcentrationModifier<DIM>::~WntConcentrationModifier()
{
}

template<unsigned DIM>
void WntConcentrationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void WntConcentrationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}


template<unsigned DIM>
void WntConcentrationModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();


    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        assert(mLengthSet);

        // Calculate height to get wnt level
        double height =  rCellPopulation.GetLocationOfCellCentre(*cell_iter)[DIM-1];

        double wnt_level = -1.0; // Test this is changed before leaving method.

        
        if ((height > 0) && (height < GetCryptLength()))
        {
            wnt_level = 1.0 - height/GetCryptLength();
        }
        else if (height <=0)
        {
            wnt_level = 1.0;
        }
        else
        {
            wnt_level = 0.0;
        }

        assert(wnt_level >= 0.0);

        // Store the cell's volume in CellData
        cell_iter->GetCellData()->SetItem("wnt_level", wnt_level);
    }
}



template<unsigned DIM>
double WntConcentrationModifier<DIM>::GetCryptLength()
{
    return mCryptLength;
}

template<unsigned DIM>
void WntConcentrationModifier<DIM>::SetCryptLength(double cryptLength)
{
    assert(cryptLength > 0.0);
    mCryptLength = cryptLength;
    mLengthSet = true;
}

template<unsigned DIM>
void WntConcentrationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
//TODO output PARAMETERS

    // Now call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class WntConcentrationModifier<1>;
template class WntConcentrationModifier<2>;
template class WntConcentrationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WntConcentrationModifier)

