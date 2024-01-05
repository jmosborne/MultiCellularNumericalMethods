/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "SimpleWntUniformCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "Debug.hpp"

SimpleWntUniformCellCycleModel::SimpleWntUniformCellCycleModel()
    : mMinCellCycleDuration(12),
      mMaxCellCycleDuration(14),
      mWntThreshold(0.5)
{
}

SimpleWntUniformCellCycleModel::SimpleWntUniformCellCycleModel(const SimpleWntUniformCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration),
     mWntThreshold(rModel.mWntThreshold)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration and the cell's proliferative type are
     * (re)set as soon as InitialiseDaughterCell() is called on the
     * new cell-cycle model.
     */
}

AbstractCellCycleModel* SimpleWntUniformCellCycleModel::CreateCellCycleModel()
{
    return new SimpleWntUniformCellCycleModel(*this);
}

void SimpleWntUniformCellCycleModel::SetCellCycleDuration()
{
    assert(mpCell != nullptr);

    double max_timestep = pow(2,-6);

    if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        double time = SimulationTime::Instance()->GetTime();
        double cell_id = mpCell->GetCellId();
        double step = 0.5+0.5*sin(1000*(cell_id+time));
        double duration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * step;
    
        mCellCycleDuration = max_timestep*ceil(duration/max_timestep);
        
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCellCycleDuration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

double SimpleWntUniformCellCycleModel::GetWntLevel() const
{
    assert(mpCell != nullptr);
    return mpCell->GetCellData()->GetItem("wnt_level");
}

bool SimpleWntUniformCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    // Set up under what level of Wnt stimulus a cell will divide
    if (!(mpCell->GetMutationState()->IsType<WildTypeCellMutationState>()))
    {
        NEVER_REACHED; // all cells are the same here see SimpleWntUniformCellCycleModel for mutations
    }

    double wnt_level = GetWntLevel();

    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt threshold
    if ((wnt_level >= mWntThreshold))
    {
    
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
    else
    {
        // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
        mCellCycleDuration = DBL_MAX;
    }
    
    return AbstractSimpleCellCycleModel::ReadyToDivide();
}

// void SimpleWntUniformCellCycleModel::UpdateCellCyclePhase()
// {
//     // Set up under what level of Wnt stimulus a cell will divide
//     if (!(mpCell->GetMutationState()->IsType<WildTypeCellMutationState>()))
//     {
//         NEVER_REACHED; // all cells are the same here see SimpleWntUniformCellCycleModel for mutations
//     }

//     double wnt_level = GetWntLevel();

//     // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt threshold
//     if (wnt_level >= mWntThreshold)
//     {
//         boost::shared_ptr<AbstractCellProperty> p_transit_type =
//             mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
//         mpCell->SetCellProliferativeType(p_transit_type);
//     }
//     else
//     {
//         // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
//         boost::shared_ptr<AbstractCellProperty> p_diff_type =
//             mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
//         mpCell->SetCellProliferativeType(p_diff_type);
//     }
//     AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
// }

void SimpleWntUniformCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool SimpleWntUniformCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleWntUniformCellCycleModel::GetWntThreshold() const
{
    return mWntThreshold;
}

void SimpleWntUniformCellCycleModel::SetWntThreshold(double wntThreshold)
{
    //assert(wntThreshold <= 1.0);
    //assert(wntThreshold >= 0.0);
    mWntThreshold = wntThreshold;
}

double SimpleWntUniformCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void SimpleWntUniformCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double SimpleWntUniformCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void SimpleWntUniformCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double SimpleWntUniformCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double SimpleWntUniformCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}


void SimpleWntUniformCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<WntThreshold>" << mWntThreshold << "</WntThreshold>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntUniformCellCycleModel)
