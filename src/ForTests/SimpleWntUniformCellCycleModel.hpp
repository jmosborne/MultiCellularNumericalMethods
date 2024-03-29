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

#ifndef SIMPLEWNTCUNIFORMELLCYCLEMODEL_HPP_
#define SIMPLEWNTCUNIFORMELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WntConcentration.hpp"

/**
 * Simple Wnt-dependent cell-cycle model.
 */
class SimpleWntUniformCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    /**
     * The minimum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 12 hours.
     */
    double mMinCellCycleDuration;
    
    /**
     * The maximum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 14 hours.
     */
    double mMaxCellCycleDuration;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mMinCellCycleDuration;
        archive & mMaxCellCycleDuration;
        archive & mWntThreshold;
    }

protected:

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntThreshold;

    /**
     * @return the Wnt level experienced by the cell.
     */
    double GetWntLevel() const;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    SimpleWntUniformCellCycleModel(const SimpleWntUniformCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractPhaseBasedCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called.
     */
    SimpleWntUniformCellCycleModel();

    /**
     * Overridden ReadyToDivideMethod to update cell types
     *
     * @return whether the cell is ready to divide (enter M phase).
     */
    virtual bool ReadyToDivide();

    /**
     * Overridden SetCellCycleDuration() method to add stochastic cell cycle times
     */
    void SetCellCycleDuration();

        /**
     * @return mMinCellCycleDuration
     */
    double GetMinCellCycleDuration();
    
    /**
     * Set mMinCellCycleDuration.
     *
     * @param minCellCycleDuration
     */
    void SetMinCellCycleDuration(double minCellCycleDuration);

    /**
     * @return mMaxCellCycleDuration
     */
    double GetMaxCellCycleDuration();

    /**
     * Set mMaxCellCycleDuration.
     *
     * @param maxCellCycleDuration
     */
    void SetMaxCellCycleDuration(double maxCellCycleDuration);

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return the new cell-cycle model
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     * @return whether cell can terminally differentiate
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * @return mWntThreshold
     */
    double GetWntThreshold() const;

    /**
     * Set mWntThreshold.
     *
     * @param wntThreshold the value of mWntThreshold
     */
    void SetWntThreshold(double wntThreshold);

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleWntUniformCellCycleModel)

#endif /*SIMPLEWNTCUNIFORMELLCYCLEMODEL_HPP_*/
