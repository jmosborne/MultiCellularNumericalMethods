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

#ifndef SIMPLECELLCENTREPOSITIONTRACKER_HPP_
#define SIMPLECELLCENTREPOSITIONTRACKER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "OutputFileHandler.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>

/**
 * A simple modifier that outputs cell IDs and the position of cell centres. Good
 * for plotting and comparing cell positions easily between simulation runs.
 * interval = number of timesteps between data outputs
 * cellInterval = what proportion of cells to track
 */
template<unsigned DIM>
class SimpleCellCentrePositionTracker : public AbstractCellBasedSimulationModifier<DIM,DIM>
{

private:
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

    // Output file stream
    out_stream outputFile;
    
    // Number of timesteps between data recordings
    int interval;

    // Spacing between which cells to track
    int cellInterval;

    // Full path to output file
    std::string outputDir;

public:


    //Constructor 
    SimpleCellCentrePositionTracker(int samplingInterval, int samplingCellInterval);

    //Destructor
    virtual ~SimpleCellCentrePositionTracker();


    //Getters
    int GetInterval() const;

    int GetCellInterval() const;

    std::string GetOutputDirectoryFull();

     void OutputPositions(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overriden UpdateAtEndOfTimeStep method
     * Specifies what to do in the simulation at the end of each timestep, in this case record data
     * if the time is appropriate.
     *
     * @param rCellPopulation reference to the cell population
     */
     void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * Overriden SetupSolve method
     * Specifies what to do in the simulation before the start of the time loop. In this case, 
     * open an output file
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
     void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);


     //Output any parameters associated with this class
     void OutputSimulationModifierParameters(out_stream& rParamsFile);

};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleCellCentrePositionTracker)

namespace boost
{
    namespace serialization
    {
        /**
        * Serialize information required to construct a SimpleCellCentrePositionTracker.
        */
        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const SimpleCellCentrePositionTracker<DIM> * t, const unsigned int file_version)
        {
            // Save data required to construct instance
            int interval = t->GetInterval();
            ar << interval;
            int cellInterval = t->GetCellInterval();
            ar << cellInterval;
        }

        /**
        * De-serialize constructor parameters and initialise a SimpleCellCentrePositionTracker.
        */
        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, SimpleCellCentrePositionTracker<DIM> * t, const unsigned int file_version)
        {
            // Retrieve data from archive required to construct new instance
            int interval;
            ar >> interval;
            int cellInterval;
            ar >> cellInterval;

            ::new(t)SimpleCellCentrePositionTracker<DIM>(interval, cellInterval);
        }
    }
} // namespace ...

#endif /*SIMPLECELLCENTREPOSITIONTRACKER_HPP_*/