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

#include "SimpleCellCentrePositionTracker.hpp"


//Constructor, initialises sampling intervals and sets output file to null
template<unsigned DIM>
SimpleCellCentrePositionTracker<DIM>::SimpleCellCentrePositionTracker(int samplingInterval, int samplingCellInterval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      outputFile(NULL),
      interval(samplingInterval),
      cellInterval(samplingCellInterval),
      outputDir(std::string())
{}

//Empty destructor 
template<unsigned DIM>
SimpleCellCentrePositionTracker<DIM>::~SimpleCellCentrePositionTracker(){}


//Getter methods
template<unsigned DIM>
int SimpleCellCentrePositionTracker<DIM>::GetInterval() const
{
  return interval;
};
template<unsigned DIM>
int SimpleCellCentrePositionTracker<DIM>::GetCellInterval() const
{
  return cellInterval;
};

template<unsigned DIM>
std::string SimpleCellCentrePositionTracker<DIM>::GetOutputDirectoryFull(){
  return outputDir;
};


//Open an output file PositionData.txt in the simulation directory
template<unsigned DIM>
void SimpleCellCentrePositionTracker<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  outputDir = rOutputFileHandler.GetOutputDirectoryFullPath();
  outputFile = rOutputFileHandler.OpenOutputFile("PositionData.txt");

  OutputPositions(rCellPopulation);
}


//At each timestep, if the time is a sampling time, loop through all cells and output data to file.
template<unsigned DIM>
void SimpleCellCentrePositionTracker<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  if(SimulationTime::Instance()->GetTimeStepsElapsed() % GetInterval() == 0)
  {
    OutputPositions(rCellPopulation);
  }

  //If the simulation is finished, close the output file.
  if(SimulationTime::Instance()->IsFinished())
  {
    outputFile->close();
  }

}

//At each timestep, if the time is a sampling time, loop through all cells and output data to file.
template<unsigned DIM>
void SimpleCellCentrePositionTracker<DIM>::OutputPositions(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
  // if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr)
  // {
    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // Get the index, radius and damping constant of this node
        unsigned id = node_iter->GetIndex();
        
        if(id % GetCellInterval() == 0)
        { 
            c_vector<double, DIM> loc = node_iter->rGetLocation();
            *outputFile << SimulationTime::Instance()->GetTime() << "\t" << id; 
            for(unsigned i=0; i<DIM; i++)
            {
                *outputFile << std::setprecision(20) << "\t" << loc[i]; 
            }
            *outputFile << "\n";
        }
      
    }
  // }
  // else
  // {

  //   for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
  //   cell_iter != rCellPopulation.End(); ++cell_iter)
  //   { 

  //     unsigned id = cell_iter->GetCellId(); 

  //     if(id % GetCellInterval() == 0){ 

  //       c_vector<double, DIM> loc = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
  //       *outputFile << SimulationTime::Instance()->GetTime() << "\t" << id; 
  //       for(int i=0; i<DIM; i++){
  //         *outputFile << std::setprecision(20) << "\t" << loc[i]; 
  //       }
  //       *outputFile << "\n";
  //     }
  //   }

    outputFile->flush();
  // }
}




//Output this class's parameters to a log file
template<unsigned DIM>
void SimpleCellCentrePositionTracker<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SamplePositionDataEveryXTimesteps>" << GetInterval() << "</SamplePositionDataEveryXTimesteps>\n";
  *rParamsFile << "\t\t\t<RecordPositionEveryXCells>" << GetCellInterval() << "</RecordPositionEveryXCells>\n";
  
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SimpleCellCentrePositionTracker<1>;
template class SimpleCellCentrePositionTracker<2>;
template class SimpleCellCentrePositionTracker<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleCellCentrePositionTracker)
