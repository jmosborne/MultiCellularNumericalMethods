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

#ifndef ADAMSMOULTONNUMERICALMETHOD_HPP_
#define ADAMSMOULTONNUMERICALMETHOD_HPP_

#include "AbstractNumericalMethod.hpp"
#include "ModifiedSimplePetscNonlinearSolver.hpp"


/**
* Implements the Adams Moulton 2 time stepping method. Very slow right now, since it uses a
* numerical Jacobian. Do not recommend.
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AdamsMoultonNumericalMethod : public AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> >(*this);
    }

    /**
    * Size of the forward Euler step applied to the initial positions, to help with convergence
    */
    double mFeStepSize;

    /**
    * A PETSc nonlinear solver, for solving the AM2 update equation
    */
    ModifiedSimplePetscNonlinearSolver* mpNonlinearSolver;

    /**
    * PETSc nonlinear solver tolerance
    */
    double mSolverTolerance;

    /**
    * A helper variable that allows PETSc solver related functions to access the current time step size dt
    */
    double mSolverDt;

public:	
	
    /**
     * Constructor
     */
	AdamsMoultonNumericalMethod();

    /**
     * Destructor
     */
	virtual ~AdamsMoultonNumericalMethod();

    /**
     * Updates all node positions
     *
     * @param dt Time step size
     */
    virtual void UpdateAllNodePositions(double dt);

    /**
     * Sets the size of the initial forward Euler step applied to help with convergence
     *
     * @param stepSize The new step size
    */
	void SetFeStepSize(double stepSize);

    /**
     * Get the size of the initial forward Euler step
     *
     * @return mFeStepSize
     */
    const double& GetFeStepSize() const;

    /**
     * Get the PETSc nonlinear solver tolerance
     *
     * @return mSolverTolerance
     */
    const double& GetSolverTolerance() const;

    /**
     * Set the PETSc nonlinear solver tolerance
     *
     * @param solverTolerance The new tolerance
     */
    void SetSolverTolerance(double solverTolerance);

    /**
     * Computes the residual for a particular guess at the new node positions
     *
     * @param currentGuess A PETSc Vec containing a guess at the new node positions
     * @param residualVector The amount by which the guessed node positions fail to satisfy the AM2 update equation
     */
    void ADAMSMOULTONComputeResidual(const Vec currentGuess, Vec residualVector);

    /*
     * Outputs any additional numerical method details to the parameters file
     *
     * @param rParamsFile Reference to the parameter output filestream
     */
    virtual void OutputNumericalMethodParameters(out_stream& rParamsFile);

};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode ADAMSMOULTON_ComputeResidual(SNES snes, Vec currentGuess, Vec residualVector, void* pContext)
{
    AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>* pStepper = (AdamsMoultonNumericalMethod<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pStepper->ADAMSMOULTONComputeResidual(currentGuess, residualVector);

    return 0;
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AdamsMoultonNumericalMethod)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an AdamsMoultonNumericalMethod.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const double stepSize = t->GetFeStepSize();
    ar << stepSize;

    const double tol = t->GetSolverTolerance();
    ar << tol;
}

/**
 * De-serialize constructor parameters and initialise a AdamsMoultonNumericalMethod.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double stepSize;
    ar >> stepSize;

    double tol;
    ar >> tol;

    // Invoke constructor to initialise instance
    ::new(t)AdamsMoultonNumericalMethod<ELEMENT_DIM,SPACE_DIM>();
    t->SetFeStepSize(stepSize);
    t->SetSolverTolerance(tol);
}
}
} // namespace

#endif /*ADAMSMOULTONNUMERICALMETHOD_HPP_*/
