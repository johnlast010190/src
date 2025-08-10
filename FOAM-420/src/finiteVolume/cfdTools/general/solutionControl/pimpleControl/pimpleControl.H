/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2021 Esi Ltd

Class
    Foam::pimpleControl

Description
    PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop.

    May also be used to for PISO-based algorithms as PISO controls are a
    sub-set of PIMPLE controls.

\*---------------------------------------------------------------------------*/

#ifndef pimpleControl_H
#define pimpleControl_H

#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

//- Declare that pimpleControl will be used
#define PIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pimpleControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleControl
:
    public solutionControl
{
    // Private member functions

        //- Disallow default bitwise copy construct
        pimpleControl(const pimpleControl&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleControl&);


protected:

    // Protected data

        // Solution controls

            //- Flag to indicate whether to solve for the flow
            bool solveFlow_;

            //- Maximum number of PIMPLE correctors
            label nCorrPIMPLE_;

            //- Maximum number of PISO correctors
            label nCorrPISO_;

            //- Current PISO corrector
            label corrPISO_;

            //- Flag to indicate whether to update density in SIMPLE
            //  rather than PISO mode
            bool SIMPLErho_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Converged flag
            bool converged_;


    // Protected Member Functions

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        ) override;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("PIMPLE");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        pimpleControl(fvMesh& mesh, const word& dictName="PIMPLE");

        //- Construct from object registry mesh and the name of control sub-dictionary
        pimpleControl(fvMesh& mesh, const objectRegistry& obr, const word& dictName="PIMPLE");

    //- Destructor
    virtual ~pimpleControl();


    // Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControls() override;

        // Access

            //- Maximum number of outer correctors
            virtual label nOuterCorr() const override
            {
                return nCorrPIMPLE();
            }

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Current PISO corrector index
            inline label corrPISO() const;

            //- Flag to indicate whether to update density in SIMPLE
            //  rather than PISO mode
            inline bool SIMPLErho() const;

            inline bool turbOnFinalIterOnly() const
            {
                return turbOnFinalIterOnly_;
            }


        // Solution control

            //- PIMPLE loop
            virtual bool loop() override;

            //- Pressure corrector loop control
            inline bool correct();

            //- Return true to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Return true for first PIMPLE (outer) iteration
            inline bool firstIter() const;

            //- Return true fore final PIMPLE (outer) iteration
            inline bool finalIter() const;

            //- Return true for final inner iteration
            inline bool finalInnerIter() const;

            //- Return true to solve for flow
            inline bool solveFlow() const;

            //- Return true to solve for turbulence
            inline bool turbCorr();

	    //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/solutionControl/pimpleControl/pimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
