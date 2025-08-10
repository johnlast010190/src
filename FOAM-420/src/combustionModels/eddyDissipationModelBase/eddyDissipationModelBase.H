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
    (c) 2016 OpenCFD Ltd

Class
    Foam::combustionModels::eddyDissipationModelBase

Description
    Standard Eddy Dissipation Model based on the assumption that the
    reaction rates are controlled by the turbulece time length scale.
    The chemical reaction rate is governed by the large-eddy mixing scale,
    (ie k/epsilon). Combustion takes place whenever the turbulence is present.


SourceFiles
    eddyDissipationModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationModelBase_H
#define eddyDissipationModelBase_H

#include "singleStepCombustion/singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationModelBase Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationModelBase
:
    public singleStepCombustion<CombThermoType, ThermoType>
{
    // Private data

        //- EDM constant
        scalar CEDC_;


    // Private Member Functions

        //- Disallow copy construct
        eddyDissipationModelBase(const eddyDissipationModelBase&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationModelBase&);


public:

    // Constructors

        //- Construct from components
        eddyDissipationModelBase
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~eddyDissipationModelBase();


    // Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Correct combustion rate
        void correct();

        //- Calculate time scale
        virtual tmp<volScalarField> timeScale() = 0;


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "eddyDissipationModelBase/eddyDissipationModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
