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
    Foam::combustionModels::eddyDissipationDiffusionModel

Group
    grpCombustionModels

Description
    Eddy dissipation model based on the principle of mixed is burnt.

    The time scale is determined by the minimum of the turbulence
    and diffusion time scales.


SourceFiles
    eddyDissipationDiffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationDiffusionModel_H
#define eddyDissipationDiffusionModel_H

#include "eddyDissipationModelBase/eddyDissipationModelBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                Class eddyDissipationDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class eddyDissipationDiffusionModel
:
    public eddyDissipationModelBase<CombThermoType, ThermoType>
{
    // Private data

        //- Diffussivity constant
        scalar Cd_;


    // Private memeber functions

        //- Disallow copy construct
        eddyDissipationDiffusionModel(const eddyDissipationDiffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationDiffusionModel&);


public:

    //- Runtime type information
    TypeName("EDM");


    // Constructors

        //- Construct from components
        eddyDissipationDiffusionModel
        (
            const word& modelType,
            const objectRegistry& obr,
            const word& combustionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~eddyDissipationDiffusionModel();


    // Member Functions

        //- Calculate time scale
        virtual tmp<volScalarField> timeScale();

         //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "eddyDissipationDiffusionModel/eddyDissipationDiffusionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
