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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::radiation::TaylorFosterAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    TaylorFoster radiation absorption and emission coefficients for continuous
    phase. The gas absorption coeff. is increased by a fraction related to
    smoke/soot according to:
    ksoot = b1*rhom[1+bT(T-2000)]
    with rhom soot mass concentration in [kg/m3] and b1 and bT constants.

SourceFiles
    TaylorFosterAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiationTaylorFosterAbsorptionEmission_H
#define radiationTaylorFosterAbsorptionEmission_H

#include "submodels/absorptionEmissionModel/constantAbsorptionEmission/constantAbsorptionEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class constantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class TaylorFosterAbsorptionEmission
:
    public constantAbsorptionEmission
{
private:
    word rhoName_;

    word sootName_;

    word TName_;

public:

    //- Runtime type information
    TypeName("TaylorFosterAbsorptionEmission");


    // Constructors

    //- Construct from components
    TaylorFosterAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~TaylorFosterAbsorptionEmission();


    // Member Functions

    //- Absorption coefficient for continuous phase
    virtual tmp<volScalarField> aCont(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
