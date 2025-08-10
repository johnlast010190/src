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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2020-2021 Esi Ltd.

Class
    Foam::decoupledEulerian::GunnKinzer

Description
    Drag model based on Gunn & Kinzer table.

SourceFiles
    GunnKinzer.C

\*---------------------------------------------------------------------------*/

#ifndef GunnKinzer_H
#define GunnKinzer_H

#include "disperseEulerian/forceModels/dragModels/dragModel/dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace decoupledEulerian
{

class phase;

/*---------------------------------------------------------------------------*\
                           Class GunnKinzer Declaration
\*---------------------------------------------------------------------------*/

class GunnKinzer
:
    public dragModel
{
public:

    //- Use temperature-dependent material properties (y/n)
    Switch autoMaterialProp_;

    //- Runtime type information
    TypeName("GunnKinzer");


    // Constructors

        //- Construct from a dictionary and a phase pair
        GunnKinzer
        (
            const dictionary& dict,
            const phase& phase,
            const bool registerObject
        );


    //- Destructor
    virtual ~GunnKinzer();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- Return the phase-intensive drag coefficient K
        //  used in the momentum equations (alphad*K*(Uc-Ud))
        virtual tmp<volScalarField> K() const;

        //- Return temperature-dependent material properties
        const dimensionedScalar mua(const dimensionedScalar& T) const;
        const dimensionedScalar rhoa(const dimensionedScalar& T) const;
        const dimensionedScalar rhop(const dimensionedScalar& T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decoupledEulerian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
