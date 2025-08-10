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
    (c) 2011 Esi Ltd.

Class
    Foam::hybridBlendingSource

Description
    hybrid blending source based on A. Travin,M. Shur, M. Strelets, and
    P. Spalart. Physical and numerical upgrades in the detached-eddy
    simulation of complex turbulent flows. In Proceedings of the
    412th Euromech Colloquium on LES and Complex Transitional and Turbulent
    Flows, Munich, Germany, 2000.

SourceFiles
    hybridBlendingSource.C

\*---------------------------------------------------------------------------*/

#ifndef hybridBlendingSource_H
#define hybridBlendingSource_H

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingSources
{

/*---------------------------------------------------------------------------*\
                   Class hybridBlendingSource Declaration
\*---------------------------------------------------------------------------*/

class hybridBlendingSource
:
    public blendingSource
{
private:

    // Private Data

        //- model coefficients
        scalar lRef_;
        scalar Uref_;
        scalar CH1_;
        scalar CH2_;
        scalar CH3_;
        scalar CDES_;
        scalar Cmu_;
        scalar OmegaLim_;
        scalar sigmaMax_;
        scalar sigmaMin_;
        Switch approxTanh_;

    // Private Member Functions

        //- return effective kinematic viscosity
        tmp<volScalarField> nuEff();

        //- return velocity
        tmp<volVectorField> U();

        //- return delta
        tmp<volScalarField> delta();

        //- Disallow default bitwise copy construct
        hybridBlendingSource(const hybridBlendingSource&);

        //- Disallow default bitwise assignment
        void operator=(const hybridBlendingSource&);

public:

    //- Runtime type information
        TypeName("hybrid");

    // Constructors

        //- Construct from components
        hybridBlendingSource(const objectRegistry& obr, const fvMesh& mesh, const dictionary& data);

    //- Destructor
    virtual ~hybridBlendingSource(){}


    // Member Functions

        //- source field calculation
        virtual tmp<surfaceScalarField> sourceField();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingSources
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
