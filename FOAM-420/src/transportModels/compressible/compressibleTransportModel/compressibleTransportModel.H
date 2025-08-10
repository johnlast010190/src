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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::compressibleTransportModel

Description
    Base-class for all transport models used by the compressible turbulence
    models.

SourceFiles
    compressibleTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTransportModel_H
#define compressibleTransportModel_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class compressibleTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleTransportModel
{
    // Private Member Functions

        //- Disallow copy construct
        compressibleTransportModel(const compressibleTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const compressibleTransportModel&);


public:

    //- Runtime type information
    TypeName("compressibleTransportModel");


    // Constructors

        //- Construct from components
        compressibleTransportModel
        ();


    //- Destructor
    virtual ~compressibleTransportModel();


    // Member Functions

        //- Return the dynamic laminar viscosity
        virtual tmp<volScalarField> mu() const = 0;

        //- Return the dynamic laminar viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const = 0;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
