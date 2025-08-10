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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::functionObjects::heatBalance

Group
    grpUtilitiesFunctionObjects

Description

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    heatBalance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_heatBalance_H
#define functionObjects_heatBalance_H

#include "fields/volFields/volFieldsFwd.H"
#include "cfdTools/general/include/fvCFD.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class heatBalance Declaration
\*---------------------------------------------------------------------------*/

class heatBalance
:
    public fvMeshFunctionObject,
    public writeFile
{

    // Private data

        //- on/off switch
        bool compressible_;

        // Optional names for required physical fields
        word phiName_;

        // Fields defining heat flux
        autoPtr<surfaceScalarField> heatFlux_;
        dimensionedScalar refTemp_;
        dimensionedScalar domainHeat_;
        dimensionedScalar domainHeatOld_;
        dimensionedScalar DpDt_;


    // Private functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Returns diffusivity from turbulence models
        tmp<volScalarField> diffusivity();

        //- Returns density
        tmp<volScalarField> density();

        //- Returns enthalpy
        tmp<volScalarField> enthalpy();

        //- Disallow default bitwise copy construct
        heatBalance(const heatBalance&);

        //- Disallow default bitwise assignment
        void operator=(const heatBalance&);


public:

    //- Runtime type information
    TypeName("heatBalance");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        heatBalance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~heatBalance();


    // Member Functions
        // Sample
        void calculate();
        dimensionedScalar calcInitialDomainHeat();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
