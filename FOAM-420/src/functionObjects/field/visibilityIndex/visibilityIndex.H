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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::functionObjects::visibilityIndex

Group
    grpFieldFunctionObjects

Description
    Calculates and writes the visibilityIndex as a volScalarField
    for simulation with smoke transport.

    VisibilityIndex = min(C/(km*smoke_density),Lmax)

    How to specify the functionObject in functions{}:

    visibility
    {
        type            visibilityIndex;
    rho             rhoEff; //rho for compressible solvers
    C               3;  //Optional
    Lmax            30; //Optional
    kappam          8700; //Optional
        executeControl  writeTime;
        writeControl    writeTime;
        libs            ( "libfieldFunctionObjects.so" );
        regions         1 ( region0 );
    }

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    visibilityIndex.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_visibilityIndex_H
#define functionObjects_visibilityIndex_H

#include "fieldExpression/fieldExpression.H"
#include "fields/volFields/volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class visibilityIndex Declaration
\*---------------------------------------------------------------------------*/
class visibilityIndex
:
    public fieldExpression
{
    // Private Member Functions

        //- Name of density field
        word rhoName_;

        // mass specific extinction coefficient [m2/kg]
        scalar kappam_;

        // C constant (object viewed across smoke)
        scalar C_;

        // maximum distance (visibility clipper)
        scalar Lmax_;

        //- Calculate the visibility index field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("visibilityIndex");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        visibilityIndex
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~visibilityIndex();

    // Member Functions

    virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
