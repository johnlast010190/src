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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Class
    Foam::IDDESDelta2

Description
    IDDESDelta2 used by the IDDES (improved low Re Spalart-Allmaras DES model)
    The min and max delta are calculated using the double distance of the min or
    max from the face centre to the cell centre.

SourceFiles
    IDDESDelta2.C

\*---------------------------------------------------------------------------*/

#ifndef IDDESDelta2Delta_H
#define IDDESDelta2Delta_H

#include "LES/LESdeltas/maxFaceDelta/maxFaceDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class IDDESDelta2 Declaration
\*---------------------------------------------------------------------------*/

class IDDESDelta2
:
    public LESdelta
{
    // Private data

        maxFaceDelta hmax_;

        scalar Cw_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        IDDESDelta2(const IDDESDelta2&);
        void operator=(const IDDESDelta2&);

        //- Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("IDDESDelta2");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        IDDESDelta2
        (
            const word& name,
            const turbulenceModel& turbulence_,
            const dictionary&
        );


    // Destructor

        ~IDDESDelta2()
        {}


    // Member Functions

        //- Read the LESdelta dictionary
        void read(const dictionary&);

        //- Return the hmax delta field
        const volScalarField& hmax() const
        {
            return hmax_;
        }

        // Correct values
        void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
