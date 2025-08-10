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
    (c) 2016 OpenCFD ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::maxDeltaxyz

Description
    Delta calculated by taking the maximum distance between the cell centre
    and any face centre.  For a regular hex cell, the computed delta will
    equate to half of the cell width; accordingly, the deltaCoeff model
    coefficient should be set to 2 for this case.

SourceFiles
    maxDeltaxyz.C

\*---------------------------------------------------------------------------*/

#ifndef maxDeltaxyzDelta_H
#define maxDeltaxyzDelta_H

#include "LES/LESdeltas/LESdelta/LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class maxDeltaxyz Declaration
\*---------------------------------------------------------------------------*/

class maxDeltaxyz
:
    public LESdelta
{
    // Private data

        //- Model coefficient
        scalar deltaCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        maxDeltaxyz(const maxDeltaxyz&);
        void operator=(const maxDeltaxyz&);


public:

    //- Runtime type information
    TypeName("maxDeltaxyz");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        maxDeltaxyz
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~maxDeltaxyz()
    {}


    // Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
