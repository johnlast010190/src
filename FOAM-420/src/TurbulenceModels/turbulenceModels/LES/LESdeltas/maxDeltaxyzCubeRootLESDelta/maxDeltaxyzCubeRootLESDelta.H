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

Class
    Foam::maxDeltaxyzCubeRootLESDelta

Description
    Maximum delta between maxDeltaxyz and cubeRootVolDelta.

SourceFiles
    maxDeltaxyzCubeRootLESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef maxDeltaxyzCubeRootLESDelta_H
#define maxDeltaxyzCubeRootLESDelta_H

#include "LES/LESdeltas/LESdelta/LESdelta.H"
#include "LES/LESdeltas/maxDeltaxyz/maxDeltaxyz.H"
#include "LES/LESdeltas/cubeRootVolDelta/cubeRootVolDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                   Class maxDeltaxyzCubeRootLESDelta Declaration
\*---------------------------------------------------------------------------*/

class maxDeltaxyzCubeRootLESDelta
:
    public LESdelta
{

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        maxDeltaxyzCubeRootLESDelta(const maxDeltaxyzCubeRootLESDelta&);
        void operator=(const maxDeltaxyzCubeRootLESDelta&);

        // Calculate the delta values
        void calcDelta();

        //- maxDeltaxyz model
        maxDeltaxyz maxDeltaxyz_;

        //- cubeRootVolDelta model
        cubeRootVolDelta cubeRootVolDelta_;


public:

    //- Runtime type information
    TypeName("maxDeltaxyzCubeRoot");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        maxDeltaxyzCubeRootLESDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~maxDeltaxyzCubeRootLESDelta()
    {}


    // Member Functions

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
