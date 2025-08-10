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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::fluidHeight

Group
    grpReportFunctionObjects

Description
    Reports maximum height of fluid interface based on direction of gravity

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    fluidHeight.C
    IOfluidHeight.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fluidHeight_H
#define functionObjects_fluidHeight_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "meshes/primitiveShapes/point/pointFieldFwd.H"
#include "containers/Lists/PackedList/PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class fluidHeight Declaration
\*---------------------------------------------------------------------------*/

class fluidHeight
:
    public fvMeshFunctionObject,
    public writeFile
{
public:

    // Public enumerations

        //- Name of VOF field, default is "alpha1"
        word vofName_;

        //- Interface threshold value, default is 0.5
        scalar tValue_;

        //- Region in which to compute fluid height
        PackedList<1> reportCells_;
        List<PackedList<1>> reportFaces_;
        List<word> regionNames_;


    // Private Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return true if the pressure field corresponds to kinematic pressure
        bool isKinematicPressure();

        //- Disallow default bitwise copy construct
        fluidHeight(const fluidHeight&);

        //- Disallow default bitwise assignment
        void operator=(const fluidHeight&);


public:

    //- Runtime type information
    TypeName("fluidHeight");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fluidHeight
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fluidHeight();


    // Member Functions

        //- Read the fluidHeight data
        virtual bool read(const dictionary&);

        //- Execute, calculate fluid height and write to screen
        virtual bool execute();

        //- Does nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
