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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

Class
    Foam::functionObjects::ventilationDrag

Group
    grpForcesFunctionObjects

Description
    Calculates the ventilation drag for specified wallVelocity patches, all
    belonging to the same reference frame.


Usage
    Example of function object specification:
    \verbatim
    ventilationDrag1
    {
        type               ventilationDrag;
        patches            ( wheels_RL_COUPE-rear-wheels-left wheels_RL_COUPE-rear-rim-left ); // list of patches of type wallVelocity belonging to the same rotational reference frame (separated by whitespace)
        Uinf               ( 38.89 0.0 0.0 );        // freestream velocity (equal to vehicle velocity)
        rhoInf             1.293;                    // freestream density (air density)
        referenceArea      1.79296625;
        functionObjectLibs ("libforces.so" );        // name of the library that has to be loaded to access the function object
    }
    \endverbatim


SourceFiles
    ventilationDrag.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ventilationDrag_H
#define functionObjects_ventilationDrag_H

#include "forceCoeffs/forceCoeffs.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    class forceCoeffs;

/*---------------------------------------------------------------------------*\
                           Class ventilationDrag Declaration
\*---------------------------------------------------------------------------*/

class ventilationDrag
:
    public forces
{
    // Public data
    public:

        //- Wheel velocity field
        word Uname_;

    // Protected data
    protected:


    // Private data
    private:

        // Free-stream conditions

        //- Velocity magnitude
        scalar magUInf_;


        // Reference scales

        //- Area
        scalar Aref_;


        // Reference frames

        //- Current reference frame pointer
        const coordinateFrame* coorFramePtr_;


    protected:

    // Protected data

        //- ventilationDrag file ptr
        autoPtr<OFstream> ventilationDragFilePtr_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct and assignment
        ventilationDrag(const ventilationDrag&);
        void operator=(const ventilationDrag&);


public:

    //- Runtime type information
    TypeName("ventilationDrag");

    // Constructors

        //- Construct for given objectRegistry and dictionary
        ventilationDrag
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    // Destructor

        virtual ~ventilationDrag();


    // Member Functions

        //- Read the parameters
        virtual bool read(const dictionary&) override;

        //- Execute
        virtual bool execute() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //