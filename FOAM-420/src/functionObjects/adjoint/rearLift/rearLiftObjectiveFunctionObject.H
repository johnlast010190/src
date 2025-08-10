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
    (c) 2019 Esi Ltd.

Authors
    2019. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

Class
    Foam::functionObjects::objectiveFunctionObject::
    rearLiftObjectiveFunctionObject

Description
    Calculates the value of rear lift.

SourceFiles
    rearLiftObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef rearLiftObjectiveFunctionObject_H
#define rearLiftObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class rearLiftObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class rearLiftObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- Direction of lift
        vector forceDirection_;

        //- Direction of moment contributing to the objective
        vector momentDirection_;

        //- Rotation centre
        //  Usually placed on the ground, in the middle of the front and
        //  rear axis.
        vector rotationCentre_;

        //- Distance between the front and rear axis
        scalar axisLength_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        rearLiftObjectiveFunctionObject
        (
            const rearLiftObjectiveFunctionObject&
        );
        void operator=(const rearLiftObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("rearLift");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        rearLiftObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        rearLiftObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~rearLiftObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the value of the cost function and write it to file
        virtual bool execute();

        //- Track and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
