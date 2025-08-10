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
    volumeVelocityLimitObjectiveFunctionObject

Description
    Calculates the volume integral of the magnitude of the velocity in areas
    where the velocity is above some threshold value.

SourceFiles
    volumeVelocityLimitObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumeVelocityLimitObjectiveFunctionObject_H
#define volumeVelocityLimitObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumeVelocityLimitObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumeVelocityLimitObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- Threshold value for the velocity
        scalar velocityLimitValue_;

    // Protected Member Functions

        //- Add cost of cell to objective value
        void addCostOfCell(scalar magU, scalar V);

        //- Add contribution from cell in the source of the adjoint momentum
        void addCellInUaSource(vector& source, vector U, scalar V) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        volumeVelocityLimitObjectiveFunctionObject
        (
            const volumeVelocityLimitObjectiveFunctionObject&
        );
        void operator=(const volumeVelocityLimitObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("volumeVelocityLimit");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        volumeVelocityLimitObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        volumeVelocityLimitObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~volumeVelocityLimitObjectiveFunctionObject();


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
