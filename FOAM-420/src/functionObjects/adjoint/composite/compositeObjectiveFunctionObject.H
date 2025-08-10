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
    compositeObjectiveFunctionObject

Description
    Combines objectiveFunctionObjects to create a new objective.

    This objective is not available as a function object yet.

SourceFiles
    compositeObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef compositeObjectiveFunctionObject_H
#define compositeObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class objective;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class compositeObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class compositeObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- PtrList of runTime selectable objective functions
        PtrList<functionObject> compositeObjectives_;

        //- normalised objective weights
        scalarField normalisedCompositeWeights_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        compositeObjectiveFunctionObject
        (
            const compositeObjectiveFunctionObject&
        );
        void operator=(const compositeObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("composite");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        compositeObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        compositeObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~compositeObjectiveFunctionObject();


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
