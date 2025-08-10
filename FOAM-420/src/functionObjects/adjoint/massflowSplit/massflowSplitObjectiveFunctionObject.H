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
    massflowSplitObjectiveFunctionObject

Description
    Function object to calculate the proportional mass flow split of outlets
    with respect to the total inlet mass flow.

SourceFiles
    massflowSplitObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef massflowSplitObjectiveFunctionObject_H
#define massflowSplitObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class massflowSplitObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class massflowSplitObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- mass flow split fraction
        scalarList massflowFraction_;

        //- Write the header of the output file in postProcessing
        virtual void writeFileHeader(Ostream&);
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        massflowSplitObjectiveFunctionObject
        (
            const massflowSplitObjectiveFunctionObject&
        );
        void operator=(const massflowSplitObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("massflowSplit");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        massflowSplitObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        massflowSplitObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~massflowSplitObjectiveFunctionObject();


    // Member Functions

        //- read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate cost function value
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
