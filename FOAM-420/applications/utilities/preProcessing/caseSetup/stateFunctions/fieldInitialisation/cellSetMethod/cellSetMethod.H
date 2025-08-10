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
    (c) 2016 Esi Ltd.

Class
    Foam::fieldInits::cellSet

Description
    cellSet initialization method.

SourceFiles
    cellSet.C

\*--------------------------------------------------------------------*/

#ifndef cellSetMethod_H
#define cellSetMethod_H

#include "fieldInitialisation/fieldInit/fieldInit.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{

/*--------------------------------------------------------------------*\
                              Class cellSetMethod
\*--------------------------------------------------------------------*/

class cellSetMethod
:
    public fieldInit
{
    // Private data

    // Private Functions


        //- other value type set functions
        template<class Type>
        bool setCellValue
        (
            const PtrList<entry>& setDicts,
            const PtrList<cellSet>& cellSets
        ) const;

public:

    //- Runtime type information
    TypeName("cellSet");


    // Constructors

        //- construct from initMethod
        cellSetMethod
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    //- Destructor
    virtual ~cellSetMethod()
    {}


    // Member Functions

        //- Correct the fieldInit model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- scalar value set function
template<>
bool cellSetMethod::setCellValue<scalar>
(
    const PtrList<entry>& setDicts,
    const PtrList<cellSet>& cellSets
) const;


} // End namespace fieldInits
} // End namespace Foam


#ifdef NoRepository
    #include "cellSetMethodTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
