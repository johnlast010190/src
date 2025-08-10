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
    uniformityObjectiveFunctionObject

Description
    Calculates uniformity on patch.

SourceFiles
    uniformityObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef uniformityObjectiveFunctionObject_H
#define uniformityObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class uniformityObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class uniformityObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    // Protected data

        //- total objective patch area
        scalar objPatchA_;

        //- coefficient for power function in objective source
        scalar sourceCoeff_;

        //- patch edge distance for near-wall objective damping
        PtrList<scalarField> patchMask_;


private:

    // Private Member Functions

        //- calculate objective weigthing mask to decrease weight of near-wall
        Xfer<PtrList<scalarField>> createPatchMask
            (const fvMesh&, const dictionary&) const;

        //- calculate total area of objective patches
        scalar objectivePatchArea() const;

        //- calculate mean normal velocity magnitude
        scalar targetVelocity() const;

        //- Disallow default bitwise copy construct and assignment
        uniformityObjectiveFunctionObject
        (
            const uniformityObjectiveFunctionObject&
        );
        void operator=(const uniformityObjectiveFunctionObject&);


public:

    //- Runtime type information
    TypeName("uniformity");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        uniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        uniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~uniformityObjectiveFunctionObject();


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
