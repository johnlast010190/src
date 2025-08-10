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
    (c) 2017-2023 Esi Ltd.

Class
    Foam::patchDistanceFunction

Description
    Class for providing disance coordinate on patch

Usage
    \table
        Property     | Description             | Required    | Default value
        distanceType | distance function       | yes         |
        origin       | position vector         | no          | (0 0 0)
        axis         | direction vector        | no          | (0 0 1)
    \endtable

    Example of the input specification - wallDistance:
    \verbatim
    <dictionary>
    {
        distanceType    wallDistance;
        origin          (0 0 0); //might be used in windProfileVelocity
        axis            (0 0 1); //might be used in windProfileVelocity
    }

    Example of the input specification - vector:
    \verbatim
    <dictionary>
    {
        distanceType    vector;
        origin          (0 0 0);
        axis            (0 0 1);
    }

    Example of the input specification - x|y|z:
    \verbatim
    <dictionary>
    {
        distanceType    z;
        origin          (0 0 0);
    }

    Example of the input specification - cylindrical radius:
    \verbatim
    <dictionary>
    {
        distanceType    radius;
        origin          (0 0 0);
        axis            (0 0 1);
    }


    \endverbatim


SourceFiles
    patchDistanceFunction.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistanceFunction_H
#define patchDistanceFunction_H

#include "db/typeInfo/typeInfo.H"
#include "fields/Fields/Field/Field.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class patchDistanceFunction Declaration
\*---------------------------------------------------------------------------*/

class patchDistanceFunction
:
    public tmp<patchDistanceFunction>::refCount
{

public:

    // Public data types
        enum distTypes
        {
            dtX,
            dtY,
            dtZ,
            dtWall,
            dtPoint,
            dtVector,
            dtRadius
        };

protected:

        static const NamedEnum<distTypes, 7> distTypeNames_;

private:

    // Private Data

        //- Reference to mesh
        const fvPatch& patch_;

        //- Distance calculation method
        distTypes distanceType_;

        //- Copy of the dict
        dictionary dictCoeffs_;

        //- Point, an optional component used with some distance methods
        mutable point origin_;

        //- Direction, optional component used with the vectorDistance method
        mutable vector axis_;

        //- Distance
        mutable autoPtr<scalarField> value_;

        //- Direction
        mutable autoPtr<vectorField> dir_;

        //- Face size in direction
        mutable autoPtr<scalarField> delta_;

        //- Optional reference frame
        mutable coordinateFrame* coorFramePtr_;

public:

    //- Runtime type information
        TypeName("patchDistanceFunction");


    // Constructors

        //- Null constructor
        patchDistanceFunction(const fvPatch&);

        //- Construct from components
        patchDistanceFunction(const fvPatch&, const dictionary&);

        //- Copy construct
        patchDistanceFunction(const patchDistanceFunction&);

        //- Construct and return a clone
        //- Construct and return a clone
        tmp<patchDistanceFunction> clone() const
        {
            return tmp<patchDistanceFunction>
            (
                new patchDistanceFunction(*this)
            );
        }


    //- Destructor
    virtual ~patchDistanceFunction()
    {}


    // Member Functions

        //- Prepare frame if the boundary is defined in the frame
        coordinateFrame* loadFrame(const dictionary& dict);

        //- Local direction mapping
        vector localToGlobalDir(vector dir) const;

        //- return distance value
        const scalarField& value() const;

        //- return direction
        const vectorField& dir() const;

        //- return face delta in direction
        const scalarField& delta() const;

        //- update data if mesh is changing
        bool update();

        inline const vector& axis() const
        {
            return axis_;
        }

        inline const vector& origin() const
        {
            return origin_;
        }


        //- Write
        void write(Ostream&os, bool writeFrame = false) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
