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
    (c) 2020 Esi Ltd.

Class
    Foam::coordinateFrameRegistry

Description

\*---------------------------------------------------------------------------*/

#ifndef coordinateFrameRegistry_H
#define coordinateFrameRegistry_H

#include "primitives/strings/lists/wordList.H"
#include "primitives/ints/lists/labelList.H"
#include "fvMesh/fvMesh.H"

namespace Foam
{

class coordinateFrameRegistry
{

protected:

    const fvMesh& mesh_;

    //- List of regIOobjects name that use this frame
    mutable DynamicList<word> registeredObjectsName_;

    //- List of regIOobjects type that use this frame
    mutable DynamicList<word> registeredObjectsType_;

    //- List of patches attached to this frame
    mutable DynamicList<label> registeredPatches_;

    mutable bool MRFFacesNeedUpdate_;

public:


    //- Constructors
    coordinateFrameRegistry(const fvMesh& mesh);

    //- Destructor

    virtual ~coordinateFrameRegistry() = default;


    const DynamicList<word>& registeredNames() const
    {
        return registeredObjectsName_;
    }

    virtual void attachObject
    (
        const word& objectName,
        const word& objectType
    ) const;

    virtual bool attachPatch(const label& patchID) const;

    virtual bool isAttachToMRF(const label& patchID) const;

    virtual bool isPatchAttachedToFrame(const label& patchID) const;

    virtual void updateFaces() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
