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

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrameRegistry/coordinateFrameRegistry.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateFrameRegistry::coordinateFrameRegistry(const fvMesh& mesh)
:
    mesh_(mesh),
    registeredObjectsName_(),
    registeredObjectsType_(),
    registeredPatches_(),
    MRFFacesNeedUpdate_(false)
{}


void Foam::coordinateFrameRegistry::attachObject
(
    const word& objectName,
    const word& objectType
) const
{
    registeredObjectsName_.append(objectName);
    registeredObjectsType_.append(objectType);
    MRFFacesNeedUpdate_ = true;
    updateFaces();
}


bool Foam::coordinateFrameRegistry::attachPatch(const label& patchi) const
{
    bool newPatch = true;
    MRFFacesNeedUpdate_ = false;
    forAll(registeredPatches_, pI)
    {
        if (registeredPatches_[pI] == patchi)
        {
            newPatch = false;
        }
    }
    if (newPatch)
    {
        registeredPatches_.append(patchi);
        MRFFacesNeedUpdate_ = true;
    }
    updateFaces();
    return MRFFacesNeedUpdate_;
}


bool Foam::coordinateFrameRegistry::isAttachToMRF
(
    const label& patchi
) const
{
    forAll(registeredObjectsName_, objI)
    {
        const word objectNamei = registeredObjectsName_[objI];
        if (mesh_.foundObject<frameSourceFaces>(objectNamei))
        {
            const frameSourceFaces& mrf =
                mesh_.lookupObject<frameSourceFaces>(objectNamei);

            if (mrf.includedFaces()[patchi].size())
            {
                return true;
            }
        }
    }
    return false;
}


bool Foam::coordinateFrameRegistry::isPatchAttachedToFrame
(
    const label& patchi
) const
{
    forAll(registeredPatches_, pI)
    {
        if (registeredPatches_[pI] == patchi)
        {
            return true;
        }
    }
    return false;
}


void Foam::coordinateFrameRegistry::updateFaces() const
{
    if (MRFFacesNeedUpdate_)
    {
        forAll(registeredObjectsName_, objI)
        {
            word objName = registeredObjectsName_[objI];
            if (mesh_.foundObject<frameSourceFaces>(objName))
            {
                frameSourceFaces& frameSrcFaces =
                    mesh_.lookupObjectRef<frameSourceFaces>(objName);
                frameSrcFaces.updateSourceFaces(registeredPatches_);
            }
        }
    }
}


// ************************************************************************* //
