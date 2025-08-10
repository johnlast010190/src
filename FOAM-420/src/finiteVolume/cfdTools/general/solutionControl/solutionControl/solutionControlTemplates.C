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
    (c) 2016-2021 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "volMesh/volMesh.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::solutionControl::storePrevIter() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    objectRegistry& obrRef = const_cast<objectRegistry&>(obr_);
    HashTable<GeoField*>
        flds(obrRef.objectRegistry::lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        const word& fName = fld.name();

        size_t prevIterField = fName.find("PrevIter");

        if
        (
            (prevIterField == word::npos)
            && (mesh_.solution().relaxField(fName)
            || (storeVars() && mesh_.solution().relaxEquation(fName)))
        )
        {
            if (debug)
            {
                Info<< algorithmName_ << ": storing previous iter for "
                    << fName << endl;
            }

            fld.storePrevIter();
        }
    }
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class Type, class ...Args>
Type& Foam::solutionControl::lookupOrCreateType
(
    fvMesh& mesh, const objectRegistry& obr, Args&&... args
)
{
    if (!obr.foundObject<Type>(solutionControl::typeName))
    {
        if (solutionControl::debug)
        {
            Pout<< "solutionControl::lookupOrCreateType(const fvMesh&, const objectRegistry&) : "
                << "constructing solution control "
                << " for region " << obr.name() << endl;
        }
        regIOobject::store(new Type(mesh, obr, std::forward<Args>(args)...));
    }

    return
        obr.lookupObjectRef<Type>
        (
            solutionControl::typeName
        );
}

// ************************************************************************* //
