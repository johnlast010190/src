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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::randomise::calcRandomised()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (foundObject<VolFieldType>(fieldName_, false))
    {
        const VolFieldType& field = lookupObject<VolFieldType>(fieldName_);

        resultName_ = fieldName_ & "Random";

        tmp<VolFieldType> rfieldt(new VolFieldType(field));
        VolFieldType& rfield = rfieldt.ref();

        Random rand(1234567);

        forAll(field, celli)
        {
            Type rndPert;
            rand.randomise01(rndPert);
            rndPert = 2.0*rndPert - pTraits<Type>::one;
            rndPert /= mag(rndPert);
            rfield[celli] += magPerturbation_*rndPert;
        }

        return store(resultName_, rfieldt);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
