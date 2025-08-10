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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "FieldNeighbourValueConstraint.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::FieldNeighbourValueConstraint<Type>::FieldNeighbourValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::FieldNeighbourValueConstraint<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        const dictionary& fixesCellsDict = coeffs_.subDict("fieldNames");

        fieldNames_.setSize(fixesCellsDict.size());

        label i = 0;
        forAllConstIter(dictionary, fixesCellsDict, iter)
        {
            fieldNames_[i] = iter().keyword();
            i++;
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::fv::FieldNeighbourValueConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInformation
        << "FieldNeighbourValueConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    // compute values by interpolation from neigghbour cells
    List<Type> fieldValues(cells_.size(), pTraits<Type>::zero);
    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    forAll(cells_, idx)
    {
        const labelList& nbrs = psi.mesh().cellCells()[cells_[idx]];
        forAll(nbrs, cI)
        {
            fieldValues[idx] += psi[nbrs[cI]]/nbrs.size();
        }
    }

    eqn.setValues(cells_, fieldValues);
}


// ************************************************************************* //
