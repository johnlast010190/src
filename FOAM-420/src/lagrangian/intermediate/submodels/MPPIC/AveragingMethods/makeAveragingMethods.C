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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/Fields/Field/Field.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "meshes/polyMesh/polyMeshTetDecomposition/polyMeshTetDecomposition.H"

#include "submodels/MPPIC/AveragingMethods/Basic/Basic.H"
#include "submodels/MPPIC/AveragingMethods/Dual/Dual.H"
#include "submodels/MPPIC/AveragingMethods/Moment/Moment.H"

// Scalar interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::scalar>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        AveragingMethod<Foam::scalar>,
        dictionary
    );
}

// Vector interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::vector>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        Foam::AveragingMethod<Foam::vector>,
        dictionary
    );
}


// Basic interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::scalar>>
    addBasicscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::vector>>
    addBasicvectorConstructorToTable_;


// Dual interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Dual<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Dual<Foam::scalar>>
    addDualscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Dual<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Dual<Foam::vector>>
    addDualvectorConstructorToTable_;


// Moment interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Moment<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Moment<Foam::scalar>>
    addMomentscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Moment<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Moment<Foam::vector>>
    addMomentvectorConstructorToTable_;


// ************************************************************************* //
