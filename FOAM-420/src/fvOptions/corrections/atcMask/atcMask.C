/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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
    © 2018 ESI Ltd.

Author
    2018. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

\*---------------------------------------------------------------------------*/

#include "atcMask.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atcMask, 0);
    addToRunTimeSelectionTable
    (
        option,
        atcMask,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atcMask::atcMask
(
    const word& name,
    const word& modelType,
    const dictionary& optionDict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, optionDict, obr),
    fieldName_(coeffs_.lookup("fieldName")),
    sensor_
    (
        sensor<vector>::New
        (
            mesh_,
            optionDict.subDict("sensor")
        )
    ),
    min_(readScalar(coeffs_.lookup("min"))),
    max_(readScalar(coeffs_.lookup("max"))),
    writeMask_(coeffs_.lookupOrDefault<Switch>("writeMask", false))
{
    fieldNames_.setSize(1, fieldName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::atcMask::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("min") >> min_;
        coeffs_.lookup("max") >> max_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::atcMask::correct(volScalarField& atcMask)
{
    volScalarField mask
    (
        IOobject
        (
            "atcMask_" + name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mask", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    );

    scalar alpha = 1/(max_ - min_);
    scalar beta = 0.5*(1 - alpha*(min_ + max_));

    tmp<scalarField> sensorField = sensor_().valueField();

    sensorField = max(sensorField, min_);
    sensorField = min(sensorField, max_);

    mask.primitiveFieldRef() = alpha*sensorField + beta;
    mask.primitiveFieldRef() = min(mask.primitiveFieldRef(), scalar(1.0));
    mask.correctBoundaryConditions();

    if (mesh_.time().outputTime())
    {
        sensor_().writeField();

        if (writeMask_ == true)
        {
            mask.write();
        }
    }

    atcMask *= scalar(1.0) - mask;
}


// ************************************************************************* //
