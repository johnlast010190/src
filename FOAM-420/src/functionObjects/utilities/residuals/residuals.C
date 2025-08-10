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
    (c) 2015 OpenCFD Ltd.
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "residuals/residuals.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "VectorN/primitives/vector4/vector4.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(residuals, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        residuals,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::residuals::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Residuals");
    writeCommented(os, "Time");

    forAll(fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];

        writeFileHeader<scalar>(os, fieldName);
        writeFileHeader<vector>(os, fieldName);
        writeFileHeader<sphericalTensor>(os, fieldName);
        writeFileHeader<symmTensor>(os, fieldName);
        writeFileHeader<tensor>(os, fieldName);
        writeFileBlockHeader<vector4>(os, fieldName);
    }

    os << endl;
}


void Foam::functionObjects::residuals::createField(const word& fieldName)
{
    if (!writeFields_)
    {
        return;
    }

    const word residualName("initialResidual_" + fieldName);

    if (!mesh_.foundObject<IOField<scalar>>(residualName))
    {
        IOField<scalar>* fieldPtr =
            new IOField<scalar>
            (
                IOobject
                (
                    residualName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Field<scalar>(mesh_.nCells(), scalar(0))
            );

        fieldPtr->store();
    }
}


void Foam::functionObjects::residuals::writeField(const word& fieldName) const
{
    const word residualName("initialResidual_" + fieldName);

    const IOField<scalar>* residualPtr =
        mesh_.lookupObjectPtr<IOField<scalar>>(residualName);

    if (residualPtr && mesh_.time().outputTime())
    {
        volScalarField residual
        (
            IOobject
            (
                residualName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0),
            zeroGradientFvPatchField<scalar>::typeName
        );

        residual.primitiveFieldRef() = *residualPtr;
        residual.correctBoundaryConditions();

        residual.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::residuals
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    fieldSet_(),
    writeFields_(false),
    initialised_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::~residuals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::residuals::read(const dictionary& dict)
{
    Log << type() << " " << name() <<  " read:" << nl;
    fvMeshFunctionObject::read(dict);

    wordList allFields(dict.lookup("fields"));
    wordHashSet uniqueFields(allFields);
    fieldSet_ = uniqueFields.toc();

    writeFields_ = dict.lookupOrDefault("writeFields", false);

    return true;
}


bool Foam::functionObjects::residuals::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    // Note: delaying initialisation until after first iteration so that
    // we can find wildcard fields
    if (!initialised_)
    {
        writeFileHeader(file());

        if (writeFields_)
        {
            forAll(fieldSet_, fieldi)
            {
                const word& fieldName = fieldSet_[fieldi];

                initialiseField<scalar>(fieldName);
                initialiseField<vector>(fieldName);
                initialiseField<sphericalTensor>(fieldName);
                initialiseField<symmTensor>(fieldName);
                initialiseField<tensor>(fieldName);
                initialiseBlockField<vector4>(fieldName);
            }
        }

        initialised_ = true;
    }

    if (Pstream::master())
    {
        writeTime(file());
    }

    forAll(fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];

        writeResidual<scalar>(fieldName);
        writeResidual<vector>(fieldName);
        writeResidual<sphericalTensor>(fieldName);
        writeResidual<symmTensor>(fieldName);
        writeResidual<tensor>(fieldName);
        writeBlockResidual<vector4>(fieldName);
    }

    if (Pstream::master())
    {
        file() << endl;
    }

    return true;
}


bool Foam::functionObjects::residuals::write()
{
    return true;
}


// ************************************************************************* //
