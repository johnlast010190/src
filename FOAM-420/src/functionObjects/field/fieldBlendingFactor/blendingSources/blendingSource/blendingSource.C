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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldBlendingFactor/blendingSources/blendingSource/blendingSource.H"
#include "finiteVolume/fvc/fvcAverage.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(blendingSource, 0);
defineRunTimeSelectionTable(blendingSource, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<surfaceScalarField> blendingSource::convertValueToBlendingFactor
(
    const surfaceScalarField& valueField,
    bool verbose
)
{
    if (verbose) Info<< ": " << type() << " ";
    label nfaces = 0;

    tmp<surfaceScalarField> tblendingField
    (
        new surfaceScalarField
        (
            IOobject
            (
                valueField.name() + "BlendingFactor",
                valueField.instance(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0",dimless, 0.0)
        )
    );
    surfaceScalarField& blendingField = tblendingField.ref();

    const scalarField& ivf = valueField;
    scalarField& ibf = blendingField;
    const labelUList& owner = mesh_.owner();

    forAll(owner, facei)
    {
        ibf[facei] = ramp_()(ivf[facei]);
        if (ibf[facei] < 1) nfaces++;
    }


    surfaceScalarField::Boundary& blendingbf = blendingField.boundaryFieldRef();
    forAll(blendingbf, patchi)
    {
        if (blendingbf[patchi].coupled())
        {
            bool singularPatch = true;
            const polyPatchList& patches = mesh_.boundaryMesh();

            if
            (
                isA<processorPolyPatch>(patches[patchi])
            )
            {
                const processorPolyPatch& procPatch
                    = refCast<const processorPolyPatch>(patches[patchi]);

                //only add lower numbered proc face to avoid duplication
                if (procPatch.myProcNo() > procPatch.neighbProcNo())
                {
                    singularPatch = false;
                }
            }

            forAll(blendingbf[patchi], facei)
            {
                scalar& bbF = blendingbf[patchi][facei];
                bbF = ramp_()(valueField.boundaryField()[patchi][facei]);
                if (bbF < 1 && singularPatch)
                {
                    nfaces++;
                }
            }
        }
    }

    if (verbose) reduce(nfaces, sumOp<label>());
    if (verbose) Info<< nfaces << flush;

    return tblendingField;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blendingSource::blendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ramp_(blendingFunction::New(dict)),
    writeBlendingField_(dict.lookupOrDefault<Switch>("write", false)),
    obr_(obr),
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<blendingSource> blendingSource::New
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& blendingSourceDict
)
{
    word sourceType = blendingSourceDict.lookup("type");

    const auto ctor = ctorTableLookup("blending source type", dictionaryConstructorTable_(), sourceType);
    return autoPtr<blendingSource>(ctor(obr, mesh, blendingSourceDict));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> blendingSource::blendingFactor()
{
    return convertValueToBlendingFactor(sourceField());
}

void  blendingSource::write()
{

    tmp<surfaceScalarField> sf(sourceField());

    tmp<surfaceScalarField> bf(convertValueToBlendingFactor(sf(), false));

    if (writeBlendingField_)
    {
        volScalarField source("average"+sf->name(), fvc::average(sf())());
        source.write();
        volScalarField blend("average"+bf->name(),fvc::average(bf())());
        blend.write();
    }
}

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
