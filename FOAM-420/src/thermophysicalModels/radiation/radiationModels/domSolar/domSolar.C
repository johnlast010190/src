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
    (c) 2010-2016 Esi Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "radiationModels/domSolar/domSolar.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/Vector/lists/vectorList.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "submodels/boundaryRadiationProperties/boundaryRadiationProperties.H"
#include "primitives/strings/wordRes/wordRes.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "primitives/strings/lists/stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(domSolar, 0);
        addToRadiationRunTimeSelectionTables(domSolar);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::domSolar::initialise(const dictionary& dict)
{
    //find and mark solar load source patches
    patchSourceCoeffs_.setSize(mesh_.boundaryMesh().size());
    forAll(patchSourceCoeffs_, patchI)
    {
        patchSourceCoeffs_.set
        (
            patchI,
            new scalarField(mesh_.boundaryMesh()[patchI].size(), 0.0)
        );
    }

    List<wordRe> sourcePatches(dict.lookup("sourcePatches"));
    wordRes spm(sourcePatches);
    const fvBoundaryMesh& bm = mesh_.boundary();

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(mesh_);

    forAll(bm, patchI)
    {
        if (spm.match(bm[patchI].name()))
        {
            if
            (
                !boundaryRadiation.radBoundaryProperties()[patchI].empty()
            )
            {
                patchSourceCoeffs_[patchI]
                    = boundaryRadiation.transmissivity(patchI);
            }
            else
            {
                patchSourceCoeffs_[patchI] = 1.0;
            }
        }
    }

    // populate directional source fields
    PtrList<entry> sources
    (
        (dict.lookupEntryPtr("sources",false,false))->stream()
    );

    IenvPtr_.setSize(sources.size());

    // count sources
    label sourcesSize = 0;
    forAllIter(PtrList<entry>, sources, iter)
    {
        if (iter().isDict())
        {
            sourcesSize++;
        }
    }
    patchLocalSourceCoeffs_.setSize(sourcesSize);

    label index = 0;
    forAllIter(PtrList<entry>, sources, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const dictionary& envRadDict = iter().dict();

        const word& name = iter().keyword();

        // use local source patches if specified
        if (envRadDict.found("localSourcePatches"))
        {
            // init
            patchLocalSourceCoeffs_.set
            (
                index,
                new FieldField<Field, scalar>(mesh_.boundaryMesh().size())
            );
            forAll(patchLocalSourceCoeffs_[index], patchI)
            {
                patchLocalSourceCoeffs_[index].set
                (
                    patchI,
                    new scalarField(mesh_.boundaryMesh()[patchI].size(), 0.0)
                );
            }

            // local source patches for i-source (index)
            List<wordRe> localSourcePatches(envRadDict.lookup("localSourcePatches"));
            wordRes lspm(localSourcePatches);

            forAll(bm, patchI)
            {
                if (lspm.match(bm[patchI].name()))
                {
                    if
                    (
                        !boundaryRadiation.radBoundaryProperties()[patchI].empty()
                    )
                    {
                        patchLocalSourceCoeffs_[index][patchI]
                            = boundaryRadiation.transmissivity(patchI);
                    }
                    else
                    {
                        patchLocalSourceCoeffs_[index][patchI] = 1.0;
                    }
                }
            }

            IenvPtr_.set
            (
                index,
                new radiantField
                (
                    word("I"+name),
                    mesh_,
                    envRadDict,
                    patchLocalSourceCoeffs_[index]
                )
            );
            index ++;

        }
        else
        {
            IenvPtr_.set
            (
                index++,
                new radiantField
                (
                    word("I"+name),
                    mesh_,
                    envRadDict,
                    patchSourceCoeffs_
                )
            );
        }
    } // end source iter

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::domSolar::domSolar(const volScalarField& T)
:
    radiationModel(typeName, T),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(true)
{
    initialise(coeffs_);
}


Foam::radiation::domSolar::domSolar
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(true)
{
    initialise(coeffs_);
}


Foam::radiation::domSolar::domSolar
(
    const dictionary& dict,
    const volScalarField& T,
    const word radWallFieldName
)
:
    radiationModel("none", T),
    IenvPtr_(),
    Qenv_
    (
        IOobject
        (
            radWallFieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    patchSourceCoeffs_(0),
    patchLocalSourceCoeffs_(),
    primary_(false)
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::domSolar::~domSolar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::domSolar::calculate()
{
    // reset environment heat flux to zero
    forAll(Qenv_.boundaryField(), bI)
    {
        const fvPatch& p = mesh_.boundary()[bI];

        if (!p.coupled())
        {
            Qenv_.boundaryFieldRef()[bI] = 0;
        }
    }

    // update radiant fields
    forAll(IenvPtr_, iI)
    {
        IenvPtr_[iI].calculate(Qenv_);
    }

    if (primary_) //apply absorptivity to scale incident radiation flux
    {
        const boundaryRadiationProperties& boundaryRadiation =
            boundaryRadiationProperties::New(mesh_);

        forAll(Qenv_.boundaryField(), patchID)
        {
            if
            (
                !boundaryRadiation.radBoundaryProperties()[patchID].empty()
            )
            {
                //only single band support
                Qenv_.boundaryFieldRef()[patchID]
                    *= boundaryRadiation.absorptivity(patchID);
            }
        }
    }
}


bool Foam::radiation::domSolar::read()
{
    initialise(coeffs_);

    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}



Foam::tmp<Foam::volScalarField> Foam::radiation::domSolar::Rp() const
{
    //this model has no volumetric contribution
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::domSolar::Ru() const
{
    //this model has no volumetric contribution
    return tmp<DimensionedField<Foam::scalar, Foam::volMesh>>
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Ru", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

// ************************************************************************* //
