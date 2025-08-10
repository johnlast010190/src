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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "waveRunup/waveRunup.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(waveRunup, 0);
    addToRunTimeSelectionTable(functionObject, waveRunup, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::waveRunup::waveRunup
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    alphaName_("alpha1"),
    patchNames_(0),
    contourLevel_(0.5),
    volHeight_(-GREAT),
    patchHeights_(0)
{
    read(dict);
    writeFileHeader(file());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::waveRunup::~waveRunup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::functionObjects::waveRunup::writeFileHeader(Ostream& os) const
{
    Log << "    Logging wave heights to file: " << fileName_
        << endl;

    writeCommented(os,"Time");
    writeDelimited(os,"height-in-volume");

    forAll(patchNames_, patchI)
    {
        writeDelimited(os,"height-on-"+patchNames_[patchI]);
    }

    os << endl;
}


bool Foam::functionObjects::waveRunup::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    //calculate
    calculate();

    //Write to screen
    //net flux
    Log << "Maximum wave height"<<nl;
    Log << "  in volume:  " << volHeight_ << nl;
    forAll(patchNames_, patchI)
    {
        Log << "  on patch " << patchNames_[patchI] << ":  ";
        Log << patchHeights_[patchI] << nl;
    }
    Log << endl;

    //write to file
    file() << obr_.time().timeName() << " ";
    file() << volHeight_ << " ";

    forAll(patchNames_, patchI)
    {
        file() << patchHeights_[patchI] << " ";
    }
    file() << endl;

    return true;
}


void Foam::functionObjects::waveRunup::calculate()
{
    const volScalarField& alpha = lookupObject<volScalarField>(alphaName_);
    const uniformDimensionedVectorField& g =
        lookupObject<uniformDimensionedVectorField>("g");
    volScalarField waveHeight
    (
        IOobject
        (
            "waveHeight",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        ),
        mesh_,
        dimensionedScalar("waveHeight", dimLength, -GREAT),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(mesh_.neighbour(), faceI)
    {
        const label own = mesh_.owner()[faceI];
        const label nei = mesh_.neighbour()[faceI];
        if
        (
            (alpha[own] < contourLevel_ && alpha[nei] >= contourLevel_) ||
            (alpha[own] >= contourLevel_ && alpha[nei] < contourLevel_)
        )
        {
            scalar r = (alpha[own]-contourLevel_)
                /stabilise(alpha[own]-alpha[nei], SMALL);
            scalar h = -(r*mesh_.C()[nei]+(1-r)*mesh_.C()[own]) & g.value();
            waveHeight[own] = max(waveHeight[own],h);
            waveHeight[nei] = max(waveHeight[nei],h);
        }
    }
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const scalarField pif( alpha.boundaryField()[patchI].patchInternalField() );
            const scalarField pnf( alpha.boundaryField()[patchI].patchNeighbourField() );
            const vectorField cif( mesh_.C().boundaryField()[patchI].patchInternalField() );
            const vectorField cnf( mesh_.C().boundaryField()[patchI].patchNeighbourField() );
            const labelUList& fc = mesh_.boundary()[patchI].faceCells();
            forAll(pif, faceI)
            {
                if
                (
                    (pif[faceI] < contourLevel_ && pnf[faceI] >= contourLevel_)
                    ||
                    (pif[faceI] >= contourLevel_ && pnf[faceI] < contourLevel_)
                )
                {
                    scalar r = (pif[faceI]-contourLevel_)
                        /stabilise(pif[faceI]-pnf[faceI], SMALL);
                    scalar h = -(r*cnf[faceI]+(1-r)*cif[faceI]) & g.value();
                    waveHeight[fc[faceI]] = max(waveHeight[fc[faceI]],h);
                }
            }
        }
    }
    waveHeight.correctBoundaryConditions();

    volHeight_ = gMax(waveHeight.internalField());
    forAll(patchNames_, patchI)
    {
        patchHeights_[patchI] =
            gMax
            (
                waveHeight.boundaryField()
                [
                    mesh_.boundary().findPatchID(patchNames_[patchI])
                ]
            );
    }

}

bool Foam::functionObjects::waveRunup::write()
{
    return true;
}


bool Foam::functionObjects::waveRunup::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    //patch names
    {
        wordList tfns = dict.lookup("patches");
        patchNames_ = tfns;
        patchHeights_.setSize(tfns.size(), -GREAT);
    }

    contourLevel_ = dict.lookupOrDefault("contourLevel", 0.5);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
