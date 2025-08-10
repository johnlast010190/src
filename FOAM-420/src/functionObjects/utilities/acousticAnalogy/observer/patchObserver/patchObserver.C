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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/observer/patchObserver/patchObserver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchObserver, 0);                                      \
    addToRunTimeSelectionTable(observer, patchObserver, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchObserver::readOrCreatePressure()
{
    if (!mesh().foundObject<volScalarField>("acPressure"))
    {
        autoPtr<volScalarField> pfieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "acPressure",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("acPr", dimensionSet(0,2,-2,0,0,0,0), 0.0)
            )
        );

        pfieldPtr.ptr()->store();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchObserver::patchObserver
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
observer(name, mesh, dict),
patchName_(dict.lookup("name")),
patchID_(mesh.boundaryMesh().findPatchID(patchName_)),
writeFile_(dict.lookupOrDefault<bool>("writeFile",false)),
writeAsField_(dict.lookupOrDefault<bool>("writeAsField",true))
{
    setPositions(gatherAndScatterPositions());
    pPrime(scalarField(size(),0.));

    if (writeAsField_)
    {
        readOrCreatePressure();
    }
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * *//

Foam::vectorField Foam::patchObserver::gatherAndScatterPositions()
{
    List<vectorField> posLoc(Pstream::nProcs());
    posLoc[Pstream::myProcNo()] = mesh().boundary()[patchID_].Cf();

    Pstream::allGatherList(posLoc);

    //tmp<vectorField> posGlobTmp(new vectorField);
    //vectorField& posGlob = posGlobTmp();

    vectorField posGlob = ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());
    //tmp<vectorField> posGlobTmp(posGlob);

    return posGlob;
    //return posGlobTmp;
}


bool Foam::patchObserver::write() const
{
    if (writeAsField_)
    {
        const volScalarField& acPc =
            mesh().lookupObject<volScalarField>("acPressure");
        volScalarField& acP = const_cast<volScalarField&>(acPc);
        const scalarField& acPrCur = pPrime();
        label offset = 0;

        for (int pI=0; pI<Pstream::nProcs(); pI++)
        {
            label newSize(0);
            if (pI==Pstream::myProcNo())
            {
                scalarField& acPrB = acP.boundaryFieldRef()[patchID_];
                newSize = acPrB.size();
                forAll(acPrB, fI)
                {
                    acPrB[fI] = acPrCur[fI+offset];
                }
            }
            reduce<label>(newSize, sumOp<label>());
            offset += newSize;
        }
    }

    return writeFile_;
}

// ************************************************************************* //
