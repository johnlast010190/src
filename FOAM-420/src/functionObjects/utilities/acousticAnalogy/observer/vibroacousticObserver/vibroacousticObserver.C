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

#include "acousticAnalogy/observer/vibroacousticObserver/vibroacousticObserver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/Fields/scalarField/scalarIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vibroacousticObserver, 0);                                      \
    addToRunTimeSelectionTable(observer, vibroacousticObserver, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vibroacousticObserver::vibroacousticObserver(
                                    const word& name,
                                    const fvMesh& mesh,
                                    const dictionary& dict
                                )
:
observer(name, mesh, dict),
localObsSize_(Pstream::nProcs(),0),
curIter_(0),
writeFile_(dict.lookupOrDefault<bool>("writeFile",false))
{
    setPositions(gatherAndScatterPositions());
    pPrime(scalarField(size(),0.));
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vectorField Foam::vibroacousticObserver::gatherAndScatterPositions()
{
    fileName dirPath;

    if (Pstream::parRun())
        dirPath = mesh().time().path()/".."/"CAA"/"processor"+std::to_string(Pstream::myProcNo());
    else
        dirPath = mesh().time().path()/"CAA";

    vectorIOField obs
    (
        IOobject
        (
            "observers",
            dirPath,
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        vectorField(0,vector::zero)
    );

    localObsSize_[Pstream::myProcNo()] = obs.size();
    reduce<labelList> (localObsSize_,sumOp<labelList>());

    List<vectorField> posLoc(Pstream::nProcs(),vectorField(0,vector::zero));
    posLoc[Pstream::myProcNo()] = obs;
    Pstream::allGatherList(posLoc);

    vectorField posGlob = ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());

    //Info<< "Observer size: " << endl;
    //forAll(posLoc, obsI)
    //{
    //    Info<< obsI << " : " << posLoc[obsI].size() << endl;
    //    Info<<       " : " << localObsSize_[obsI] << endl;
    //}

    return posGlob;
}

//- Return true to write on FO file
bool Foam::vibroacousticObserver::write() const
{
    fileName dirPath;

    if (Pstream::parRun())
        dirPath = fileName("..")/".."/"CAA"/"processor"+std::to_string(Pstream::myProcNo());
    else
        dirPath = fileName("..")/"CAA";

    scalarIOField acP
    (
        IOobject
        (
            "pac_"+std::to_string(curIter_++),
            mesh().time().system(),
            dirPath,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalarField(0,0.)
    );

    acP.setSize(localObsSize_[Pstream::myProcNo()]);

    const scalarField& acPrCur = pPrime();
    label offset=0;
    for (int pI=0; pI<Pstream::nProcs(); pI++)
    {
        label newSize(0);

        if (pI==Pstream::myProcNo())
        {
            newSize = acP.size();
            forAll(acP, pI)
            {
                acP[pI] = acPrCur[pI+offset];
            }

            if (newSize>0) acP.write();
        }

        reduce<label>(newSize, sumOp<label>());
        offset += newSize;
    }

    return writeFile_;
}

// ************************************************************************* //
