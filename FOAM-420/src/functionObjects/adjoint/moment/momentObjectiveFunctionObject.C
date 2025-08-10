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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "moment/momentObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "sources/derived/MRFSource/MRFSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(momentObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        momentObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::momentObjectiveFunctionObject::
momentObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    zoneNames_
    (
        objectiveDict.lookupOrDefault<wordList>
        (
            "volumeZonesToOptimize",
            wordList(0)
        )
    )
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::momentObjectiveFunctionObject::
momentObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    momentObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::momentObjectiveFunctionObject::
~momentObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::momentObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    zoneNames_ = dict.lookupOrDefault<wordList>
        (
            "volumeZonesToOptimize",
            wordList(0)
        );

    return true;
}


bool Foam::functionObjects::momentObjectiveFunctionObject::execute()
{
    tmp<volScalarField> pp = objectiveFunctionObject::P();
    const volScalarField::Boundary& Pw = pp().boundaryField();

    tmp<volSymmTensorField> tau = devRhoReff();

    const volSymmTensorField::Boundary& tauw = tau().boundaryField();
    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    vector force = vector::zero;
    scalar power = 0.0;
    fv::options& fvOpt(fv::options::New(this->mesh_));
    for (int fvI=0; fvI < (fvOpt.optionList::size()); fvI++)
    {
        fv::option& fvOptI = fvOpt.optionList::operator[](fvI);
        if (fvOptI.isActive()) //TODO: should also check if it belongs in the region
        {
            if (fvOptI.type() == "MRFSource")
            {
                fv::MRFSource& mrfS = refCast<fv::MRFSource>(fvOptI);

                label zoneNamesN = zoneNames_.size();
                bool zoneOpt(false);
                if (zoneNamesN == 0)
                {
                    zoneOpt=true;
                }
                forAll(zoneNames_,zNI)
                {
                    if (zoneNames_[zNI] == mrfS.cellSetName())
                    {
                        zoneOpt = true;
                    }
                }
                if (zoneOpt)
                {
                    const vector& axisPatch = mrfS.axis();
                    const scalar& omegaPatch = mag(mrfS.Omega());
                    const vector& centerOfRotation = mrfS.CofR();
                    const labelListList& incFaces =
                        mrfS.getFrameSourceFaces().includedFaces();
                    vector omegaV = axisPatch*omegaPatch;
                    forAll(mesh_.boundary(), pI)
                    {
                        tmp<vectorField> normal = mesh_.boundary()[pI].nf();
                        vectorField pf( Sfb[pI]*Pw[pI] );
                        vectorField vf
                        (
                            Sfb[pI]^((Sfb[pI]/mesh_.magSf().boundaryField()[pI])
                            & tauw[pI])
                        );
                        force += sum(pf + vf);
                        forAll(incFaces[pI],i)
                        {
                            label pF = incFaces[pI][i];
                            vector Cp_n = mesh_.C().boundaryField()[pI][pF];
                            scalar Cp_naxis = Cp_n & axisPatch;
                            vector localVector
                            (
                                (
                                    Cp_n - (Cp_naxis*axisPatch)
                                  - centerOfRotation
                                )
                               ^omegaV
                            );
                            power +=  ((pf[pF]+vf[pF]) & localVector);
                        }
                    }
                }
            }
        }
    }

    reduce(force, sumOp<vector>());
    reduce(power, sumOp<scalar>());
    objectiveValue_ = power;

    Info<< type() << " " << name() << " execute:" << nl
        << "Moment = " << objectiveValue_ << " [Nm]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::momentObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
