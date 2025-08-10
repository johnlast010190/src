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

#include "uniformity/uniformityObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "patchBoundaryDist/patchBoundaryDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(uniformityObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        uniformityObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * *Private Member Functions * * * * * * * * * * * //

Xfer<PtrList<scalarField>>
Foam::functionObjects::uniformityObjectiveFunctionObject::
createPatchMask
(
    const fvMesh& mesh,
    const dictionary& dict
) const
{
    PtrList<scalarField> patchMask(mesh.boundary().size());

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, pI)
    {
        patchMask.set(pI, new scalarField(patches[pI].size(), 1.0));
    }

    if (dict.found("patchMask"))
    {
        scalar smoothDist
            = readScalar(dict.lookup("patchMask"));

        scalar C1 = -1/sqr(smoothDist);
        scalar C2 = 2/smoothDist;

        forAll(objectivePatch_, pI)
        {
            if
            (
                (
                    patches[pI].type() == "outlet" ||
                    patches[pI].patch().physicalType() == "outlet" ||
                    patches[pI].type() == "cyclicAMI"
                ) &&
                (
                    objectivePatch_[pI]
                )
            )
            {
                patchBoundaryDist pbd(mesh, pI);
                patchMask.set(pI, new scalarField(pbd*(C1*pbd+C2)));
                forAll(patchMask[pI], fI)
                {
                    if (pbd[fI] > smoothDist)
                    {
                        patchMask[pI][fI] = 1;
                    }
                }
            }
        }
    }

    return patchMask.xfer();
}


scalar
Foam::functionObjects::uniformityObjectiveFunctionObject::
objectivePatchArea() const
{
    scalar outletArea = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet"
             || p.patch().physicalType() == "outlet"
             || p.type() == "cyclicAMI"
            )
         && (
                objectivePatch_[patchI]
            )
        )
        {
            outletArea += sum(p.magSf());
        }
    }

    reduce(outletArea, sumOp<scalar>());

    return outletArea;
}


scalar
Foam::functionObjects::uniformityObjectiveFunctionObject::
targetVelocity() const
{
    scalar Vd = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet" ||
                p.patch().physicalType() == "outlet" ||
                p.type() == "cyclicAMI"
            ) &&
            (
                objectivePatch_[patchI]
            )
        )
        {
            Vd += sum(phi().boundaryField()[patchI]);
        }
    }

    reduce(Vd, sumOp<scalar>());

    Vd /= objPatchA_;

    //to shore things up if the initial outlet velocity is negative or zero
    Vd = sign(Vd)*(max(mag(Vd), SMALL));

    return Vd;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::uniformityObjectiveFunctionObject::
uniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    objPatchA_(objectivePatchArea()),
    sourceCoeff_(objectiveDict.lookupOrDefault<scalar>("powerScale", 1)),
    patchMask_(createPatchMask(mesh_, objectiveDict))
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::uniformityObjectiveFunctionObject::
uniformityObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    uniformityObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::uniformityObjectiveFunctionObject::
~uniformityObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::uniformityObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::uniformityObjectiveFunctionObject::execute()
{
    scalar Vd = targetVelocity();
    objectiveValue_ = 0;
    scalar uniformity = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet" ||
                p.patch().physicalType() == "outlet" ||
                p.type() == "cyclicAMI"
            ) &&
            (
                objectivePatch_[patchI]
            )
        )
        {
            scalarField Udev
            (
                mag(phi().boundaryField()[patchI]
                /p.magSf()/Vd - 1.0)
            );

            uniformity += sum
            (
                p.magSf()*Udev
            );

            scalarField Vdev2( magSqr(U().boundaryField()[patchI] - Vd*p.nf()) );

            objectiveValue_
                += sum(patchMask_[patchI]*0.5*p.magSf()*(sourceCoeff_*Vdev2
                + (1-sourceCoeff_)*sqr(Vdev2)));
        }
    }

    reduce(uniformity, sumOp<scalar>());
    uniformity = 1 - 0.5*uniformity/objPatchA_;

    reduce(objectiveValue_, sumOp<scalar>());

    Info<< type() << " " << name() << " execute:" << nl
        << "Outlet velocity SD = " << objectiveValue_
        << " [m/s]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::uniformityObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
