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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "forceCoeffs/forceCoeffs.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());

        if (nBin_ > 1)
        {
            CmBinFilePtr_ = createFile("CmBin");
            writeBinHeader("Moment coefficient bins", CmBinFilePtr_());
            CdBinFilePtr_ = createFile("CdBin");
            writeBinHeader("Drag coefficient bins", CdBinFilePtr_());
            ClBinFilePtr_ = createFile("ClBin");
            writeBinHeader("Lift coefficient bins", ClBinFilePtr_());
        }
    }
}


void Foam::functionObjects::forceCoeffs::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, "Force coefficients");
    writeHeaderValue(os, "liftDir", liftDir_);
    writeHeaderValue(os, "dragDir", dragDir_);
    writeHeaderValue(os, "pitchAxis", pitchAxis_);
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", csys().origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeDelimited(os, "Cm");
    writeDelimited(os, "Cd");
    writeDelimited(os, "Cl");
    writeDelimited(os, "Cl(f)");
    writeDelimited(os, "Cl(r)");
    os  << endl;
}


void Foam::functionObjects::forceCoeffs::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointi)
    {
        binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
        os << tab << binPoints[pointi].x();
    }
    os << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointi)
    {
        os << tab << binPoints[pointi].y();
    }
    os << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointi)
    {
        os << tab << binPoints[pointi].z();
    }
    os << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; j++)
    {
        word jn(Foam::name(j) + ':');
        writeDelimited(os, jn + "total");
        writeDelimited(os, jn + "pressure");
        writeDelimited(os, jn + "viscous");

        if (porosity_)
        {
            writeDelimited(os, jn + "porous");
        }
    }

    os  << endl;
}


void Foam::functionObjects::forceCoeffs::writeIntegratedData
(
    const word& title,
    const List<Field<scalar>>& coeff
) const
{
    if (!log)
    {
        return;
    }

    const scalar pressure = sum(coeff[0]);
    const scalar viscous = sum(coeff[1]);
    const scalar porous = sum(coeff[2]);
    const scalar total = pressure + viscous + porous;

    Info<< "        " << title << "       : " << total << token::TAB
        << "("
        << "pressure: " << pressure << token::TAB
        << "viscous: " << viscous;

    if (porosity_)
    {
        Info<< token::TAB << "porous: " << porous;
    }

    Info<< ")" << endl;
}


void Foam::functionObjects::forceCoeffs::writeBinData
(
    const List<Field<scalar>>& coeffs,
    Ostream& os
) const
{
    writeTime(os);

    for (label bini = 0; bini < nBin_; bini++)
    {
        const scalar total =
            coeffs[0][bini] + coeffs[1][bini] + coeffs[2][bini];

        os  << tab << total << tab << coeffs[0][bini]
            << tab << coeffs[1][bini];

        if (porosity_)
        {
            os  << tab << coeffs[2][bini];
        }
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict),
    liftDir_(Zero),
    dragDir_(Zero),
    pitchAxis_(Zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0),
    liftDragFormat_(false),
    coeffFilePtr_(),
    CmBinFilePtr_(),
    CdBinFilePtr_(),
    ClBinFilePtr_()
{
    forceCoeffs::read(dict);
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    // Free stream velocity
    vector Uinf(Zero);
    if (dict.found("magUInf"))
    {
        magUInf_ = dict.lookup<scalar>("magUInf");
    }
    else if (dict.found("Uinf"))
    {
        Uinf = dict.lookup("Uinf");
        magUInf_ = mag(Uinf);

        if (magUInf_ < SMALL)
        {
            FatalErrorInFunction
                << "mag(Uinf) must be greater than zero."
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Either Uinf or magUInf have to be defined."
            << exit(FatalError);
    }

    // Lift direction
    const wordList liftKeys({"liftDirection", "liftDir"});
    liftDir_ = normalised(dict.lookupOrDefault<vector>(liftKeys, Zero));
    if (mag(liftDir_) < VSMALL)
    {
        FatalErrorInFunction
            << "Either " << liftKeys[0] << " or " << liftKeys[1]
            << " must be defined and magnitude bigger than zero."
            << exit(FatalError);
    }

    // Drag direction
    const wordList dragKeys({"dragDirection", "dragDir"});
    dragDir_ = normalised(dict.lookupOrDefault<vector>(dragKeys, Zero));
    if (!dict.found(dragKeys[0]) && !dict.found(dragKeys[1]))
    {
        dragDir_ = normalised(Uinf);
    }
    if (mag(dragDir_) < VSMALL)
    {
        FatalErrorInFunction
            << "Drag direction must have magnitude bigger than zero."
            << exit(FatalError);
    }

    // Pitch axis
    if (dict.found("pitchAxis"))
    {
        pitchAxis_ = dict.lookup<vector>("pitchAxis");
    }
    else
    {
        pitchAxis_ = -(dragDir_ ^ liftDir_);
        if (mag(pitchAxis_) < VSMALL)
        {
            WarningInFunction
                << "Drag and lift directions are the same." << endl;
        }
    }

    pitchAxis_ = normalised(pitchAxis_);

    // If case is compressible we must read rhoInf (store in rhoRef_) to
    // calculate the reference dynamic pressure
    // - note: for incompressible, rhoRef_ is already initialised
    if (rhoName_ != "rhoInf")
    {
        dict.lookup("rhoInf") >> rhoRef_;
    }

    // Reference length
    lRef_ =
        dict.lookupOrDefault<scalar>
        (
            wordList({"referenceLength", "lRef"}),
            Zero
        );
    if (lRef_ < VSMALL)
    {
        FatalErrorInFunction
            << "Either referenceLength or lRef has to be defined "
            << "and be grater than zero."
            << exit(FatalError);
    }

    // Reference area
    Aref_ =
        dict.lookupOrDefault<scalar>
        (
            wordList({"referenceArea", "Aref"}),
            Zero
        );
    if (Aref_ < VSMALL)
    {
        FatalErrorInFunction
            << "Either referenceArea or Aref has to be defined "
            << "and be grater than zero."
            << exit(FatalError);
    }

    if (writeFields_)
    {
        if (!mesh_.foundObject<volVectorField>(fieldName("forceCoeff")))
        {
            auto forceCoeffPtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName("forceCoeff"),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("0", dimless, Zero)
                )
            );

            mesh_.objectRegistry::store(forceCoeffPtr);
        }

        if (!mesh_.foundObject<volVectorField>(fieldName("momentCoeff")))
        {
            volVectorField* momentCoeffPtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName("momentCoeff"),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("0", dimless, Zero)
                )
            );

            mesh_.objectRegistry::store(momentCoeffPtr);
        }
    }

    return true;
}


bool Foam::functionObjects::forceCoeffs::execute()
{
    forces::calcForcesMoment();

    createFiles();

    const scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

    // Storage for pressure, viscous and porous contributions to coeffs
    List<Field<scalar>> momentCoeffs(3);
    List<Field<scalar>> dragCoeffs(3);
    List<Field<scalar>> liftCoeffs(3);
    forAll(liftCoeffs, i)
    {
        momentCoeffs[i].setSize(nBin_);
        dragCoeffs[i].setSize(nBin_);
        liftCoeffs[i].setSize(nBin_);
    }

    // Calculate coefficients
    scalar CmTot = 0;
    scalar CdTot = 0;
    scalar ClTot = 0;
    forAll(liftCoeffs, i)
    {
        momentCoeffs[i] = (moment_[i] & pitchAxis())/(Aref_*pDyn*lRef_);
        dragCoeffs[i] = (force_[i] & dragDir())/(Aref_*pDyn);
        liftCoeffs[i] = (force_[i] & liftDir())/(Aref_*pDyn);

        CmTot += sum(momentCoeffs[i]);
        CdTot += sum(dragCoeffs[i]);
        ClTot += sum(liftCoeffs[i]);
    }

    scalar ClfTot = ClTot/2.0 + CmTot;
    scalar ClrTot = ClTot/2.0 - CmTot;

    Log << type() << " " << name() << " execute:" << nl
        << "    Coefficients" << nl;

    writeIntegratedData("Cm", momentCoeffs);
    writeIntegratedData("Cd", dragCoeffs);
    writeIntegratedData("Cl", liftCoeffs);

    Log << "        Cl(f)    : " << ClfTot << nl
        << "        Cl(r)    : " << ClrTot << nl
        << endl;

    if (writeToFile())
    {
        writeTime(coeffFilePtr_());
        coeffFilePtr_()
            << tab << CmTot << tab  << CdTot
            << tab << ClTot << tab << ClfTot << tab << ClrTot << endl;


        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                forAll(liftCoeffs, i)
                {
                    for (label bini = 1; bini < nBin_; bini++)
                    {
                        liftCoeffs[i][bini] += liftCoeffs[i][bini-1];
                        dragCoeffs[i][bini] += dragCoeffs[i][bini-1];
                        momentCoeffs[i][bini] += momentCoeffs[i][bini-1];
                    }
                }
            }

            writeBinData(dragCoeffs, CdBinFilePtr_());
            writeBinData(liftCoeffs, ClBinFilePtr_());
            writeBinData(momentCoeffs, CmBinFilePtr_());
        }
    }

    // Write state/results information
    {
        setResult("Cm", CmTot);
        setResult("Cd", CdTot);
        setResult("Cl", ClTot);
        setResult("Cl(f)", ClfTot);
        setResult("Cl(r)", ClrTot);
    }

    if (writeFields_)
    {
        const volVectorField& force =
            lookupObject<volVectorField>(fieldName("force"));

        const volVectorField& moment =
            lookupObject<volVectorField>(fieldName("moment"));

        volVectorField& forceCoeff =
            lookupObjectRef<volVectorField>(fieldName("forceCoeff"));

        volVectorField& momentCoeff =
            lookupObjectRef<volVectorField>(fieldName("momentCoeff"));

        dimensionedScalar f0("f0", dimForce, Aref_*pDyn);
        dimensionedScalar m0("m0", dimForce*dimLength, Aref_*lRef_*pDyn);

        forceCoeff.forceAssign(force/f0);
        momentCoeff.forceAssign(moment/m0);
    }

    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{
    if (writeFields_)
    {
        const volVectorField& forceCoeff =
            lookupObject<volVectorField>(fieldName("forceCoeff"));

        const volVectorField& momentCoeff =
            lookupObject<volVectorField>(fieldName("momentCoeff"));

        forceCoeff.write();
        momentCoeff.write();
    }

    return true;
}


// ************************************************************************* //
