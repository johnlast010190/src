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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2015 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "forces/forces.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "coordinate/systems/cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);
    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::getp()
{

    const volScalarField& p = lookupObject<volScalarField>(pName_);

    if (p.dimensions() == dimensionSet(1, -1, -2, 0, 0, 0, 0))
    {
        return p - dimensionedScalar("Pref", dimPressure, pRef_);
    }
    else
    {
        dimensionedScalar dimRho
        (
            "rho",
            dimensionSet(1, -3, 0, 0, 0, 0, 0),
            rhoRef_
        );

        return dimRho*p - dimensionedScalar("Pref", dimPressure, pRef_);
    }
}


Foam::word Foam::functionObjects::forces::fieldName(const word& name) const
{
    return this->name() + "_" + name;
}


void Foam::functionObjects::forces::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !forceFilePtr_.valid())
    {
        if (combineFiles_)
        {
            forceFilePtr_ = createFile("forces");
            writeIntegratedHeader();
        }
        else
        {
            forceFilePtr_ = createFile("force");
            writeIntegratedHeader("Force", forceFilePtr_());
            momentFilePtr_ = createFile("moment");
            writeIntegratedHeader("Moment", momentFilePtr_());
        }

        if (nBin_ > 1)
        {
            if (combineFiles_)
            {
                forceBinFilePtr_ = createFile("forces_bins");
                writeBinHeader();
            }
            else
            {
                forceBinFilePtr_ = createFile("forceBin");
                writeBinHeader("Force", forceBinFilePtr_());
                momentBinFilePtr_ = createFile("momentBin");
                writeBinHeader("Moment", momentBinFilePtr_());
            }
        }
    }
}


void Foam::functionObjects::forces::writeIntegratedHeader()
{
    Ostream& os = forceFilePtr_();
    writeHeaderValue(os, "CofR", csys().origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeDelimited(os, "forces(total,pressure,viscous,porous)");
    writeDelimited(os, "moment(total,pressure,viscous,porous)");

    os  << endl;
}


void Foam::functionObjects::forces::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", csys().origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeDelimited(os, "(total_x total_y total_z)");
    writeDelimited(os, "(pressure_x pressure_y pressure_z)");
    writeDelimited(os, "(viscous_x viscous_y viscous_z)");

    if (porosity_)
    {
        writeDelimited(os, "(porous_x porous_y porous_z)");
    }

    os  << endl;
}


void Foam::functionObjects::forces::writeBinHeader()
{
    Ostream& os = forceBinFilePtr_();

    os << "# bins      : " << nBin_ << nl
       << "# start     : " << binMin_ << nl
       << "# delta     : " << binDx_ << nl
       << "# direction : " << binDir() << nl
       << "# Time";

    for (label j = 0; j < nBin_; j++)
    {
        const word jn('[' + Foam::name(j) + ']');

        os << tab
           << "forces" << jn << "(total,pressure,viscous,porous) "
           << "moment" << jn << "(total,pressure,viscous,porous)";
    }

    os << endl;
}



void Foam::functionObjects::forces::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header + " bins");
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir());

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointi)
    {
        binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir();
        os  << tab << binPoints[pointi].x();
    }
    os  << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].y();
    }
    os  << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].z();
    }
    os  << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; j++)
    {
        const word jn(Foam::name(j) + ':');
        os  << tab << jn << "(total_x total_y total_z)"
            << tab << jn << "(pressure_x pressure_y pressure_z)"
            << tab << jn << "(viscous_x viscous_y viscous_z)";

        if (porosity_)
        {
            os  << tab << jn << "(porous_x porous_y porous_z)";
        }
    }

    os << endl;
}



void Foam::functionObjects::forces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_ << " or p:" << pName_
                << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_
                << exit(FatalError);
        }
    }

    initialiseBins();

    initialised_ = true;
}


void Foam::functionObjects::forces::initialiseBins()
{
    if (nBin_ > 1)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Determine extents of patches
        binMin_ = GREAT;
        scalar binMax = -GREAT;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            const polyPatch& pp = pbm[patchi];
            scalarField d(pp.faceCentres() & binDir());
            binMin_ = min(min(d), binMin_);
            binMax = max(max(d), binMax);
        }

        // Include porosity
        if (porosity_)
        {
            const HashTable<const porosityModel*> models =
                obr_.lookupClass<porosityModel>();

            const scalarField dd(mesh_.C() & binDir());

            forAllConstIter(HashTable<const porosityModel*>, models, iter)
            {
                const porosityModel& pm = *iter();
                const labelList& cellZoneIDs = pm.cellZoneIDs();

                forAll(cellZoneIDs, i)
                {
                    label zonei = cellZoneIDs[i];
                    const cellZone& cZone = mesh_.cellZones()[zonei];

                    bool calculateForce = checkZone(cZone);

                    if (calculateForce)
                    {
                        const scalarField d(dd, cZone);
                        binMin_ = min(min(d), binMin_);
                        binMax = max(max(d), binMax);
                    }
                }
            }
        }

        reduce(binMin_, minOp<scalar>());
        reduce(binMax, maxOp<scalar>());

        // Slightly boost binMax so that region of interest is fully
        // within bounds
        binMax = 1.0001*(binMax - binMin_) + binMin_;

        binDx_ = (binMax - binMin_)/scalar(nBin_);

        // Create the bin points used for writing
        binPoints_.setSize(nBin_);
        forAll(binPoints_, i)
        {
            binPoints_[i] = (i + 0.5)*binDir()*binDx_;
        }
    }

    // Allocate storage for forces and moments
    forAll(force_, i)
    {
        force_[i].setSize(nBin_, vector::zero);
        moment_[i].setSize(nBin_, vector::zero);
    }
}


void Foam::functionObjects::forces::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;
    moment_[2] = Zero;

    if (writeFields_)
    {
        volVectorField& force =
            const_cast<volVectorField&>
            (
                lookupObject<volVectorField>(fieldName("force"))
            );

        force.forceAssign(dimensionedVector("0", force.dimensions(), Zero));

        volVectorField& moment =
            const_cast<volVectorField&>
            (
                lookupObject<volVectorField>(fieldName("moment"))
            );

        moment.forceAssign(dimensionedVector("0", moment.dimensions(), Zero));
    }
}


bool Foam::functionObjects::forces::checkZone(const cellZone& cZone)
{
    bool calculateForce = true;
    if (porosityZones_.size())
    {
        const word& cZoneName = cZone.name();
        calculateForce = false;
        forAll(porosityZones_, zI)
        {
            if (porosityZones_[zI] == cZoneName)
            {
                calculateForce = true;
            }
        }
    }
    return calculateForce;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
            lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu("nu", transportProperties.lookup("nu"));

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);
        ::abort();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::functionObjects::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::functionObjects::forces::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fT);
        moment_[2][0] += sum(Md^fP);
    }
    else
    {
        scalarField dd((d & binDir()) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), force_[0].size() - 1);

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fP[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fT[i];
            moment_[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::functionObjects::forces::addToFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        lookupObjectRef<volVectorField>(fieldName("force"));

    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT + fP;

    volVectorField& moment =
        lookupObjectRef<volVectorField>(fieldName("moment"));

    vectorField& pm = moment.boundaryFieldRef()[patchi];
    pm += Md;
}


void Foam::functionObjects::forces::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        lookupObjectRef<volVectorField>(fieldName("force"));

    volVectorField& moment =
        lookupObjectRef<volVectorField>(fieldName("moment"));

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
        moment[celli] += Md[i];
    }
}


void Foam::functionObjects::forces::writeIntegratedForceMoment()
{
    vectorField forceN(force_[0]);
    vectorField forceT(force_[1]);
    vectorField forceP(force_[2]);
    vectorField momentN(moment_[0]);
    vectorField momentT(moment_[1]);
    vectorField momentP(moment_[2]);

    if (localSystem_)
    {
        forceN = csys().localVector(force_[0]);
        forceT = csys().localVector(force_[1]);
        forceP = csys().localVector(force_[2]);
        momentN = csys().localVector(moment_[0]);
        momentT = csys().localVector(moment_[1]);
        momentP = csys().localVector(moment_[2]);
    }

    vector forcepressure = sum(forceN);
    vector forceviscous = sum(forceT);
    vector forceporous = sum(forceP);
    vector forcetotal = forcepressure + forceviscous + forceporous;

    vector momentpressure = sum(momentN);
    vector momentviscous = sum(momentT);
    vector momentporous = sum(momentP);
    vector momenttotal = momentpressure + momentviscous + momentporous;

    Log << "    forces(total,pressure,viscous,porous) = ("
        << forcetotal << ","
        << forcepressure << ","
        << forceviscous << ","
        << forceporous << ")" << nl
        << "    moment(total,pressure,viscous,porous) = ("
        << momenttotal << ","
        << momentpressure << ","
        << momentviscous << ","
        << momentporous << ")"
        << nl;

    if (writeToFile())
    {
        Ostream& os = forceFilePtr_();
        writeTime(os);
        os << tab << "("
        << forcetotal << ","
        << forcepressure << ","
        << forceviscous << ","
        << forceporous
        << ") "
        << "("
        << momenttotal << ","
        << momentpressure << ","
        << momentviscous << ","
        << momentporous
        << ")";
        os << endl;
    }
}


void Foam::functionObjects::forces::writeIntegratedForceMoment
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2,
    autoPtr<OFstream>& osPtr
) const
{
    vector pressure = sum(fm0);
    vector viscous = sum(fm1);
    vector porous = sum(fm2);
    vector total = pressure + viscous + porous;

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << total << nl
        << "        Pressure : " << pressure << nl
        << "        Viscous  : " << viscous << nl;

    if (porosity_)
    {
        Log << "        Porous   : " << porous << nl;
    }

    if (writeToFile())
    {
        Ostream& os = osPtr();

        writeTime(os);

        os  << tab << total
            << tab << pressure
            << tab << viscous;

        if (porosity_)
        {
            os  << tab << porous;
        }

        os  << endl;
    }
}


void Foam::functionObjects::forces::writeForces()
{
    Log << type() << " " << name() << " write:" << nl;

    if (combineFiles_)
    {
        writeIntegratedForceMoment();
    }
    else
    {
        if (localSystem_)
        {
            writeIntegratedForceMoment
            (
                "forces",
                csys().localVector(force_[0]),
                csys().localVector(force_[1]),
                csys().localVector(force_[2]),
                forceFilePtr_
            );

            writeIntegratedForceMoment
            (
                "moments",
                csys().localVector(moment_[0]),
                csys().localVector(moment_[1]),
                csys().localVector(moment_[2]),
                momentFilePtr_
            );
        }
        else
        {
            writeIntegratedForceMoment
            (
                "forces",
                force_[0],
                force_[1],
                force_[2],
                forceFilePtr_
            );

            writeIntegratedForceMoment
            (
                "moments",
                moment_[0],
                moment_[1],
                moment_[2],
                momentFilePtr_
            );
        }
    }

    Log << endl;
}


void Foam::functionObjects::forces::writeBinnedForceMoment()
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    List<vectorField> f(force_);
    List<vectorField> m(moment_);

    if (localSystem_)
    {
        forAll(f, compi)
        {
            f[compi] = csys().localVector(force_[compi]);
            m[compi] = csys().localVector(moment_[compi]);
        }
    }

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            forAll(f, compi)
            {
                f[compi][i] += f[compi][i - 1];
                m[compi][i] += m[compi][i - 1];
            }
        }
    }

    Ostream& os = forceBinFilePtr_();

    writeTime(os);

    forAll(f[0], i)
    {
        vector totalforce = f[0][i] + f[1][i] + f[2][i];
        vector totalmoment = m[0][i] + m[1][i] + m[2][i];

        os  << tab
            << "(" << totalforce << "," << f[0][i] << ","
            << f[1][i] << "," << f[2][i] << ") "
            << "(" << totalmoment << "," <<m[0][i] << ","
            << m[1][i] << "," << m[2][i] << ")";
    }

    os  << nl;
}


void Foam::functionObjects::forces::writeBinnedForceMoment
(
    const List<vectorField>& fm,
    autoPtr<OFstream>& osPtr
) const
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    List<vectorField> f(fm);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];
        }
    }

    Ostream& os = osPtr();

    writeTime(os);

    forAll(f[0], i)
    {
        vector total = f[0][i] + f[1][i] + f[2][i];

        os  << tab << total
            << tab << f[0][i]
            << tab << f[1][i];

        if (porosity_)
        {
            os  << tab << f[2][i];
        }
    }

    os  << nl;
}


void Foam::functionObjects::forces::writeBins()
{
    if (combineFiles_)
    {
        writeBinnedForceMoment();
    }
    else
    {
        List<vectorField> f(force_);
        List<vectorField> m(moment_);

        if (localSystem_)
        {
            forAll(f, compi)
            {
                f[compi] = csys().localVector(force_[compi]);
                m[compi] = csys().localVector(moment_[compi]);
            }
        }

        writeBinnedForceMoment(f, forceBinFilePtr_);
        writeBinnedForceMoment(m, momentBinFilePtr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool baseClass
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    coorFramePtr_(nullptr),
    definedInFrame_(false),
    localSystem_(false),
    porosity_(false),
    porosityZones_(wordList::null()),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false),
    combineFiles_(false)
{
    if (!baseClass)
    {
        forces::read(dict);
        Log << endl;
    }
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool baseClass
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(obr_, name),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    coorFramePtr_(nullptr),
    definedInFrame_(false),
    localSystem_(false),
    porosity_(false),
    porosityZones_(wordList::null()),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false),
    combineFiles_(false)
{
    if (!baseClass)
    {
        forces::read(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forces::~forces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::coordinateSystem& Foam::functionObjects::forces::csys() const
{
    if (coorFramePtr_)
    {
        return coorFramePtr_->coorSys();
    }
    return coordSys_;
}


Foam::vector Foam::functionObjects::forces::binDir() const
{
    if (coorFramePtr_ && definedInFrame_)
    {
        return coorFramePtr_->coorSys().globalVector(binDir_);
    }
    return binDir_;
}


bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    //set patches used by forces
    patchSet_.clear();

    if (dict.found("patches"))
    {
        patchSet_ =
            pbm.patchSet(wordReList(dict.lookup("patches")), false, true);
    }
    //patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional entries U and p
        const wordList pNames({"p", "pName"});
        pName_ = dict.lookupOrDefault<word>(pNames, "p");

        const wordList UNames({"U", "UName"});
        UName_ = dict.lookupOrDefault<word>(UNames, "U");

        const wordList rhoNames({"rho", "rhoName"});
        rhoName_ = dict.lookupOrDefault<word>(rhoNames, "rho");

        if (dict.found("rhoInf") && !dict.found("rho") && !dict.found("rhoName"))
        {
            rhoName_ = "rhoInf";
        }

        if (dict.found("rhoName") && dict.found("rho"))
        {
            WarningInFunction
                << "\n    Since both entries 'rho' and 'rhoName' "
                << " were found in forces function object, '"
                << rhoNames[0] << "' value from 'rho' entry "
                << "is being used for forces calculation\n" << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            dict.lookup("rhoInf") >> rhoRef_;
        }

        // Reference pressure, 0 by default
        const wordList pRefNames({"pRef", "Pref"});
        pRef_ = dict.lookupOrDefault<scalar>(pRefNames, 0.0);
    }

    if (dict.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
        definedInFrame_ =
            dict.lookupOrDefault<Switch>("definedInFrame", false);
        localSystem_ = true;
        if (dict.found("CofR") || dict.found("referencePoint"))
        {
            WarningInFunction
                << "Reference frame is switched on. "
                << "\"CofR\" and \"referencePoint\" have no effect."
                << endl;
        }
    }
    else
    {
        // Centre of rotation for moment calculations
        // specified directly, from coordinate system, or implicitly (0 0 0)
        coordSys_.clear();
        if (!dict.readIfPresent<point>("CofR", coordSys_.origin()))
        {
            if (!dict.readIfPresent<point>("referencePoint", coordSys_.origin()))
            {
                // The 'coordinateSystem' sub-dictionary is optional,
                // but enforce use of a cartesian system.
                if (dict.found(coordinateSystem::typeName_()))
                {
                    // New() for access to indirect (global) coordinate system
                    coordSys_ =
                        coordinateSystem::New
                        (
                            obr_, dict, coordinateSystem::typeName_()
                        );
                }
                else
                {
                    coordSys_ = coordSystem::cartesian(dict);
                }

                localSystem_ = true;
            }
        }
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
        porosityZones_ = dict.lookupOrDefault<wordList>
        (
            "porosityZones", wordList::null()
        );
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));

        // check for new format of liftDrag (uses "nBins")
        if (binDict.found("nBin"))
        {
            binDict.lookup("nBin") >> nBin_;

            if (nBin_ < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Number of bins (nBin) must be zero or greater"
                    << exit(FatalIOError);
            }
            else if (nBin_ == 0)
            {
                // Case of no bins equates to a single bin to collect all data
                nBin_ = 1;
            }
            else
            {
                binDict.lookup("cumulative") >> binCumulative_;
                binDict.lookup("direction") >> binDir_;
                binDir_ /= mag(binDir_);
            }
        }
    }

    writeFields_ = dict.lookupOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        if (!mesh_.foundObject<volVectorField>(fieldName("force")))
        {
            volVectorField* forcePtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName("force"),
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("0", dimForce, Zero)
                )
            );

            mesh_.objectRegistry::store(forcePtr);
        }

        if (!mesh_.foundObject<volVectorField>(fieldName("moment")))
        {
            volVectorField* momentPtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName("moment"),
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("0", dimForce*dimLength, Zero)
                )
            );

            mesh_.objectRegistry::store(momentPtr);
        }
    }

    dict.readIfPresent("combineFiles", combineFiles_);

    return true;
}


void Foam::functionObjects::forces::calcForcesMoment()
{
    if (coorFramePtr_ && definedInFrame_)
    {
        initialiseBins();
    }

    initialise();

    resetFields();

    if (directForceDensity_)
    {
        const volVectorField& fD = lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - csys().origin()
            );

            scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            // Porous force
            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }
    else
    {
        tmp<volScalarField> p = getp();

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary& devRhoReffb =
            tdevRhoReff().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - csys().origin()
            );

            vectorField fN(Sfb[patchi]*p->boundaryField()[patchi]);

            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // Non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zonei = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - csys().origin());

                const vectorField fDummy(Md.size(), Zero);

                bool calculateForce = checkZone(cZone);

                if (calculateForce)
                {
                    addToFields(cZone, Md, fDummy, fDummy, fP);
                    applyBins(Md, fDummy, fDummy, fP, d);
                }
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::functionObjects::forces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


Foam::vector Foam::functionObjects::forces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]);
}


bool Foam::functionObjects::forces::execute()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        writeBins();

        Log << endl;
    }

    // Write state/results information
    setResult("normalForce", sum(force_[0]));
    setResult("tangentialForce", sum(force_[1]));
    setResult("porousForce", sum(force_[2]));

    setResult("normalMoment", sum(moment_[0]));
    setResult("tangentialMoment", sum(moment_[1]));
    setResult("porousMoment", sum(moment_[2]));

    return true;
}


bool Foam::functionObjects::forces::write()
{
    if (writeFields_)
    {
        lookupObject<volVectorField>(fieldName("force")).write();
        lookupObject<volVectorField>(fieldName("moment")).write();
    }

    return true;
}


// ************************************************************************* //
