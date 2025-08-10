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
    (c) 2010-2017 Esi Ltd.
    (c) 2012-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/porosityModel/porosityModel/porosityModelList.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModelList::porosityModelList
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    bool readFromFvOptions
)
:
    PtrList<porosityModel>(),
    obr_(obr),
    mesh_(mesh)
{
    reset(dict, readFromFvOptions);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModelList::~porosityModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::porosityModelList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No porosity models active" << endl;
    }

    return a;
}


void Foam::porosityModelList::reset
(
    const dictionary& dict,
    bool readFromFvOptions
)
{

    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            if (readFromFvOptions)
            {
                word modelType = iter().dict().lookup("type");

                if (modelType == "explicitPorositySource")
                {
                    count++;
                }
            }
            else
            {
                count++;
            }
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();

            if (readFromFvOptions)
            {
                word modelType = iter().dict().lookup("type");

                const dictionary& modelDict
                    = iter().dict().subDict(modelType+"Coeffs");

                bool active
                    = iter().dict().lookupOrDefault<bool>("active", true);

                if
                (
                    modelType == "explicitPorositySource"
                 && selectionMode(modelDict)
                 && active
                )
                {
                    word cellZoneName = modelDict.lookup("cellZone");

                    this->set
                    (
                        i++,
                        porosityModel::New(name, obr_, mesh_, modelDict, cellZoneName)
                    );
                }
            }
            else
            {
                const dictionary& modelDict = iter().dict();

                this->set
                (
                    i++,
                    porosityModel::New(name, obr_, mesh_, modelDict)
                );
            }
        }
    }
}


bool Foam::porosityModelList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        porosityModel& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::porosityModelList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::porosityModelList::addResistance
(
    fvVectorMatrix& UEqn
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}

void Foam::porosityModelList::addResistance
(
    fvBlockMatrix<vector>& UEqn
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}

void Foam::porosityModelList::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, rho, mu);
    }
}

void Foam::porosityModelList::addResistance
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, rho, mu);
    }
}

void Foam::porosityModelList::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, AU, correctAUprocBC);
    }
}

void Foam::porosityModelList::addAdjointResistance
(
    fvVectorMatrix& UaEqn,
    const volVectorField& Uprimal
)
{
    forAll(*this, i)
    {
        this->operator[](i).addAdjointResistance(UaEqn, Uprimal);
    }
}


void Foam::porosityModelList::addAdjointResistance
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
)
{
    forAll(*this, i)
    {
        this->operator[](i).addAdjointResistance(UaEqn, rho, mu, Uprimal);
    }
}


void Foam::porosityModelList::addAdjointResistance
(
    fvBlockMatrix<vector>& UEqn,
    const volVectorField& Ua
)
{
    forAll(*this, i)
    {
        this->operator[](i).addAdjointResistance(UEqn, Ua);
    }
}


void Foam::porosityModelList::addAdjointResistance
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Ua
)
{
    forAll(*this, i)
    {
        this->operator[](i).addAdjointResistance(UEqn, rho, mu, Ua);
    }
}


void Foam::porosityModelList::addForceMoment
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& porousForce,
    vectorField& porousMoment
)
{
    forAll(*this, i)
    {
        porousForce += this->operator[](i).force(U, rho, mu);
    }
    porousMoment += (porousForce ^ mesh_.C());
}


Foam::wordList Foam::porosityModelList::names() const
{
    const PtrList<porosityModel>& pm = *this;

    wordList lst(pm.size());

    forAll(pm, pmI)
    {
        lst[pmI] = pm[pmI].name();
    }

    return lst;
}


Foam::label Foam::porosityModelList::findPorosityModelID
(
    const word& name
) const
{
    const PtrList<porosityModel>& pm = *this;

    forAll(pm, pmI)
    {
        if (pm[pmI].name() == name)
        {
            return pmI;
        }
    }

    // not found
    return -1;
}


bool Foam::porosityModelList::selectionMode(const dictionary& dict)
{
    word selectionMode = dict.lookup("selectionMode");

    if (selectionMode == "cellZone")
    {
        return true;
    }

    FatalErrorInFunction
        << "The Porous region must be specified as a cellZone.  Current "
        << "selection mode is " << selectionMode
        << exit(FatalError);

        return false;
}




const Foam::porosityModel& Foam::porosityModelList::operator[]
(
    const word& name
) const
{
    const label pmI = findPorosityModelID(name);

    if (pmI < 0)
    {
        FatalErrorInFunction
            << "Porous model named " << name << " not found." << nl
            << "Available porous model names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](pmI);
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const porosityModelList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
