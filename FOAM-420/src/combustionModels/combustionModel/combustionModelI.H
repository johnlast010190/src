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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::combustionModel::mesh() const
{
    return mesh_;
}


inline const Foam::compressibleTurbulenceModel&
Foam::combustionModel::turbulence() const
{
    if (!turbulencePtr_)
    {
        if
        (
            obr_.foundObject<compressibleTurbulenceModel>
            (
                IOobject::groupName(turbulenceModel::propertiesName, phaseName_)
            )
        )
        {
            turbulencePtr_ =
                obr_.lookupObjectRefPtr<compressibleTurbulenceModel>
                (
                    IOobject::groupName(turbulenceModel::propertiesName, phaseName_)
                );
        }
        else
        {
            FatalErrorInFunction
                << "turbulencePtr_ is empty. Please use "
                << "combustionModel::setTurbulence "
                << "(compressibleTurbulenceModel&)"
                << abort(FatalError);
        }
    }

    return *turbulencePtr_;
}


inline const Foam::volScalarField& Foam::combustionModel::rho() const
{
    return turbulence().rho();
}


inline Foam::tmp<Foam::surfaceScalarField> Foam::combustionModel::phi() const
{
    return turbulence().alphaRhoPhi();
}


inline const Foam::Switch& Foam::combustionModel::active() const
{
    return active_;
}


inline void Foam::combustionModel::setTurbulence
(
    compressibleTurbulenceModel& turbModel
)
{
    turbulencePtr_ = &turbModel;
}


inline const Foam::dictionary& Foam::combustionModel::coeffs() const
{
    return coeffs_;
}

// ************************************************************************* //
