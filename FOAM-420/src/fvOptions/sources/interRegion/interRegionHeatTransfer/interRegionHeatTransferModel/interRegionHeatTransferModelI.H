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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::word&
Foam::fv::interRegionHeatTransferModel::nbrRegionName() const
{
    return nbrRegionName_;
}


inline const Foam::meshToMesh&
Foam::fv::interRegionHeatTransferModel::meshInterp() const
{
    if (!meshInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "Interpolation object not set"
            << abort(FatalError);
    }

    return meshInterpPtr_();
}


inline const Foam::volScalarField&
Foam::fv::interRegionHeatTransferModel::htc() const
{
    return  htc_;
}


inline const Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel() const
{
    if (nbrModel_ == nullptr)
    {
        FatalErrorInFunction
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


inline Foam::fv::interRegionHeatTransferModel&
Foam::fv::interRegionHeatTransferModel::nbrModel()
{
    if (nbrModel_ == nullptr)
    {
        FatalErrorInFunction
            << "Neighbour model not set"
            << abort(FatalError);
    }

    return *nbrModel_;
}


// ************************************************************************* //
