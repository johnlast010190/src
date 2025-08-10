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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::foamCoupledControl::nCorr() const
{
    return nCorr_;
}


inline Foam::label Foam::foamCoupledControl::corr() const
{
    return corr_;
}


inline Foam::label& Foam::foamCoupledControl::corr()
{
    return corr_;
}


inline Foam::label Foam::foamCoupledControl::stabilityMode() const
{
    return stabilityMode_;
}


inline Foam::Switch Foam::foamCoupledControl::choiCorrection() const
{
    return choiCorrection_;
}


inline Foam::scalar Foam::foamCoupledControl::procRelax() const
{
    return procRelax_;
}


inline Foam::scalar Foam::foamCoupledControl::meshQualityRelax() const
{
    return meshQualityRelax_;
}


inline Foam::scalar Foam::foamCoupledControl::nonOrthoThreshold() const
{
    return nonOrthoThreshold_;
}


inline Foam::scalar Foam::foamCoupledControl::skewnessThreshold() const
{
    return skewnessThreshold_;
}


inline Foam::scalar Foam::foamCoupledControl::localCellProcRelax() const
{
    return localCellProcRelax_;
}


inline bool Foam::foamCoupledControl::correct()
{
    setFirstIterFlag();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corr = " << corr_ << endl;
    }

    if (corr_ <= nCorr_)
    {
        return true;
    }
    else
    {
        corr_ = 0;

        setFirstIterFlag();

        return false;
    }
}


inline bool Foam::foamCoupledControl::storeInitialResiduals() const
{
    // Start from second PIMPLE iteration
    return corr_ == 2;
}


inline bool Foam::foamCoupledControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::foamCoupledControl::finalIter() const
{
    return converged_ || (corr_ == nCorr_);
}


inline bool Foam::foamCoupledControl::solveFlow() const
{
    return solveFlow_;
}


inline bool Foam::foamCoupledControl::turbCorr()
{
    if (turbOnFinalIterOnly_)
    {
        if (finalIter())
        {
            setFirstIterFlag(true, true);
        }
    }
    else
    {
        if (firstIter())
        {
            setFirstIterFlag(true, true);
        }
        else
        {
            setFirstIterFlag();
        }
    }

    return !turbOnFinalIterOnly_ || finalIter();
}


// ************************************************************************* //
