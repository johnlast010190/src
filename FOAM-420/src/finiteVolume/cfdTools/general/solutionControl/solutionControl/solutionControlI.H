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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary& Foam::solutionControl::dict() const
{
    const dictionary* dictPtr =
        mesh_.solution().dict().subDictPtr(solutionControl::typeName);
    return
        dictPtr
      ? dictPtr->optionalSubDict(this->type()+"Coeffs")
      : mesh_.solution().dict().subDict(algorithmName_);
}


inline Foam::label Foam::solutionControl::corr() const
{
    return corr_;
}


inline Foam::label Foam::solutionControl::corrNonOrtho() const
{
    return corrNonOrtho_;
}


inline Foam::label Foam::solutionControl::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline bool Foam::solutionControl::finalNonOrthogonalIter() const
{
    return corrNonOrtho_ == nNonOrthCorr_ + 1;
}


inline bool Foam::solutionControl::momentumPredictor() const
{
    return momentumPredictor_;
}


inline bool Foam::solutionControl::transonic() const
{
    return transonic_;
}


inline bool Foam::solutionControl::consistent() const
{
    return consistent_;
}


inline bool Foam::solutionControl::modifiedMomentumInterp() const
{
    return modifiedMomentumInterp_;
}

inline bool Foam::solutionControl::adjointOptimization() const
{
    return adjointOptimization_;
}

inline bool Foam::solutionControl::frozenFlow() const
{
    return frozenFlow_;
}


inline bool Foam::solutionControl::correctNonOrthogonal()
{
    setFirstIterFlag();

    corrNonOrtho_++;

    if (debug)
    {
        Info<< algorithmName_ << " correctNonOrthogonal: corrNonOrtho = "
            << corrNonOrtho_ << endl;
    }

    if (corrNonOrtho_ <= nNonOrthCorr_ + 1)
    {
        return true;
    }
    else
    {
        corrNonOrtho_ = 0;

        setFirstIterFlag();

        return false;
    }
}


// ************************************************************************* //
