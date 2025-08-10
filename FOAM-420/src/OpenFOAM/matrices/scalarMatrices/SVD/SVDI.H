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

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T Foam::SVD::sign(const T& a, const T& b)
{
    //return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
    return b >= 0 ? a : -a;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::scalarRectangularMatrix& Foam::SVD::U() const
{
    return U_;
}


inline const Foam::scalarRectangularMatrix& Foam::SVD::V() const
{
    return V_;
}


inline const Foam::scalarDiagonalMatrix& Foam::SVD::S() const
{
    return S_;
}


inline bool Foam::SVD::converged() const
{
    return converged_;
}


inline Foam::label Foam::SVD::nZeros() const
{
    return nZeros_;
}


inline Foam::scalar Foam::SVD::minNonZeroS() const
{
    scalar minS = S_[0];
    for (label i=1; i<S_.size(); i++)
    {
        scalar s = S_[i];
        if (s > VSMALL && s < minS) minS = s;
    }
    return minS;
}


// ************************************************************************* //
