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

Author
    Georgios Karpouzas

Class
    Foam::Tensor4thRank

Description
    4th rank tensor.In Einstein's notation: t_{ijkn}

SourceFiles
    Tensor4thRankI.H

\*---------------------------------------------------------------------------*/

#ifndef Tensor4thRankTemplate_H
#define Tensor4thRankTemplate_H

#include "primitives/VectorSpace/VectorSpace.H"
#include "primitives/Tensor/Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class Tensor4thRank Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor4thRank
:
    public VectorSpace<Tensor4thRank<Cmpt>, Cmpt, 81>
{

public:

    // Member constants

        enum
        {
            rank = 4
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const Tensor4thRank zero;
        static const Tensor4thRank one;
        static const Tensor4thRank max;
        static const Tensor4thRank min;


        //- Component labeling enumeration
        enum components
        {
            XXXX, XXXY, XXXZ, XXYX, XXYY, XXYZ, XXZX, XXZY, XXZZ,
            XYXX, XYXY, XYXZ, XYYX, XYYY, XYYZ, XYZX, XYZY, XYZZ,
            XZXX, XZXY, XZXZ, XZYX, XZYY, XZYZ, XZZX, XZZY, XZZZ,

            YXXX, YXXY, YXXZ, YXYX, YXYY, YXYZ, YXZX, YXZY, YXZZ,
            YYXX, YYXY, YYXZ, YYYX, YYYY, YYYZ, YYZX, YYZY, YYZZ,
            YZXX, YZXY, YZXZ, YZYX, YZYY, YZYZ, YZZX, YZZY, YZZZ,

            ZXXX, ZXXY, ZXXZ, ZXYX, ZXYY, ZXYZ, ZXZX, ZXZY, ZXZZ,
            ZYXX, ZYXY, ZYXZ, ZYYX, ZYYY, ZYYZ, ZYZX, ZYZY, ZYZZ,
            ZZXX, ZZXY, ZZXZ, ZZYX, ZZYY, ZZYZ, ZZZX, ZZZY, ZZZZ
        };


    // Constructors

        //- Construct null
        inline Tensor4thRank();

        //- Construct given VectorSpace
        inline Tensor4thRank
        (
            const VectorSpace<Tensor4thRank<Cmpt>, Cmpt, 81>&
        );

        inline Tensor4thRank(const Tensor<Cmpt>&);

        //- Construct given the nine components
        inline Tensor4thRank
        (
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,

            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,

            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt,
            const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt, const Cmpt
        );

        //- Construct from Istream
        Tensor4thRank(Istream&);


    // Member Functions

        // Access

        inline void operator=(const Tensor<Cmpt>&);
};


template<class Cmpt>
class typeOfRank<Cmpt, 4>
{
public:

    typedef Tensor4thRank<Cmpt> type;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "primitives/Tensor4thRank/Tensor4thRankI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
