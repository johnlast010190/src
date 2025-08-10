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
    (c) 2017-2021 OpenCFD Ltd.

Class
    Foam::kahipDecomp

Description
    Domain decomposition using KaHIP
    http://algo2.iti.kit.edu/documents/kahip/

    When run in parallel will collect the entire graph on to the master,
    decompose and send back.

    Coefficients dictionary: \a kahipCoeffs, \a coeffs.

    \verbatim
    numberOfSubdomains   N;
    method               kahip;

    kahipCoeffs
    {
        config          fast;
        imbalance       0.01;
    }
    \endverbatim

    Method coefficients:
    \table
        Property  | Description                          | Required | Default
        config    | fast / eco / strong                  | no       | fast
        imbalance | imbalance on cells between domains   | no       | 0.01
        seed      | initial value for random number generator | no  | 0
    \endtable

SourceFiles
    kahipDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef kahipDecomp_H
#define kahipDecomp_H

#include "metisLikeDecomp/metisLikeDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class kahipDecomp Declaration
\*---------------------------------------------------------------------------*/

class kahipDecomp
:
    public metisLikeDecomp
{
protected:

    // Protected Member Functions

        //- Call kahip with options from dictionary.
        virtual label decomposeSerial
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cellWeights,
            List<int>& decomp
        ) const;


        //- No copy construct
        kahipDecomp(const kahipDecomp&) = delete;

        //- No copy assignment
        void operator=(const kahipDecomp&) = delete;


public:

    //- The predefined KaHIP configuration types
    enum configs
    {
        FAST, //!< default
        ECO,
        STRONG,
        FASTSOCIAL,
        ECOSOCIAL,
        STRONGSOCIAL
    };


    //- The selection names for predefined KaHIP configurations
    static const NamedEnum<configs, 6> configNames;


    //- Runtime type information
    TypeName("kahip");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit kahipDecomp(const dictionary& decompDict);


    //- Destructor
    virtual ~kahipDecomp() = default;


    // Member Functions

        virtual bool parallelAware() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
