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
    (c) 2017 Esi Ltd.

Class
    Foam::layerPointStack

Description
    Object containing information for layer point stacks

SourceFiles
    layerPointStack.C

\*---------------------------------------------------------------------------*/

#ifndef layerPointStack_H
#define layerPointStack_H

#include "cfdTools/general/include/fvCFD.H"
#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class layerPointStack;
Istream& operator>>(Istream& is, septernion&);
Ostream& operator<<(Ostream& os, const septernion& C);

class layerPointStack
{
private:

    bool baseStack_;

    List<vector> stackPoints_;

    labelList pointStackLabels_;

    label syncAddressing_;

    label layerNu_;

public:

    // Constructors

        //- Construct null
        layerPointStack();

        // Construct from components
        layerPointStack
        (
            const bool& baseStack,
            const List<vector>& stackPoints,
            const labelList& pointStackLabels,
            const label& layerNu
        );

        //- Construct from Istream.
        layerPointStack(Istream&);

    //Destructor
    ~layerPointStack(){};

    // Member functions

        //Edit

            void setAddress(const label& add);
            void increment(const vector& p);
            void substitute(const List<vector>& newStack);
            void merge(const List<layerPointStack>& lPSList);
            static void merge(layerPointStack& ps1, const layerPointStack& ps2);

        //Access

            label size() const{return stackPoints().size();}
            const List<vector>& stackPoints() const {return stackPoints_;}
            const labelList& pointStackLabels() const {return pointStackLabels_;}
            const label& address() const{return syncAddressing_;}
            const bool& baseStack() const{return baseStack_;}
            const label& layerNu() const{return layerNu_;}
            label getIndex(const vector& point) const;

    // Member operators

    // IOstream operators
    friend Istream& operator>>(Istream& is, layerPointStack&);
    friend Ostream& operator<<(Ostream& os, const layerPointStack&);

};

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline layerPointStack transform(const tensor&, const layerPointStack i)
{
    return i;
}

inline bool operator==(const layerPointStack& ps1, const layerPointStack& ps2)
{
    return
    (
        (ps1.stackPoints()==ps2.stackPoints())
     && (ps1.pointStackLabels()==ps2.pointStackLabels())
     && (ps1.baseStack()==ps2.baseStack())
     && (ps1.address()==ps2.address())
    );
}

inline bool operator!=(const layerPointStack& ps1, const layerPointStack& ps2)
{
    return !operator==(ps1, ps2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
