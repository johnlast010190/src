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
    (c) 2011-2021 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::transformer

Description
    Vector-tensor class used to perform translations, rotations and scaling
    operations in 3D space.

SourceFiles
    transformerI.H
    transformer.C
    transformerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef transformer_H
#define transformer_H

#include "primitives/Tensor/tensor/tensor.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "fields/Fields/transformField/transformField.H"
#include "VectorN/primitives/vector4/vector4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class transformer;
Istream& operator>>(Istream& is, transformer&);
Ostream& operator<<(Ostream& os, const transformer& C);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the inverse of the given transformer
inline transformer inv(const transformer& tr);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const transformer& tr1, const transformer& tr2);

inline bool operator!=(const transformer& tr1, const transformer& tr2);

inline transformer operator+(const transformer& tr, const vector& t);

inline transformer operator+(const vector& t, const transformer& tr);

inline transformer operator-(const transformer& tr, const vector& t);

inline transformer operator&(const transformer& tr1, const transformer& tr2);


/*---------------------------------------------------------------------------*\
                         Class transformer Declaration
\*---------------------------------------------------------------------------*/

class transformer
{
    // Private data

        //- Translation vector
        vector t_;

        //- True if the translation vector is non-zero
        bool translates_;

        //- Transformation tensor
        tensor T_;

        //- True if the transformation tensor has component vectors of
        //  differing magnitudes
        bool scales_;

        //- True if the transformation tensor has off-diagonal terms
        bool rotates_;


    // Private constructors

        //- Construct given a translation vector and transformation tensor
        inline transformer
        (
            const vector& t,
            const bool translates,
            const tensor& T,
            const bool scales,
            const bool rotates
        );


public:

    // Static Data Members

        static const char* const typeName;

        static const transformer zero;

        static const transformer I;

        //- Null transformer
        static const transformer null;


    // Static Member Functions

        //- Construct a pure translation transformer
        inline static transformer translation(const vector& t);

        //- Construct a pure scaling transformer
        inline static transformer scaling(const tensor& T);

        //- Construct a pure rotation transformer
        inline static transformer rotation(const tensor& T);


    // Constructors

        //- Construct null (i.e., no transformation)
        inline transformer();

        //- Construct from Istream
        explicit transformer(Istream&);


    // Member Functions

        // Access

            //- Return the translation vector
            inline const vector& t() const;

            //- Return true if the transformer performs pure translation
            //  (i.e. the translation vector is non-zero and the transformation
            //  tensor is I)
            inline bool translates() const;

            //- Return the transformation tensor
            inline const tensor& T() const;

            //- Return the inverse transformation tensor
            inline tensor invT() const;

            //- Return true if the transformer performs pure scaling
            //  (i.e. the transformation tensor is diagonal)
            inline bool scales() const;

            //- Return true if the transformer performs pure rotation
            //  (i.e. the transformation tensor is orthogonal)
            inline bool rotates() const;

            //- Return true if the transformer transforms a type
            //  (i.e. scales or rotates)
            inline bool transforms() const;

            //- Return true if the transformer transforms the given type
            //  (i.e. scales or rotates)
            template<typename Type>
            inline bool transforms() const;

            //- Return true if the transformer transforms a point
            // (i.e. translates or scales or rotates)
            inline bool transformsPosition() const;


        // Transform

            //- Transform the given position
            inline vector transformPosition(const vector& v) const;

            //- Transform the given pointField
            void transformPosition(pointField&, const pointField&) const;

            //- Transform the given pointField
            tmp<pointField> transformPosition(const pointField&) const;

            //- Inverse transform the given position
            inline vector invTransformPosition(const vector& v) const;

            //- Inverse transform the given pointField
            void invTransformPosition(pointField&, const pointField&) const;

            //- Inverse transform the given pointField
            tmp<pointField> invTransformPosition(const pointField&) const;

            //- Transform the given type
            template<class Type>
            Type transform(const Type&) const;

            //- Transform the given field
            template<class Type>
            void transform(Field<Type>&, const Field<Type>&) const;

            //- Transform the given field
            template<class Type>
            tmp<Field<Type>> transform(const Field<Type>&) const;

            //- Transform the given field
            template<class Type>
            tmp<Field<Type>> transform(const tmp<Field<Type>>&) const;

            //- Transform the given container
            template<class Type, template<class> class Container>
            void transformList(Container<Type>&) const;

            //- Transform the given vector field, with specialisation
            //  for vector4 type
            template<class Type>
            Type transformVector(const Type& vec) const
            {
                return vec;
            }

            //- Inverse transform the given type
            template<class Type>
            Type invTransform(const Type&) const;

            //- Inverse transform the given field
            template<class Type>
            void invTransform(Field<Type>&, const Field<Type>&) const;

            //- Inverse transform the given field
            template<class Type>
            tmp<Field<Type>> invTransform(const Field<Type>&) const;

            //- Inverse transform the given field
            template<class Type>
            tmp<Field<Type>> invTransform(const tmp<Field<Type>>&) const;

            //- Inverse transform the given container
            template<class Type, template<class> class Container>
            void invTransformList(Container<Type>&) const;


    // Global Functions

        //- Return the inverse of the given transformer
        friend inline transformer inv(const transformer& tr);


    // Global Operators

        friend inline bool operator==
        (
            const transformer& tr1,
            const transformer& tr2
        );

        friend inline bool operator!=
        (
            const transformer& tr1,
            const transformer& tr2
        );

        friend inline transformer operator&
        (
            const transformer& tr1,
            const transformer& tr2
        );


    // IOstream Operators

        friend Istream& operator>>(Istream& is, transformer&);

        friend Ostream& operator<<(Ostream& os, const transformer&);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return a string representation of a transformer
word name(const transformer&);

//- Data associated with transformer type are contiguous
template<>
constexpr bool contiguous<transformer>() {return true;}

// Template specialisations

template<>
tmp<Field<bool>> transformer::transform(const Field<bool>&) const;

template<>
tmp<Field<bool>> transformer::transform(const tmp<Field<bool>>&) const;

template<>
tmp<Field<label>> transformer::transform(const Field<label>&) const;

template<>
tmp<Field<label>> transformer::transform(const tmp<Field<label>>&) const;

template<>
tmp<Field<scalar>> transformer::transform(const Field<scalar>&)
const;

template<>
tmp<Field<scalar>> transformer::transform(const tmp<Field<scalar>>&)
const;

template<>
vector transformer::transformVector(const vector&)
const;

template<>
vector4 transformer::transformVector(const vector4&)
const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/transform/transformer/transformerI.H"

#ifdef NoRepository
    #include "primitives/transform/transformer/transformerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
