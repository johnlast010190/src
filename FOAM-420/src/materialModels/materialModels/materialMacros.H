#ifndef materialMacros_H
#define materialMacros_H

#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define field functions
#define defineFieldFunctionsTwo(Type,funcName,field1Ptr,field2Ptr)       \
    virtual inline bool is##funcName##Const() const override             \
    {                                                                    \
        return field1Ptr->isConst() && field2Ptr->isConst();             \
    }                                                                    \
                                                                         \
    virtual inline Type funcName##Cell(const label celli) const override \
    {                                                                    \
        return funcName                                                  \
        (                                                                \
            field1Ptr->operator[](celli),                                \
            field2Ptr->operator[](celli)                                 \
        );                                                               \
    }                                                                    \
                                                                         \
    /* Calculate internal field */                                       \
    virtual tmp<Type##Field> funcName##Internal() const override         \
    {                                                                    \
        tmp<scalarField> tfield1(field1Ptr->primitiveField());           \
        scalarField& field1 = tfield1.ref();                             \
        tmp<scalarField> tfield2(field2Ptr->primitiveField());           \
        scalarField& field2 = tfield2.ref();                             \
        tmp<Type##Field> tInternalField(new Type##Field(mesh().nCells()));\
        Type##Field& fieldCells = tInternalField.ref();                  \
        forAll(fieldCells, celli)                                        \
        {                                                                \
            fieldCells[celli] =                                          \
                funcName(field1[celli], field2[celli]);                  \
        }                                                                \
        return tInternalField;                                           \
    }                                                                    \
                                                                         \
    /* Get patch fields references */                                    \
    virtual tmp<Type##Field> funcName##Patch                             \
    (                                                                    \
        const label patchi                                               \
    ) const override                                                     \
    {                                                                    \
        tmp<scalarField> tfield1(field1Ptr->boundaryField()[patchi]);    \
        scalarField& field1 = tfield1.ref();                             \
        tmp<scalarField> tfield2(field2Ptr->boundaryField()[patchi]);    \
        scalarField& field2 = tfield2.ref();                             \
        tmp<Type##Field> tPatchField(new Type##Field(field1.size()));    \
        Type##Field& patchField = tPatchField.ref();                     \
        forAll(patchField, facei)                                        \
        {                                                                \
            patchField[facei] = funcName(field1[facei], field2[facei]);  \
        }                                                                \
        return tPatchField;                                              \
    }                                                                    \
                                                                         \
    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>             \
    funcName##Geometric() const override                                 \
    {                                                                    \
        typedef GeometricField<Type, fvPatchField, volMesh> volField;    \
        tmp<volField> tField                                             \
        (                                                                \
            new volField                                                 \
            (                                                            \
                IOobject                                                 \
                (                                                        \
                    #funcName,                                           \
                    mesh().time().timeName(),                            \
                    obr_,                                                \
                    IOobject::NO_READ,                                   \
                    IOobject::NO_WRITE,                                  \
                    false                                                \
                ),                                                       \
                mesh(),                                                  \
                dimensioned<Type>(#funcName, funcName##Model::dimensions(), Zero)\
            )                                                            \
        );                                                               \
        volField& field = tField.ref();                                  \
                                                                         \
        /* Calculate internal field */                                   \
        field.primitiveFieldRef() = funcName##Internal();                \
                                                                         \
        /* Calculate bounary patches */                                  \
        forAll(field.boundaryField(), patchi)                            \
        {                                                                \
            field.boundaryFieldRef()[patchi].forceAssign(funcName##Patch(patchi)); \
        }                                                                \
        return tField;                                                   \
    }


#define defineFieldFunctionsOne(Type,funcName,fieldPtr)                  \
    virtual inline bool is##funcName##Const() const override             \
    {                                                                    \
        return fieldPtr->isConst();                                                    \
    }                                                                    \
                                                                         \
    virtual inline Type funcName##Cell(const label celli) const override \
    {                                                                    \
        return funcName(fieldPtr->operator[](celli));                    \
    }                                                                    \
                                                                         \
    /* Calculate internal field */                                       \
    virtual tmp<Type##Field> funcName##Internal() const override         \
    {                                                                    \
        tmp<scalarField> tfield(fieldPtr->primitiveField());             \
        scalarField& field = tfield.ref();                               \
        tmp<Type##Field> tInternalField(new Type##Field(field.size()));  \
        Type##Field& fieldCells = tInternalField.ref();                  \
        forAll(fieldCells, celli)                                        \
        {                                                                \
            fieldCells[celli] = funcName(field[celli]);                  \
        }                                                                \
        return tInternalField;                                           \
    }                                                                    \
                                                                         \
    /* Get patch fields references */                                    \
    virtual tmp<Type##Field> funcName##Patch                             \
    (                                                                    \
        const label patchi                                               \
    ) const override                                                     \
    {                                                                    \
        tmp<scalarField> tfield(fieldPtr->boundaryField()[patchi]);      \
        scalarField& field = tfield.ref();                               \
        tmp<Type##Field> tPatchField(new Type##Field(field.size()));     \
        Type##Field& patchField = tPatchField.ref();                     \
        forAll(patchField, facei)                                        \
        {                                                                \
            patchField[facei] = funcName(field[facei]);                  \
        }                                                                \
        return tPatchField;                                              \
    }                                                                    \
                                                                         \
    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>             \
    funcName##Geometric() const override                                 \
    {                                                                    \
        typedef GeometricField<Type, fvPatchField, volMesh> volField;    \
        tmp<volField> tField                                             \
        (                                                                \
            new volField                                                 \
            (                                                            \
                IOobject                                                 \
                (                                                        \
                    #funcName,                                           \
                    mesh().time().timeName(),                            \
                    obr_,                                                \
                    IOobject::NO_READ,                                   \
                    IOobject::NO_WRITE,                                  \
                    false                                                \
                ),                                                       \
                mesh(),                                                  \
                dimensioned<Type>(#funcName, funcName##Model::dimensions(), Zero)\
            )                                                            \
        );                                                               \
        volField& field = tField.ref();                                  \
                                                                         \
        /* Calculate internal field */                                   \
         field.primitiveFieldRef() = funcName##Internal();               \
                                                                         \
        /* Calculate bounary patches */                                  \
        forAll(field.boundaryField(), patchi)                            \
        {                                                                \
            field.boundaryFieldRef()[patchi].forceAssign(funcName##Patch(patchi));  \
        }                                                                \
        return tField;                                                   \
    }


#define defineFieldFunctionsZero(Type,funcName)                          \
    virtual inline bool is##funcName##Const() const override             \
    {                                                                    \
        return true;                                                     \
    }                                                                    \
                                                                         \
    virtual inline Type funcName##Cell(const label celli) const override \
    {                                                                    \
        return funcName();                                               \
    }                                                                    \
                                                                         \
    /* Calculate internal field */                                       \
    virtual tmp<Type##Field> funcName##Internal() const override         \
    {                                                                    \
        tmp<Type##Field> tInternalField(new Type##Field(mesh().nCells()));\
        Type##Field& fieldCells = tInternalField.ref();                  \
        forAll(fieldCells, celli)                                        \
        {                                                                \
            fieldCells[celli] = funcName();                              \
        }                                                                \
        return tInternalField;                                           \
    }                                                                    \
                                                                         \
    /* Get patch fields references */                                    \
    virtual tmp<Type##Field> funcName##Patch                             \
    (                                                                    \
        const label patchi                                               \
    ) const override                                                     \
    {                                                                    \
        tmp<Type##Field> tPatchField(new Type##Field(mesh().boundary()[patchi].size()));\
        Type##Field& patchField = tPatchField.ref();                     \
        forAll(patchField, facei)                                        \
        {                                                                \
            patchField[facei] = funcName();                              \
        }                                                                \
        return tPatchField;                                              \
    }                                                                    \
                                                                         \
    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>             \
    funcName##Geometric() const override                                 \
    {                                                                    \
        typedef GeometricField<Type, fvPatchField, volMesh> volField;    \
        tmp<volField> tField                                             \
        (                                                                \
            new volField                                                 \
            (                                                            \
                IOobject                                                 \
                (                                                        \
                    #funcName,                                           \
                    mesh().time().timeName(),                            \
                    obr_,                                                \
                    IOobject::NO_READ,                                   \
                    IOobject::NO_WRITE,                                  \
                    false                                                \
                ),                                                       \
                mesh(),                                                  \
                dimensioned<Type>(#funcName, funcName##Model::dimensions(), Zero)\
            )                                                            \
        );                                                               \
        volField& field = tField.ref();                                  \
                                                                         \
        /* Calculate internal field */                                   \
        Type##Field& fieldCells = field.primitiveFieldRef();             \
        forAll(fieldCells, celli)                                        \
        {                                                                \
            fieldCells[celli] = funcName();                              \
        }                                                                \
                                                                         \
        /* Calculate bounary patches */                                  \
        forAll(field.boundaryField(), patchi)                            \
        {                                                                \
            Type##Field& patchField = field.boundaryFieldRef()[patchi];  \
            forAll(patchField, facei)                                    \
            {                                                            \
                patchField[facei] = funcName();                          \
            }                                                            \
        }                                                                \
        return tField;                                                   \
    }


#define declareTypeFunctionTwo(returnType, funcName, Type1, Type2)       \
    virtual inline returnType funcName                                   \
    (                                                                    \
        const Type1 val1,                                                \
        const Type2 val2                                                 \
    ) const;

#define declareTypeFunctionOne(returnType, funcName, Type1)              \
    virtual inline returnType funcName(const Type1 val) const;

#define declareTypeFunctionZero(returnType, funcName)                    \
    virtual inline returnType funcName() const;

#define defineScalarFieldFunctionsNone(funcName)                         \
    inline scalar funcName() const { return 0.0; }                       \
    defineFieldFunctionsZero(scalar,funcName)

#define defineVectorFieldFunctionsNone(funcName)                         \
    inline vector funcName() const { return vector::zero; }              \
    defineFieldFunctionsZero(vector,funcName)

#define defineTensorFieldFunctionsNone(funcName)                         \
    inline tensor funcName() const { return tensor::zero; }              \
    defineFieldFunctionsZero(tensor,funcName)

#define defineScalarFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)      \
    declareTypeFunctionTwo(scalar, funcName, scalar, scalar)             \
    defineFieldFunctionsTwo(scalar,funcName,field1Ptr,field2Ptr)

#define defineVectorFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)      \
    declareTypeFunctionTwo(vector, funcName, scalar, scalar)             \
    defineFieldFunctionsTwo(vector,funcName,field1Ptr,field2Ptr)

#define defineTensorFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)      \
    declareTypeFunctionTwo(tensor, funcName, scalar, scalar)             \
    defineFieldFunctionsTwo(tensor,funcName,field1Ptr,field2Ptr)

#define defineScalarFieldFunctionsOne(funcName,fieldPtr)                 \
    declareTypeFunctionOne(scalar, funcName, scalar)                     \
    defineFieldFunctionsOne(scalar,funcName,fieldPtr)

#define defineVectorFieldFunctionsOne(funcName,fieldPtr)                 \
    declareTypeFunctionOne(vector, funcName, scalar)                     \
    defineFieldFunctionsOne(vector,funcName,fieldPtr)

#define defineTensorFieldFunctionsOne(funcName,fieldPtr)                 \
    declareTypeFunctionOne(tensor, funcName, scalar)                     \
    defineFieldFunctionsOne(tensor,funcName,fieldPtr)

#define defineScalarFieldFunctionsZero(funcName)                         \
    declareTypeFunctionZero(scalar, funcName)                            \
    defineFieldFunctionsZero(scalar,funcName)

#define defineVectorFieldFunctionsZero(funcName)                         \
    declareTypeFunctionZero(vector, funcName)                            \
    defineFieldFunctionsZero(vector,funcName)

#define defineTensorFieldFunctionsZero(funcName)                         \
    declareTypeFunctionZero(tensor, funcName)                            \
    defineFieldFunctionsZero(tensor,funcName)

//- Declare functions
#define declareVirtualFieldFunctions(Type,funcName,virtual)              \
    virtual inline Type funcName##Cell(const label celli) const override;\
                                                                         \
    virtual tmp<Type##Field> funcName##Internal() const override;        \
                                                                         \
    virtual tmp<Type##Field> funcName##Patch                             \
    (                                                                    \
        const label patchi                                               \
    ) const override;                                                    \
                                                                         \
    virtual tmp<GeometricField<Type, fvPatchField, volMesh>>             \
    funcName##Geometric() const override;

//- Declare functions
#define declareFieldFunctions(Type,funcName)                             \
    declareVirtualFieldFunctions(Type,funcName,)


#define declareScalarFieldFunctions(funcName)                            \
    declareVirtualFieldFunctions(scalar,funcName,)

#define declareVectorFieldFunctions(funcName)                            \
    declareVirtualFieldFunctions(vector,funcName,)

#define declareTensorFieldFunctions(funcName)                            \
    declareVirtualFieldFunctions(tensor,funcName,)


#define defineMaterialModel(funcName,type,PatchField,GeoMesh)            \
    namespace Foam                                                       \
    {                                                                    \
    class funcName##Model                                                \
    :                                                                    \
        public baseModels<type,PatchField,GeoMesh>                       \
    {                                                                    \
    public:                                                              \
        TypeName(#funcName);                                             \
        virtual ~funcName##Model(){delete boundary_;};                   \
        static const dimensionSet modelDims;                             \
        static const bool isMixture;                                     \
        class Boundary                                                   \
        :                                                                \
            public baseModels<type,PatchField,GeoMesh>::Boundary         \
        {                                                                \
        public:                                                          \
            Boundary(const funcName##Model& base):funcName##Model_(base) \
            {};                                                          \
            const funcName##Model& funcName##Model_;                     \
            virtual ~Boundary(){};                                       \
            virtual tmp<Field<type>> operator[](const label patchi) const\
            {                                                            \
                return funcName##Model_.funcName##Patch(patchi);         \
            }                                                            \
        };                                                               \
        typedef type ModelType;                                          \
        static const word tempName;                                      \
        Boundary* boundary_;                                             \
        funcName##Model(){boundary_ = new Boundary(*this);};             \
        virtual bool is##funcName##Const() const = 0;                    \
        virtual tmp<GeometricField<type,PatchField,GeoMesh>>             \
        funcName##Geometric() const = 0;                                 \
        virtual tmp<type##Field> funcName##Internal() const = 0;         \
        virtual tmp<type##Field> funcName##Patch                         \
        (                                                                \
            label const patchi                                           \
        ) const = 0;                                                     \
        virtual type funcName##Cell(label const celli) const = 0;        \
        virtual bool isConst() const override                            \
        {                                                                \
            return is##funcName##Const();                                \
        }                                                                \
        virtual tmp<type##Field> primitiveField() const override         \
        {                                                                \
            return funcName##Internal();                                 \
        }                                                                \
        virtual baseModels<type,PatchField,GeoMesh>::Boundary&           \
        boundaryField() const override                                   \
        {                                                                \
            return *boundary_;                                           \
        }                                                                \
        virtual const dimensionSet& dimensions() const override          \
        {                                                                \
            return modelDims;                                            \
        }                                                                \
        virtual type operator[](const label celli) const override        \
        {                                                                \
            return funcName##Cell(celli);                                \
        }                                                                \
        virtual tmp<GeometricField<type,PatchField,GeoMesh>>             \
        operator()() const override                                      \
        {                                                                \
            return funcName##Geometric();                                \
        }                                                                \
        virtual void updateTable(const word&) override = 0;              \
        virtual const depList& dep() const override = 0;                 \
        virtual const word funcType() const override                     \
        {                                                                \
            return #funcName;                                            \
        }                                                                \
        virtual const word modelType() const override                    \
        {                                                                \
            return #type;                                                \
        }                                                                \
        virtual bool read() override = 0;                                \
    };                                                                   \
    }


#define declareMaterialModelDims(funcName,dim,type,PatchField,GeoMesh,mix)  \
    defineTypeNameAndDebug(funcName##Model, 0);                         \
    const dimensionSet funcName##Model::modelDims = dim;                \
    const bool funcName##Model::isMixture = mix;                \
    baseModels<type,PatchField,GeoMesh>::addToTable<funcName##Model> test##funcName##_;\
    const word funcName##Model::tempName = #type;


#define defineScalarMaterialModel(funcName)                              \
    defineMaterialModel(funcName,scalar,fvPatchField,volMesh)

#define defineVectorMaterialModel(funcName)                              \
    defineMaterialModel(funcName,vector,fvPatchField,volMesh)

#define defineTensorMaterialModel(funcName)                              \
    defineMaterialModel(funcName,tensor,fvPatchField,volMesh)


#define declareMaterialModelScalarMixture(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,scalar,fvPatchField,volMesh,true)

#define declareMaterialModelVectorMixture(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,vector,fvPatchField,volMesh,true)

#define declareMaterialModelTensorMixture(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,tensor,fvPatchField,volMesh,true)

#define declareMaterialModelScalarDims(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,scalar,fvPatchField,volMesh,false)

#define declareMaterialModelVectorDims(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,vector,fvPatchField,volMesh,false)

#define declareMaterialModelTensorDims(funcName,dims)                    \
    declareMaterialModelDims(funcName,dims,tensor,fvPatchField,volMesh,false)

#define makeMatThermo(BaseThermo,MixtureType,MixtureThermo,matTypeName)       \
    typedef MixtureThermo<BaseThermo, MixtureType>                            \
    MixtureThermo##BaseThermo##MixtureType;                                   \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        MixtureThermo##BaseThermo##MixtureType,                               \
        #matTypeName,                                                         \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        basicThermo,                                                          \
        MixtureThermo##BaseThermo##MixtureType,                               \
        objectRegistry                                                        \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        fluidThermo,                                                          \
        MixtureThermo##BaseThermo##MixtureType,                               \
        objectRegistry                                                        \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        BaseThermo,                                                           \
        MixtureThermo##BaseThermo##MixtureType,                               \
        objectRegistry                                                        \
    );


#define makeMatSolidThermo(BaseThermo,MixtureThermo,MixtureType,matTypeName)  \
    typedef matHeThermo<BaseThermo, MixtureType> matHeThermo##MixtureType;    \
    typedef MixtureThermo<BaseThermo, MixtureType>                            \
        MixtureThermo##BaseThermo##MixtureType;                               \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        MixtureThermo##BaseThermo##MixtureType,                               \
        #matTypeName,                                                         \
        0                                                                     \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        basicThermo,                                                          \
        MixtureThermo##BaseThermo##MixtureType,                               \
        objectRegistry                                                        \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        BaseThermo,                                                           \
        MixtureThermo##BaseThermo##MixtureType,                               \
        objectRegistry                                                        \
    );                                                                        \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        BaseThermo,                                                           \
        MixtureThermo##BaseThermo##MixtureType,                               \
        dictionary                                                            \
    );




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //