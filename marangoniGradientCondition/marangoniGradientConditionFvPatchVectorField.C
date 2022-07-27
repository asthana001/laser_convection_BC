/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "marangoniGradientConditionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"
#include "symmTransformField.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::marangoniGradientConditionFvPatchVectorField::
marangoniGradientConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
	dSigDT_(Zero),
    rhoValue_(1.0),
    rhoName_("rho")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symm(sqr(this->patch().nf()));
}


Foam::marangoniGradientConditionFvPatchVectorField::
marangoniGradientConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    dSigDT_(dict.lookupOrDefault<scalar>("dSigDT",0.0)),
    rhoValue_(dict.lookupOrDefault<scalar>("rho",-1)),
    rhoName_(dict.lookupOrDefault<word>("rhoName","none"))
{
    fvPatchField<vector>::operator=(patchInternalField());
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symm(sqr(this->patch().nf()));
}

Foam::marangoniGradientConditionFvPatchVectorField::
marangoniGradientConditionFvPatchVectorField
(
    const marangoniGradientConditionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    dSigDT_(ptf.dSigDT_),
    rhoValue_(ptf.rhoValue_),
    rhoName_(ptf.rhoName_)
{
	refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symm(sqr(this->patch().nf()));
}


Foam::marangoniGradientConditionFvPatchVectorField::
marangoniGradientConditionFvPatchVectorField
(
    const marangoniGradientConditionFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    dSigDT_(ptf.dSigDT_),
    rhoValue_(ptf.rhoValue_),
    rhoName_(ptf.rhoName_)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symm(sqr(this->patch().nf()));
}


Foam::marangoniGradientConditionFvPatchVectorField::
marangoniGradientConditionFvPatchVectorField
(
    const marangoniGradientConditionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    dSigDT_(ptf.dSigDT_),
    rhoValue_(ptf.rhoValue_),
    rhoName_(ptf.rhoName_)
{
	refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = symm(sqr(this->patch().nf()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::marangoniGradientConditionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    getValues();
    
    vectorField normalValue(transform(valueFraction(), refValue()));
    
    vectorField gradValue
        (this->patchInternalField() + refGrad()/this->patch().deltaCoeffs());

    vectorField transformGradValue
        (transform(I - valueFraction(), gradValue));
        
//    Info<<"updateCoeffs transform= "<<transformGradValue<<nl;
//    Info<<"updateCoeffs normal= "<<normalValue<<nl; 
    operator==
    (
    	normalValue + transformGradValue
    ); 

    directionMixedFvPatchVectorField::updateCoeffs();
//    directionMixedFvPatchVectorField::evaluate();
}


void Foam::marangoniGradientConditionFvPatchVectorField::
evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    getValues();
    
    vectorField normalValue(transform(valueFraction(), refValue()));
    
    vectorField gradValue
        (this->patchInternalField() + refGrad()/this->patch().deltaCoeffs());

    vectorField transformGradValue
        (transform(I - valueFraction(), gradValue));
    vectorField::operator=
    (
        normalValue + transformGradValue
    );
    transformFvPatchField::evaluate();
}


Foam::tmp<Foam::vectorField>
Foam::marangoniGradientConditionFvPatchVectorField::snGrad() const
{
    const vectorField pif(this->patchInternalField());

    tmp<vectorField> normalValue = transform(valueFraction(), refValue());

    tmp<vectorField> gradValue = pif + refGrad()/this->patch().deltaCoeffs();

    tmp<vectorField> transformGradValue =
        transform(I - valueFraction(), gradValue);

    return
        (normalValue + transformGradValue - pif)*
        this->patch().deltaCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::marangoniGradientConditionFvPatchVectorField::snGradTransformDiag() const
{
    vectorField diag(valueFraction().size());
    diag.replace
    (
        vector::X,
        sqrt(mag(valueFraction().component(symmTensor::XX)))
    );
    diag.replace
    (
        vector::Y,
        sqrt(mag(valueFraction().component(symmTensor::YY)))
    );
    diag.replace
    (
        vector::Z,
        sqrt(mag(valueFraction().component(symmTensor::ZZ)))
    );

    return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}


void Foam::marangoniGradientConditionFvPatchVectorField::
write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("dSigDT") << dSigDT_ << token::END_STATEMENT << nl;

    //- No need to write both rhoName and value
    if(rhoName_!="none")
        os.writeKeyword("rhoName") << rhoName_ <<token::END_STATEMENT << nl;
    if(rhoValue_!=-1)
        os.writeKeyword("rho") << rhoValue_ << token::END_STATEMENT << nl;
}


void Foam::marangoniGradientConditionFvPatchVectorField::getValues()
{

	const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    if(!db().foundObject<scalarField>("T"))
    {
        Info<< "This is a thermophysical condition. "
            << "Use temperature condition also" << nl;
        refGrad() = Zero;
        refValue() = Zero;
        valueFraction() = Zero;
        return;
    }

    const volScalarField& T = this->internalField().mesh().lookupObject<volScalarField>("T");
    const tmp<volVectorField> tgradT(fvc::grad(T.oldTime()));
    const vectorField& tGrad = tgradT().boundaryField()[this->patch().index()];
    vectorField nHat (this->patch().nf());
    vectorField tGradPlane(transform(I-sqr(nHat),tGrad));
    scalarField nuEff(turbModel.nuEff(patch().index()));
    vectorField tau1;
    if(db().foundObject<volScalarField>(rhoName_))
    {
            const scalarField& rhop = this->patch().lookupPatchField<volScalarField,scalar>(rhoName_);
            tau1 = (dSigDT_/(rhop*nuEff))*tGradPlane;
            rhoValue_=-1;
    }

    else
    {
        rhoName_="none";
        if(rhoValue_<0)
            FatalErrorInFunction
                << "Did not find registered density field " << rhoName_
                << " or no constant density 'rho' specified"
                << exit(FatalError);
        else
        {
            tau1 = (dSigDT_/(rhoValue_*nuEff))*tGradPlane;
        }
    }


    refGrad() = tau1;
    refValue() = Zero;
    valueFraction() = symm(sqr(nHat));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        marangoniGradientConditionFvPatchVectorField
    );
}

// ************************************************************************* //
