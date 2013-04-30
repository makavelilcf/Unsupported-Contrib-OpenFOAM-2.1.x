/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "convectiveHeatTransferTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
	/// Stefan-Bolzmann constant
	static const dimensionedScalar SigmaSB
	(
		"SigmaSB",
		dimensionSet(1, 0, -3, -4, 0, 0, 0),
		scalar(5.67037321e-08)
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::convectiveHeatTransferTemperatureFvPatchScalarField::convectiveHeatTransferTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    lambdaEffName_("undefinedlambdaEff"),
    T_ext_Name_("undefinedText"),
    htcName_("undefinedHTC"),
    T_air_Name_("undefinedTair"),
    Emissivity_Name_("undefinedEmissivity")
{}


Foam::convectiveHeatTransferTemperatureFvPatchScalarField::convectiveHeatTransferTemperatureFvPatchScalarField
(
    const convectiveHeatTransferTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    lambdaEffName_(ptf.lambdaEffName_),
    T_ext_Name_(ptf.T_ext_Name_),
    htcName_(ptf.htcName_),
    T_air_Name_(ptf.T_air_Name_),
    Emissivity_Name_(ptf.Emissivity_Name_)
{}


Foam::convectiveHeatTransferTemperatureFvPatchScalarField::convectiveHeatTransferTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, dict),
    lambdaEffName_(dict.lookup("lambda")),
    T_ext_Name_(dict.lookup("T_ext")),
    htcName_(dict.lookup("HTC")),
    T_air_Name_(dict.lookup("T_air")),
    Emissivity_Name_(dict.lookup("e"))
{
    evaluate();
}


Foam::convectiveHeatTransferTemperatureFvPatchScalarField::convectiveHeatTransferTemperatureFvPatchScalarField
(
    const convectiveHeatTransferTemperatureFvPatchScalarField& ptf
)
:
    fvPatchField<scalar>(ptf),
    lambdaEffName_(ptf.lambdaEffName_),
    T_ext_Name_(ptf.T_ext_Name_),
    htcName_(ptf.htcName_),
    T_air_Name_(ptf.T_air_Name_),
    Emissivity_Name_(ptf.Emissivity_Name_)
{}


Foam::convectiveHeatTransferTemperatureFvPatchScalarField::convectiveHeatTransferTemperatureFvPatchScalarField
(
    const convectiveHeatTransferTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    lambdaEffName_(ptf.lambdaEffName_),
    T_ext_Name_(ptf.T_ext_Name_),
    htcName_(ptf.htcName_),
    T_air_Name_(ptf.T_air_Name_),
    Emissivity_Name_(ptf.Emissivity_Name_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::convectiveHeatTransferTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}


void Foam::convectiveHeatTransferTemperatureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
}



void Foam::convectiveHeatTransferTemperatureFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);

    const fvPatchField<scalar> & T_ext_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(T_ext_Name_);

    const scalarField A_ = HTC_/(this->patch().deltaCoeffs()*lambda_);

    const fvPatchField<scalar> & T_air_ =
		this->patch().lookupPatchField<volScalarField, scalar>(T_air_Name_);

    const fvPatchField<scalar> & emissivity_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(Emissivity_Name_);

    const scalarField B_ = this->patch().deltaCoeffs()*lambda_;

    const scalarField & T_w_ = *this;

    const scalarField qRad_ =
        emissivity_
       *SigmaSB.value()
       *(
	        Foam::pow(T_w_, scalar(4))
	      - Foam::pow(T_air_, scalar(4))
        );

    Field<scalar>::operator =
    (
        (
            A_*T_ext_
          + this->patchInternalField()
          - qRad_/B_
        )
        /
        (1.0 + A_)
    );

    fvPatchField<scalar>::evaluate();
}

Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveHeatTransferTemperatureFvPatchScalarField::snGrad() const
{
    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);

    const fvPatchField<scalar> & T_ext_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(T_ext_Name_);

    const scalarField A_ = HTC_/(this->patch().deltaCoeffs()*lambda_);

    const fvPatchField<scalar> & T_air_ =
		this->patch().lookupPatchField<volScalarField, scalar>(T_air_Name_);

    const fvPatchField<scalar> & emissivity_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(Emissivity_Name_);

    const scalarField B_ = this->patch().deltaCoeffs()*lambda_;

    const scalarField & T_w_ = *this;

    /// Radiative heat flux
    const scalarField qRad_ =
		emissivity_
	   *SigmaSB.value()
	   *(
			Foam::pow(T_w_, scalar(4))
    	  - Foam::pow(T_air_, scalar(4))
		);

    return
    (
		(
			A_*(T_ext_ - this->patchInternalField())
		  - qRad_/B_
		)
        /
        (1.0 + A_)
        *
        this->patch().deltaCoeffs()
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveHeatTransferTemperatureFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);

    const scalarField A_ = lambda_*this->patch().deltaCoeffs();

    return (A_/(A_ + HTC_));
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveHeatTransferTemperatureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);


    const fvPatchField<scalar> & T_ext_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(T_ext_Name_);

    const scalarField A_ = HTC_/(this->patch().deltaCoeffs()*lambda_);

    const fvPatchField<scalar> & T_air_ =
		this->patch().lookupPatchField<volScalarField, scalar>(T_air_Name_);

    const fvPatchField<scalar> & emissivity_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(Emissivity_Name_);

    const scalarField B_ = this->patch().deltaCoeffs()*lambda_;

    const scalarField & T_w_ = *this;

    /// Radiative heat flux
    const scalarField qRad_ =
		emissivity_
	   *SigmaSB.value()
	   *(
			Foam::pow(T_w_, scalar(4))
    	  - Foam::pow(T_air_, scalar(4))
		);

    return
    ( 
        (
			A_*T_ext_
		  - qRad_/B_
        )
       /(1 + A_)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveHeatTransferTemperatureFvPatchScalarField::gradientInternalCoeffs() const
{
    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);

    const scalarField A_ = HTC_/(this->patch().deltaCoeffs()*lambda_);

    return
    ( 
      - A_*this->patch().deltaCoeffs()
       /(1 + A_)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> > Foam::convectiveHeatTransferTemperatureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const fvPatchField<scalar> & lambda_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(lambdaEffName_);

    const fvPatchField<scalar> & HTC_ =
		this->patch().lookupPatchField<volScalarField, scalar>(htcName_);

    const fvPatchField<scalar> & T_ext_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(T_ext_Name_);

    const scalarField A_ = HTC_/(this->patch().deltaCoeffs()*lambda_);

    const fvPatchField<scalar> & T_air_ =
		this->patch().lookupPatchField<volScalarField, scalar>(T_air_Name_);

    const fvPatchField<scalar> & emissivity_ =
    	this->patch().lookupPatchField<volScalarField, scalar>(Emissivity_Name_);

    const scalarField B_ = this->patch().deltaCoeffs()*lambda_;

    const scalarField & T_w_ = *this;

    /// Radiative heat flux
    const scalarField qRad_ =
		emissivity_
	   *SigmaSB.value()
	   *(
			Foam::pow(T_w_, scalar(4))
    	  - Foam::pow(T_air_, scalar(4))
		);

    return
    (
		(
			A_*T_ext_
		  - qRad_/B_
		)
       /(1 + A_)
       *this->patch().deltaCoeffs()
    );
}


void Foam::convectiveHeatTransferTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("lambda") << lambdaEffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("T_ext") << T_ext_Name_ << token::END_STATEMENT << nl;
    os.writeKeyword("HTC") << htcName_ << token::END_STATEMENT << nl;
    os.writeKeyword("T_air") << T_air_Name_ << token::END_STATEMENT << nl;
    os.writeKeyword("e") << Emissivity_Name_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        convectiveHeatTransferTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
