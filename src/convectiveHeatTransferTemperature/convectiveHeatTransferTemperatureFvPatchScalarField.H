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

Class
    Foam::convectiveHeatTransferTemperatureFvPatchScalarField

Description
    Foam::convectiveHeatTransferTemperatureFvPatchScalarField

SourceFiles
    convectiveHeatTransferTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveHeatTransferTemperatureFvPatchScalarField_H
#define convectiveHeatTransferTemperatureFvPatchScalarField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class convectiveHeatTransferTemperatureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class convectiveHeatTransferTemperatureFvPatchScalarField
:
    public fvPatchScalarField
{
    // Private data

        //- Name of effective thermal conductivity field
        word lambdaEffName_;

        /// External temperature for the convective heat transfer
        word T_ext_Name_;

        /// Heat transfer coefficient
        word htcName_;

        /// External temperature for the radiative heat transfer
        word T_air_Name_;

        /// Surface emissivity
        word Emissivity_Name_;

public:

    //- Runtime type information
    TypeName("convRadHeatTrans");


    // Constructors

        //- Construct from patch and internal field
        convectiveHeatTransferTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveHeatTransferTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given convectiveHeatTransferTemperatureFvPatchScalarField onto a new patch
        convectiveHeatTransferTemperatureFvPatchScalarField
        (
            const convectiveHeatTransferTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        convectiveHeatTransferTemperatureFvPatchScalarField
        (
            const convectiveHeatTransferTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new convectiveHeatTransferTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        convectiveHeatTransferTemperatureFvPatchScalarField
        (
            const convectiveHeatTransferTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new convectiveHeatTransferTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }


        // Return defining fields
/*
            virtual Field<Type>& refValue()
            {
                return T_ext_;
            }

            virtual const Field<Type>& refValue() const
            {
                return T_ext_;
            }

            virtual Field<Type>& refGrad()
            {
                return HTC_;
            }

            virtual const Field<Type>& refGrad() const
            {
                return HTC_;
            }
*/

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<scalar> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<scalar>&) {}

        virtual void operator=(const fvPatchScalarField&) {}
        virtual void operator+=(const fvPatchScalarField&) {}
        virtual void operator-=(const fvPatchScalarField&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<scalar>&) {}
        virtual void operator-=(const Field<scalar>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const scalar&) {}
        virtual void operator+=(const scalar&) {}
        virtual void operator-=(const scalar&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
