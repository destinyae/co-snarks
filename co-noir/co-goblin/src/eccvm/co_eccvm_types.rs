use crate::eccvm::co_ecc_op_queue::{
    CoECCOpQueue, CoEccvmOpsTable, CoUltraEccOpsTable, CoUltraOp, CoVMOperation,
};
use crate::eccvm::co_ecc_op_queue::{MSMRow, ScalarMul};
use ark_ec::AffineRepr;
use ark_ec::CurveGroup;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::Zero;
use co_acvm::mpc::NoirWitnessExtensionProtocol;
use co_builder::flavours::eccvm_flavour::ECCVMFlavour;
use co_builder::prelude::NUM_DISABLED_ROWS_IN_SUMCHECK;
use co_builder::prelude::NUM_TRANSLATION_EVALUATIONS;
use co_builder::prelude::Polynomial;
use co_builder::prelude::offset_generator_scaled;
use co_builder::{
    HonkProofResult,
    prelude::{HonkCurve, ProverCrs},
};
use co_ultrahonk::prelude::Polynomials;
use co_ultrahonk::prelude::SharedSmallSubgroupIPAProver;
use co_ultrahonk::prelude::SharedUnivariate;
use co_ultrahonk::prelude::SharedUnivariateTrait;
use common::CoUtils;
use common::shared_polynomial::SharedPolynomial;
use common::{
    mpc::NoirUltraHonkProver,
    transcript::{Transcript, TranscriptFieldType, TranscriptHasher},
};
use goblin::{NUM_WNAF_DIGIT_BITS, NUM_WNAF_DIGITS_PER_SCALAR};
use goblin::{POINT_TABLE_SIZE, WNAF_DIGITS_PER_ROW};
use mpc_core::MpcState;
use mpc_net::Network;
use num_bigint::BigUint;
use std::marker::PhantomData;

#[derive(Default)]
pub(crate) struct SharedTranslationData<T: NoirUltraHonkProver<P>, P: CurveGroup> {
    // M(X) whose Lagrange coefficients are given by (m_0 || m_1 || ... || m_{NUM_TRANSLATION_EVALUATIONS-1} || 0 || ... || 0)
    pub(crate) concatenated_polynomial_lagrange: SharedPolynomial<T, P>,

    // M(X) + Z_H(X) * R(X), where R(X) is a random polynomial of length = WITNESS_MASKING_TERM_LENGTH
    pub(crate) masked_concatenated_polynomial: SharedPolynomial<T, P>,
    // Interpolation domain {1, g, \ldots, g^{SUBGROUP_SIZE - 1}} required for Lagrange interpolation
    pub(crate) interpolation_domain: Vec<P::ScalarField>,
}

impl<T: NoirUltraHonkProver<P>, P: HonkCurve<TranscriptFieldType>> SharedTranslationData<T, P> {
    pub(crate) fn new(interpolation_domain: Vec<P::ScalarField>) -> Self {
        Self {
            concatenated_polynomial_lagrange: SharedPolynomial::new_zero(P::SUBGROUP_SIZE),
            masked_concatenated_polynomial: SharedPolynomial::new_zero(P::SUBGROUP_SIZE * 2),
            interpolation_domain,
        }
    }
    pub(crate) fn construct_translation_data<
        H: TranscriptHasher<TranscriptFieldType>,
        N: Network,
    >(
        transcript_polynomials: &[&Vec<<T as NoirUltraHonkProver<P>>::ArithmeticShare>],
        transcript: &mut Transcript<TranscriptFieldType, H>,
        crs: &ProverCrs<P>,
        net: &N,
        state: &mut T::State,
    ) -> HonkProofResult<Self> {
        // Create interpolation domain required for Lagrange interpolation
        let mut interpolation_domain = vec![P::ScalarField::one(); P::SUBGROUP_SIZE];
        let subgroup_generator = P::get_subgroup_generator();
        for idx in 1..P::SUBGROUP_SIZE {
            interpolation_domain[idx] = interpolation_domain[idx - 1] * subgroup_generator;
        }

        let mut translation_data = Self::new(interpolation_domain);

        // Concatenate the last entries of the `translation_polynomials`.

        translation_data.compute_concatenated_polynomials(transcript_polynomials, net, state)?;

        // Commit to M(X) + Z_H(X)*R(X), where R is a random polynomial of WITNESS_MASKING_TERM_LENGTH.
        let commitment = CoUtils::commit::<T, P>(
            translation_data.masked_concatenated_polynomial.as_ref(),
            crs,
        );
        let open = T::open_point(commitment, net, state)?;
        transcript.send_point_to_verifier::<P>(
            "Translation:concatenated_masking_term_commitment".to_string(),
            open.into(),
        );

        Ok(translation_data)
    }

    pub fn compute_small_ipa_prover<H: TranscriptHasher<TranscriptFieldType>, N: Network>(
        &mut self,
        evaluation_challenge_x: P::ScalarField,
        batching_challenge_v: P::ScalarField,
        transcript: &mut Transcript<TranscriptFieldType, H>,
        net: &N,
        state: &mut T::State,
    ) -> HonkProofResult<SharedSmallSubgroupIPAProver<T, P>> {
        let mut small_ipa_prover = SharedSmallSubgroupIPAProver::<T, P> {
            interpolation_domain: self.interpolation_domain.to_owned(),
            concatenated_polynomial: self.masked_concatenated_polynomial.to_owned(),
            libra_concatenated_lagrange_form: self.concatenated_polynomial_lagrange.to_owned(),
            challenge_polynomial: Polynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::SUBGROUP_SIZE,
            ),
            challenge_polynomial_lagrange: Polynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::SUBGROUP_SIZE,
            ),
            grand_sum_polynomial_unmasked: SharedPolynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::SUBGROUP_SIZE,
            ),
            grand_sum_polynomial: SharedPolynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::MASKED_GRAND_SUM_LENGTH,
            ),
            grand_sum_lagrange_coeffs: vec![
                T::ArithmeticShare::default();
                SharedSmallSubgroupIPAProver::<T, P>::SUBGROUP_SIZE
            ],
            grand_sum_identity_polynomial: SharedPolynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::GRAND_SUM_IDENTITY_LENGTH,
            ),
            grand_sum_identity_quotient: SharedPolynomial::new_zero(
                SharedSmallSubgroupIPAProver::<T, P>::QUOTIENT_LENGTH,
            ),
            claimed_inner_product: P::ScalarField::zero(),
            prefix_label: "Translator:".to_string(),
            phantom_data: PhantomData,
        };

        small_ipa_prover
            .compute_eccvm_challenge_polynomial(evaluation_challenge_x, batching_challenge_v);

        let mut claimed_inner_product = T::ArithmeticShare::default();
        for idx in 0..P::SUBGROUP_SIZE {
            let tmp = T::mul_with_public(
                small_ipa_prover.challenge_polynomial_lagrange[idx],
                self.concatenated_polynomial_lagrange[idx],
            );
            T::add_assign(&mut claimed_inner_product, tmp);
        }
        let claimed_inner_product = T::open_many(&[claimed_inner_product], net, state)?[0];
        transcript.send_fr_to_verifier::<P>(
            "Translation:masking_term_eval".to_string(),
            claimed_inner_product,
        );
        small_ipa_prover.claimed_inner_product = claimed_inner_product;

        Ok(small_ipa_prover)
    }

    fn compute_concatenated_polynomials<N: Network>(
        &mut self,
        transcript_polynomials: &[&Vec<<T as NoirUltraHonkProver<P>>::ArithmeticShare>],
        net: &N,
        state: &mut T::State,
    ) -> HonkProofResult<()> {
        const WITNESS_MASKING_TERM_LENGTH: usize = 2;
        let circuit_size = transcript_polynomials[0].len();

        let mut coeffs_lagrange_subgroup = vec![T::ArithmeticShare::default(); P::SUBGROUP_SIZE];

        // Extract the Lagrange coefficients of the concatenated masking term from the transcript polynomials
        for poly_idx in 0..NUM_TRANSLATION_EVALUATIONS {
            for idx in 0..NUM_DISABLED_ROWS_IN_SUMCHECK {
                let idx_to_populate = poly_idx * NUM_DISABLED_ROWS_IN_SUMCHECK + idx;
                coeffs_lagrange_subgroup[idx_to_populate as usize] = transcript_polynomials
                    [poly_idx as usize]
                    [circuit_size - NUM_DISABLED_ROWS_IN_SUMCHECK as usize + idx as usize];
            }
        }
        self.concatenated_polynomial_lagrange = SharedPolynomial::new(coeffs_lagrange_subgroup);

        // Generate the masking term
        let masking_scalars =
            SharedUnivariate::<T, P, WITNESS_MASKING_TERM_LENGTH>::get_random(net, state)?;

        // Compute monomial coefficients of the concatenated polynomial
        let concatenated_monomial_form_unmasked = SharedPolynomial::<T, P>::interpolate_from_evals(
            &self.interpolation_domain,
            &self.concatenated_polynomial_lagrange.coefficients,
            P::SUBGROUP_SIZE,
        );

        self.masked_concatenated_polynomial =
            SharedPolynomial::new_zero(P::SUBGROUP_SIZE + WITNESS_MASKING_TERM_LENGTH);
        for idx in 0..P::SUBGROUP_SIZE {
            self.masked_concatenated_polynomial[idx] = concatenated_monomial_form_unmasked[idx];
        }

        // Mask the polynomial in monomial form.
        for idx in 0..WITNESS_MASKING_TERM_LENGTH {
            self.masked_concatenated_polynomial[idx] = T::sub(
                self.masked_concatenated_polynomial[idx],
                masking_scalars.evaluations_as_ref()[idx],
            );

            T::add_assign(
                &mut self.masked_concatenated_polynomial[P::SUBGROUP_SIZE + idx],
                masking_scalars.evaluations_as_ref()[idx],
            );
        }
        Ok(())
    }
}

struct CoVMState<C: HonkCurve<TranscriptFieldType>, T: NoirWitnessExtensionProtocol<C::BaseField>> {
    pc: T::AcvmType,
    count: T::AcvmType,
    accumulator: T::AcvmPoint<C>,
    msm_accumulator: T::AcvmPoint<C>,
    is_accumulator_empty: T::AcvmType, //bool
}

impl<C: HonkCurve<TranscriptFieldType>, T: NoirWitnessExtensionProtocol<C::BaseField>> Clone
    for CoVMState<C, T>
{
    fn clone(&self) -> Self {
        Self {
            pc: self.pc,
            count: self.count,
            accumulator: self.accumulator,
            msm_accumulator: self.msm_accumulator,
            is_accumulator_empty: self.is_accumulator_empty,
        }
    }
}
impl<C: HonkCurve<TranscriptFieldType>, T: NoirWitnessExtensionProtocol<C::BaseField>>
    CoVMState<C, T>
{
    fn new() -> Self {
        Self {
            pc: T::AcvmType::default(),
            count: T::AcvmType::default(),
            accumulator: T::AcvmPoint::<C>::default(),
            msm_accumulator: T::AcvmPoint::from(offset_generator_scaled::<C>().into()),
            is_accumulator_empty: T::AcvmType::from(C::BaseField::one()), //true
        }
    }

    fn process_mul(
        entry: &CoVMOperation<T, C>,
        updated_state: &mut CoVMState<C, T>,
        state: &CoVMState<C, T>,
        driver: &mut T,
    ) {
        // TODO FLORIN
        let p = entry.base_point;
        let r = state.msm_accumulator;
        //TACEO TODO Can we batch these scalar muls?
        let mul = driver.scalar_mul_many(&[p], &[entry.mul_scalar_full])[0];
        updated_state.msm_accumulator = driver.add_points(r, mul);
    }

    fn process_add(
        entry: &CoVMOperation<T, C>,
        updated_state: &mut CoVMState<C, T>,
        old_state: &CoVMState<C, T>,
        is_accumulator_empty: T::OtherAcvmType<C>,
        driver: &mut T,
    ) -> eyre::Result<()> {
        // TODO FLORIN
        let mul = driver.mul_with_public_other(-C::ScalarField::one(), is_accumulator_empty);
        let inv = driver.add_other(T::OtherAcvmType::from(C::ScalarField::one()), mul);
        let other = driver.add_points(old_state.accumulator, entry.base_point);
        let mul = driver.scalar_mul_many(&[entry.base_point, other], &[is_accumulator_empty, inv]);
        let result = driver.add_points(mul[0], mul[1]);
        updated_state.accumulator = result;

        updated_state.is_accumulator_empty =
            driver.point_is_zero_many(&[updated_state.accumulator])?[0];
        Ok(())
    }

    // TODO FLORIN: explain what is going on here
    fn process_msm_transition(
        row: &mut CoTranscriptRow<C, T>,
        updated_state: &mut CoVMState<C, T>,
        old_state: &CoVMState<C, T>,
        is_accumulator_empty: T::OtherAcvmType<C>,
        msm_transition_is_zero: T::OtherAcvmType<C>,
        driver: &mut T,
    ) -> eyre::Result<()> {
        //TODO FLORIN
        let mul = driver.mul_with_public_other(-C::ScalarField::one(), is_accumulator_empty);
        let inv = driver.add_other(T::OtherAcvmType::from(C::ScalarField::one()), mul);
        let if_value = driver.add_points(
            updated_state.msm_accumulator,
            T::AcvmPoint::from(-offset_generator_scaled::<C>().into()),
        );
        let mut else_value =
            driver.add_points(old_state.accumulator, updated_state.msm_accumulator);
        else_value = driver.add_points(
            else_value,
            T::AcvmPoint::from(-offset_generator_scaled::<C>().into()),
        );
        let mul = driver.scalar_mul_many(&[if_value, else_value], &[is_accumulator_empty, inv]);
        let result = driver.add_points(mul[0], mul[1]);
        let mul = driver.scalar_mul_many(
            &[driver.sub_points(result, updated_state.accumulator)],
            &[msm_transition_is_zero],
        )[0];
        updated_state.accumulator = driver.add_points(updated_state.accumulator, mul);

        let msm_output = driver.sub_points(
            updated_state.msm_accumulator,
            T::AcvmPoint::from(offset_generator_scaled::<C>().into()),
        );
        let is_zero = driver.point_is_zero_many(&[msm_output, updated_state.accumulator])?;

        //TODO FLORIN MAYBE WE NEED TO MULTIPLY THIS STILL WITH msm_transition_is_zero
        updated_state.is_accumulator_empty = is_zero[1];

        //TACEO TODO: Batch this is_zero check with others
        row.transcript_msm_infinity = is_zero[0];
        Ok(())
    }

    fn populate_transcript_row(
        row: &mut CoTranscriptRow<C, T>,
        base_point_infinity: T::AcvmType,
        entry: &CoVMOperation<T, C>,
        state: &CoVMState<C, T>,
        msm_transition: T::AcvmType,
        driver: &mut T,
    ) -> eyre::Result<()> {
        row.accumulator_empty = state.is_accumulator_empty;
        row.q_add = entry.op_code.add;
        row.q_mul = entry.op_code.mul;
        row.q_eq = entry.op_code.eq;
        row.q_reset_accumulator = entry.op_code.reset;
        row.msm_transition = msm_transition;
        row.pc = state.pc;
        row.msm_count = state.count;
        // row.msm_count_zero_at_transition =
        // (state.count + num_muls == 0) && entry.op_code.mul && next_not_msm; // We do this already outside of the function
        //TACEO TODO Batch this function
        let points = driver.pointshare_to_field_shares(entry.base_point)?;

        if entry.op_code.add || entry.op_code.mul || entry.op_code.eq {
            let mut inv = driver.mul_with_public(-C::BaseField::one(), base_point_infinity);
            driver.add_assign_with_public(C::BaseField::one(), &mut inv);
            let mul = driver.mul_many(&[points.0, points.1], &[inv, inv])?;
            row.base_x = mul[0].to_owned();
            row.base_y = mul[1].to_owned();
        }
        // row.base_x = if (entry.op_code.add || entry.op_code.mul || entry.op_code.eq)
        //     && !base_point_infinity
        // {
        //     entry
        //         .base_point
        //         .x()
        //         .expect("Base point x should not be zero")
        // } else {
        //     C::BaseField::zero()
        // };
        // row.base_y = if (entry.op_code.add || entry.op_code.mul || entry.op_code.eq)
        //     && !base_point_infinity
        // {
        //     entry
        //         .base_point
        //         .y()
        //         .expect("Base point y should not be zero")
        // } else {
        //     C::BaseField::zero()
        // };
        row.base_infinity = if entry.op_code.add || entry.op_code.mul || entry.op_code.eq {
            base_point_infinity
        } else {
            T::AcvmType::from(C::BaseField::zero())
        };
        row.z1 = if entry.op_code.mul {
            entry.z1
        } else {
            T::AcvmType::default()
        };
        row.z2 = if entry.op_code.mul {
            entry.z2
        } else {
            T::AcvmType::default()
        };
        // row.z1_zero = entry.z1.is_zero(); // We do this already outside of the function
        // row.z2_zero = entry.z2.is_zero(); // We do this already outside of the function
        row.opcode = entry.op_code.value();
        Ok(())
    }
}

fn add_affine_coordinates_to_transcript<
    C: HonkCurve<TranscriptFieldType>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
>(
    transcript_state: &mut [CoTranscriptRow<C, T>],
    accumulator_trace: &[T::AcvmPoint<C>],
    msm_accumulator_trace: &[T::AcvmPoint<C>],
    intermediate_accumulator_trace: &[T::AcvmPoint<C>],
    driver: &mut T,
) -> eyre::Result<()> {
    let (xs, ys, _) = driver.pointshare_to_field_shares_many(
        &[
            accumulator_trace,
            msm_accumulator_trace,
            intermediate_accumulator_trace,
        ]
        .concat(),
    )?;

    let len_acc = accumulator_trace.len();
    let len_msm = msm_accumulator_trace.len();
    let (acc_xs, rest) = xs.split_at(len_acc);
    let (msm_xs, int_xs) = rest.split_at(len_msm);
    let (acc_ys, rest) = ys.split_at(len_acc);
    let (msm_ys, int_ys) = rest.split_at(len_msm);

    for i in 0..accumulator_trace.len() {
        let row = &mut transcript_state[i + 1];
        row.accumulator_x = acc_xs[i];
        row.accumulator_y = acc_ys[i];
        row.msm_output_x = msm_xs[i];
        row.msm_output_y = msm_ys[i];
        row.transcript_msm_intermediate_x = int_xs[i];
        row.transcript_msm_intermediate_y = int_ys[i];
    }
    Ok(())
}

#[expect(clippy::too_many_arguments)]
fn compute_inverse_trace_coordinates<
    C: HonkCurve<TranscriptFieldType>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
>(
    msm_transition: T::AcvmType,
    row: &CoTranscriptRow<C, T>,
    intermediate_accumulator_trace_x: T::AcvmType,
    intermediate_accumulator_trace_y: T::AcvmType,
    transcript_msm_x_inverse_trace: &mut T::AcvmType,
    msm_accumulator_trace_x: T::AcvmType,
    msm_accumulator_trace_infinity: T::AcvmType,
    accumulator_trace_x: T::AcvmType,
    accumulator_trace_y: T::AcvmType,
    inverse_trace_x: &mut T::AcvmType,
    inverse_trace_y: &mut T::AcvmType,
    driver: &mut T,
) -> eyre::Result<()> {
    // let msm_output_infinity = intermediate_accumulator_trace.is_zero();
    let row_msm_infinity = row.transcript_msm_infinity;
    //TODO FLORIN: can do this over scalarfield also and remove some of these functions
    let mul = driver.mul_with_public(-C::BaseField::one(), row_msm_infinity);
    let inv_row_msm_infinity = driver.add(T::AcvmType::from(C::BaseField::one()), mul);

    let mul = driver.mul_with_public(-C::BaseField::one(), msm_accumulator_trace_infinity);
    let inv_accumulator_trace_infinity = driver.add(T::AcvmType::from(C::BaseField::one()), mul);

    let bb_infinity_default =
        driver.mul_with_public(C::get_bb_infinity_default(), msm_accumulator_trace_infinity);

    let mul = driver.mul_many(
        &[msm_accumulator_trace_x],
        &[inv_accumulator_trace_infinity],
    )?;
    let mut result = driver.add(mul[0].to_owned(), bb_infinity_default);
    driver.add_assign_with_public(
        -offset_generator_scaled::<C>()
            .x()
            .expect("Offset generator x-coordinate should not be zero"),
        &mut result,
    );
    let mul = driver.mul_with_public(-C::BaseField::one(), msm_transition);
    let inv_msm_transition = driver.add(T::AcvmType::from(C::BaseField::one()), mul);
    let mul = driver.mul_many(
        &[
            result,
            msm_transition,
            msm_transition,
            inv_msm_transition,
            inv_msm_transition,
        ],
        &[
            inv_row_msm_infinity,
            intermediate_accumulator_trace_x,
            intermediate_accumulator_trace_y,
            row.base_x,
            row.base_y,
        ],
    )?;

    *transcript_msm_x_inverse_trace = driver.mul(result, inv_row_msm_infinity)?;

    let res_x = driver.add(mul[1], mul[3]);
    let res_y = driver.add(mul[2], mul[4]);

    let (lhsx, lhsy) = (res_x, res_y);

    let (rhsx, rhsy) = (accumulator_trace_x, accumulator_trace_y);

    *inverse_trace_x = driver.add(lhsx, rhsx); //lhsx - rhsx;
    *inverse_trace_y = driver.add(lhsy, rhsy); //lhsy - rhsy;

    Ok(())
}

struct CoTranscriptRow<
    C: HonkCurve<TranscriptFieldType>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
> {
    transcript_msm_infinity: T::AcvmType, //bool
    accumulator_empty: T::AcvmType,
    q_add: bool,
    q_mul: bool,
    q_eq: bool,
    q_reset_accumulator: bool,
    msm_transition: T::AcvmType,
    pc: T::AcvmType,
    msm_count: T::AcvmType,
    msm_count_zero_at_transition: T::AcvmType,
    base_x: T::AcvmType,
    base_y: T::AcvmType,
    base_infinity: T::AcvmType,
    z1: T::AcvmType,
    z2: T::AcvmType,
    z1_zero: T::AcvmType,
    z2_zero: T::AcvmType,
    opcode: u32,

    accumulator_x: T::AcvmType,
    accumulator_y: T::AcvmType,
    msm_output_x: T::AcvmType,
    msm_output_y: T::AcvmType,
    transcript_msm_intermediate_x: T::AcvmType,
    transcript_msm_intermediate_y: T::AcvmType,

    transcript_add_x_equal: T::AcvmType,
    transcript_add_y_equal: T::AcvmType,

    base_x_inverse: T::AcvmType,
    base_y_inverse: T::AcvmType,
    transcript_add_lambda: T::AcvmType,
    transcript_msm_x_inverse: T::AcvmType,
    msm_count_at_transition_inverse: T::AcvmType,
}

impl<C: HonkCurve<TranscriptFieldType>, T: NoirWitnessExtensionProtocol<C::BaseField>> Default
    for CoTranscriptRow<C, T>
{
    fn default() -> Self {
        Self {
            transcript_msm_infinity: T::AcvmType::default(),
            accumulator_empty: T::AcvmType::default(),
            q_add: false,
            q_mul: false,
            q_eq: false,
            q_reset_accumulator: false,
            msm_transition: T::AcvmType::default(),
            pc: T::AcvmType::default(),
            msm_count: T::AcvmType::default(),
            msm_count_zero_at_transition: T::AcvmType::default(),
            base_x: T::AcvmType::default(),
            base_y: T::AcvmType::default(),
            base_infinity: T::AcvmType::default(),
            z1: T::AcvmType::default(),
            z2: T::AcvmType::default(),
            z1_zero: T::AcvmType::default(),
            z2_zero: T::AcvmType::default(),
            opcode: 0,
            accumulator_x: T::AcvmType::default(),
            accumulator_y: T::AcvmType::default(),
            msm_output_x: T::AcvmType::default(),
            msm_output_y: T::AcvmType::default(),
            transcript_msm_intermediate_x: T::AcvmType::default(),
            transcript_msm_intermediate_y: T::AcvmType::default(),
            transcript_add_x_equal: T::AcvmType::default(),
            transcript_add_y_equal: T::AcvmType::default(),
            base_x_inverse: T::AcvmType::default(),
            base_y_inverse: T::AcvmType::default(),
            transcript_add_lambda: T::AcvmType::default(),
            transcript_msm_x_inverse: T::AcvmType::default(),
            msm_count_at_transition_inverse: T::AcvmType::default(),
        }
    }
}

fn finalize_transcript<
    C: HonkCurve<TranscriptFieldType>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
>(
    updated_state: &CoVMState<C, T>,
    driver: &mut T,
) -> eyre::Result<CoTranscriptRow<C, T>>
where
    <C as CurveGroup>::BaseField: PrimeField,
{
    let mut final_row = CoTranscriptRow::<C, T>::default();

    let (result_x, result_y, _) = driver.pointshare_to_field_shares(updated_state.accumulator)?; //TODO FLORIN: batch this outside?

    final_row.accumulator_x = result_x;
    final_row.accumulator_y = result_y;

    final_row.pc = updated_state.pc;
    final_row.accumulator_empty = updated_state.is_accumulator_empty;
    Ok(final_row)
}

fn compute_rows<
    C: HonkCurve<TranscriptFieldType>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
>(
    vm_operations: &[CoVMOperation<T, C>],
    total_number_of_muls: T::AcvmType,
    driver: &mut T,
) -> eyre::Result<Vec<CoTranscriptRow<C, T>>> {
    // TODO FLORIN: REPLACE WITH CMUXES where possible

    let num_vm_entries = vm_operations.len();
    // The transcript contains an extra zero row at the beginning and the accumulated state at the end
    let transcript_size = num_vm_entries + 2;
    let mut transcript_state = Vec::with_capacity(transcript_size);

    // These vectors track quantities that we need to invert.
    // We fill these vectors and then perform batch inversions to amortize the cost of FF inverts
    let mut inverse_trace_x = vec![T::AcvmType::default(); num_vm_entries];
    let mut inverse_trace_y = vec![T::AcvmType::default(); num_vm_entries];
    let mut transcript_msm_x_inverse_trace = vec![T::AcvmType::default(); num_vm_entries];
    let mut msm_count_at_transition_inverse_trace = vec![T::AcvmType::default(); num_vm_entries];

    let mut msm_accumulator_trace: Vec<_> = vec![T::AcvmPoint::<C>::default(); num_vm_entries];
    let mut accumulator_trace: Vec<_> = vec![T::AcvmPoint::<C>::default(); num_vm_entries];
    let mut intermediate_accumulator_trace: Vec<_> =
        vec![T::AcvmPoint::<C>::default(); num_vm_entries];

    let mut state = CoVMState::<C, T> {
        pc: total_number_of_muls,
        count: T::AcvmType::default(),
        accumulator: T::AcvmPoint::<C>::default(),
        msm_accumulator: T::AcvmPoint::<C>::from(offset_generator_scaled::<C>().into()),
        is_accumulator_empty: T::AcvmType::from(C::BaseField::one()), //true
    };

    let mut updated_state = CoVMState::<C, T>::new();

    // add an empty row. 1st row all zeroes because of our shiftable polynomials
    transcript_state.push(CoTranscriptRow::<C, T>::default());

    // during the first iteration over the ECCOpQueue, the operations are being performed using Jacobian
    // coordinates and the base point coordinates are recorded in the transcript. at the same time, the transcript
    // logic is being populated

    let mut base_points = Vec::with_capacity(num_vm_entries);
    let mut entry_z1 = Vec::with_capacity(num_vm_entries);
    let mut entry_z2 = Vec::with_capacity(num_vm_entries);
    for entry in vm_operations.iter() {
        entry_z1.push(entry.z1);
        entry_z2.push(entry.z2);
        base_points.push(entry.base_point);
    }

    let mut is_zero_results = driver.is_zero_many(&[entry_z1, entry_z2].concat())?;
    let (z1_zero_results_slice, z2_zero_results_slice) =
        is_zero_results.split_at_mut(base_points.len());
    let (z1_is_zero_unchanged, z2_is_zero_unchanged) = (
        z1_zero_results_slice.to_vec(),
        z2_zero_results_slice.to_vec(),
    );
    driver.scale_many_in_place(z1_zero_results_slice, -C::BaseField::one());
    driver.add_scalar_in_place(z1_zero_results_slice, C::BaseField::one());
    driver.scale_many_in_place(z2_zero_results_slice, -C::BaseField::one());
    driver.add_scalar_in_place(z2_zero_results_slice, C::BaseField::one());
    let num_mul_partial = driver.add_many(z1_zero_results_slice, z2_zero_results_slice);

    let base_points_is_zero = driver.point_is_zero_many(&base_points)?;
    let mut base_points_is_zero_modified = base_points_is_zero.clone();
    driver.scale_many_in_place(&mut base_points_is_zero_modified, -C::BaseField::one());
    driver.add_scalar_in_place(&mut base_points_is_zero_modified, C::BaseField::one());
    let num_mul = driver.mul_many(&num_mul_partial, &base_points_is_zero_modified)?;

    for i in 0..num_vm_entries {
        let mut row = CoTranscriptRow::<C, T>::default();
        let entry = &vm_operations[i];
        updated_state = state.clone();

        let is_mul: bool = entry.op_code.mul;
        let is_add: bool = entry.op_code.add;
        // let z1_zero: bool = if is_mul { entry.z1.is_zero() } else { true };
        // let z2_zero: bool = if is_mul { entry.z2.is_zero() } else { true };

        // let base_point_infinity = entry.base_point.is_zero();
        let num_muls = num_mul[i];

        updated_state.pc = driver.add(state.pc, num_muls);

        if entry.op_code.reset {
            updated_state.is_accumulator_empty = T::AcvmType::from(C::BaseField::one()); //true;
            updated_state.accumulator = T::AcvmPoint::<C>::default();
            updated_state.msm_accumulator =
                T::AcvmPoint::from(offset_generator_scaled::<C>().into());
        }

        let last_row = i == (num_vm_entries - 1);

        // msm transition = current row is doing a lookup to validate output = msm output
        // i.e. next row is not part of MSM and current row is part of MSM
        //   or next row is irrelevant and current row is a straight MUL
        let next_not_msm = last_row || !vm_operations[i + 1].op_code.mul;

        //     // we reset the count in updated state if we are not accumulating and not doing an msm
        let mut msm_transition_public = true;
        let mut is_zero = T::AcvmType::default();
        // is_mul && next_not_msm && (state.count + num_muls > 0);
        if !(is_mul && next_not_msm) {
            msm_transition_public = false;
        } else {
            let add = driver.add(state.count, num_muls);
            is_zero = driver.is_zero_many(&[add])?[0]; //TODO FLORIN BATCH WITH OTHER IS_ZERO
        }
        let msm_transition = driver.mul_with_public(
            C::BaseField::from(entry.op_code.mul && next_not_msm),
            is_zero,
        );
        row.msm_count_zero_at_transition = msm_transition; // This happens in bb inside populate_transcript_row, for simplicity we do it here 
        // we want state.count + num_muls > 0, hence we invert the is_zero result
        is_zero = driver.add(T::AcvmType::from(C::BaseField::one()), is_zero);

        // determine ongoing msm and update the respective counter
        let current_ongoing_msm = is_mul && !next_not_msm;

        updated_state.count = if current_ongoing_msm {
            driver.add(state.count, num_muls)
        } else {
            T::AcvmType::default()
        };

        if is_mul {
            CoVMState::<C, T>::process_mul(entry, &mut updated_state, &state, driver);
        }

        let old_state_accumulator_is_zero = driver.is_zero_many(&[state.is_accumulator_empty])?[0];

        if msm_transition_public {
            CoVMState::<C, T>::process_msm_transition(
                &mut row,
                &mut updated_state,
                &state,
                T::convert_fields(&[old_state_accumulator_is_zero])?[0],
                T::convert_fields(&[is_zero])?[0],
                driver,
            )?; //TODO NEED TO MULTYIPLY/CORRECT THIS WITH THE IS_ZEROCHECK
        } else {
            msm_accumulator_trace[i] = T::AcvmPoint::<C>::default();
            intermediate_accumulator_trace[i] = T::AcvmPoint::<C>::default();
        }

        if is_add {
            CoVMState::<C, T>::process_add(
                entry,
                &mut updated_state,
                &state,
                T::convert_fields(&[old_state_accumulator_is_zero])?[0],
                driver,
            )?;
        }

        row.z1_zero = z1_is_zero_unchanged[i]; // We do this already outside of the function
        row.z2_zero = z2_is_zero_unchanged[i]; // We do this already outside of the function

        // populate the first group of TranscriptRow entries
        CoVMState::<C, T>::populate_transcript_row(
            &mut row,
            base_points_is_zero[i],
            entry,
            &state,
            msm_transition,
            driver,
        )?;

        msm_count_at_transition_inverse_trace[i] = driver.add(state.count, num_muls);

        //      update the accumulators
        accumulator_trace[i] = state.accumulator;
        let msm_transition_as_scalarfield = T::convert_fields(&[msm_transition])?[0];
        let mul = driver.scalar_mul_many(
            &[
                updated_state.msm_accumulator,
                driver.add_points(
                    updated_state.msm_accumulator,
                    T::AcvmPoint::from(-offset_generator_scaled::<C>().into()),
                ),
            ],
            &[msm_transition_as_scalarfield, msm_transition_as_scalarfield],
        );
        msm_accumulator_trace[i] = mul[0];
        intermediate_accumulator_trace[i] = mul[1];

        state = updated_state.clone();

        if is_mul && next_not_msm {
            state.msm_accumulator = T::AcvmPoint::from(offset_generator_scaled::<C>().into());
        }
        transcript_state.push(row);
    }

    // add required affine coordinates to the transcript

    let accumulator_trace_len = accumulator_trace.len();
    let msm_accumulator_trace_len = msm_accumulator_trace.len();
    let (xs, ys, inf) = driver.pointshare_to_field_shares_many(
        &[
            accumulator_trace.clone(),
            msm_accumulator_trace.clone(),
            intermediate_accumulator_trace.clone(),
        ]
        .concat(),
    )?;
    //TODO FLORIN: check sizes
    let len_acc = accumulator_trace_len;
    let len_msm = msm_accumulator_trace_len;
    let (acc_xs, rest) = xs.split_at(len_acc);
    let (msm_xs, int_xs) = rest.split_at(len_msm);
    let (acc_ys, rest) = ys.split_at(len_acc);
    let (msm_ys, int_ys) = rest.split_at(len_msm);
    let (acc_inf, rest) = inf.split_at(len_acc);
    let (msm_inf, _int_inf) = rest.split_at(len_msm);

    for i in 0..accumulator_trace_len {
        let row = &mut transcript_state[i + 1];
        row.accumulator_x = acc_xs[i];
        row.accumulator_y = acc_ys[i];
        row.msm_output_x = msm_xs[i];
        row.msm_output_y = msm_ys[i];
        row.transcript_msm_intermediate_x = int_xs[i];
        row.transcript_msm_intermediate_y = int_ys[i];
    }

    // process the slopes when adding points or results of MSMs. to increase efficiency, we use batch inversion
    // after the loop

    let mut is_zero_vm_point = Vec::new();
    for i in 0..accumulator_trace_len {
        let entry = &vm_operations[i];
        if entry.op_code.add {
            is_zero_vm_point.push(vm_operations[i].base_point);
        } else {
            is_zero_vm_point.push(intermediate_accumulator_trace[i]);
        }
    }
    //TODO FLORIN: BATCH ALL THESE as much as possible
    let (vm_points_x, vm_points_y, vm_points_inf) =
        driver.pointshare_to_field_shares_many(&is_zero_vm_point)?; //TODO FLORIN BATCH THIS INTO THE ABOVE CALL
    let vm_x_squared = driver.mul_many(&vm_points_x, &vm_points_x)?; //TODO FLORIN BATCH THIS INTO THE ABOVE CALL
    let vm_x_squared_times_3 = driver.scale_many(&vm_x_squared, C::BaseField::from(3u32));
    let vm_y_doubled = driver.add_many(&vm_points_y, &vm_points_y);
    let acc_x_minus_vm_x = driver.sub_many(acc_xs, &vm_points_x);
    let acc_y_minus_vm_y = driver.sub_many(acc_ys, &vm_points_y);

    let vm_inf_and_acc_inf = driver.mul_many(&vm_points_inf, acc_inf)?;
    let scale = driver.scale_many(&vm_inf_and_acc_inf, -C::BaseField::one());
    let inv_vm_inf_and_acc_inf = driver.add_scalar(&scale, C::BaseField::one());
    let to_cmp_x = driver.sub_many(&vm_points_x, acc_xs);
    let to_cmp_y = driver.sub_many(&vm_points_y, acc_ys);
    let is_zero_transcript_add = driver.is_zero_many(&[to_cmp_x, to_cmp_y].concat())?;
    let (is_zero_transcript_add_x, is_zero_transcript_add_y) =
        is_zero_transcript_add.split_at(is_zero_transcript_add.len() / 2);

    let transcript_add_values = driver.add_many(
        &is_zero_transcript_add[..num_vm_entries],
        &vm_inf_and_acc_inf[..num_vm_entries],
    );
    let transcript_add_values = driver.is_zero_many(&transcript_add_values)?; //TODO FLORIN: is there a better way to do this? 
    let (transcript_add_x_equal, transcript_add_y_equal) =
        transcript_add_values.split_at(transcript_add_values.len() / 2);
    let mul = driver.mul_many(is_zero_transcript_add_x, is_zero_transcript_add_y)?; //(accumulator_x == vm_x) && (accumulator_y == vm_y
    let next_mul = driver.mul_many(&mul, &inv_vm_inf_and_acc_inf)?; //(accumulator_x == vm_x) && (accumulator_y == vm_y) && !vm_infinity && !accumulator_infinity
    let scale = driver.scale_many(transcript_add_x_equal, -C::BaseField::one());
    let inv_transcript_add_x_equal = driver.add_scalar(&scale, C::BaseField::one());
    let else_mul = driver.mul_many(&vm_inf_and_acc_inf, &inv_transcript_add_x_equal)?; //(vm_infinity && accumulator_infinity) && !((accumulator_x == vm_x) && (accumulator_y == vm_y))
    let lambda_denom_1 = driver.mul_many(&vm_y_doubled, &next_mul)?; //vm_y + vm_y if (accumulator_x == vm_x) && (accumulator_y == vm_y) && !vm_infinity && !accumulator_infinity
    let lambda_denom_2 = driver.mul_many(&acc_x_minus_vm_x, &else_mul)?; //accumulator_x - vm_x if (vm_infinity && accumulator_infinity) && !((accumulator_x == vm_x) && (accumulator_y == vm_y))
    let lambda_num_1 = driver.mul_many(&vm_x_squared_times_3, &next_mul)?; //vm_x * vm_x * 3 if (accumulator_x == vm_x) && (accumulator_y == vm_y) && !vm_infinity && !accumulator_infinity
    let lambda_num_2 = driver.mul_many(&acc_y_minus_vm_y, &else_mul)?; //accumulator_y - vm_y if (vm_infinity && accumulator_infinity) && !((accumulator_x == vm_x) && (accumulator_y == vm_y))
    let mut add_lambda_numerator = driver.add_many(&lambda_num_1, &lambda_num_2);
    let mut add_lambda_denominator = driver.add_many(&lambda_denom_1, &lambda_denom_2);
    for i in 0..accumulator_trace_len {
        let row = &mut transcript_state[i + 1];
        let msm_transition = row.msm_transition;

        // compute the differences between point coordinates
        // TACEO Note: We compute everything and then multiply it by the msm_transition is_zero afterwards
        compute_inverse_trace_coordinates::<C, T>(
            msm_transition,
            row,
            int_xs[i],
            int_ys[i],
            &mut transcript_msm_x_inverse_trace[i],
            msm_xs[i],
            msm_inf[i],
            acc_xs[i],
            acc_ys[i],
            &mut inverse_trace_x[i],
            &mut inverse_trace_y[i],
            driver,
        )?;

        row.transcript_add_x_equal = transcript_add_x_equal[i]; //(vm_x == accumulator_x) || (vm_infinity && accumulator_infinity);
        row.transcript_add_y_equal = transcript_add_y_equal[i]; //(vm_y == accumulator_y) || (vm_infinity && accumulator_infinity);

        // compute the numerators and denominators of slopes between the points
    }
    let mut tmp_transcript_add_x_equal = Vec::new(); //TODO FLORIN
    let mut tmp_transcript_add_y_equal = Vec::new(); //TODO FLORIN
    let mut tmp_add_lambda_numerator = Vec::new(); //TODO FLORIN
    let mut tmp_add_lambda_denominator = Vec::new(); //TODO FLORIN
    let mut tmp_inverse_trace_x = Vec::new(); //TODO FLORIN
    let mut tmp_inverse_trace_y = Vec::new(); //TODO FLORIN
    let mut tmp_msm_transition = Vec::new(); //TODO FLORIN
    let mut indices = Vec::new(); //TODO FLORIN
    for i in 0..accumulator_trace_len {
        let row = &transcript_state[i + 1];
        if vm_operations[i].op_code.add {
            continue;
        } else {
            tmp_transcript_add_x_equal.push(row.transcript_add_x_equal);
            tmp_transcript_add_y_equal.push(row.transcript_add_y_equal);
            tmp_add_lambda_numerator.push(add_lambda_numerator[i]);
            tmp_add_lambda_denominator.push(add_lambda_denominator[i]);
            tmp_inverse_trace_x.push(inverse_trace_x[i]);
            tmp_inverse_trace_y.push(inverse_trace_y[i]);
            tmp_msm_transition.push(row.msm_transition);
            indices.push(i);
        }
    }
    let mul = driver.mul_many(
        &[
            tmp_transcript_add_x_equal,
            tmp_transcript_add_y_equal,
            tmp_add_lambda_numerator,
            tmp_add_lambda_denominator,
            tmp_inverse_trace_x,
            tmp_inverse_trace_y,
        ]
        .concat(),
        &[
            tmp_msm_transition.clone(),
            tmp_msm_transition.clone(),
            tmp_msm_transition.clone(),
            tmp_msm_transition.clone(),
            tmp_msm_transition.clone(),
            tmp_msm_transition.clone(),
        ]
        .concat(),
    )?; //TODO FLORIN is it possible to make this nicer
    for (j, i) in indices.iter().enumerate() {
        transcript_state[i + 1].transcript_add_x_equal = mul[j];
        transcript_state[i + 1].transcript_add_y_equal = mul[j + indices.len()];
        add_lambda_numerator[*i] = mul[j + 2 * indices.len()];
        add_lambda_denominator[*i] = mul[j + 3 * indices.len()];
        inverse_trace_x[*i] = mul[j + 4 * indices.len()];
        inverse_trace_y[*i] = mul[j + 5 * indices.len()];
    }

    // Perform all required inversions at once
    //TODO FLORIN: INVERT THESE:
    // ark_ff::batch_inversion(&mut inverse_trace_x);
    // ark_ff::batch_inversion(&mut inverse_trace_y);
    // ark_ff::batch_inversion(&mut transcript_msm_x_inverse_trace);
    // ark_ff::batch_inversion(&mut add_lambda_denominator);
    // ark_ff::batch_inversion(&mut msm_count_at_transition_inverse_trace);

    // Populate the fields of the transcript row containing inverted scalars
    let mul = driver.mul_many(&add_lambda_numerator, &add_lambda_denominator)?;
    for i in 0..num_vm_entries {
        let row = &mut transcript_state[i + 1];
        row.base_x_inverse = inverse_trace_x[i];
        row.base_y_inverse = inverse_trace_y[i];
        row.transcript_msm_x_inverse = transcript_msm_x_inverse_trace[i];
        row.transcript_add_lambda = mul[i];
        row.msm_count_at_transition_inverse = msm_count_at_transition_inverse_trace[i];
    }

    // process the final row containing the result of the sequence of group ops in ECCOpQueue
    let final_row = finalize_transcript(&updated_state, driver)?;
    transcript_state.push(final_row);

    Ok(transcript_state)
}
#[derive(Debug)]
struct PointTablePrecomputationRow<
    C: CurveGroup<BaseField: PrimeField>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
> {
    s1: T::AcvmType,
    s2: T::AcvmType,
    s3: T::AcvmType,
    s4: T::AcvmType,
    s5: T::AcvmType,
    s6: T::AcvmType,
    s7: T::AcvmType,
    s8: T::AcvmType,
    skew: T::AcvmType,
    point_transition: bool,
    pc: u32,
    round: u32,
    scalar_sum: T::AcvmType,
    precompute_accumulator: T::AcvmPoint<C>,
    precompute_double: T::AcvmPoint<C>,
}

impl<C: CurveGroup<BaseField: PrimeField>, T: NoirWitnessExtensionProtocol<C::BaseField>> Default
    for PointTablePrecomputationRow<C, T>
{
    fn default() -> Self {
        Self {
            s1: T::AcvmType::default(),
            s2: T::AcvmType::default(),
            s3: T::AcvmType::default(),
            s4: T::AcvmType::default(),
            s5: T::AcvmType::default(),
            s6: T::AcvmType::default(),
            s7: T::AcvmType::default(),
            s8: T::AcvmType::default(),
            skew: T::AcvmType::default(),
            point_transition: false,
            pc: 0,
            round: 0,
            scalar_sum: T::AcvmType::default(),
            precompute_accumulator: T::AcvmPoint::<C>::default(),
            precompute_double: T::AcvmPoint::<C>::default(),
        }
    }
}
impl<C: CurveGroup<BaseField: PrimeField>, T: NoirWitnessExtensionProtocol<C::BaseField>> Clone
    for PointTablePrecomputationRow<C, T>
{
    fn clone(&self) -> Self {
        Self {
            s1: self.s1,
            s2: self.s2,
            s3: self.s3,
            s4: self.s4,
            s5: self.s5,
            s6: self.s6,
            s7: self.s7,
            s8: self.s8,
            skew: self.skew,
            point_transition: self.point_transition,
            pc: self.pc,
            round: self.round,
            scalar_sum: self.scalar_sum.clone(),
            precompute_accumulator: self.precompute_accumulator,
            precompute_double: self.precompute_double,
        }
    }
}

impl<C: HonkCurve<TranscriptFieldType>, T: NoirWitnessExtensionProtocol<C::BaseField>>
    PointTablePrecomputationRow<C, T>
{
    fn compute_rows(
        msms: &[ScalarMul<T, C>],
        driver: &mut T,
    ) -> eyre::Result<Vec<PointTablePrecomputationRow<C, T>>> {
        let num_rows_per_scalar = NUM_WNAF_DIGITS_PER_SCALAR / WNAF_DIGITS_PER_ROW;
        let num_precompute_rows = num_rows_per_scalar * msms.len() + 1;
        let mut precompute_state =
            vec![PointTablePrecomputationRow::<C, T>::default(); num_precompute_rows];

        // Start with an empty row (shiftable polynomials must have 0 as the first coefficient)
        precompute_state[0] = PointTablePrecomputationRow::<C, T>::default();

        // current impl doesn't work if not 4
        assert_eq!(WNAF_DIGITS_PER_ROW, 4);
        let mut wnaf_digits = Vec::with_capacity(msms.len() * NUM_WNAF_DIGITS_PER_SCALAR);
        for msm in msms.iter() {
            wnaf_digits.extend(msm.wnaf_digits);
        }

        let precomputed = driver.compute_rows(&wnaf_digits)?;
        let mut index = 0;
        for (j, entry) in msms.iter().enumerate() {
            let mut scalar_sum = T::AcvmType::default();

            for i in 0..num_rows_per_scalar {
                let mut row = PointTablePrecomputationRow {
                    s1: precomputed[index].0[0],
                    s2: precomputed[index].0[1],
                    s3: precomputed[index].0[2],
                    s4: precomputed[index].0[3],
                    s5: precomputed[index].0[4],
                    s6: precomputed[index].0[5],
                    s7: precomputed[index].0[6],
                    s8: precomputed[index].0[7],
                    ..Default::default()
                };

                // TODO FLORIN: Maybe do this already in the garbled circuit

                // row.s1 = slice0base2 >> 2;
                // row.s2 = slice0base2 & 3;
                // row.s3 = slice1base2 >> 2;
                // row.s4 = slice1base2 & 3;
                // row.s5 = slice2base2 >> 2;
                // row.s6 = slice2base2 & 3;
                // row.s7 = slice3base2 >> 2;
                // row.s8 = slice3base2 & 3;

                // let slice0 = slices[i * WNAF_DIGITS_PER_ROW];
                // let slice1 = slices[i * WNAF_DIGITS_PER_ROW + 1];
                // let slice2 = slices[i * WNAF_DIGITS_PER_ROW + 2];
                // let slice3 = slices[i * WNAF_DIGITS_PER_ROW + 3];

                // let slice0base2 = (slice0 + 15) / 2;
                // let slice1base2 = (slice1 + 15) / 2;
                // let slice2base2 = (slice2 + 15) / 2;
                // let slice3base2 = (slice3 + 15) / 2;

                /*    let mut slice0base2 = slice0; // + 15) / 2;
                driver.add_assign_with_public(C::BaseField::from(15), &mut slice0base2);
                slice0base2 = driver.mul_with_public(
                    C::BaseField::from(2)
                        .inverse()
                        .expect("2 should be invertible"),
                    slice0base2,
                );
                let mut slice1base2 = slice1; // + 15) / 2;
                driver.add_assign_with_public(C::BaseField::from(15), &mut slice1base2);
                slice1base2 = driver.mul_with_public(
                    C::BaseField::from(2)
                        .inverse()
                        .expect("2 should be invertible"),
                    slice1base2,
                );
                let mut slice2base2 = slice2; // + 15) / 2;
                driver.add_assign_with_public(C::BaseField::from(15), &mut slice2base2);
                slice2base2 = driver.mul_with_public(
                    C::BaseField::from(2)
                        .inverse()
                        .expect("2 should be invertible"),
                    slice2base2,
                );
                let mut slice3base2 = slice3; // + 15) / 2;
                driver.add_assign_with_public(C::BaseField::from(15), &mut slice3base2);
                slice3base2 = driver.mul_with_public(
                    C::BaseField::from(2)
                        .inverse()
                        .expect("2 should be invertible"),
                    slice3base2,
                ); */

                // Convert into 2-bit chunks

                let last_row = i == num_rows_per_scalar - 1;
                row.skew = if last_row {
                    entry.wnaf_skew
                } else {
                    T::AcvmType::default()
                };
                row.scalar_sum = scalar_sum;

                // Ensure slice1 is positive for the first row of each scalar sum
                let row_chunk = precomputed[index].1; //slice3 + (slice2 << 4) + (slice1 << 8) + (slice0 << 12);
                let chunk_negative = precomputed[index].2;
                let truthy = driver.mul_with_public(-C::BaseField::one(), row_chunk);
                let summand = driver.cmux(chunk_negative, truthy, row_chunk)?;

                let factor = 1 << (NUM_WNAF_DIGIT_BITS * WNAF_DIGITS_PER_ROW);
                scalar_sum = driver.mul_with_public(C::BaseField::from(factor), scalar_sum);
                // scalar_sum <<= NUM_WNAF_DIGIT_BITS * WNAF_DIGITS_PER_ROW;
                driver.add_assign(&mut scalar_sum, summand);

                row.round = i as u32;
                row.point_transition = last_row;
                row.pc = entry.pc;

                // We don't do this assert here
                // if last_row {
                //     assert_eq!(
                //         scalar_sum.clone() - BigUint::from(entry.wnaf_skew as u64),
                //         entry.scalar
                //     );
                // }

                row.precompute_double = entry.precomputed_table[POINT_TABLE_SIZE].to_owned();
                // fill accumulator in reverse order i.e. first row = 15[P], then 13[P], ..., 1[P]
                row.precompute_accumulator =
                    entry.precomputed_table[POINT_TABLE_SIZE - 1 - i].to_owned();
                precompute_state[j * num_rows_per_scalar + i + 1] = row;
                index += 1;
            }
        }
        // precompute_state
        todo!("Only after scalarmul is implemented")
    }
}

// pub fn construct_from_builder<
//     C: HonkCurve<TranscriptFieldType>,
//     T: NoirUltraHonkProver<C::CycleGroup>,
//     A: NoirUltraHonkProver<C>,
//     N: Network,
// >(
//     op_queue: &mut CoECCOpQueue<T, C::CycleGroup>,
//     net: &N,
//     state_: &mut T::State,
// ) -> eyre::Result<Polynomials<A::ArithmeticShare, C::ScalarField, ECCVMFlavour>>
// where
//     A::ArithmeticShare: From<T::ArithmeticShare>,
// {
//     let eccvm_ops = op_queue.get_eccvm_ops().to_vec();
//     let number_of_muls = op_queue.get_number_of_muls();
//     let transcript_rows = compute_rows::<C::CycleGroup, T, N>(&eccvm_ops, number_of_muls)
//         .expect("Failed to compute transcript rows");
//     let msms = op_queue.get_msms();
//     let point_table_rows = PointTablePrecomputationRow::<C::CycleGroup, T>::compute_rows(
//         &msms.iter().flat_map(|msm| msm.clone()).collect::<Vec<_>>(),
//     );
//     let result = MSMRow::<C::CycleGroup, T>::compute_rows_msms(
//         &msms,
//         number_of_muls,
//         op_queue.get_num_msm_rows(),
//         net,
//         state_,
//     );

//     todo!()
// }
