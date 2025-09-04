use ark_ec::AffineRepr;
use ark_ec::CurveGroup;
use ark_ff::Field;
use ark_ff::{One, PrimeField, Zero};
use co_acvm::mpc::NoirWitnessExtensionProtocol;
use co_builder::TranscriptFieldType;
use co_builder::prelude::HonkCurve;
use co_builder::prelude::offset_generator;
use common::{mpc::NoirUltraHonkProver, shared_polynomial::SharedPolynomial};
use goblin::{
    ADDITIONS_PER_ROW, NUM_WNAF_DIGIT_BITS, NUM_WNAF_DIGITS_PER_SCALAR, POINT_TABLE_SIZE,
    WNAF_MASK,
    prelude::{EccOpCode, EccOpsTable},
};
use mpc_core::MpcState;
use mpc_net::Network;
use num_bigint::BigUint;
use std::array;

pub(crate) const TABLE_WIDTH: usize = 4; // dictated by the number of wires in the Ultra arithmetization
pub(crate) const NUM_ROWS_PER_OP: usize = 2; // A single ECC op is split across two width-4 rows

pub(crate) type CoEccvmOpsTable<T, C> = EccOpsTable<CoVMOperation<T, C>>;

pub(crate) struct CoUltraEccOpsTable<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub(crate) table: EccOpsTable<CoUltraOp<T, C>>,
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>>
    CoUltraEccOpsTable<T, C>
{
    pub fn ultra_table_size(&self) -> usize {
        self.table.size() * NUM_ROWS_PER_OP
    }

    pub fn current_ultra_subtable_size(&self) -> usize {
        self.table.get()[0].len() * NUM_ROWS_PER_OP
    }

    pub fn previous_ultra_table_size(&self) -> usize {
        self.ultra_table_size() - self.current_ultra_subtable_size()
    }

    pub fn create_new_subtable(&mut self, size_hint: usize) {
        self.table.create_new_subtable(size_hint);
    }

    // pub fn construct_table_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     let poly_size = self.ultra_table_size();
    //     let subtable_start_idx = 0; // include all subtables
    //     let subtable_end_idx = self.table.num_subtables();

    //     self.construct_column_polynomials_from_subtables(
    //         poly_size,
    //         subtable_start_idx,
    //         subtable_end_idx,
    //         id,
    //     )
    // }

    // pub fn construct_previous_table_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     let poly_size = self.previous_ultra_table_size();
    //     let subtable_start_idx = 1; // exclude the 0th subtable
    //     let subtable_end_idx = self.table.num_subtables();

    //     self.construct_column_polynomials_from_subtables(
    //         poly_size,
    //         subtable_start_idx,
    //         subtable_end_idx,
    //         id,
    //     )
    // }

    // pub fn construct_current_ultra_ops_subtable_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     let poly_size = self.current_ultra_subtable_size();
    //     let subtable_start_idx = 0;
    //     let subtable_end_idx = 1; // include only the 0th subtable

    //     self.construct_column_polynomials_from_subtables(
    //         poly_size,
    //         subtable_start_idx,
    //         subtable_end_idx,
    //         id,
    //     )
    // }

    // fn construct_column_polynomials_from_subtables(
    //     &self,
    //     poly_size: usize,
    //     subtable_start_idx: usize,
    //     subtable_end_idx: usize,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     let mut column_polynomials: [SharedPolynomial<T, C>; TABLE_WIDTH] =
    //         array::from_fn(|_| SharedPolynomial::new_zero(poly_size));

    //     let mut i = 0;
    //     for subtable_idx in subtable_start_idx..subtable_end_idx {
    //         let subtable = &self.table.get()[subtable_idx];
    //         for op in subtable {
    //             column_polynomials[0][i] =
    //                 T::promote_to_trivial_share(id, C::ScalarField::from(op.op_code.value()));
    //             column_polynomials[1][i] = op.x_lo;
    //             column_polynomials[2][i] = op.x_hi;
    //             column_polynomials[3][i] = op.y_lo;
    //             i += 1;
    //             column_polynomials[0][i] = T::ArithmeticShare::default(); // only the first 'op' field is utilized
    //             column_polynomials[1][i] = op.y_hi;
    //             column_polynomials[2][i] = op.z_1;
    //             column_polynomials[3][i] = op.z_2;
    //             i += 1;
    //         }
    //     }
    //     column_polynomials
    // }
}

#[derive(Default)]
pub struct CoVMOperation<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub op_code: EccOpCode,
    pub base_point: T::AcvmPoint<C>,
    pub z1: T::AcvmType, //TODO FLORIN: I think this does not have to be a binary share (It is a uint256 in bb)
    pub z2: T::AcvmType, //TODO FLORIN: I think this does not have to be a binary share (It is a uint256 in bb)
    pub mul_scalar_full: T::OtherAcvmType<C>,
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Clone
    for CoVMOperation<T, C>
{
    fn clone(&self) -> Self {
        Self {
            op_code: self.op_code.clone(),
            base_point: self.base_point,
            z1: self.z1,
            z2: self.z2,
            mul_scalar_full: self.mul_scalar_full,
        }
    }
}

pub struct CoUltraOp<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub op_code: EccOpCode,
    pub x_lo: T::AcvmType,
    pub x_hi: T::AcvmType,
    pub y_lo: T::AcvmType,
    pub y_hi: T::AcvmType,
    pub z_1: T::AcvmType,
    pub z_2: T::AcvmType,
    pub return_is_infinity: T::AcvmType,
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Clone
    for CoUltraOp<T, C>
{
    fn clone(&self) -> Self {
        Self {
            op_code: self.op_code.clone(),
            x_lo: self.x_lo,
            x_hi: self.x_hi,
            y_lo: self.y_lo,
            y_hi: self.y_hi,
            z_1: self.z_1,
            z_2: self.z_2,
            return_is_infinity: self.return_is_infinity,
        }
    }
}

#[derive(Debug)]
pub struct CoEccvmRowTracker<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub cached_num_muls: T::AcvmType,
    pub cached_active_msm_count: T::AcvmType,
    pub num_transcript_rows: u32,
    pub num_precompute_table_rows: T::AcvmType,
    pub num_msm_rows: T::AcvmType,
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Default
    for CoEccvmRowTracker<T, C>
{
    fn default() -> Self {
        Self {
            cached_num_muls: T::AcvmType::default(),
            cached_active_msm_count: T::AcvmType::default(),
            num_transcript_rows: 0,
            num_precompute_table_rows: T::AcvmType::default(),
            num_msm_rows: T::AcvmType::default(),
        }
    }
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>>
    CoEccvmRowTracker<T, C>
{
    pub fn get_number_of_muls(&self, driver: &mut T) -> T::AcvmType {
        driver.add(
            self.cached_num_muls.to_owned(),
            self.cached_active_msm_count.to_owned(),
        )
    }

    pub fn num_eccvm_msm_rows(
        msm_size: T::ArithmeticShare,
        // id: <T::State as MpcState>::PartyID,
    ) -> T::ArithmeticShare {
        // let rows_per_wnaf_digit = (msm_size / ADDITIONS_PER_ROW)
        //     + if msm_size % ADDITIONS_PER_ROW != 0 {
        //         1
        //     } else {
        //         0
        //     };
        // let num_rows_for_all_rounds = (NUM_WNAF_DIGITS_PER_SCALAR + 1) * rows_per_wnaf_digit;
        // let num_double_rounds = NUM_WNAF_DIGITS_PER_SCALAR - 1;
        // T::add_with_public(
        //     C::ScalarField::from(num_double_rounds as u32),
        //     num_rows_for_all_rounds,
        //     id,
        // )
        todo!()
    }

    pub fn num_eccvm_msm_rows_public(msm_size: usize) -> u32 {
        let rows_per_wnaf_digit = (msm_size / ADDITIONS_PER_ROW)
            + if msm_size % ADDITIONS_PER_ROW != 0 {
                1
            } else {
                0
            };
        let num_rows_for_all_rounds = (NUM_WNAF_DIGITS_PER_SCALAR + 1) * rows_per_wnaf_digit;
        let num_double_rounds = NUM_WNAF_DIGITS_PER_SCALAR - 1;
        (num_rows_for_all_rounds + num_double_rounds) as u32
    }

    pub fn get_num_msm_rows(&self) -> T::AcvmType {
        // let mut msm_rows = self.num_msm_rows as usize + 2;
        // if self.cached_active_msm_count > 0 {
        //     msm_rows += Self::num_eccvm_msm_rows(self.cached_active_msm_count as usize) as usize;
        // }
        // msm_rows
        todo!()
    }

    pub fn get_num_rows(&self) -> T::AcvmType {
        // let transcript_rows = self.num_transcript_rows as usize + 2;
        // let mut msm_rows = self.num_msm_rows as usize + 2;
        // let mut precompute_rows = self.num_precompute_table_rows as usize + 1;
        // if self.cached_active_msm_count > 0 {
        //     msm_rows += Self::num_eccvm_msm_rows(self.cached_active_msm_count as usize) as usize;
        //     precompute_rows += Self::get_precompute_table_row_count_for_single_msm(
        //         self.cached_active_msm_count as usize,
        //     ) as usize;
        // }
        // std::cmp::max(transcript_rows, std::cmp::max(msm_rows, precompute_rows))
        todo!()
    }

    pub fn get_precompute_table_row_count_for_single_msm(msm_count: T::AcvmType) -> T::AcvmType {
        // let num_precompute_rows_per_scalar = NUM_WNAF_DIGITS_PER_SCALAR / WNAF_DIGITS_PER_ROW;
        // (msm_count * num_precompute_rows_per_scalar) as u32
        todo!()
    }
}

pub struct CoECCOpQueue<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub(crate) eccvm_ops_table: CoEccvmOpsTable<T, C>,
    pub(crate) ultra_ops_table: CoUltraEccOpsTable<T, C>,
    pub(crate) accumulator: T::AcvmPoint<C>,
    pub(crate) eccvm_ops_reconstructed: Vec<CoVMOperation<T, C>>,
    pub ultra_ops_reconstructed: Vec<CoUltraOp<T, C>>,
    pub(crate) eccvm_row_tracker: CoEccvmRowTracker<T, C>,
}

pub struct test<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>>
{
    pub(crate) accumulator: T::AcvmPoint<C>,
    pub test: T::AcvmType,
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>>
    CoECCOpQueue<T, C>
{
    // Initialize a new subtable of ECCVM ops and Ultra ops corresponding to an individual circuit
    pub fn initialize_new_subtable(&mut self) {
        self.eccvm_ops_table.create_new_subtable(0);
        self.ultra_ops_table.create_new_subtable(0);
    }

    // // Construct polynomials corresponding to the columns of the full aggregate ultra ecc ops table
    // pub fn construct_ultra_ops_table_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     self.ultra_ops_table.construct_table_columns(id)
    // }

    // // Construct polys corresponding to the columns of the aggregate ultra ops table, excluding the most recent subtable
    // pub fn construct_previous_ultra_ops_table_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     self.ultra_ops_table.construct_previous_table_columns(id)
    // }

    // // Construct polynomials corresponding to the columns of the current subtable of ultra ecc ops
    // pub fn construct_current_ultra_ops_subtable_columns(
    //     &self,
    //     id: <T::State as MpcState>::PartyID,
    // ) -> [SharedPolynomial<T, C>; TABLE_WIDTH] {
    //     self.ultra_ops_table
    //         .construct_current_ultra_ops_subtable_columns(id)
    // }
    // Reconstruct the full table of eccvm ops in contiguous memory from the independent subtables
    pub fn construct_full_eccvm_ops_table(&mut self) {
        self.eccvm_ops_reconstructed = self.eccvm_ops_table.get_reconstructed();
    }

    // Reconstruct the full table of ultra ops in contiguous memory from the independent subtables
    pub fn construct_full_ultra_ops_table(&mut self) {
        self.ultra_ops_reconstructed = self.ultra_ops_table.table.get_reconstructed();
    }

    pub fn get_ultra_ops_table_num_rows(&self) -> usize {
        self.ultra_ops_table.ultra_table_size()
    }

    pub fn get_current_ultra_ops_subtable_num_rows(&self) -> usize {
        self.ultra_ops_table.current_ultra_subtable_size()
    }
    // Get the full table of ECCVM ops in contiguous memory; construct it if it has not been constructed already
    pub fn get_eccvm_ops(&mut self) -> &Vec<CoVMOperation<T, C>> {
        if self.eccvm_ops_reconstructed.is_empty() {
            self.construct_full_eccvm_ops_table();
        }
        &self.eccvm_ops_reconstructed
    }

    /**
     * @brief Get the number of rows in the 'msm' column section, for all msms in the circuit
     */
    pub fn get_num_msm_rows(&self) -> T::AcvmType {
        self.eccvm_row_tracker.get_num_msm_rows()
    }

    /**
     * @brief Get the number of rows for the current ECCVM circuit
     */
    pub fn get_num_rows(&self) -> T::AcvmType {
        self.eccvm_row_tracker.get_num_rows()
    }

    /**
     * @brief get number of muls for the current ECCVM circuit
     */
    pub fn get_number_of_muls(&self, driver: &mut T) -> T::AcvmType {
        self.eccvm_row_tracker.get_number_of_muls(driver)
    }

    /**
     * @brief A fuzzing only method for setting eccvm ops directly
     *
     */
    pub fn set_eccvm_ops_for_fuzzing(&mut self, eccvm_ops_in: Vec<CoVMOperation<T, C>>) {
        self.eccvm_ops_reconstructed = eccvm_ops_in;
    }

    pub fn get_accumulator(&self) -> &T::AcvmPoint<C> {
        &self.accumulator
    }
}

#[expect(dead_code)]
//TODO FLORIN: Rename?
pub(crate) struct ScalarMul<
    T: NoirWitnessExtensionProtocol<C::BaseField>,
    C: CurveGroup<BaseField: PrimeField>,
> {
    pub(crate) pc: T::AcvmType,
    pub(crate) scalar: T::AcvmType,
    pub(crate) base_point: T::AcvmPoint<C>,
    pub(crate) wnaf_digits: [T::AcvmType; NUM_WNAF_DIGITS_PER_SCALAR],
    pub(crate) wnaf_skew: bool,
    // size bumped by 1 to record base_point.dbl()
    pub(crate) precomputed_table: [T::AcvmPoint<C>; POINT_TABLE_SIZE + 1],
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Default
    for ScalarMul<T, C>
{
    fn default() -> Self {
        Self {
            pc: T::AcvmType::default(),
            scalar: T::AcvmType::default(),
            base_point: T::AcvmPoint::<C>::default(),
            wnaf_digits: [T::AcvmType::default(); NUM_WNAF_DIGITS_PER_SCALAR],
            wnaf_skew: false,
            precomputed_table: [T::AcvmPoint::<C>::default(); POINT_TABLE_SIZE + 1],
        }
    }
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Clone
    for ScalarMul<T, C>
{
    fn clone(&self) -> Self {
        Self {
            pc: self.pc,
            scalar: self.scalar,
            base_point: self.base_point,
            wnaf_digits: self.wnaf_digits,
            wnaf_skew: self.wnaf_skew,
            precomputed_table: self.precomputed_table,
        }
    }
}

pub(crate) type Msm<C, T> = Vec<ScalarMul<T, C>>;

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: HonkCurve<TranscriptFieldType>>
    CoECCOpQueue<T, C>
{
    pub(crate) fn get_msms<N: Network>(&mut self, driver: &mut T) -> eyre::Result<Vec<Msm<C, T>>> {
        let num_muls = self.get_number_of_muls(driver);

        let compute_precomputed_table =
            |base_point: T::AcvmPoint<C>| -> [T::AcvmPoint<C>; POINT_TABLE_SIZE + 1] {
                let d2 = driver.scalar_mul_public_scalar(base_point, C::ScalarField::from(2u32));
                let mut table = [T::AcvmPoint::default(); POINT_TABLE_SIZE + 1];
                table[POINT_TABLE_SIZE] = d2.into();
                table[POINT_TABLE_SIZE / 2] = base_point;

                for i in 1..(POINT_TABLE_SIZE / 2) {
                    table[i + POINT_TABLE_SIZE / 2] =
                        driver.add_points(table[i + POINT_TABLE_SIZE / 2 - 1], d2);
                }

                for i in 0..(POINT_TABLE_SIZE / 2) {
                    table[i] = driver.scalar_mul_public_scalar(
                        table[POINT_TABLE_SIZE - 1 - i],
                        -C::ScalarField::one(),
                    );
                }

                //TODO FLORIN: Make this nicer
                let mut result = [T::AcvmPoint::default(); POINT_TABLE_SIZE + 1];
                for (i, point) in table.iter().enumerate() {
                    result[i] = *point;
                }
                result
            };

        let compute_wnaf_digits = |mut scalar: T::AcvmType| -> [i32; NUM_WNAF_DIGITS_PER_SCALAR] {
            let mut output = [0; NUM_WNAF_DIGITS_PER_SCALAR];
            let mut previous_slice = 0;

            for i in 0..NUM_WNAF_DIGITS_PER_SCALAR {
                let raw_slice = &scalar & BigUint::from(WNAF_MASK);
                let is_even = (&raw_slice & BigUint::one()) == BigUint::zero();
                let mut wnaf_slice = if let Some(&digit) = raw_slice.to_u32_digits().first() {
                    digit as i32
                } else {
                    0
                };

                if i == 0 && is_even {
                    wnaf_slice += 1;
                } else if is_even {
                    const BORROW_CONSTANT: i32 = 1 << NUM_WNAF_DIGIT_BITS;
                    previous_slice -= BORROW_CONSTANT;
                    wnaf_slice += 1;
                }

                if i > 0 {
                    output[NUM_WNAF_DIGITS_PER_SCALAR - i] = previous_slice;
                }
                previous_slice = wnaf_slice;

                scalar >>= NUM_WNAF_DIGIT_BITS;
            }

            assert!(scalar.is_zero());
            output[0] = previous_slice;

            output
        };

        let mut msm_count = T::AcvmType::default();
        let mut active_mul_count = T::AcvmType::default();
        let mut msm_opqueue_index = Vec::new();
        let mut msm_mul_index = Vec::new();
        let mut msm_sizes = Vec::new();

        let eccvm_ops = self.get_eccvm_ops();
        //TODO FLORIN: Do only the ones needed
        let mut z1s = Vec::with_capacity(eccvm_ops.len());
        let mut z2s = Vec::with_capacity(eccvm_ops.len());
        let mut base_points = Vec::with_capacity(eccvm_ops.len());
        for op in eccvm_ops.iter() {
            z1s.push(op.z1);
            z2s.push(op.z2);
            base_points.push(op.base_point);
        }
        //TODO FLORIN Optimize this
        let is_zeros = driver.is_zero_many(&[z1s.as_slice(), z2s.as_slice()].concat())?;
        let is_zero_z1s = &is_zeros[0..z1s.len()];
        let is_zero_z2s = &is_zeros[z1s.len()..];
        let scale = driver.scale_many(is_zero_z1s, -C::BaseField::one());
        let inv_is_zero_z1s = driver.add_scalar(&scale, C::BaseField::one());
        let scale = driver.scale_many(is_zero_z2s, -C::BaseField::one());
        let inv_is_zero_z2s = driver.add_scalar(&scale, C::BaseField::one());
        let is_zero_base_points = driver.point_is_zero_many(&base_points)?;
        let scale = driver.scale_many(&is_zero_base_points, -C::BaseField::one());
        let inv_is_zero_base_points = driver.add_scalar(&scale, C::BaseField::one());
        let added = driver.add_many(&inv_is_zero_z1s, &inv_is_zero_z2s);
        let mul = driver.mul_many(&added, &inv_is_zero_base_points)?;
        let op_indices: Vec<C::BaseField> = (0..eccvm_ops.len())
            .map(|i| C::BaseField::from(i as u32))
            .collect();
        let op_indices = driver.mul_with_public_many(&op_indices, &mul);

        //TACEO TODO: This has to be optimized
        for (op_idx, op) in eccvm_ops.iter().enumerate() {
            if op.op_code.mul {
                // if (op.z1 != BigUint::zero() || op.z2 != BigUint::zero())
                //     && !op.base_point.is_zero()
                // {
                //     msm_mul_index.push((msm_count, active_mul_count));
                // }
                let mul = driver.mul_many(&[msm_count, active_mul_count], &[mul[op_idx]; 2])?;
                msm_mul_index.push((mul[0], mul[1]));
                msm_opqueue_index.push(op_indices[op_idx]);
                // active_mul_count +=
                //     (op.z1 != BigUint::zero()) as usize + (op.z2 != BigUint::zero()) as usize;
                driver.add_assign(&mut active_mul_count, mul[op_idx]);
            } else {
                //if active_mul_count > 0 {
                let is_zero = driver.is_zero_many(&[active_mul_count])?[0];
                let mut inv_is_zero = driver.mul_with_public(-C::BaseField::one(), is_zero);
                driver.add_assign_with_public(C::BaseField::one(), &mut inv_is_zero);
                let mul = driver.mul_many(&[active_mul_count], &[inv_is_zero])?;

                msm_sizes.push(mul[0]);

                driver.add_assign(&mut msm_count, inv_is_zero);
                // msm_count += 1;
                active_mul_count = is_zero;
            }
        }

        let active_mul_count_is_zero = driver.is_zero_many(&[active_mul_count])?[0];
        let mul = driver.mul_with_public(-C::BaseField::one(), active_mul_count_is_zero);
        let inv_is_zero = driver.add(T::AcvmType::from(C::BaseField::one()), mul);
        let mul = driver.mul_many(&[active_mul_count], &[inv_is_zero])?;
        if eccvm_ops.last().is_some_and(|op| op.op_code.mul) {
            msm_sizes.push(mul[0]);
            driver.add_assign(&mut msm_count, inv_is_zero);
        }

        let mut result: Vec<Msm<C, T>> = Vec::with_capacity(msm_count);
        for size in &msm_sizes {
            result.push(vec![ScalarMul::default(); *size]);
        }

        msm_opqueue_index
            .iter()
            .enumerate()
            .for_each(|(i, &op_idx)| {
                let op = &eccvm_ops[op_idx];
                let (msm_index, mut mul_index) = msm_mul_index[i];

                if op.z1 != BigUint::zero() && !op.base_point.is_zero() {
                    result[msm_index][mul_index] = ScalarMul {
                        pc: 0,
                        scalar: op.z1.clone(),
                        base_point: op.base_point,
                        wnaf_digits: compute_wnaf_digits(op.z1.clone()),
                        wnaf_skew: (op.z1.clone() & BigUint::from(1u32)) == BigUint::zero(),
                        precomputed_table: compute_precomputed_table(op.base_point),
                    };
                    mul_index += 1;
                }

                if op.z2 != BigUint::zero() && !op.base_point.is_zero() {
                    let endo_point = C::g1_affine_from_xy(
                        op.base_point.x().expect("BasePoint should not be zero")
                            * C::get_cube_root_of_unity(),
                        -op.base_point.y().expect("BasePoint should not be zero"),
                    );
                    result[msm_index][mul_index] = ScalarMul {
                        pc: 0,
                        scalar: op.z2.clone(),
                        base_point: endo_point,
                        wnaf_digits: compute_wnaf_digits(op.z2.clone()),
                        wnaf_skew: (op.z2.clone() & BigUint::from(1u32)) == BigUint::zero(),
                        precomputed_table: compute_precomputed_table(endo_point),
                    };
                }
            });

        let mut pc = num_muls;
        for msm in &mut result {
            for mul in msm {
                mul.pc = pc;
                pc -= 1;
            }
        }

        Ok(result)
    }
}
pub(crate) struct AddState<
    C: CurveGroup<BaseField: PrimeField>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
> {
    pub add: bool,
    pub slice: T::AcvmType,
    pub point: T::AcvmPoint<C>,
    pub lambda: T::AcvmType,
    pub collision_inverse: T::AcvmType,
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Default
    for AddState<C, T>
{
    fn default() -> Self {
        Self {
            add: false,
            slice: T::AcvmType::default(),
            point: T::AcvmPoint::<C>::default(),
            lambda: T::AcvmType::default(),
            collision_inverse: T::AcvmType::default(),
        }
    }
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Clone
    for AddState<C, T>
{
    fn clone(&self) -> Self {
        Self {
            add: self.add,
            slice: self.slice,
            point: self.point,
            lambda: self.lambda,
            collision_inverse: self.collision_inverse,
        }
    }
}

pub(crate) struct MSMRow<
    C: CurveGroup<BaseField: PrimeField>,
    T: NoirWitnessExtensionProtocol<C::BaseField>,
> {
    // Counter over all half-length scalar muls used to compute the required MSMs
    pub(crate) pc: T::AcvmType,
    // The number of points that will be scaled and summed
    pub(crate) msm_size: u32,
    pub(crate) msm_count: u32,
    pub(crate) msm_round: u32,
    pub(crate) msm_transition: bool,
    pub(crate) q_add: bool,
    pub(crate) q_double: bool,
    pub(crate) q_skew: bool,
    pub(crate) add_state: [AddState<C, T>; 4],
    pub(crate) accumulator_x: T::AcvmType,
    pub(crate) accumulator_y: T::AcvmType,
    phantom: std::marker::PhantomData<T>,
}
impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Default
    for MSMRow<C, T>
{
    fn default() -> Self {
        Self {
            pc: T::AcvmType::default(),
            msm_size: 0,
            msm_count: 0,
            msm_round: 0,
            msm_transition: false,
            q_add: false,
            q_double: false,
            q_skew: false,
            add_state: array::from_fn(|_| AddState::default()),
            accumulator_x: T::AcvmType::default(),
            accumulator_y: T::AcvmType::default(),
            phantom: std::marker::PhantomData,
        }
    }
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: CurveGroup<BaseField: PrimeField>> Clone
    for MSMRow<C, T>
{
    fn clone(&self) -> Self {
        Self {
            pc: self.pc,
            msm_size: self.msm_size,
            msm_count: self.msm_count,
            msm_round: self.msm_round,
            msm_transition: self.msm_transition,
            q_add: self.q_add,
            q_double: self.q_double,
            q_skew: self.q_skew,
            add_state: self.add_state.clone(),
            accumulator_x: self.accumulator_x,
            accumulator_y: self.accumulator_y,
            phantom: std::marker::PhantomData,
        }
    }
}

impl<T: NoirWitnessExtensionProtocol<C::BaseField>, C: HonkCurve<TranscriptFieldType>>
    MSMRow<C, T>
{
    pub(crate) fn compute_rows_msms<N: Network>(
        msms: &[Msm<C, T>],
        total_number_of_muls: T::AcvmType,
        num_msm_rows: usize,
        driver: &mut T,
    ) -> eyre::Result<(Vec<Self>, [Vec<T::AcvmType>; 2])> {
        let num_rows_in_read_counts_table = driver.mul_with_public(
            C::BaseField::from((POINT_TABLE_SIZE / 2) as u32),
            total_number_of_muls,
        );

        // (total_number_of_muls as usize) * (POINT_TABLE_SIZE / 2);
        let mut point_table_read_counts = [Vec::new(), Vec::new()];

        let mut update_read_count = |point_idx: T::AcvmType,
                                     slice: T::AcvmType,
                                     driver: &mut T|
         -> eyre::Result<()> {
            let row_index_offset = driver.mul_with_public(C::BaseField::from(8), point_idx);
            let digit_is_negative = driver.lt(slice, T::AcvmType::default())?;
            let mul = driver.mul_with_public(-C::BaseField::one(), digit_is_negative);
            let inv_digit_is_negative = driver.add(T::AcvmType::from(C::BaseField::one()), mul);
            let mut relative_row_idx = slice; //((slice + 15) / 2) as usize;
            driver.add_assign_with_public(C::BaseField::from(15), &mut relative_row_idx);
            relative_row_idx = driver.mul_with_public(
                C::BaseField::from(2)
                    .inverse()
                    .expect("2 should have an inverse..."),
                relative_row_idx,
            );
            //TODO FLORIN DONT INITIALIZE THIS EVERY TIME
            let mut lut_1 = T::init_lut_by_acvm_type(driver, point_table_read_counts[0].clone());
            let mut lut_2 = T::init_lut_by_acvm_type(driver, point_table_read_counts[1].clone());
            let first_index = driver.add(row_index_offset, relative_row_idx);
            let mut second_index = driver.sub(row_index_offset, relative_row_idx);
            driver.add_assign_with_public(C::BaseField::from(15), &mut second_index);
            let mut first_value = driver.read_lut_by_acvm_type(first_index, &lut_1)?;
            let mut second_value = driver.read_lut_by_acvm_type(second_index, &lut_2)?;
            driver.add_assign(&mut first_value, digit_is_negative);
            driver.add_assign(&mut second_value, inv_digit_is_negative);
            driver.write_lut_by_acvm_type(first_index, first_value, &mut lut_1)?;
            driver.write_lut_by_acvm_type(second_index, second_value, &mut lut_2)?;
            point_table_read_counts[0] = T::get_shared_lut(&lut_1)?.to_vec();
            point_table_read_counts[1] = T::get_shared_lut(&lut_2)?.to_vec();

            Ok(())
        };

        let mut msm_row_counts = Vec::with_capacity(msms.len() + 1);
        msm_row_counts.push(1);

        let mut pc_values = Vec::with_capacity(msms.len() + 1);
        pc_values.push(total_number_of_muls);

        for msm in msms {
            let num_rows_required = CoEccvmRowTracker::<T, C>::num_eccvm_msm_rows_public(msm.len());
            msm_row_counts.push(
                msm_row_counts
                    .last()
                    .expect("msm_row_counts should not be empty")
                    + num_rows_required as usize,
            );
            let add = driver.add(
                *pc_values.last().expect("pc_values should not be empty"),
                T::AcvmType::from(-C::BaseField::from(msm.len() as u32)),
            );
            pc_values.push(add);
        }

        let mut msm_rows = vec![MSMRow::default(); num_msm_rows];
        msm_rows[0] = MSMRow::default();

        for (msm_idx, msm) in msms.iter().enumerate() {
            for digit_idx in 0..NUM_WNAF_DIGITS_PER_SCALAR {
                let pc = pc_values[msm_idx];
                let msm_size = msm.len();
                let num_rows_per_digit = (msm_size / ADDITIONS_PER_ROW)
                    + if msm_size % ADDITIONS_PER_ROW != 0 {
                        1
                    } else {
                        0
                    };

                for relative_row_idx in 0..num_rows_per_digit {
                    let num_points_in_row = if (relative_row_idx + 1) * ADDITIONS_PER_ROW > msm_size
                    {
                        msm_size % ADDITIONS_PER_ROW
                    } else {
                        ADDITIONS_PER_ROW
                    };
                    let offset = relative_row_idx * ADDITIONS_PER_ROW;

                    for relative_point_idx in 0..ADDITIONS_PER_ROW {
                        let point_idx = offset + relative_point_idx;
                        let add = num_points_in_row > relative_point_idx;
                        if add {
                            let slice = msm[point_idx].wnaf_digits[digit_idx];
                            update_read_count(
                                driver.add(
                                    total_number_of_muls,
                                    driver.sub(
                                        T::AcvmType::from(C::BaseField::from(point_idx as u32)),
                                        pc,
                                    ),
                                ),
                                slice,
                                driver,
                            );
                        }
                    }
                }

                if digit_idx == NUM_WNAF_DIGITS_PER_SCALAR - 1 {
                    for row_idx in 0..num_rows_per_digit {
                        let num_points_in_row = if (row_idx + 1) * ADDITIONS_PER_ROW > msm_size {
                            msm_size % ADDITIONS_PER_ROW
                        } else {
                            ADDITIONS_PER_ROW
                        };
                        let offset = row_idx * ADDITIONS_PER_ROW;

                        for relative_point_idx in 0..ADDITIONS_PER_ROW {
                            let add = num_points_in_row > relative_point_idx;
                            let point_idx = offset + relative_point_idx;
                            if add {
                                let slice = if msm[point_idx].wnaf_skew { -1 } else { -15 };
                                let sub = driver.sub(total_number_of_muls, pc);
                                update_read_count(
                                    driver.add(
                                        sub,
                                        T::AcvmType::from(C::BaseField::from(point_idx as u32)),
                                    ),
                                    T::AcvmType::from(C::BaseField::from(slice)),
                                    driver,
                                );
                            }
                        }
                    }
                }
            }
        }

        // The execution trace data for the MSM columns requires knowledge of intermediate values from *affine* point
        // addition. The naive solution to compute this data requires 2 field inversions per in-circuit group addition
        // evaluation. This is bad! To avoid this, we split the witness computation algorithm into 3 steps.
        //   Step 1: compute the execution trace group operations in *projective* coordinates
        //   Step 2: use batch inversion trick to convert all points into affine coordinates
        //   Step 3: populate the full execution trace, including the intermediate values from affine group operations
        // This section sets up the data structures we need to store all intermediate ECC operations in projective form
        let num_point_adds_and_doubles = (num_msm_rows - 2) * 4;
        let num_accumulators = num_msm_rows - 1;
        // In what fallows, either p1 + p2 = p3, or p1.dbl() = p3
        // We create 1 vector to store the entire point trace. We split into multiple containers using std::span
        // (we want 1 vector object to more efficiently batch normalize points)
        const NUM_POINTS_IN_ADDITION_RELATION: usize = 3;
        let num_points_to_normalize =
            (num_point_adds_and_doubles * NUM_POINTS_IN_ADDITION_RELATION) + num_accumulators;
        let mut p1_trace = vec![T::AcvmPoint::<C>::default(); num_point_adds_and_doubles];
        let mut p2_trace = vec![T::AcvmPoint::<C>::default(); num_point_adds_and_doubles];
        let mut p3_trace = vec![T::AcvmPoint::<C>::default(); num_point_adds_and_doubles];
        // operation_trace records whether an entry in the p1/p2/p3 trace represents a point addition or doubling
        let mut operation_trace = vec![false; num_point_adds_and_doubles];
        // accumulator_trace tracks the value of the ECCVM accumulator for each row
        let mut accumulator_trace = vec![T::AcvmPoint::<C>::default(); num_accumulators];

        // we start the accumulator at the offset generator point. This ensures we can support an MSM that produces a
        let offset_generator = T::AcvmPoint::from(offset_generator::<C>().into());
        accumulator_trace[0] = offset_generator;

        // AZTEC TODO(https://github.com/AztecProtocol/barretenberg/issues/973): Reinstate multitreading?
        // populate point trace, and the components of the MSM execution trace that do not relate to affine point
        // operations
        for msm_idx in 0..msms.len() {
            let mut accumulator = offset_generator;
            let msm = &msms[msm_idx];
            let mut msm_row_index = msm_row_counts[msm_idx];
            let msm_size = msm.len();
            let num_rows_per_digit = (msm_size / ADDITIONS_PER_ROW)
                + if msm_size % ADDITIONS_PER_ROW != 0 {
                    1
                } else {
                    0
                };
            let mut trace_index = (msm_row_counts[msm_idx] - 1) * 4;

            for digit_idx in 0..NUM_WNAF_DIGITS_PER_SCALAR {
                let pc = pc_values[msm_idx];
                for row_idx in 0..num_rows_per_digit {
                    let num_points_in_row = if (row_idx + 1) * ADDITIONS_PER_ROW > msm_size {
                        msm_size % ADDITIONS_PER_ROW
                    } else {
                        ADDITIONS_PER_ROW
                    };
                    let row = &mut msm_rows[msm_row_index];
                    let offset = row_idx * ADDITIONS_PER_ROW;
                    row.msm_transition = (digit_idx == 0) && (row_idx == 0);

                    for point_idx in 0..ADDITIONS_PER_ROW {
                        let add_state = &mut row.add_state[point_idx];
                        add_state.add = num_points_in_row > point_idx;
                        let slice = if add_state.add {
                            msm[offset + point_idx].wnaf_digits[digit_idx]
                        } else {
                            T::AcvmType::default()
                        };
                        // In the MSM columns in the ECCVM circuit, we can add up to 4 points per row.
                        // if `row.add_state[point_idx].add = true`, this indicates that we want to add the
                        // `point_idx`'th point in the MSM columns into the MSM accumulator.
                        // `add_state.slice` = A 4-bit WNAF slice of the scalar multiplier associated with the point we are adding
                        // (the specific slice chosen depends on the value of msm_round).
                        // (WNAF = windowed-non-adjacent-form. Value range is `-15, -13, ..., 15`).
                        // If `add_state.add = true`, we want `add_state.slice` to be the *compressed*
                        // form of the WNAF slice value. (compressed = no gaps in the value range. i.e. -15,
                        // -13, ..., 15 maps to 0, ..., 15).
                        add_state.slice = if add_state.add {
                            let mut tmp = slice; //((slice + 15) / 2) as usize;
                            driver.add_assign_with_public(C::BaseField::from(15), &mut tmp);
                            tmp = driver.mul_with_public(
                                C::BaseField::from(2)
                                    .inverse()
                                    .expect("2 should have an inverse..."),
                                tmp,
                            );
                            tmp
                        } else {
                            T::AcvmType::default()
                        };
                        add_state.point = if add_state.add {
                            let lut = driver.init_lut_by_acvm_point(
                                msm[offset + point_idx].precomputed_table.to_vec(),
                            );
                            driver.read_lut_by_acvm_point(add_state.slice, &lut)?
                        } else {
                            T::AcvmPoint::<C>::default()
                        };

                        let p1 = accumulator;
                        let p2 = add_state.point;
                        accumulator = if add_state.add {
                            driver.add_points(accumulator, add_state.point)
                        } else {
                            p1
                        };
                        p1_trace[trace_index] = p1;
                        p2_trace[trace_index] = p2;
                        p3_trace[trace_index] = accumulator;
                        operation_trace[trace_index] = false;
                        trace_index += 1;
                    }
                    accumulator_trace[msm_row_index] = accumulator;
                    row.q_add = true;
                    row.q_double = false;
                    row.q_skew = false;
                    row.msm_round = digit_idx as u32;
                    row.msm_size = msm_size as u32;
                    row.msm_count = offset as u32;
                    row.pc = pc;
                    msm_row_index += 1;
                }
                // doubling
                if digit_idx < NUM_WNAF_DIGITS_PER_SCALAR - 1 {
                    let row = &mut msm_rows[msm_row_index];
                    row.msm_transition = false;
                    row.msm_round = (digit_idx + 1) as u32;
                    row.msm_size = msm_size as u32;
                    row.msm_count = 0_u32;
                    row.q_add = false;
                    row.q_double = true;
                    row.q_skew = false;
                    for point_idx in 0..ADDITIONS_PER_ROW {
                        let add_state = &mut row.add_state[point_idx];
                        add_state.add = false;
                        add_state.slice = T::AcvmType::default();
                        add_state.point = T::AcvmPoint::default();
                        add_state.collision_inverse = T::AcvmType::default();

                        p1_trace[trace_index] = accumulator;
                        p2_trace[trace_index] = accumulator;
                        accumulator = driver.add_points(accumulator, accumulator);
                        p3_trace[trace_index] = accumulator;
                        operation_trace[trace_index] = true;
                        trace_index += 1;
                    }
                    accumulator_trace[msm_row_index] = accumulator;
                    msm_row_index += 1;
                } else {
                    for row_idx in 0..num_rows_per_digit {
                        let row = &mut msm_rows[msm_row_index];

                        let num_points_in_row = if (row_idx + 1) * ADDITIONS_PER_ROW > msm_size {
                            msm_size % ADDITIONS_PER_ROW
                        } else {
                            ADDITIONS_PER_ROW
                        };
                        let offset = row_idx * ADDITIONS_PER_ROW;
                        row.msm_transition = false;
                        for point_idx in 0..ADDITIONS_PER_ROW {
                            let add_state = &mut row.add_state[point_idx];
                            add_state.add = num_points_in_row > point_idx;
                            add_state.slice = if add_state.add {
                                if msm[offset + point_idx].wnaf_skew {
                                    T::AcvmType::from(C::BaseField::from(7))
                                } else {
                                    T::AcvmType::default()
                                }
                            } else {
                                T::AcvmType::default()
                            };

                            add_state.point = if add_state.add {
                                // msm[offset + point_idx].precomputed_table[add_state.slice as usize]
                                let lut = driver.init_lut_by_acvm_point(
                                    msm[offset + point_idx].precomputed_table.to_vec(),
                                );
                                driver.read_lut_by_acvm_point(add_state.slice, &lut)?
                            } else {
                                T::AcvmPoint::<C>::default()
                            };
                            let add_predicate = if add_state.add {
                                msm[offset + point_idx].wnaf_skew
                            } else {
                                false
                            };
                            let p1 = accumulator;
                            accumulator = if add_predicate {
                                driver.add_points(accumulator, add_state.point)
                            } else {
                                accumulator
                            };
                            p1_trace[trace_index] = p1;
                            p2_trace[trace_index] = add_state.point;
                            p3_trace[trace_index] = accumulator;
                            operation_trace[trace_index] = false;
                            trace_index += 1;
                        }
                        row.q_add = false;
                        row.q_double = false;
                        row.q_skew = true;
                        row.msm_round = (digit_idx + 1) as u32;
                        row.msm_size = msm_size as u32;
                        row.msm_count = offset as u32;
                        row.pc = pc;
                        accumulator_trace[msm_row_index] = accumulator;
                        msm_row_index += 1;
                    }
                }
            }
        }

        // Normalize the points in the point trace
        let mut points_to_normalize = Vec::with_capacity(num_points_to_normalize);
        points_to_normalize.extend_from_slice(&p1_trace);
        points_to_normalize.extend_from_slice(&p2_trace);
        points_to_normalize.extend_from_slice(&p3_trace);
        points_to_normalize.extend_from_slice(&accumulator_trace);

        let p1_trace = &points_to_normalize[0..num_point_adds_and_doubles];
        let p2_trace =
            &points_to_normalize[num_point_adds_and_doubles..num_point_adds_and_doubles * 2];
        let accumulator_trace =
            &points_to_normalize[num_point_adds_and_doubles * 3..num_points_to_normalize];

        // inverse_trace is used to compute the value of the `collision_inverse` column in the ECCVM.
        let mut inverse_trace = Vec::with_capacity(num_point_adds_and_doubles);
        for operation_idx in 0..num_point_adds_and_doubles {
            //TODO FLORIN: BATCH THIS
            let (tmp1_x, tmp1_y, _) = driver.pointshare_to_field_shares(p1_trace[operation_idx])?;
            let (tmp2_x, _, _) = driver.pointshare_to_field_shares(p1_trace[operation_idx])?;
            if operation_trace[operation_idx] {
                inverse_trace.push(driver.add(tmp1_x, tmp1_y));
            } else {
                inverse_trace.push(driver.sub(tmp2_x, tmp1_x));
            }
        }

        // TODO FLORIN INVERT THIS
        // ark_ff::batch_inversion(&mut inverse_trace);

        // complete the computation of the ECCVM execution trace, by adding the affine intermediate point data
        // i.e. row.accumulator_x, row.accumulator_y, row.add_state[0...3].collision_inverse,
        // row.add_state[0...3].lambda
        for msm_idx in 0..msms.len() {
            let msm = &msms[msm_idx];
            let mut trace_index = (msm_row_counts[msm_idx] - 1) * ADDITIONS_PER_ROW;
            let mut msm_row_index = msm_row_counts[msm_idx];
            // 1st MSM row will have accumulator equal to the previous MSM output
            // (or point at infinity for 1st MSM)
            let mut accumulator_index = msm_row_counts[msm_idx] - 1;
            let msm_size = msm.len();
            let num_rows_per_digit = (msm_size / ADDITIONS_PER_ROW)
                + (if msm_size % ADDITIONS_PER_ROW != 0 {
                    1
                } else {
                    0
                });

            for digit_idx in 0..NUM_WNAF_DIGITS_PER_SCALAR {
                for _ in 0..num_rows_per_digit {
                    let row = &mut msm_rows[msm_row_index];
                    let normalized_accumulator = &accumulator_trace[accumulator_index];
                    //TODO FLORIN: BATCH THIS
                    let (normalized_accumulator_x, normalized_accumulator_y, _) =
                        driver.pointshare_to_field_shares(*normalized_accumulator)?;
                    row.accumulator_x = normalized_accumulator_x;
                    row.accumulator_y = normalized_accumulator_y;
                    for point_idx in 0..ADDITIONS_PER_ROW {
                        let add_state = &mut row.add_state[point_idx];
                        let inverse = &inverse_trace[trace_index];
                        let p1 = &p1_trace[trace_index];
                        let p2 = &p2_trace[trace_index];
                        add_state.collision_inverse = if add_state.add {
                            *inverse
                        } else {
                            T::AcvmType::default()
                        };
                        add_state.lambda = if add_state.add {
                            //TODO FLORIN: BATCH THIS
                            let (p1_x, p1_y, _) = driver.pointshare_to_field_shares(*p1)?;
                            let (p2_x, p2_y, _) = driver.pointshare_to_field_shares(*p2)?;
                            let sub = driver.sub(p2_y, p1_y);
                            driver.mul(sub, *inverse)?
                        } else {
                            T::AcvmType::default()
                        };
                        trace_index += 1;
                    }
                    accumulator_index += 1;
                    msm_row_index += 1;
                }

                if digit_idx < NUM_WNAF_DIGITS_PER_SCALAR - 1 {
                    let row = &mut msm_rows[msm_row_index];
                    let normalized_accumulator = &accumulator_trace[accumulator_index];
                    let (normalized_accumulator_x, normalized_accumulator_y, _) =
                        driver.pointshare_to_field_shares(*normalized_accumulator)?;
                    let acc_x = normalized_accumulator_x;
                    let acc_y = normalized_accumulator_y;
                    row.accumulator_x = acc_x;
                    row.accumulator_y = acc_y;
                    for point_idx in 0..ADDITIONS_PER_ROW {
                        let add_state = &mut row.add_state[point_idx];
                        add_state.collision_inverse = T::AcvmType::default();
                        // TODO FLORIN: BATCH THIS
                        let (p1_x, _, _) =
                            driver.pointshare_to_field_shares(p1_trace[trace_index])?;
                        let dx = &p1_x;
                        let inverse = &inverse_trace[trace_index];
                        // TODO FLORIN: BATCH THIS
                        let three_dx = driver.mul_with_public(C::BaseField::from(3), *dx);
                        let three_dx_dx = driver.mul(three_dx, *dx)?;
                        add_state.lambda = driver.mul(three_dx_dx, *inverse)?; //((*dx + dx + dx) * dx) * inverse;
                        trace_index += 1;
                    }
                    accumulator_index += 1;
                    msm_row_index += 1;
                } else {
                    for row_idx in 0..num_rows_per_digit {
                        let row = &mut msm_rows[msm_row_index];
                        let normalized_accumulator = &accumulator_trace[accumulator_index];
                        let offset = row_idx * ADDITIONS_PER_ROW;
                        // TODO FLORIN: BATCH THIS
                        let (normalized_accumulator_x, normalized_accumulator_y, _) =
                            driver.pointshare_to_field_shares(*normalized_accumulator)?;
                        row.accumulator_x = normalized_accumulator_x;
                        row.accumulator_y = normalized_accumulator_y;
                        for point_idx in 0..ADDITIONS_PER_ROW {
                            let add_state = &mut row.add_state[point_idx];
                            let add_predicate = if add_state.add {
                                msm[offset + point_idx].wnaf_skew
                            } else {
                                false
                            };

                            let inverse = &inverse_trace[trace_index];
                            let p1 = &p1_trace[trace_index];
                            let p2 = &p2_trace[trace_index];
                            add_state.collision_inverse = if add_predicate {
                                *inverse
                            } else {
                                T::AcvmType::default()
                            };
                            add_state.lambda = if add_predicate {
                                //TODO FLORIN: BATCH THIS
                                let (_, p1_y, _) = driver.pointshare_to_field_shares(*p1)?;
                                let (_, p2_y, _) = driver.pointshare_to_field_shares(*p2)?;
                                let sub = driver.sub(p2_y, p1_y);
                                driver.mul(sub, *inverse)?
                            } else {
                                T::AcvmType::default()
                            };
                            trace_index += 1;
                        }
                        accumulator_index += 1;
                        msm_row_index += 1;
                    }
                }
            }
        }

        // populate the final row in the MSM execution trace.
        // we always require 1 extra row at the end of the trace, because the accumulator x/y coordinates for row `i`
        // are present at row `i+1`
        let final_accumulator = accumulator_trace
            .last()
            .expect("Should have at least one accumulator");
        let final_row = &mut msm_rows.last_mut().expect("Should have at least one row");
        final_row.pc = *pc_values.last().expect("Should have at least one pc value");
        final_row.msm_transition = true;
        let (final_x, final_y, _) = driver.pointshare_to_field_shares(*final_accumulator)?;
        final_row.accumulator_x = final_x;
        final_row.accumulator_y = final_y;
        final_row.msm_size = 0;
        final_row.msm_count = 0;
        final_row.q_add = false;
        final_row.q_double = false;
        final_row.q_skew = false;
        final_row.add_state = [
            AddState {
                add: false,
                slice: T::AcvmType::default(),
                point: T::AcvmPoint::<C>::default(),
                lambda: T::AcvmType::default(),
                collision_inverse: T::AcvmType::default(),
            },
            AddState {
                add: false,
                slice: T::AcvmType::default(),
                point: T::AcvmPoint::<C>::default(),
                lambda: T::AcvmType::default(),
                collision_inverse: T::AcvmType::default(),
            },
            AddState {
                add: false,
                slice: T::AcvmType::default(),
                point: T::AcvmPoint::<C>::default(),
                lambda: T::AcvmType::default(),
                collision_inverse: T::AcvmType::default(),
            },
            AddState {
                add: false,
                slice: T::AcvmType::default(),
                point: T::AcvmPoint::<C>::default(),
                lambda: T::AcvmType::default(),
                collision_inverse: T::AcvmType::default(),
            },
        ];

        Ok((msm_rows, point_table_read_counts))
    }
}
