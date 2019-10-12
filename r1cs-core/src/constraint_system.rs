use std::marker::PhantomData;
use algebra::Field;

use crate::{Index, Variable, LinearCombination, SynthesisError};

/// Represents a constraint system that can allocate new variables and enforce new 
/// constraints between these variables.
pub trait ConstraintSystem<F: Field>: Sized {
    /// Represents the type of the "root" of this constraint system
    /// so that nested namespaces can minimize indirection.
    type Root: ConstraintSystem<F>;

    /// Return the "one" input variable
    fn one() -> Variable {
        Variable::new_unchecked(Index::Input(0))
    }

    /// Allocate a private variable in the constraint system. The provided
    /// function is used to determine the assignment of the variable. The
    /// given `annotation` function is invoked in testing contexts in order
    /// to derive a unique name for this variable in the current namespace.
    fn alloc<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Allocate a public variable in the constraint system. The provided
    /// function is used to determine the assignment of the variable.
    fn alloc_input<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Enforce that `A` * `B` = `C`. The `annotation` function is invoked in
    /// testing contexts in order to derive a unique name for the constraint
    /// in the current namespace.
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LB: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LC: FnOnce(LinearCombination<F>) -> LinearCombination<F>;

    /// Create a new (sub)namespace and enter into it. Not intended
    /// for downstream use; use `namespace` instead.
    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR;

    /// Exit out of the existing namespace. Not intended for
    /// downstream use; use `namespace` instead.
    fn pop_namespace(&mut self);

    /// Gets the "root" constraint system, bypassing the namespacing.
    /// Not intended for downstream use; use `namespace` instead.
    fn get_root(&mut self) -> &mut Self::Root;

    /// Begin a namespace for this constraint system.
    fn ns<'a, NR, N>(&'a mut self, name_fn: N) -> Namespace<'a, F, Self::Root>
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        self.get_root().push_namespace(name_fn);

        Namespace(self.get_root(), PhantomData)
    }

    /// Output the number of constraints in the system.
    fn num_constraints(&self) -> usize;

    /// Output the number of variables in the system.
    fn num_variables(&self) -> usize;

    /// Output the maximum number of non-zero entries across the A, B, and C
    /// matrices in the system.
    fn num_non_zero(&self) -> usize;
}

/// This is a "namespaced" constraint system which borrows a constraint system
/// (pushing a namespace context) and, when dropped, pops out of the namespace
/// context.
pub struct Namespace<'a, F: Field, CS: ConstraintSystem<F>>(&'a mut CS, PhantomData<F>);

/// Computations are expressed in terms of rank-1 constraint systems (R1CS).
/// The `generate_constraints` method is called to generate constraints for
/// both CRS generation and for proving.
pub trait ConstraintSynthesizer<F: Field> {
    /// Drives generation of new constraints inside `CS`.
    fn generate_constraints<CS: ConstraintSystem<F>>(self, cs: &mut CS) -> Result<(), SynthesisError>;
}


impl<F: Field, CS: ConstraintSystem<F>> ConstraintSystem<F> for Namespace<'_, F, CS> {
    type Root = CS::Root;

    #[inline]
    fn one() -> Variable {
        CS::one()
    }

    #[inline]
    fn alloc<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.0.alloc(annotation, f)
    }

    #[inline]
    fn alloc_input<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.0.alloc_input(annotation, f)
    }

    #[inline]
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LB: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LC: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
    {
        self.0.enforce(annotation, a, b, c)
    }

    // Downstream users who use `namespace` will never interact with these
    // functions and they will never be invoked because the namespace is
    // never a root constraint system.

    #[inline]
    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        panic!("only the root's push_namespace should be called");
    }

    #[inline]
    fn pop_namespace(&mut self) {
        panic!("only the root's pop_namespace should be called");
    }

    #[inline]
    fn get_root(&mut self) -> &mut Self::Root {
        self.0.get_root()
    }

    #[inline]
    fn num_constraints(&self) -> usize {
        self.0.num_constraints()
    }

    #[inline]
    fn num_variables(&self) -> usize {
        self.0.num_variables()
    }

    #[inline]
    fn num_non_zero(&self) -> usize {
        self.0.num_non_zero()
    }
}

impl<F: Field, CS: ConstraintSystem<F>> Drop for Namespace<'_, F, CS> {
    #[inline]
    fn drop(&mut self) {
        self.get_root().pop_namespace()
    }
}

/// Convenience implementation of ConstraintSystem<F> for mutable references to
/// constraint systems.
impl<F: Field, CS: ConstraintSystem<F>> ConstraintSystem<F> for &mut CS {
    type Root = CS::Root;

    #[inline]
    fn one() -> Variable {
        CS::one()
    }

    #[inline]
    fn alloc<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc(annotation, f)
    }

    #[inline]
    fn alloc_input<FN, A, AR>(&mut self, annotation: A, f: FN) -> Result<Variable, SynthesisError>
    where
        FN: FnOnce() -> Result<F, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc_input(annotation, f)
    }

    #[inline]
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LB: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
        LC: FnOnce(LinearCombination<F>) -> LinearCombination<F>,
    {
        (**self).enforce(annotation, a, b, c)
    }

    #[inline]
    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        (**self).push_namespace(name_fn)
    }

    #[inline]
    fn pop_namespace(&mut self) {
        (**self).pop_namespace()
    }

    #[inline]
    fn get_root(&mut self) -> &mut Self::Root {
        (**self).get_root()
    }

    #[inline]
    fn num_constraints(&self) -> usize {
        (**self).num_constraints()
    }

    #[inline]
    fn num_variables(&self) -> usize {
        (**self).num_variables()
    }

    #[inline]
    fn num_non_zero(&self) -> usize {
        (**self).num_non_zero()
    }
}

/// A "standard" implementation of the `ConstraintSystem` trait that should
/// suffice for most use cases.
pub struct StandardConstraintSystem<F: Field> {
    /// The mode in which the constraint system is operating. `self` can either 
    /// be in setup mode (i.e., `self.mode == Mode::Setup`) or in proving mode 
    /// (i.e., `self.mode == Mode::Prove`. If we are in proving mode, then we 
    /// have the additional option of whether or not to construct the A, B, and
    /// C matrices (see below).
    pub mode: Mode,
    /// The number of variables that are "public inputs" to the constraint system.
    pub num_input_variables: usize,
    /// The number of variables that are "private inputs" to the constraint system.
    pub num_witness_variables: usize,
    /// The number of constraints in the constraint system.
    pub num_constraints: usize,
    /// The number of non_zero entries in the A matrix.
    pub a_num_non_zero: usize,
    /// The number of non_zero entries in the B matrix.
    pub b_num_non_zero: usize,
    /// The number of non_zero entries in the C matrix.
    pub c_num_non_zero: usize,
    /// The A matrix. This is empty when `self.mode == Mode::Prove { construct_matrices = false }`.
    pub a: Vec<Vec<(F, Index)>>,
    /// The B matrix. This is empty when `self.mode == Mode::Prove { construct_matrices = false }`.
    pub b: Vec<Vec<(F, Index)>>,
    /// The C matrix. This is empty when `self.mode == Mode::Prove { construct_matrices = false }`.
    pub c: Vec<Vec<(F, Index)>>,
    /// Assignments to the public input variables. This is empty if `self.mode == Mode::Setup`.
    pub input_assignment: Vec<F>,
    /// Assignments to the private input variables. This is empty if `self.mode == Mode::Setup`.
    pub witness_assignment: Vec<F>,
}

/// Defines the mode of operation of a `StandardConstraintSystem`.
#[derive(Eq, PartialEq)]
pub enum Mode {
    /// Indicates to the `StandardConstraintSystem` that it should only generate 
    /// constraint matrices and not populate the variable assignments.
    Setup,
    /// Indicates to the `StandardConstraintSystem` that it populate the variable 
    /// assignments. If additionally `construct_matrices == true`, then generate
    /// the matrices as in the `Setup` case.
    Prove { 
        /// If `construct_matrices == true`, then generate
        /// the matrices as in the `Setup` case.
        construct_matrices: bool 
    },
}

impl<F: Field> StandardConstraintSystem<F> {
    #[inline]
    fn make_row(l: &LinearCombination<F>) -> Vec<(F, Index)> {
        l.as_ref()
            .iter()
            .map(|(var, coeff)| (*coeff, var.get_unchecked()))
            .collect()
    }

    /// Construct an ampty `StandardConstraintSystem`.
    pub fn new() -> Self {
        Self {
            num_input_variables: 1,
            num_witness_variables: 0,
            num_constraints: 0,
            a_num_non_zero: 0,
            b_num_non_zero: 0,
            c_num_non_zero: 0,
            a: Vec::new(),
            b: Vec::new(),
            c: Vec::new(),
            input_assignment: Vec::new(),
            witness_assignment: Vec::new(),
            mode: Mode::Setup,
        }
    }

    /// Check whether `self.mode == Mode::Setup`.
    pub fn in_setup_mode(&self) -> bool {
        self.mode == Mode::Setup
    }

    /// Check whether or not `self` will construct matrices.
    pub fn should_construct_matrices(&self) -> bool {
        match self.mode {
            Mode::Setup => true,
            Mode::Prove { construct_matrices } => construct_matrices,
        }
    }
}

impl<ConstraintF: Field> ConstraintSystem<ConstraintF> for StandardConstraintSystem<ConstraintF> {
    type Root = Self;

    #[inline]
    fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't invoke the
        // function for obtaining one.

        let index = self.num_witness_variables;
        self.num_witness_variables += 1;

        if !self.in_setup_mode() {
            self.witness_assignment.push(f()?);
        }
        Ok(Variable::new_unchecked(Index::Aux(index)))
    }

    #[inline]
    fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't invoke the
        // function for obtaining one.

        let index = self.num_input_variables;
        self.num_input_variables += 1;

        if !self.in_setup_mode() {
            self.input_assignment.push(f()?);
        }
        Ok(Variable::new_unchecked(Index::Input(index)))
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LB: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LC: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
    {
        if self.should_construct_matrices() {
            let a_row = Self::make_row(&a(LinearCombination::zero()));
            self.a_num_non_zero += a_row.len();
            self.a.push(a_row);

            let b_row = Self::make_row(&b(LinearCombination::zero()));
            self.b_num_non_zero += b_row.len();
            self.b.push(b_row);

            let c_row = Self::make_row(&c(LinearCombination::zero()));
            self.c_num_non_zero += c_row.len();
            self.c.push(c_row);
        }
        self.num_constraints += 1;
    }

    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn pop_namespace(&mut self) {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }

    fn num_constraints(&self) -> usize {
        self.num_constraints
    }

    fn num_variables(&self) -> usize {
        self.num_input_variables + self.num_witness_variables
    }

    fn num_non_zero(&self) -> usize {
        *[self.a_num_non_zero, self.b_num_non_zero, self.c_num_non_zero]
            .iter()
            .max()
            .expect("iterator is not empty")
    }
}
