//! Internal algorithm modules for Bessel function computation.
//!
//! These modules implement the core numerical routines from Amos Algorithm 644.
//! All functions are `pub(crate)` — they are not part of the public API.
//!
//! # Return value convention
//!
//! Most internal routines return `Result<usize, Error>` where the `usize`
//! value (`nz`) indicates the number of trailing output components set to zero
//! due to underflow. The upper-interface functions (`zbesj`, `zbesy`, etc.)
//! translate these into the public `Result<(usize, Accuracy), Error>`.
//!
//! Some routines (notably `zbknu`) return `nz` as a raw `i32` internally:
//! - `nz >= 0`: number of underflowed trailing components
//! - `nz == -1`: overflow detected
//! - `nz == -2`: convergence failure
//!
//! These sentinel values are converted to the appropriate `Error` variant
//! before reaching the public API.
//!
//! # Clippy suppressions
//!
//! Most algorithm modules carry module-level `#![allow(...)]` attributes:
//!
//! - `clippy::excessive_precision` / `clippy::approx_constant` — Fortran DATA
//!   constants are transcribed at full f64 precision for 1:1 traceability against
//!   zbsubs.f. Letting clippy round them would silently break verification.
//! - `clippy::too_many_arguments` — internal functions mirror Fortran subroutine
//!   signatures (often 8–12 parameters). Restructuring would obscure the
//!   correspondence with the reference implementation.

pub(crate) mod constants;

// K function core path
pub(crate) mod bknu;
pub(crate) mod gamln;
pub(crate) mod kscl;
pub(crate) mod s1s2;
pub(crate) mod shch;
pub(crate) mod uchk;

// H function + Wronskian
pub(crate) mod rati;
pub(crate) mod wrsk;

// I, J, Y + analytic continuation
pub(crate) mod acon;
pub(crate) mod asyi;
pub(crate) mod binu;
pub(crate) mod mlri;
pub(crate) mod seri;
pub(crate) mod uoik;

// Region 1 uniform asymptotic
pub(crate) mod uni1;
pub(crate) mod unik;
pub(crate) mod unk1;

// Airy analytic continuation helper
pub(crate) mod acai;

// Region 2 uniform asymptotic + dispatchers
pub(crate) mod buni;
pub(crate) mod bunk;
pub(crate) mod unhj;
pub(crate) mod uni2;
pub(crate) mod unk2;
