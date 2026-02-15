//! Internal algorithm modules for Bessel function computation.
//!
//! These modules implement the core numerical routines from Amos Algorithm 644.
//! All functions are `pub(crate)` — they are not part of the public API.

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
