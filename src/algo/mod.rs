//! Internal algorithm modules for Bessel function computation.
//!
//! These modules implement the core numerical routines from Amos Algorithm 644.
//! All functions are `pub(crate)` — they are not part of the public API.

// Phase 2: K function core path
pub(crate) mod bknu;
pub(crate) mod gamln;
pub(crate) mod kscl;
pub(crate) mod s1s2;
pub(crate) mod shch;
pub(crate) mod uchk;

// Phase 3: H function + Wronskian
pub(crate) mod rati;
pub(crate) mod wrsk;

// Phase 4: I, J, Y + analytic continuation
pub(crate) mod acon;
pub(crate) mod asyi;
pub(crate) mod binu;
pub(crate) mod mlri;
pub(crate) mod seri;
pub(crate) mod uoik;

// Phase 5a: Region 1 uniform asymptotic
pub(crate) mod uni1;
pub(crate) mod unik;
pub(crate) mod unk1;
