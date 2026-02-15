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
