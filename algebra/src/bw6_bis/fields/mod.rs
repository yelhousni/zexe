#[cfg(any(feature = "bw6_bis", feature = "edwards_bw6_bis"))]
pub mod fr;
#[cfg(any(feature = "bw6_bis", feature = "edwards_bw6_bis"))]
pub use self::fr::*;

#[cfg(feature = "bw6_bis")]
pub mod fq;
#[cfg(feature = "bw6_bis")]
pub use self::fq::*;

#[cfg(feature = "bw6_bis")]
pub mod fq3;
#[cfg(feature = "bw6_bis")]
pub use self::fq3::*;

#[cfg(feature = "bw6_bis")]
pub mod fq6;
#[cfg(feature = "bw6_bis")]
pub use self::fq6::*;

#[cfg(all(feature = "bw6_bis", test))]
mod tests;
