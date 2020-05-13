#[cfg(any(feature = "sw6_bis", feature = "edwards_sw6_bis"))]
pub mod fr;
#[cfg(any(feature = "sw6_bis", feature = "edwards_sw6_bis"))]
pub use self::fr::*;

#[cfg(feature = "sw6_bis")]
pub mod fq;
#[cfg(feature = "sw6_bis")]
pub use self::fq::*;

#[cfg(feature = "sw6_bis")]
pub mod fq3;
#[cfg(feature = "sw6_bis")]
pub use self::fq3::*;

#[cfg(feature = "sw6_bis")]
pub mod fq6;
#[cfg(feature = "sw6_bis")]
pub use self::fq6::*;

#[cfg(all(feature = "sw6_bis", test))]
mod tests;
