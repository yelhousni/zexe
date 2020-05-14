#[cfg(feature = "bw6_bis")]
mod curves;
mod fields;

#[cfg(feature = "bw6_bis")]
pub use curves::*;
pub use fields::*;
