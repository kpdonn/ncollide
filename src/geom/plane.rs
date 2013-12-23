//!
//! Support mapping based Plane geometry.
//!

use nalgebra::na::AlgebraicVec;
use nalgebra::na;

/**
 * Implicit description of a plane.
 *
 *   - `V`: type of the plane normal.
 */
#[deriving(Eq, ToStr, Clone, Encodable, Decodable)]
pub struct Plane<N, V> {
    /// The plane normal.
    normal: V
}

impl<N: Algebraic, V: AlgebraicVec<N>> Plane<N, V> {
    /// Builds a new plane from its center and its normal.
    #[inline]
    pub fn new(normal: V) -> Plane<N, V> {
        unsafe { Plane::new_normalized(na::normalize(&normal)) }
    }

    /// Builds a new plane from its center and its normal.
    #[inline]
    pub unsafe fn new_normalized(normal: V) -> Plane<N, V> {
        Plane {
            normal: normal
        }
    }
}

impl<N, V: Clone> Plane<N, V> {
    /// The plane normal.
    #[inline]
    pub fn normal(&self) -> V {
        self.normal.clone()
    }
}
