//!
//! Support mapping based Cylinder geometry.
//!

use nalgebra::na::Cast;

/// Implicit description of a cylinder geometry with its principal axis aligned with the `x` axis.
#[deriving(Eq, ToStr, Clone, Encodable, Decodable)]
pub struct Cylinder<N> {
    priv half_height: N,
    priv radius:      N,
    priv margin:      N
}

impl<N: Signed + Cast<f32>> Cylinder<N> {
    /// Creates a new cylinder.
    ///
    /// # Arguments:
    ///     * `half_height` - the half length of the cylinder along the `x` axis.
    ///     * `radius` - the length of the cylinder along all other axis.
    pub fn new(half_height: N, radius: N) -> Cylinder<N> {
        Cylinder::new_with_margin(half_height, radius, Cast::from(0.04))
    }

    /// Creates a new cylinder.
    ///
    /// # Arguments:
    ///     * `half_height` - the half length of the cylinder along the `x` axis.
    ///     * `radius` - the length of the cylinder along all other axis.
    pub fn new_with_margin(half_height: N, radius: N, margin: N) -> Cylinder<N> {
        assert!(half_height.is_positive() && radius.is_positive());

        Cylinder {
            half_height: half_height,
            radius:      radius,
            margin:      margin
        }
    }
}

impl<N: Clone> Cylinder<N> {
    /// The cylinder half length along the `x` axis.
    pub fn half_height(&self) -> N {
        self.half_height.clone()
    }

    /// The radius of the cylinder along all but the `x` axis.
    pub fn radius(&self) -> N {
        self.radius.clone()
    }

    /// Size of the margin around the cylinder.
    pub fn margin(&self) -> N {
        self.margin.clone()
    }
}
