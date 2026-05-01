use crate::geometry::Transformation;
use crate::geometry::convex_hull;
use crate::geometry::fail_fast::{piers, pole};
use crate::geometry::geo_traits::{Transformable, TransformableFrom};
use crate::geometry::primitives::Circle;
use crate::geometry::primitives::Edge;
use crate::geometry::primitives::Point;
use crate::geometry::primitives::SPolygon;
use serde::{Deserialize, Serialize};

use anyhow::Result;

#[derive(Clone, Debug)]
/// Surrogate representation of a [`SPolygon`] - a 'light-weight' representation that
/// is fully contained in the original [`SPolygon`].
/// Used for *fail-fast* collision detection.
pub struct SPSurrogate {
    /// Set of [poles](pole::generate_surrogate_poles)
    pub poles: Vec<Circle>,
    /// Set of [piers](piers::generate_piers)
    pub piers: Vec<Edge>,
    /// Indices of the vertices in the [`SPolygon`] that form the convex hull
    pub convex_hull_indices: Vec<usize>,
    /// Points that form a hull approximation for the [`SPolygon`], including extrema of arc segments.
    pub convex_hull_points: Vec<Point>,
    /// The area of the convex hull of the [`SPolygon`].
    pub convex_hull_area: f64,
    /// The configuration used to generate the surrogate
    pub config: SPSurrogateConfig,
}

impl SPSurrogate {
    /// Creates a new [`SPSurrogate`] from a [`SPolygon`] and a configuration.
    /// Expensive operations are performed here!
    pub fn new(simple_poly: &SPolygon, config: SPSurrogateConfig) -> Result<Self> {
        let convex_hull_indices = convex_hull::convex_hull_indices(simple_poly);
        let convex_hull_points = convex_hull::convex_hull_points(simple_poly);
        let convex_hull_area = SPolygon::calculate_area(&convex_hull_points).max(simple_poly.area);
        let poles = pole::generate_surrogate_poles(simple_poly, &config.n_pole_limits)?;
        let n_ff_poles = usize::min(config.n_ff_poles, poles.len());
        let relevant_poles_for_piers = &poles[0..n_ff_poles]; //poi + all poles that will be checked during fail fast are relevant for piers
        let piers =
            piers::generate_piers(simple_poly, config.n_ff_piers, relevant_poles_for_piers)?;

        Ok(Self {
            convex_hull_indices,
            convex_hull_points,
            poles,
            piers,
            convex_hull_area,
            config,
        })
    }

    pub fn ff_poles(&self) -> &[Circle] {
        &self.poles[0..self.config.n_ff_poles]
    }

    pub fn ff_piers(&self) -> &[Edge] {
        &self.piers
    }
}

impl Transformable for SPSurrogate {
    fn transform(&mut self, t: &Transformation) -> &mut Self {
        //destructuring pattern used to ensure that the code is updated accordingly when the struct changes
        let Self {
            convex_hull_indices: _,
            convex_hull_points,
            poles,
            piers,
            convex_hull_area: _,
            config: _,
        } = self;

        //transform poles
        convex_hull_points.iter_mut().for_each(|p| {
            p.transform(t);
        });

        poles.iter_mut().for_each(|c| {
            c.transform(t);
        });

        //transform piers
        piers.iter_mut().for_each(|p| {
            p.transform(t);
        });

        self
    }
}

impl TransformableFrom for SPSurrogate {
    fn transform_from(&mut self, reference: &Self, t: &Transformation) -> &mut Self {
        debug_assert!(self.poles.len() == reference.poles.len());
        debug_assert!(self.piers.len() == reference.piers.len());

        //destructuring pattern used to ensure that the code is updated accordingly when the struct changes
        let Self {
            convex_hull_indices: _,
            convex_hull_points,
            poles,
            piers,
            convex_hull_area: _,
            config: _,
        } = self;

        for (p, ref_p) in convex_hull_points
            .iter_mut()
            .zip(reference.convex_hull_points.iter())
        {
            p.transform_from(ref_p, t);
        }

        for (pole, ref_pole) in poles.iter_mut().zip(reference.poles.iter()) {
            pole.transform_from(ref_pole, t);
        }

        for (pier, ref_pier) in piers.iter_mut().zip(reference.piers.iter()) {
            pier.transform_from(ref_pier, t);
        }

        self
    }
}

/// maximum number of definable pole limits, increase if needed
const N_POLE_LIMITS: usize = 3;

/// Configuration of the [`SPSurrogate`](crate::geometry::fail_fast::SPSurrogate) generation
#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq)]
pub struct SPSurrogateConfig {
    ///Limits on the number of poles to be generated at different coverage levels.
    ///For example: [(100, 0.0), (20, 0.75), (10, 0.90)]:
    ///While the coverage is below 75% the generation will stop at 100 poles.
    ///If 75% coverage with 20 or more poles the generation will stop.
    ///If 90% coverage with 10 or more poles the generation will stop.
    pub n_pole_limits: [(usize, f64); N_POLE_LIMITS],
    ///Number of poles to test during fail-fast
    pub n_ff_poles: usize,
    ///number of piers to test during fail-fast
    pub n_ff_piers: usize,
}

impl SPSurrogateConfig {
    pub fn none() -> Self {
        Self {
            n_pole_limits: [(0, 0.0); N_POLE_LIMITS],
            n_ff_poles: 0,
            n_ff_piers: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn rounded_top_shape() -> SPolygon {
        SPolygon::new_with_bulges(
            vec![
                Point(-1.0, 0.0),
                Point(1.0, 0.0),
                Point(1.0, 2.0),
                Point(-1.0, 2.0),
            ],
            vec![0.0, 0.0, 1.0, 0.0],
        )
        .unwrap()
    }

    #[test]
    fn surrogate_hull_area_accounts_for_curved_boundary() {
        let shape = rounded_top_shape();
        let surrogate = SPSurrogate::new(&shape, SPSurrogateConfig::none()).unwrap();

        assert!(
            surrogate
                .convex_hull_points
                .iter()
                .any(|point| point.0.abs() < 1.0e-3 && (point.1 - 3.0).abs() < 1.0e-3)
        );
        assert!(surrogate.convex_hull_area >= shape.area);
        assert!(surrogate.convex_hull_area > 4.0 + PI / 2.0 - 1.0e-3);
    }
}
