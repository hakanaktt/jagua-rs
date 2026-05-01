use std::hash::{Hash, Hasher};

use crate::geometry::Transformation;
use crate::geometry::geo_traits::{CollidesWith, DistanceTo, Transformable, TransformableFrom};

/// A Point in 2D space with x and y coordinates
#[derive(Debug, Clone, PartialEq, Copy)]
pub struct Point(pub f64, pub f64);

impl Transformable for Point {
    fn transform(&mut self, t: &Transformation) -> &mut Self {
        let Point(x, y) = self;
        (*x, *y) = TRANSFORM_FORMULA(*x, *y, t);
        self
    }
}

impl TransformableFrom for Point {
    fn transform_from(&mut self, reference: &Self, t: &Transformation) -> &mut Self {
        let Point(x, y) = self;
        (*x, *y) = TRANSFORM_FORMULA(reference.0, reference.1, t);
        self
    }
}

const TRANSFORM_FORMULA: fn(f64, f64, &Transformation) -> (f64, f64) = |x, y, t| -> (f64, f64) {
    let m = t.matrix();
    let t_x = m[0][0].into_inner() * x + m[0][1].into_inner() * y + m[0][2].into_inner() * 1.0;
    let t_y = m[1][0].into_inner() * x + m[1][1].into_inner() * y + m[1][2].into_inner() * 1.0;
    (t_x, t_y)
};

impl Point {
    pub fn x(&self) -> f64 {
        self.0
    }

    pub fn y(&self) -> f64 {
        self.1
    }
}

impl DistanceTo<Point> for Point {
    #[inline(always)]
    fn distance_to(&self, other: &Point) -> f64 {
        ((self.0 - other.0).powi(2) + (self.1 - other.1).powi(2)).sqrt()
    }

    #[inline(always)]
    fn sq_distance_to(&self, other: &Point) -> f64 {
        (self.0 - other.0).powi(2) + (self.1 - other.1).powi(2)
    }
}

impl Eq for Point {}

impl Hash for Point {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let x = self.0.to_bits();
        let y = self.1.to_bits();
        x.hash(state);
        y.hash(state);
    }
}

impl From<Point> for (f64, f64) {
    fn from(p: Point) -> Self {
        (p.0, p.1)
    }
}

impl From<(f64, f64)> for Point {
    fn from((x, y): (f64, f64)) -> Self {
        Point(x, y)
    }
}

impl<T> CollidesWith<T> for Point
where
    T: CollidesWith<Point>,
{
    fn collides_with(&self, other: &T) -> bool {
        other.collides_with(self)
    }
}
