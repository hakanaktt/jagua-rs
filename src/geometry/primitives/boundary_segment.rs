use crate::geometry::Transformation;
use crate::geometry::geo_enums::GeoPosition;
use crate::geometry::geo_traits::{
    CollidesWith, DistanceTo, SeparationDistance, Transformable, TransformableFrom,
};
use crate::geometry::primitives::{Arc, Circle, Edge, Point, Rect};

/// Boundary segment of a closed shape: either a straight line or a circular arc.
#[derive(Clone, Debug, PartialEq, Copy)]
pub enum BoundarySegment {
    Line(Edge),
    Arc(Arc),
}

impl BoundarySegment {
    pub fn start(&self) -> Point {
        match self {
            BoundarySegment::Line(edge) => edge.start,
            BoundarySegment::Arc(arc) => arc.start,
        }
    }

    pub fn end(&self) -> Point {
        match self {
            BoundarySegment::Line(edge) => edge.end,
            BoundarySegment::Arc(arc) => arc.end,
        }
    }

    pub fn bbox(&self) -> Rect {
        match self {
            BoundarySegment::Line(edge) => edge.bbox(),
            BoundarySegment::Arc(arc) => arc.bbox,
        }
    }

    pub fn collides_at(&self, edge: &Edge) -> Vec<Point> {
        match self {
            BoundarySegment::Line(line) => line.collides_at(edge).into_iter().collect(),
            BoundarySegment::Arc(arc) => arc.collides_at(edge),
        }
    }
}

impl From<Edge> for BoundarySegment {
    fn from(edge: Edge) -> Self {
        BoundarySegment::Line(edge)
    }
}

impl From<Arc> for BoundarySegment {
    fn from(arc: Arc) -> Self {
        BoundarySegment::Arc(arc)
    }
}

impl Transformable for BoundarySegment {
    fn transform(&mut self, t: &Transformation) -> &mut Self {
        match self {
            BoundarySegment::Line(edge) => {
                edge.transform(t);
            }
            BoundarySegment::Arc(arc) => {
                arc.transform(t);
            }
        }
        self
    }
}

impl TransformableFrom for BoundarySegment {
    fn transform_from(&mut self, reference: &Self, t: &Transformation) -> &mut Self {
        match reference {
            BoundarySegment::Line(reference_edge) => match self {
                BoundarySegment::Line(edge) => {
                    edge.transform_from(reference_edge, t);
                }
                _ => {
                    *self = BoundarySegment::Line(*reference_edge);
                    self.transform(t);
                }
            },
            BoundarySegment::Arc(reference_arc) => match self {
                BoundarySegment::Arc(arc) => {
                    arc.transform_from(reference_arc, t);
                }
                _ => {
                    *self = BoundarySegment::Arc(*reference_arc);
                    self.transform(t);
                }
            },
        }
        self
    }
}

impl CollidesWith<Point> for BoundarySegment {
    fn collides_with(&self, point: &Point) -> bool {
        match self {
            BoundarySegment::Line(edge) => edge.sq_distance_to(point) == 0.0,
            BoundarySegment::Arc(arc) => arc.collides_with(point),
        }
    }
}

impl CollidesWith<Edge> for BoundarySegment {
    fn collides_with(&self, edge: &Edge) -> bool {
        match self {
            BoundarySegment::Line(line) => line.collides_with(edge),
            BoundarySegment::Arc(arc) => arc.collides_with(edge),
        }
    }
}

impl CollidesWith<Arc> for BoundarySegment {
    fn collides_with(&self, arc: &Arc) -> bool {
        match self {
            BoundarySegment::Line(edge) => edge.collides_with(arc),
            BoundarySegment::Arc(self_arc) => self_arc.collides_with(arc),
        }
    }
}

impl CollidesWith<Rect> for BoundarySegment {
    fn collides_with(&self, rect: &Rect) -> bool {
        match self {
            BoundarySegment::Line(edge) => edge.collides_with(rect),
            BoundarySegment::Arc(arc) => arc.collides_with(rect),
        }
    }
}

impl CollidesWith<Circle> for BoundarySegment {
    fn collides_with(&self, circle: &Circle) -> bool {
        match self {
            BoundarySegment::Line(edge) => circle.collides_with(edge),
            BoundarySegment::Arc(arc) => arc.collides_with(circle),
        }
    }
}

impl CollidesWith<BoundarySegment> for BoundarySegment {
    fn collides_with(&self, other: &BoundarySegment) -> bool {
        match other {
            BoundarySegment::Line(edge) => self.collides_with(edge),
            BoundarySegment::Arc(arc) => self.collides_with(arc),
        }
    }
}

impl CollidesWith<BoundarySegment> for Edge {
    fn collides_with(&self, segment: &BoundarySegment) -> bool {
        segment.collides_with(self)
    }
}

impl CollidesWith<BoundarySegment> for Arc {
    fn collides_with(&self, segment: &BoundarySegment) -> bool {
        match segment {
            BoundarySegment::Line(edge) => self.collides_with(edge),
            BoundarySegment::Arc(arc) => self.collides_with(arc),
        }
    }
}

impl CollidesWith<BoundarySegment> for Rect {
    fn collides_with(&self, segment: &BoundarySegment) -> bool {
        segment.collides_with(self)
    }
}

impl CollidesWith<BoundarySegment> for Circle {
    fn collides_with(&self, segment: &BoundarySegment) -> bool {
        segment.collides_with(self)
    }
}

impl DistanceTo<Point> for BoundarySegment {
    fn distance_to(&self, point: &Point) -> f64 {
        self.sq_distance_to(point).sqrt()
    }

    fn sq_distance_to(&self, point: &Point) -> f64 {
        match self {
            BoundarySegment::Line(edge) => edge.sq_distance_to(point),
            BoundarySegment::Arc(arc) => arc.sq_distance_to(point),
        }
    }
}

impl SeparationDistance<Point> for BoundarySegment {
    fn separation_distance(&self, point: &Point) -> (GeoPosition, f64) {
        let (position, sq_distance) = self.sq_separation_distance(point);
        (position, sq_distance.sqrt())
    }

    fn sq_separation_distance(&self, point: &Point) -> (GeoPosition, f64) {
        match self {
            BoundarySegment::Line(edge) => {
                let sq_distance = edge.sq_distance_to(point);
                if sq_distance == 0.0 {
                    (GeoPosition::Interior, sq_distance)
                } else {
                    (GeoPosition::Exterior, sq_distance)
                }
            }
            BoundarySegment::Arc(arc) => arc.sq_separation_distance(point),
        }
    }
}
