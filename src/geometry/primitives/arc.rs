use crate::geometry::Transformation;
use crate::geometry::geo_enums::GeoPosition;
use crate::geometry::geo_traits::{
    CollidesWith, DistanceTo, SeparationDistance, Transformable, TransformableFrom,
};
use crate::geometry::primitives::{Circle, Edge, Point, Rect};
use anyhow::{Result, ensure};
use std::f32::consts::{FRAC_PI_2, PI, TAU};

const ANGLE_EPSILON: f32 = 1.0e-5;
const DISTANCE_EPSILON: f32 = 1.0e-4;

/// Circular arc segment between two [`Point`]s.
#[derive(Clone, Debug, PartialEq, Copy)]
pub struct Arc {
    pub center: Point,
    pub radius: f32,
    pub start: Point,
    pub end: Point,
    /// Signed sweep angle in radians. Positive values follow the counter-clockwise direction.
    pub sweep: f32,
    pub bbox: Rect,
}

impl Arc {
    pub fn try_new(
        center: Point,
        radius: f32,
        start: Point,
        end: Point,
        sweep: f32,
    ) -> Result<Self> {
        ensure!(
            radius.is_finite() && radius > 0.0,
            "invalid arc radius: {radius}",
        );
        ensure!(
            center.0.is_finite()
                && center.1.is_finite()
                && start.0.is_finite()
                && start.1.is_finite()
                && end.0.is_finite()
                && end.1.is_finite(),
            "invalid arc point: center={center:?}, start={start:?}, end={end:?}",
        );
        ensure!(
            sweep.is_finite() && sweep.abs() > ANGLE_EPSILON && sweep.abs() <= TAU + ANGLE_EPSILON,
            "invalid arc sweep: {sweep}",
        );
        ensure!(
            (start.distance_to(&center) - radius).abs()
                <= DISTANCE_EPSILON.max(radius * DISTANCE_EPSILON),
            "arc start point is not on the circle: {start:?}",
        );
        ensure!(
            (end.distance_to(&center) - radius).abs()
                <= DISTANCE_EPSILON.max(radius * DISTANCE_EPSILON),
            "arc end point is not on the circle: {end:?}",
        );

        let bbox = calculate_bounding_box(center, radius, start, end, sweep);
        Ok(Self {
            center,
            radius,
            start,
            end,
            sweep,
            bbox,
        })
    }

    pub fn try_from_bulge(start: Point, end: Point, bulge: f32) -> Result<Self> {
        ensure!(
            bulge.is_finite() && bulge != 0.0,
            "invalid arc bulge: {bulge}"
        );
        ensure!(start != end, "arc endpoints must be distinct: {start:?}");

        let chord = Edge::try_new(start, end)?;
        let chord_len = chord.length();
        let abs_bulge = bulge.abs();
        let radius = chord_len * (abs_bulge.powi(2) + 1.0) / (4.0 * abs_bulge);
        let sagitta = abs_bulge * chord_len / 2.0;
        let center_offset = radius - sagitta;
        let chord_dx = end.0 - start.0;
        let chord_dy = end.1 - start.1;

        let mut offset_x = -center_offset * chord_dy / chord_len;
        let mut offset_y = center_offset * chord_dx / chord_len;
        if bulge < 0.0 {
            offset_x = -offset_x;
            offset_y = -offset_y;
        }

        let center = Point(
            start.0 + chord_dx / 2.0 + offset_x,
            start.1 + chord_dy / 2.0 + offset_y,
        );
        let sweep = 4.0 * bulge.atan();

        Self::try_new(center, radius, start, end, sweep)
    }

    #[cfg(feature = "curves")]
    pub fn try_from_cavalier(
        start: cavalier_contours::polyline::PlineVertex<f32>,
        end: cavalier_contours::polyline::PlineVertex<f32>,
    ) -> Result<Self> {
        Self::try_from_bulge(Point(start.x, start.y), Point(end.x, end.y), start.bulge)
    }

    #[cfg(feature = "curves")]
    pub fn to_cavalier_pair(
        &self,
    ) -> (
        cavalier_contours::polyline::PlineVertex<f32>,
        cavalier_contours::polyline::PlineVertex<f32>,
    ) {
        (
            cavalier_contours::polyline::PlineVertex::new(self.start.0, self.start.1, self.bulge()),
            cavalier_contours::polyline::PlineVertex::new(self.end.0, self.end.1, 0.0),
        )
    }

    pub fn bulge(&self) -> f32 {
        (self.sweep / 4.0).tan()
    }

    pub fn start_angle(&self) -> f32 {
        angle_of(self.center, self.start)
    }

    pub fn end_angle(&self) -> f32 {
        angle_of(self.center, self.end)
    }

    pub fn point_at_angle(&self, angle: f32) -> Point {
        point_on_circle(self.center, self.radius, angle)
    }

    pub fn cardinal_extrema(&self) -> Vec<Point> {
        [0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2]
            .into_iter()
            .filter(|angle| self.contains_angle(*angle))
            .map(|angle| self.point_at_angle(angle))
            .collect()
    }

    pub fn contains_angle(&self, angle: f32) -> bool {
        angle_in_sweep(angle, self.start_angle(), self.sweep)
    }

    pub fn contains_point_projection(&self, point: &Point) -> bool {
        if point.sq_distance_to(&self.center) <= DISTANCE_EPSILON.powi(2) {
            return true;
        }
        self.contains_angle(angle_of(self.center, *point))
    }

    pub fn closest_point_on_arc(&self, point: &Point) -> Point {
        if point.sq_distance_to(&self.center) <= DISTANCE_EPSILON.powi(2) {
            return self.start;
        }

        let projected_angle = angle_of(self.center, *point);
        if self.contains_angle(projected_angle) {
            self.point_at_angle(projected_angle)
        } else if point.sq_distance_to(&self.start) <= point.sq_distance_to(&self.end) {
            self.start
        } else {
            self.end
        }
    }

    pub fn collides_at(&self, edge: &Edge) -> Vec<Point> {
        let Point(x1, y1) = edge.start;
        let Point(x2, y2) = edge.end;
        let Point(cx, cy) = self.center;
        let dx = x2 - x1;
        let dy = y2 - y1;
        let fx = x1 - cx;
        let fy = y1 - cy;

        let a = dx * dx + dy * dy;
        let b = 2.0 * (fx * dx + fy * dy);
        let c = fx * fx + fy * fy - self.radius.powi(2);
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < -DISTANCE_EPSILON {
            return vec![];
        }

        let mut intersections = Vec::new();
        if discriminant.abs() <= DISTANCE_EPSILON {
            self.push_edge_intersection(edge, -b / (2.0 * a), &mut intersections);
        } else {
            let sqrt_discriminant = discriminant.sqrt();
            self.push_edge_intersection(
                edge,
                (-b - sqrt_discriminant) / (2.0 * a),
                &mut intersections,
            );
            self.push_edge_intersection(
                edge,
                (-b + sqrt_discriminant) / (2.0 * a),
                &mut intersections,
            );
        }
        intersections
    }

    fn push_edge_intersection(&self, edge: &Edge, t: f32, intersections: &mut Vec<Point>) {
        if !(-DISTANCE_EPSILON..=1.0 + DISTANCE_EPSILON).contains(&t) {
            return;
        }
        let clamped_t = t.clamp(0.0, 1.0);
        let point = Point(
            edge.start.0 + (edge.end.0 - edge.start.0) * clamped_t,
            edge.start.1 + (edge.end.1 - edge.start.1) * clamped_t,
        );
        if self.collides_with(&point)
            && intersections
                .iter()
                .all(|other: &Point| other.sq_distance_to(&point) > DISTANCE_EPSILON.powi(2))
        {
            intersections.push(point);
        }
    }
}

impl Transformable for Arc {
    fn transform(&mut self, t: &Transformation) -> &mut Self {
        self.center.transform(t);
        self.start.transform(t);
        self.end.transform(t);
        self.bbox =
            calculate_bounding_box(self.center, self.radius, self.start, self.end, self.sweep);
        self
    }
}

impl TransformableFrom for Arc {
    fn transform_from(&mut self, reference: &Self, t: &Transformation) -> &mut Self {
        self.center.transform_from(&reference.center, t);
        self.start.transform_from(&reference.start, t);
        self.end.transform_from(&reference.end, t);
        self.radius = reference.radius;
        self.sweep = reference.sweep;
        self.bbox =
            calculate_bounding_box(self.center, self.radius, self.start, self.end, self.sweep);
        self
    }
}

impl CollidesWith<Point> for Arc {
    fn collides_with(&self, point: &Point) -> bool {
        (point.distance_to(&self.center) - self.radius).abs()
            <= DISTANCE_EPSILON.max(self.radius * DISTANCE_EPSILON)
            && self.contains_point_projection(point)
    }
}

impl CollidesWith<Edge> for Arc {
    fn collides_with(&self, edge: &Edge) -> bool {
        if !self.bbox.collides_with(&edge.bbox()) {
            return false;
        }
        !self.collides_at(edge).is_empty()
    }
}

impl CollidesWith<Arc> for Arc {
    fn collides_with(&self, other: &Arc) -> bool {
        if !self.bbox.collides_with(&other.bbox) {
            return false;
        }

        let center_distance = self.center.distance_to(&other.center);
        if center_distance <= DISTANCE_EPSILON {
            if (self.radius - other.radius).abs() > DISTANCE_EPSILON {
                return false;
            }
            return self.contains_angle(other.start_angle())
                || self.contains_angle(other.end_angle())
                || other.contains_angle(self.start_angle())
                || other.contains_angle(self.end_angle());
        }

        if center_distance > self.radius + other.radius + DISTANCE_EPSILON
            || center_distance < (self.radius - other.radius).abs() - DISTANCE_EPSILON
        {
            return false;
        }

        let a = (self.radius.powi(2) - other.radius.powi(2) + center_distance.powi(2))
            / (2.0 * center_distance);
        let h_sq = self.radius.powi(2) - a.powi(2);
        if h_sq < -DISTANCE_EPSILON {
            return false;
        }

        let Point(cx1, cy1) = self.center;
        let Point(cx2, cy2) = other.center;
        let ux = (cx2 - cx1) / center_distance;
        let uy = (cy2 - cy1) / center_distance;
        let base = Point(cx1 + a * ux, cy1 + a * uy);

        if h_sq.abs() <= DISTANCE_EPSILON {
            return self.collides_with(&base) && other.collides_with(&base);
        }

        let h = h_sq.sqrt();
        let p1 = Point(base.0 - uy * h, base.1 + ux * h);
        let p2 = Point(base.0 + uy * h, base.1 - ux * h);
        (self.collides_with(&p1) && other.collides_with(&p1))
            || (self.collides_with(&p2) && other.collides_with(&p2))
    }
}

impl CollidesWith<Rect> for Arc {
    fn collides_with(&self, rect: &Rect) -> bool {
        if !self.bbox.collides_with(rect) {
            return false;
        }
        if rect.collides_with(&self.start) || rect.collides_with(&self.end) {
            return true;
        }
        rect.edges().iter().any(|edge| self.collides_with(edge))
    }
}

impl CollidesWith<Circle> for Arc {
    fn collides_with(&self, circle: &Circle) -> bool {
        self.distance_to(&circle.center) <= circle.radius
    }
}

impl CollidesWith<Arc> for Edge {
    fn collides_with(&self, arc: &Arc) -> bool {
        arc.collides_with(self)
    }
}

impl CollidesWith<Arc> for Rect {
    fn collides_with(&self, arc: &Arc) -> bool {
        arc.collides_with(self)
    }
}

impl CollidesWith<Arc> for Circle {
    fn collides_with(&self, arc: &Arc) -> bool {
        arc.collides_with(self)
    }
}

impl DistanceTo<Point> for Arc {
    fn distance_to(&self, point: &Point) -> f32 {
        self.sq_distance_to(point).sqrt()
    }

    fn sq_distance_to(&self, point: &Point) -> f32 {
        let closest = self.closest_point_on_arc(point);
        closest.sq_distance_to(point)
    }
}

impl SeparationDistance<Point> for Arc {
    fn separation_distance(&self, point: &Point) -> (GeoPosition, f32) {
        let (position, sq_distance) = self.sq_separation_distance(point);
        (position, sq_distance.sqrt())
    }

    fn sq_separation_distance(&self, point: &Point) -> (GeoPosition, f32) {
        let sq_distance = self.sq_distance_to(point);
        if self.collides_with(point) {
            (GeoPosition::Interior, sq_distance)
        } else {
            (GeoPosition::Exterior, sq_distance)
        }
    }
}

fn calculate_bounding_box(
    center: Point,
    radius: f32,
    start: Point,
    end: Point,
    sweep: f32,
) -> Rect {
    let mut points = vec![start, end];
    for angle in [0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2] {
        if angle_in_sweep(angle, angle_of(center, start), sweep) {
            points.push(point_on_circle(center, radius, angle));
        }
    }

    let (mut x_min, mut y_min, mut x_max, mut y_max) = (
        f32::INFINITY,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NEG_INFINITY,
    );
    for point in points {
        x_min = x_min.min(point.0);
        y_min = y_min.min(point.1);
        x_max = x_max.max(point.0);
        y_max = y_max.max(point.1);
    }

    Rect {
        x_min,
        y_min,
        x_max,
        y_max,
    }
}

fn point_on_circle(center: Point, radius: f32, angle: f32) -> Point {
    Point(
        center.0 + radius * angle.cos(),
        center.1 + radius * angle.sin(),
    )
}

fn angle_of(center: Point, point: Point) -> f32 {
    normalize_angle((point.1 - center.1).atan2(point.0 - center.0))
}

fn angle_in_sweep(angle: f32, start_angle: f32, sweep: f32) -> bool {
    if sweep.abs() >= TAU - ANGLE_EPSILON {
        return true;
    }

    if sweep >= 0.0 {
        normalize_angle(angle - start_angle) <= sweep + ANGLE_EPSILON
    } else {
        normalize_angle(start_angle - angle) <= -sweep + ANGLE_EPSILON
    }
}

fn normalize_angle(angle: f32) -> f32 {
    let normalized = angle % TAU;
    if normalized < 0.0 {
        normalized + TAU
    } else {
        normalized
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tessellated_edges(arc: &Arc, n_segments: usize) -> Vec<Edge> {
        (0..n_segments)
            .map(|i| {
                let t0 = i as f32 / n_segments as f32;
                let t1 = (i + 1) as f32 / n_segments as f32;
                Edge::try_new(
                    arc.point_at_angle(arc.start_angle() + arc.sweep * t0),
                    arc.point_at_angle(arc.start_angle() + arc.sweep * t1),
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn arc_from_bulge_matches_expected_half_circle() {
        let arc = Arc::try_from_bulge(Point(0.0, 0.0), Point(1.0, 0.0), 1.0).unwrap();
        assert!((arc.center.0 - 0.5).abs() < 1.0e-5);
        assert!(arc.center.1.abs() < 1.0e-5);
        assert!((arc.radius - 0.5).abs() < 1.0e-5);
        assert!((arc.sweep - PI).abs() < 1.0e-5);
        assert!((arc.bulge() - 1.0).abs() < 1.0e-5);
    }

    #[test]
    fn arc_bbox_includes_cardinal_extrema() {
        let arc = Arc::try_from_bulge(Point(0.0, 0.0), Point(1.0, 0.0), 1.0).unwrap();
        assert!((arc.bbox.x_min - 0.0).abs() < 1.0e-5);
        assert!((arc.bbox.x_max - 1.0).abs() < 1.0e-5);
        assert!((arc.bbox.y_min + 0.5).abs() < 1.0e-5);
        assert!(arc.bbox.y_max.abs() < 1.0e-5);
    }

    #[test]
    fn arc_edge_collision_detects_crossing_and_miss() {
        let arc = Arc::try_from_bulge(Point(0.0, 0.0), Point(1.0, 0.0), 1.0).unwrap();
        let crossing = Edge::try_new(Point(0.5, -1.0), Point(0.5, 1.0)).unwrap();
        let miss = Edge::try_new(Point(0.5, 0.1), Point(0.5, 1.0)).unwrap();
        assert!(arc.collides_with(&crossing));
        assert!(!arc.collides_with(&miss));
    }

    #[test]
    fn arc_rect_collision_matches_tessellation_oracle() {
        let arc = Arc::try_from_bulge(Point(0.0, 0.0), Point(2.0, 0.0), 0.75).unwrap();
        let rects = [
            Rect::try_new(0.9, -0.7, 1.1, -0.4).unwrap(),
            Rect::try_new(0.9, 0.2, 1.1, 0.5).unwrap(),
            Rect::try_new(-0.2, -0.1, 0.2, 0.2).unwrap(),
            Rect::try_new(1.8, -0.1, 2.2, 0.2).unwrap(),
        ];
        let oracle_edges = tessellated_edges(&arc, 256);

        for rect in rects {
            let oracle = oracle_edges.iter().any(|edge| rect.collides_with(edge));
            assert_eq!(arc.collides_with(&rect), oracle, "rect={rect:?}");
        }
    }

    #[test]
    fn arc_arc_collision_detects_circle_intersection_on_both_sweeps() {
        let a1 = Arc::try_from_bulge(Point(0.0, 0.0), Point(2.0, 0.0), 1.0).unwrap();
        let a2 = Arc::try_from_bulge(Point(1.0, -1.0), Point(1.0, 1.0), 1.0).unwrap();
        assert!(a1.collides_with(&a2));
    }

    #[test]
    fn arc_transform_preserves_radius_and_sweep() {
        let mut arc = Arc::try_from_bulge(Point(0.0, 0.0), Point(1.0, 0.0), 1.0).unwrap();
        arc.transform(&Transformation::from_translation((2.0, 3.0)).rotate(FRAC_PI_2));
        assert!((arc.radius - 0.5).abs() < 1.0e-5);
        assert!((arc.sweep - PI).abs() < 1.0e-5);
        assert!(arc.bbox.collides_with(&arc.start));
        assert!(arc.bbox.collides_with(&arc.end));
    }
}
