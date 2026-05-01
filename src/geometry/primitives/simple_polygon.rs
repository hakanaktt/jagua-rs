use std::borrow::Borrow;

use itertools::Itertools;
use ordered_float::{NotNan, OrderedFloat};

use crate::geometry::Transformation;
use crate::geometry::convex_hull::convex_hull_from_points;
use crate::geometry::fail_fast::{SPSurrogate, SPSurrogateConfig, compute_pole};
use crate::geometry::geo_enums::GeoPosition;
use crate::geometry::geo_traits::{
    CollidesWith, DistanceTo, SeparationDistance, Transformable, TransformableFrom,
};
use crate::geometry::primitives::Arc;
use crate::geometry::primitives::BoundarySegment;
use crate::geometry::primitives::Circle;
use crate::geometry::primitives::Edge;
use crate::geometry::primitives::Point;
use crate::geometry::primitives::Rect;
use crate::util::FPA;
use anyhow::{Result, bail};

const BULGE_EPSILON: f64 = 1.0e-6;
const AREA_EPSILON: f64 = 1.0e-6;
const POINT_EPSILON: f64 = 1.0e-5;

/// A Simple Polygon is a polygon that does not intersect itself and contains no holes.
/// It is a closed shape with a finite number of vertices and edges.
/// [read more](https://en.wikipedia.org/wiki/Simple_polygon)
#[derive(Clone, Debug)]
pub struct SPolygon {
    /// Set of points that form the polygon
    pub vertices: Vec<Point>,
    /// Bulge values for each boundary segment. `bulges[i]` describes the segment from `vertices[i]` to the next vertex.
    pub bulges: Vec<f64>,
    /// Bounding box
    pub bbox: Rect,
    /// Area of its interior
    pub area: f64,
    /// Maximum distance between any two points in the polygon
    pub diameter: f64,
    /// [Pole of inaccessibility](https://en.wikipedia.org/wiki/Pole_of_inaccessibility) represented as a circle
    pub poi: Circle,
    /// Optional surrogate representation of the polygon (subset of the original)
    pub surrogate: Option<SPSurrogate>,
}

impl SPolygon {
    /// Create a new simple polygon from a set of points, expensive operations are performed here! Use [Self::clone()] or [Self::transform()] to avoid recomputation.
    pub fn new(points: Vec<Point>) -> Result<Self> {
        let bulges = vec![0.0; points.len()];
        Self::new_with_bulges(points, bulges)
    }

    /// Create a new simple polygon from vertices and one bulge value per outgoing boundary segment.
    pub fn new_with_bulges(mut points: Vec<Point>, mut bulges: Vec<f64>) -> Result<Self> {
        if points.len() < 3 {
            bail!("Simple polygon must have at least 3 points: {points:?}");
        }
        if bulges.len() != points.len() {
            bail!(
                "Simple polygon must have one bulge per vertex: {} points, {} bulges",
                points.len(),
                bulges.len()
            );
        }
        if bulges.iter().any(|bulge| !bulge.is_finite()) {
            bail!("Simple polygon contains invalid bulge values: {bulges:?}");
        }
        if points.iter().unique().count() != points.len() {
            bail!("Simple polygon should not contain duplicate points: {points:?}");
        }

        let area = match SPolygon::calculate_boundary_area(&points, &bulges)? {
            area if area.abs() <= AREA_EPSILON => bail!("Simple polygon has no area: {points:?}"),
            area if area < 0.0 => {
                //edges should always be ordered counterclockwise (positive area)
                reverse_boundary(&mut points, &mut bulges);
                -area
            }
            area => area,
        };

        let diameter = SPolygon::calculate_boundary_diameter(&points, &bulges)?;
        let bbox = SPolygon::generate_boundary_bounding_box(&points, &bulges)?;
        let poi = SPolygon::calculate_poi_with_bulges(&points, &bulges, diameter)?;

        Ok(SPolygon {
            vertices: points,
            bulges,
            bbox,
            area,
            diameter,
            poi,
            surrogate: None,
        })
    }

    pub fn generate_surrogate(&mut self, config: SPSurrogateConfig) -> Result<()> {
        //regenerate the surrogate if it is not present or if the config has changed
        match &self.surrogate {
            Some(surrogate) if surrogate.config == config => {}
            _ => self.surrogate = Some(SPSurrogate::new(self, config)?),
        }
        Ok(())
    }

    pub fn vertex(&self, i: usize) -> Point {
        self.vertices[i]
    }

    pub fn edge(&self, i: usize) -> Edge {
        assert!(i < self.n_vertices(), "index out of bounds");
        let j = if i == self.n_vertices() - 1 { 0 } else { i + 1 };
        Edge {
            start: self.vertices[i],
            end: self.vertices[j],
        }
    }

    pub fn segment(&self, i: usize) -> BoundarySegment {
        assert!(i < self.n_segments(), "index out of bounds");
        let edge = self.edge(i);
        let bulge = self.bulges[i];
        if bulge.abs() <= BULGE_EPSILON {
            BoundarySegment::Line(edge)
        } else {
            BoundarySegment::Arc(
                Arc::try_from_bulge(edge.start, edge.end, bulge)
                    .expect("valid SPolygon bulge should define a valid arc"),
            )
        }
    }

    pub fn edge_iter(&self) -> impl Iterator<Item = Edge> + '_ {
        (0..self.n_vertices()).map(move |i| self.edge(i))
    }

    pub fn segment_iter(&self) -> impl Iterator<Item = BoundarySegment> + '_ {
        (0..self.n_segments()).map(move |i| self.segment(i))
    }

    pub fn tessellated_edge_iter(&self, tolerance: f64) -> impl Iterator<Item = Edge> + '_ {
        self.segment_iter()
            .flat_map(move |segment| tessellate_segment(segment, tolerance).into_iter())
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn n_segments(&self) -> usize {
        self.bulges.len()
    }

    pub fn has_arcs(&self) -> bool {
        self.bulges.iter().any(|bulge| bulge.abs() > BULGE_EPSILON)
    }

    pub fn boundary_points(&self) -> Vec<Point> {
        boundary_points_from_parts(&self.vertices, &self.bulges)
            .expect("valid SPolygon boundary should produce boundary points")
    }

    pub fn surrogate(&self) -> &SPSurrogate {
        self.surrogate.as_ref().expect("surrogate not generated")
    }

    pub fn calculate_diameter(points: Vec<Point>) -> f64 {
        //The two points furthest apart must be part of the convex hull
        let ch = convex_hull_from_points(points);

        //go through all pairs of points and find the pair with the largest distance
        let sq_diam = ch
            .iter()
            .tuple_combinations()
            .map(|(p1, p2)| p1.sq_distance_to(p2))
            .max_by_key(|sq_d| NotNan::new(*sq_d).unwrap())
            .expect("convex hull is empty");

        sq_diam.sqrt()
    }

    pub fn calculate_boundary_diameter(points: &[Point], bulges: &[f64]) -> Result<f64> {
        Ok(SPolygon::calculate_diameter(boundary_points_from_parts(
            points, bulges,
        )?))
    }

    pub fn generate_bounding_box(points: &[Point]) -> Rect {
        let (mut x_min, mut y_min) = (f64::MAX, f64::MAX);
        let (mut x_max, mut y_max) = (f64::MIN, f64::MIN);

        for point in points.iter() {
            x_min = x_min.min(point.0);
            y_min = y_min.min(point.1);
            x_max = x_max.max(point.0);
            y_max = y_max.max(point.1);
        }
        Rect::try_new(x_min, y_min, x_max, y_max).unwrap()
    }

    pub fn generate_boundary_bounding_box(points: &[Point], bulges: &[f64]) -> Result<Rect> {
        Ok(SPolygon::generate_bounding_box(
            &boundary_points_from_parts(points, bulges)?,
        ))
    }

    //https://en.wikipedia.org/wiki/Shoelace_formula
    //counterclockwise = positive area, clockwise = negative area
    pub fn calculate_area(points: &[Point]) -> f64 {
        let mut sigma: f64 = 0.0;
        for i in 0..points.len() {
            //next point
            let j = (i + 1) % points.len();

            let (x_i, y_i) = points[i].into();
            let (x_j, y_j) = points[j].into();

            sigma += (y_i + y_j) * (x_i - x_j)
        }

        0.5 * sigma
    }

    pub fn calculate_boundary_area(points: &[Point], bulges: &[f64]) -> Result<f64> {
        let chord_area = SPolygon::calculate_area(points);
        let segment_area = arc_segments(points, bulges)?
            .into_iter()
            .map(|arc| 0.5 * arc.radius.powi(2) * (arc.sweep - arc.sweep.sin()))
            .sum::<f64>();
        Ok(chord_area + segment_area)
    }

    pub fn calculate_poi(points: &[Point], diameter: f64) -> Result<Circle> {
        let bulges = vec![0.0; points.len()];
        Self::calculate_poi_with_bulges(points, &bulges, diameter)
    }

    pub fn calculate_poi_with_bulges(
        points: &[Point],
        bulges: &[f64],
        diameter: f64,
    ) -> Result<Circle> {
        //need to make a dummy simple polygon, because the pole generation algorithm
        //relies on many of the methods provided by the simple polygon struct
        let dummy_sp = {
            let bbox = SPolygon::generate_boundary_bounding_box(points, bulges)?;
            let area = SPolygon::calculate_boundary_area(points, bulges)?;
            let dummy_poi = Circle::try_new(Point(f64::MAX, f64::MAX), f64::MAX).unwrap();

            SPolygon {
                vertices: points.to_vec(),
                bulges: bulges.to_vec(),
                bbox,
                area,
                diameter,
                poi: dummy_poi,
                surrogate: None,
            }
        };

        compute_pole(&dummy_sp, &[])
    }

    pub fn centroid(&self) -> Point {
        calculate_boundary_centroid(&self.vertices, &self.bulges, self.area)
    }
}

impl Transformable for SPolygon {
    fn transform(&mut self, t: &Transformation) -> &mut Self {
        //destructuring pattern to ensure that the code is updated when the struct changes
        let SPolygon {
            vertices: points,
            bulges,
            bbox,
            area: _,
            diameter: _,
            poi,
            surrogate,
        } = self;

        //transform all points of the simple poly
        points.iter_mut().for_each(|p| {
            p.transform(t);
        });

        poi.transform(t);

        //transform the surrogate
        if let Some(surrogate) = surrogate.as_mut() {
            surrogate.transform(t);
        }

        //regenerate bounding box
        *bbox = SPolygon::generate_boundary_bounding_box(points, bulges)
            .expect("transformed polygon should have a valid bounding box");

        self
    }
}

impl TransformableFrom for SPolygon {
    fn transform_from(&mut self, reference: &Self, t: &Transformation) -> &mut Self {
        //destructuring pattern to ensure that the code is updated when the struct changes
        let SPolygon {
            vertices: points,
            bulges,
            bbox,
            area: _,
            diameter: _,
            poi,
            surrogate,
        } = self;

        for (p, ref_p) in points.iter_mut().zip(&reference.vertices) {
            p.transform_from(ref_p, t);
        }
        *bulges = reference.bulges.clone();

        poi.transform_from(&reference.poi, t);

        //transform the surrogate
        if let Some(surrogate) = surrogate.as_mut() {
            surrogate.transform_from(reference.surrogate(), t);
        }
        //regenerate bounding box
        *bbox = SPolygon::generate_boundary_bounding_box(points, bulges)
            .expect("transformed polygon should have a valid bounding box");

        self
    }
}

impl CollidesWith<Point> for SPolygon {
    fn collides_with(&self, point: &Point) -> bool {
        //based on the ray casting algorithm: https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
        match self.bbox.collides_with(point) {
            false => false,
            true => {
                let ray = horizontal_ray_from_point(*point, self.bbox);
                let mut intersections = Vec::new();
                for segment in self.segment_iter() {
                    if segment.collides_with(point) {
                        return true;
                    }

                    match segment {
                        BoundarySegment::Line(edge) => {
                            if let Some(intersection) = ray_edge_crossing(point, &edge) {
                                push_unique_intersection(&mut intersections, intersection);
                            }
                        }
                        BoundarySegment::Arc(arc) => {
                            for intersection in ray_arc_crossings(point, &ray, &arc) {
                                push_unique_intersection(&mut intersections, intersection);
                            }
                        }
                    }
                }
                intersections.len() % 2 == 1
            }
        }
    }
}

impl DistanceTo<Point> for SPolygon {
    fn distance_to(&self, point: &Point) -> f64 {
        self.sq_distance_to(point).sqrt()
    }
    fn sq_distance_to(&self, point: &Point) -> f64 {
        match self.collides_with(point) {
            true => 0.0,
            false => self
                .segment_iter()
                .map(|segment| segment.sq_distance_to(point))
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap(),
        }
    }
}

impl SeparationDistance<Point> for SPolygon {
    fn separation_distance(&self, point: &Point) -> (GeoPosition, f64) {
        let (position, sq_distance) = self.sq_separation_distance(point);
        (position, sq_distance.sqrt())
    }

    fn sq_separation_distance(&self, point: &Point) -> (GeoPosition, f64) {
        let distance_to_closest_edge = self
            .segment_iter()
            .map(|segment| segment.sq_distance_to(point))
            .min_by_key(|sq_d| OrderedFloat(*sq_d))
            .unwrap();

        match self.collides_with(point) {
            true => (GeoPosition::Interior, distance_to_closest_edge),
            false => (GeoPosition::Exterior, distance_to_closest_edge),
        }
    }
}

fn boundary_points_from_parts(points: &[Point], bulges: &[f64]) -> Result<Vec<Point>> {
    let mut boundary_points = points.to_vec();
    for arc in arc_segments(points, bulges)? {
        boundary_points.extend(arc.cardinal_extrema());
    }
    Ok(boundary_points)
}

fn arc_segments(points: &[Point], bulges: &[f64]) -> Result<Vec<Arc>> {
    validate_boundary_parts(points, bulges)?;
    let mut arcs = Vec::new();
    for i in 0..points.len() {
        if bulges[i].abs() > BULGE_EPSILON {
            let j = next_index(i, points.len());
            arcs.push(Arc::try_from_bulge(points[i], points[j], bulges[i])?);
        }
    }
    Ok(arcs)
}

fn validate_boundary_parts(points: &[Point], bulges: &[f64]) -> Result<()> {
    if points.len() != bulges.len() {
        bail!(
            "Simple polygon must have one bulge per vertex: {} points, {} bulges",
            points.len(),
            bulges.len()
        );
    }
    Ok(())
}

fn reverse_boundary(points: &mut Vec<Point>, bulges: &mut Vec<f64>) {
    let n_points = points.len();
    let old_bulges = bulges.clone();
    points.reverse();
    for (i, bulge) in bulges.iter_mut().enumerate() {
        *bulge = -old_bulges[(n_points + n_points - 2 - i) % n_points];
    }
}

fn next_index(i: usize, len: usize) -> usize {
    if i == len - 1 { 0 } else { i + 1 }
}

fn calculate_boundary_centroid(points: &[Point], bulges: &[f64], area: f64) -> Point {
    let chord_area = SPolygon::calculate_area(points);
    let chord_centroid = calculate_chord_centroid(points, chord_area);
    let mut moment_x = chord_centroid.0 * chord_area;
    let mut moment_y = chord_centroid.1 * chord_area;

    for arc in arc_segments(points, bulges).expect("valid SPolygon boundary") {
        let segment_area = 0.5 * arc.radius.powi(2) * (arc.sweep - arc.sweep.sin());
        if segment_area.abs() <= AREA_EPSILON {
            continue;
        }

        let sweep_abs = arc.sweep.abs();
        let denominator = 3.0 * (sweep_abs - sweep_abs.sin());
        if denominator.abs() <= AREA_EPSILON {
            continue;
        }
        let centroid_offset = 4.0 * arc.radius * (sweep_abs / 2.0).sin().powi(3) / denominator;
        let mid_angle = arc.start_angle() + arc.sweep / 2.0;
        let segment_centroid = Point(
            arc.center.0 + centroid_offset * mid_angle.cos(),
            arc.center.1 + centroid_offset * mid_angle.sin(),
        );
        moment_x += segment_centroid.0 * segment_area;
        moment_y += segment_centroid.1 * segment_area;
    }

    Point(moment_x / area, moment_y / area)
}

fn calculate_chord_centroid(points: &[Point], area: f64) -> Point {
    if area.abs() <= AREA_EPSILON {
        return points[0];
    }

    let mut c_x = 0.0;
    let mut c_y = 0.0;
    for i in 0..points.len() {
        let j = next_index(i, points.len());
        let Point(x_i, y_i) = points[i];
        let Point(x_j, y_j) = points[j];
        let cross = x_i * y_j - x_j * y_i;
        c_x += (x_i + x_j) * cross;
        c_y += (y_i + y_j) * cross;
    }
    Point(c_x / (6.0 * area), c_y / (6.0 * area))
}

fn tessellate_segment(segment: BoundarySegment, tolerance: f64) -> Vec<Edge> {
    match segment {
        BoundarySegment::Line(edge) => vec![edge],
        BoundarySegment::Arc(arc) => {
            let tolerance = tolerance.max(1.0e-4).min(arc.radius);
            let max_sweep = 2.0 * (1.0 - tolerance / arc.radius).acos();
            let n_segments = ((arc.sweep.abs() / max_sweep).ceil() as usize).max(1);
            (0..n_segments)
                .map(|i| {
                    let t0 = i as f64 / n_segments as f64;
                    let t1 = (i + 1) as f64 / n_segments as f64;
                    Edge::try_new(
                        arc.point_at_angle(arc.start_angle() + arc.sweep * t0),
                        arc.point_at_angle(arc.start_angle() + arc.sweep * t1),
                    )
                    .expect("arc tessellation should not create degenerate edges")
                })
                .collect()
        }
    }
}

fn horizontal_ray_from_point(point: Point, bbox: Rect) -> Edge {
    Edge {
        start: point,
        end: Point(bbox.x_max + bbox.width(), point.1),
    }
}

fn ray_edge_crossing(point: &Point, edge: &Edge) -> Option<Point> {
    let Point(_, point_y) = *point;
    let Point(start_x, start_y) = edge.start;
    let Point(end_x, end_y) = edge.end;

    if (start_y > point_y) == (end_y > point_y) {
        return None;
    }

    let t = (point_y - start_y) / (end_y - start_y);
    let intersection = Point(start_x + t * (end_x - start_x), point_y);
    (FPA(intersection.0) > FPA(point.0)).then_some(intersection)
}

fn ray_arc_crossings(point: &Point, ray: &Edge, arc: &Arc) -> Vec<Point> {
    arc.collides_at(ray)
        .into_iter()
        .filter(|intersection| FPA(intersection.0) > FPA(point.0))
        .filter(|intersection| !is_horizontal_tangent(arc, intersection))
        .collect()
}

fn is_horizontal_tangent(arc: &Arc, point: &Point) -> bool {
    (point.0 - arc.center.0).abs() <= POINT_EPSILON.max(arc.radius * POINT_EPSILON)
}

fn push_unique_intersection(intersections: &mut Vec<Point>, point: Point) {
    if intersections
        .iter()
        .all(|existing| existing.sq_distance_to(&point) > POINT_EPSILON.powi(2))
    {
        intersections.push(point);
    }
}

#[cfg(feature = "curves")]
impl SPolygon {
    pub fn to_cavalier_vertices(&self) -> Vec<cavalier_contours::polyline::PlineVertex<f64>> {
        self.vertices
            .iter()
            .zip(&self.bulges)
            .map(|(point, bulge)| {
                cavalier_contours::polyline::PlineVertex::new(point.0, point.1, *bulge)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-3,
            "actual={actual}, expected={expected}"
        );
    }

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
    fn straight_polygon_segments_remain_chord_edges() {
        let shape = SPolygon::new(vec![
            Point(0.0, 0.0),
            Point(2.0, 0.0),
            Point(2.0, 1.0),
            Point(0.0, 1.0),
        ])
        .unwrap();

        assert_eq!(shape.n_vertices(), 4);
        assert_eq!(shape.n_segments(), 4);
        assert!(!shape.has_arcs());
        assert!(shape.bulges.iter().all(|bulge| *bulge == 0.0));
        assert!(matches!(shape.segment(0), BoundarySegment::Line(_)));
        assert_eq!(shape.edge_iter().count(), shape.segment_iter().count());
        assert_close(shape.area, 2.0);
        assert_close(shape.centroid().0, 1.0);
        assert_close(shape.centroid().1, 0.5);
    }

    #[test]
    fn bulged_polygon_derived_fields_include_arc_extrema() {
        let shape = rounded_top_shape();

        assert!(shape.has_arcs());
        assert!(matches!(shape.segment(2), BoundarySegment::Arc(_)));
        assert_close(shape.area, 4.0 + PI / 2.0);
        assert_close(shape.bbox.x_min, -1.0);
        assert_close(shape.bbox.x_max, 1.0);
        assert_close(shape.bbox.y_min, 0.0);
        assert_close(shape.bbox.y_max, 3.0);
        assert_close(shape.diameter, 10.0_f64.sqrt());

        let centroid = shape.centroid();
        let semicircle_area = PI / 2.0;
        let semicircle_centroid_y = 2.0 + 4.0 / (3.0 * PI);
        let expected_centroid_y =
            (4.0 * 1.0 + semicircle_area * semicircle_centroid_y) / (4.0 + semicircle_area);
        assert_close(centroid.0, 0.0);
        assert_close(centroid.1, expected_centroid_y);

        assert!(shape.collides_with(&Point(0.0, 2.5)));
        assert!(!shape.collides_with(&Point(0.0, 3.2)));
    }

    #[test]
    fn bulged_polygon_point_inclusion_ignores_tangent_ray_touch() {
        let shape = rounded_top_shape();

        assert!(shape.collides_with(&Point(0.0, 3.0)));
        assert!(!shape.collides_with(&Point(-0.5, 3.0)));
    }

    #[test]
    fn bulged_polygon_tessellation_keeps_segment_order() {
        let shape = rounded_top_shape();
        let tessellated = shape.tessellated_edge_iter(0.01).collect_vec();

        assert!(tessellated.len() > shape.n_segments());
        assert_eq!(tessellated.first().unwrap().start, shape.vertex(0));
        assert_eq!(tessellated.last().unwrap().end, shape.vertex(0));
    }

    #[test]
    fn clockwise_bulged_polygon_reverses_vertices_and_flips_bulges() {
        let shape = SPolygon::new_with_bulges(
            vec![
                Point(-1.0, 2.0),
                Point(1.0, 2.0),
                Point(1.0, 0.0),
                Point(-1.0, 0.0),
            ],
            vec![-1.0, 0.0, 0.0, 0.0],
        )
        .unwrap();

        assert_eq!(shape.vertices, rounded_top_shape().vertices);
        assert_eq!(shape.bulges, vec![0.0, 0.0, 1.0, 0.0]);
        assert_close(shape.area, 4.0 + PI / 2.0);
    }
}

impl<T> From<T> for SPolygon
where
    T: Borrow<Rect>,
{
    fn from(r: T) -> Self {
        let r = r.borrow();
        SPolygon::new(vec![
            (r.x_min, r.y_min).into(),
            (r.x_max, r.y_min).into(),
            (r.x_max, r.y_max).into(),
            (r.x_min, r.y_max).into(),
        ])
        .unwrap()
    }
}
