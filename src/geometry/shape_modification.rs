use itertools::Itertools;
use log::{debug, error, info, warn};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

use crate::geometry::geo_traits::{CollidesWith, DistanceTo};
use crate::geometry::primitives::Edge;
use crate::geometry::primitives::Point;
use crate::geometry::primitives::SPolygon;

use crate::io::ext_repr::ExtSPolygon;
use crate::io::import;
use anyhow::{Result, bail};

#[cfg(feature = "curves")]
use cavalier_contours::polyline::{PlineSource, Polyline};

#[cfg(feature = "curves")]
const CAVALIER_MAX_BULGE: f64 = 1.0;
#[cfg(feature = "curves")]
const CAVALIER_BULGE_EPSILON: f64 = 1.0e-6;

/// Whether to strictly inflate or deflate when making any modifications to shape.
/// Depends on the [`position`](crate::collision_detection::hazards::HazardEntity::scope) of the [`HazardEntity`](crate::collision_detection::hazards::HazardEntity) that the shape represents.
#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShapeModifyMode {
    /// Modify the shape to be strictly larger than the original (superset).
    Inflate,
    /// Modify the shape to be strictly smaller than the original (subset).
    Deflate,
}

#[derive(Serialize, Deserialize, Clone, Copy, Debug, Default, PartialEq)]
pub struct ShapeModifyConfig {
    /// Maximum deviation of the simplified polygon with respect to the original polygon area as a ratio.
    /// If undefined, no simplification is performed.
    /// See [`simplify_shape`]
    pub simplify_tolerance: Option<f64>,
    /// Offset by which to inflate or deflate the polygon.
    /// If undefined, no offset is applied.
    /// See [`offset_shape`]
    pub offset: Option<f64>,
    /// Definition for narrow concavities that can be closed by a straight edge.
    /// Defined as a tuple of (max_distance_ratio, max_area_ratio) where:
    /// - max_distance_ratio: maximum distance between two vertices of a polygon to consider it a narrow concavity, defined as a fraction of the item's diameter.
    /// - max_area_ratio: maximum area of the sub-shape formed by the vertices between the two vertices, defined as a fraction of the item's area.
    ///
    /// If undefined, no narrow concavities will be closed.
    /// See [`close_narrow_concavities`]
    pub narrow_concavity_cutoff: Option<(f64, f64)>,
}

/// Simplifies a [`SPolygon`] by reducing the number of edges.
///
/// The simplified shape will either be a subset or a superset of the original shape, depending on the [`ShapeModifyMode`].
/// The procedure sequentially eliminates edges until either the change in area (ratio)
/// exceeds `max_area_delta` or the number of edges < 4.
pub fn simplify_shape(
    shape: &SPolygon,
    mode: ShapeModifyMode,
    max_area_change_ratio: f64,
) -> SPolygon {
    if shape.has_arcs() {
        warn!(
            "[PS] skipping simplification for arc-bearing shape; arc-aware simplification is deferred"
        );
        return shape.clone();
    }

    let original_area = shape.area;

    let mut ref_points = shape.vertices.clone();

    for _ in 0..shape.n_vertices() {
        let n_points = ref_points.len() as isize;
        if n_points < 4 {
            //can't simplify further
            break;
        }

        let mut corners = (0..n_points)
            .map(|i| {
                let i_prev = (i - 1).rem_euclid(n_points);
                let i_next = (i + 1).rem_euclid(n_points);
                Corner(i_prev as usize, i as usize, i_next as usize)
            })
            .collect_vec();

        if mode == ShapeModifyMode::Deflate {
            //default mode is to inflate, so we need to reverse the order of the corners and flip the corners for deflate mode
            //reverse the order of the corners
            corners.reverse();
            //reverse each corner
            corners.iter_mut().for_each(|c| c.flip());
        }

        let mut candidates = vec![];

        let mut prev_corner = corners.last().expect("corners is empty");
        let mut prev_corner_type = CornerType::from(prev_corner.to_points(&ref_points));

        //Go over all corners and generate candidates
        for corner in corners.iter() {
            let corner_type = CornerType::from(corner.to_points(&ref_points));

            //Generate a removal candidate (or not)
            match (&corner_type, &prev_corner_type) {
                (CornerType::Concave, _) => candidates.push(Candidate::Concave(*corner)),
                (CornerType::Collinear, _) => candidates.push(Candidate::Collinear(*corner)),
                (CornerType::Convex, CornerType::Convex) => {
                    candidates.push(Candidate::ConvexConvex(*prev_corner, *corner))
                }
                (_, _) => {}
            };
            (prev_corner, prev_corner_type) = (corner, corner_type);
        }

        //search the candidate with the smallest change in area that is valid
        let best_candidate = candidates
            .iter()
            .sorted_by_cached_key(|c| {
                OrderedFloat(calculate_area_delta(&ref_points, c).unwrap_or(f64::INFINITY))
            })
            .find(|c| candidate_is_valid(&ref_points, c));

        //if it is within the area change constraints, execute the candidate
        if let Some(best_candidate) = best_candidate {
            let new_shape = execute_candidate(&ref_points, best_candidate);
            let new_shape_area = SPolygon::calculate_area(&new_shape);
            let area_delta = (new_shape_area - original_area).abs() / original_area;
            if area_delta <= max_area_change_ratio {
                debug!(
                    "[PS] executed {:?} simplification causing {:.2}% area change",
                    best_candidate,
                    area_delta * 100.0
                );
                ref_points = new_shape;
            } else {
                break; //area change too significant
            }
        } else {
            break; //no candidate found
        }
    }

    //Convert it back to a simple polygon
    let simpl_shape = SPolygon::new(ref_points).unwrap();

    if simpl_shape.n_vertices() < shape.n_vertices() {
        info!(
            "[PS] simplified from {} to {} edges with {:.3}% area difference",
            shape.n_vertices(),
            simpl_shape.n_vertices(),
            (simpl_shape.area - shape.area) / shape.area * 100.0
        );
    } else {
        info!("[PS] no simplification possible within area change constraints");
    }

    simpl_shape
}

fn calculate_area_delta(shape: &[Point], candidate: &Candidate) -> Result<f64, InvalidCandidate> {
    //calculate the difference in area of the shape if the candidate were to be executed
    let area = match candidate {
        Candidate::Collinear(_) => 0.0,
        Candidate::Concave(c) => {
            //Triangle formed by i_prev, i and i_next will correspond to the change area
            let Point(x0, y0) = shape[c.0];
            let Point(x1, y1) = shape[c.1];
            let Point(x2, y2) = shape[c.2];

            let area = (x0 * y1 + x1 * y2 + x2 * y0 - x0 * y2 - x1 * y0 - x2 * y1) / 2.0;

            area.abs()
        }
        Candidate::ConvexConvex(c1, c2) => {
            let replacing_vertex = replacing_vertex_convex_convex_candidate(shape, (*c1, *c2))?;

            //the triangle formed by corner c1, c2, and replacing vertex will correspond to the change in area
            let Point(x0, y0) = shape[c1.1];
            let Point(x1, y1) = replacing_vertex;
            let Point(x2, y2) = shape[c2.1];

            let area = (x0 * y1 + x1 * y2 + x2 * y0 - x0 * y2 - x1 * y0 - x2 * y1) / 2.0;

            area.abs()
        }
    };
    Ok(area)
}

fn candidate_is_valid(shape: &[Point], candidate: &Candidate) -> bool {
    //ensure the removal/replacement does not create any self intersections
    match candidate {
        Candidate::Collinear(_) => true,
        Candidate::Concave(c) => {
            let new_edge = Edge::try_new(shape[c.0], shape[c.2]).unwrap();
            let affected_points = [shape[c.0], shape[c.1], shape[c.2]];

            //check for self-intersections
            edge_iter(shape)
                .filter(|l| !affected_points.contains(&l.start))
                .filter(|l| !affected_points.contains(&l.end))
                .all(|l| !l.collides_with(&new_edge))
        }
        Candidate::ConvexConvex(c1, c2) => {
            match replacing_vertex_convex_convex_candidate(shape, (*c1, *c2)) {
                Err(_) => false,
                Ok(new_vertex) => {
                    let new_edge_1 = Edge::try_new(shape[c1.0], new_vertex).unwrap();
                    let new_edge_2 = Edge::try_new(new_vertex, shape[c2.2]).unwrap();

                    let affected_points = [shape[c1.1], shape[c1.0], shape[c2.1], shape[c2.2]];

                    //check for self-intersections
                    edge_iter(shape)
                        .filter(|l| !affected_points.contains(&l.start))
                        .filter(|l| !affected_points.contains(&l.end))
                        .all(|l| !l.collides_with(&new_edge_1) && !l.collides_with(&new_edge_2))
                }
            }
        }
    }
}

fn edge_iter(points: &[Point]) -> impl Iterator<Item = Edge> + '_ {
    let n_points = points.len();
    (0..n_points).map(move |i| {
        let j = (i + 1) % n_points;
        Edge::try_new(points[i], points[j]).unwrap()
    })
}

fn execute_candidate(shape: &[Point], candidate: &Candidate) -> Vec<Point> {
    let mut points = shape.iter().cloned().collect_vec();
    match candidate {
        Candidate::Collinear(c) | Candidate::Concave(c) => {
            points.remove(c.1);
        }
        Candidate::ConvexConvex(c1, c2) => {
            let replacing_vertex = replacing_vertex_convex_convex_candidate(shape, (*c1, *c2))
                .expect("invalid candidate cannot be executed");
            points.remove(c1.1);
            let other_index = if c1.1 < c2.1 { c2.1 - 1 } else { c2.1 };
            points.remove(other_index);
            points.insert(other_index, replacing_vertex);
        }
    }
    points
}

fn replacing_vertex_convex_convex_candidate(
    shape: &[Point],
    (c1, c2): (Corner, Corner),
) -> Result<Point, InvalidCandidate> {
    assert_eq!(c1.2, c2.1, "non-consecutive corners {c1:?},{c2:?}");
    assert_eq!(c1.1, c2.0, "non-consecutive corners {c1:?},{c2:?}");

    let edge_prev = Edge::try_new(shape[c1.0], shape[c1.1]).unwrap();
    let edge_next = Edge::try_new(shape[c2.2], shape[c2.1]).unwrap();

    calculate_intersection_in_front(&edge_prev, &edge_next).ok_or(InvalidCandidate)
}

fn calculate_intersection_in_front(l1: &Edge, l2: &Edge) -> Option<Point> {
    //Calculates the intersection point between l1 and l2 if both were extended in front to infinity.

    //https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment
    //vector 1 = [(x1,y1),(x2,y2)[ and vector 2 = [(x3,y3),(x4,y4)[
    let Point(x1, y1) = l1.start;
    let Point(x2, y2) = l1.end;
    let Point(x3, y3) = l2.start;
    let Point(x4, y4) = l2.end;

    //used formula is slightly different to the one on wikipedia. The orientation of the line segments are flipped
    //We consider an intersection if t == ]0,1] && u == ]0,1]

    let t_nom = (x2 - x4) * (y4 - y3) - (y2 - y4) * (x4 - x3);
    let t_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);

    let u_nom = (x2 - x4) * (y2 - y1) - (y2 - y4) * (x2 - x1);
    let u_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);

    let t = if t_denom != 0.0 { t_nom / t_denom } else { 0.0 };

    let u = if u_denom != 0.0 { u_nom / u_denom } else { 0.0 };

    if t < 0.0 && u < 0.0 {
        //intersection is in front both vectors
        Some(Point(x2 + t * (x1 - x2), y2 + t * (y1 - y2)))
    } else {
        //no intersection (parallel or not in front)
        None
    }
}

#[derive(Debug, Clone)]
struct InvalidCandidate;

#[derive(Clone, Debug, PartialEq)]
enum Candidate {
    Concave(Corner),
    ConvexConvex(Corner, Corner),
    Collinear(Corner),
}

#[derive(Clone, Copy, Debug, PartialEq)]
///Corner is defined as the left hand side of points 0-1-2
struct Corner(pub usize, pub usize, pub usize);

impl Corner {
    pub fn flip(&mut self) {
        std::mem::swap(&mut self.0, &mut self.2);
    }

    pub fn to_points(self, points: &[Point]) -> [Point; 3] {
        [points[self.0], points[self.1], points[self.2]]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum CornerType {
    Concave,
    Convex,
    Collinear,
}

impl CornerType {
    pub fn from([p1, p2, p3]: [Point; 3]) -> Self {
        //returns the corner type on the left-hand side p1->p2->p3
        //From: https://algorithmtutor.com/Computational-Geometry/Determining-if-two-consecutive-segments-turn-left-or-right/

        let p1p2 = (p2.0 - p1.0, p2.1 - p1.1);
        let p1p3 = (p3.0 - p1.0, p3.1 - p1.1);
        let cross_prod = p1p2.0 * p1p3.1 - p1p2.1 * p1p3.0;

        //a positive cross product indicates that p2p3 turns to the left with respect to p1p2
        match cross_prod.partial_cmp(&0.0).expect("cross product is NaN") {
            Ordering::Less => CornerType::Concave,
            Ordering::Equal => CornerType::Collinear,
            Ordering::Greater => CornerType::Convex,
        }
    }
}

/// Offsets a [`SPolygon`] by a certain `distance` either inwards or outwards depending on the [`ShapeModifyMode`].
/// Straight polygons use the [`geo_offset`](https://crates.io/crates/geo_offset) crate; arc-bearing polygons use Cavalier Contours when the `curves` feature is enabled.
pub fn offset_shape(sp: &SPolygon, mode: ShapeModifyMode, distance: f64) -> Result<SPolygon> {
    if sp.has_arcs() {
        return offset_curved_shape(sp, mode, distance);
    }

    let offset = match mode {
        ShapeModifyMode::Deflate => -distance,
        ShapeModifyMode::Inflate => distance,
    };

    // Convert the SPolygon to a geo_types::Polygon
    let geo_poly = geo_types::Polygon::new(
        sp.vertices
            .iter()
            .map(|p| (p.0 as f64, p.1 as f64))
            .collect(),
        vec![],
    );

    // Create the offset polygon
    let geo_poly_offsets = geo_buffer::buffer_polygon_rounded(&geo_poly, offset as f64).0;

    let geo_poly_offset = match geo_poly_offsets.len() {
        0 => bail!("Offset resulted in an empty polygon"),
        1 => &geo_poly_offsets[0],
        _ => {
            // If there are multiple polygons, we take the first one.
            // This can happen if the offset creates multiple disconnected parts.
            warn!("Offset resulted in multiple polygons, taking the first one.");
            &geo_poly_offsets[0]
        }
    };

    // Convert back to internal representation (by using the import function)
    let ext_s_polygon = ExtSPolygon(
        geo_poly_offset
            .exterior()
            .points()
            .map(|p| (p.x() as f64, p.y() as f64))
            .collect_vec(),
    );

    import::import_simple_polygon(&ext_s_polygon)
}

#[cfg(feature = "curves")]
fn offset_curved_shape(sp: &SPolygon, mode: ShapeModifyMode, distance: f64) -> Result<SPolygon> {
    let offset = match mode {
        ShapeModifyMode::Deflate => distance,
        ShapeModifyMode::Inflate => -distance,
    };

    let polyline = cavalier_polyline_from_spolygon(sp)?;
    let offset_polylines = polyline.parallel_offset(offset);
    let offset_polyline = offset_polylines
        .iter()
        .max_by_key(|polyline| OrderedFloat(polyline.area().abs()))
        .ok_or_else(|| anyhow::anyhow!("Offset resulted in an empty polygon"))?;

    if offset_polylines.len() > 1 {
        warn!("Offset resulted in multiple curved polygons, taking the largest one.");
    }

    spolygon_from_cavalier_polyline(offset_polyline)
}

#[cfg(not(feature = "curves"))]
fn offset_curved_shape(_sp: &SPolygon, _mode: ShapeModifyMode, _distance: f64) -> Result<SPolygon> {
    bail!("Offsetting arc-bearing shapes requires the `curves` feature")
}

#[cfg(feature = "curves")]
fn cavalier_polyline_from_spolygon(sp: &SPolygon) -> Result<Polyline<f64>> {
    if sp
        .bulges
        .iter()
        .any(|bulge| bulge.abs() > CAVALIER_MAX_BULGE + CAVALIER_BULGE_EPSILON)
    {
        bail!(
            "Cavalier offsets require arcs with |bulge| <= 1; split larger arcs before offsetting"
        )
    }

    Ok(Polyline {
        vertex_data: sp.to_cavalier_vertices(),
        is_closed: true,
        userdata: Vec::new(),
    })
}

#[cfg(feature = "curves")]
fn spolygon_from_cavalier_polyline(polyline: &Polyline<f64>) -> Result<SPolygon> {
    let mut vertices = polyline
        .vertex_data
        .iter()
        .map(|vertex| Point(vertex.x, vertex.y))
        .collect_vec();
    let mut bulges = polyline
        .vertex_data
        .iter()
        .map(|vertex| vertex.bulge)
        .collect_vec();

    if vertices.len() > 1 && vertices.first() == vertices.last() {
        vertices.pop();
        bulges.pop();
    }

    SPolygon::new_with_bulges(vertices, bulges)
}

/// Closes narrow concavities in a [`SPolygon`] by replacing them with a straight edge, eliminating the vertices in between.
pub fn close_narrow_concavities(
    orig_shape: &SPolygon,
    mode: ShapeModifyMode,
    (cutoff_distance_ratio, cutoff_area_ratio): (f64, f64),
) -> SPolygon {
    if orig_shape.has_arcs() {
        warn!(
            "[PS] skipping narrow-concavity closing for arc-bearing shape; arc-aware simplification is deferred"
        );
        return orig_shape.clone();
    }

    let mut n_concav_closed = 0;
    let mut shape = orig_shape.clone();

    for _ in 0..shape.n_vertices() {
        let n_points = shape.n_vertices();

        let calc_vert_elim = |i, j| {
            if j > i {
                j - i - 1
            } else {
                n_points - i + j - 1
            }
        };

        let mut best_candidate = None;
        for i in 0..n_points {
            for j in 0..n_points {
                if i == j || (i + 1) % n_points == j || (j + 1) % n_points == i {
                    continue; //skip adjacent points
                }
                //Simulate the replacing edge
                let c_edge = Edge::try_new(shape.vertex(i), shape.vertex(j))
                    .expect("invalid edge in string candidate")
                    .scale(0.9999); //slightly shrink the edge to avoid self-intersections

                if c_edge.length() > cutoff_distance_ratio * shape.diameter {
                    //If the edge is too long, skip it
                    continue;
                }

                if mode == ShapeModifyMode::Inflate
                    && (shape.collides_with(&c_edge.start) || shape.collides_with(&c_edge.end))
                {
                    //If we are only allowed to inflate the shape and any end point is inside the shape, skip it
                    continue;
                } else if mode == ShapeModifyMode::Deflate
                    && !(shape.collides_with(&c_edge.start) && shape.collides_with(&c_edge.end))
                {
                    //If we are only allowed to deflate the shape and both end points are not inside the shape, skip it
                    continue;
                }

                if shape.edge_iter().any(|e| e.collides_with(&c_edge)) {
                    //If the edge collides with any edge of the shape, reject always
                    continue;
                }
                //the eliminated vertices should form a negative area (in inflation mode) or positive area (in deflation mode)
                let sub_shape_area = {
                    let sub_shape_points = if j > i {
                        shape.vertices[i..j].to_vec()
                    } else {
                        [&shape.vertices[i..], &shape.vertices[..j]].concat()
                    };
                    SPolygon::calculate_area(&sub_shape_points)
                };
                if sub_shape_area >= 0.0 {
                    //if the area is not negative, skip it
                    continue;
                }
                if sub_shape_area.abs() > cutoff_area_ratio * shape.area {
                    //if the area is too large, skip it
                    continue;
                }

                //Valid candidate found...
                match best_candidate {
                    None => {
                        //first candidate found
                        best_candidate = Some((i, j));
                    }
                    Some((best_i, best_j)) => {
                        //check the number of points that would be removed
                        if calc_vert_elim(i, j) > calc_vert_elim(best_i, best_j) {
                            best_candidate = Some((i, j));
                        }
                    }
                }
            }
        }
        if let Some((i, j)) = best_candidate {
            let mut ref_points = shape.vertices.clone();
            let start = i as isize + 1;
            let end = j as isize - 1;
            debug!(
                "[PS] closing concavity between points (idx: {}, {:?}) and (idx: {}, {:?}) with edge length {:.3} ({} vertices eliminated)",
                i,
                shape.vertex(i),
                j,
                shape.vertex(j),
                Edge::try_new(shape.vertex(i), shape.vertex(j))
                    .expect("invalid edge in string candidate")
                    .length(),
                calc_vert_elim(i, j)
            );
            if start <= end {
                // if j does not wrap around the shape
                ref_points.drain((start as usize)..=(end as usize));
            } else {
                // if j wraps around the shape
                if (start as usize) < n_points {
                    //remove from `start` to back
                    ref_points.drain(start as usize..);
                }
                if end >= 0 {
                    //remove from front to `end`
                    ref_points.drain(0..=(end as usize));
                }
            }
            shape = SPolygon::new(ref_points).expect("invalid shape after closing concavity");
            n_concav_closed += 1;
        } else {
            //no more candidates found, break the loop
            break;
        }
    }

    if n_concav_closed > 0 {
        info!(
            "[PS] [EXPERIMENTAL] closed {} concavities closer than {:.3}% of diameter and less than {:.3}% of area, reducing vertices from {} to {}",
            n_concav_closed,
            cutoff_distance_ratio * 100.0,
            cutoff_area_ratio * 100.0,
            orig_shape.n_vertices(),
            shape.n_vertices()
        );
    }

    shape
}

pub fn shape_modification_valid(orig: &SPolygon, simpl: &SPolygon, mode: ShapeModifyMode) -> bool {
    //make sure each point of the original shape is either in the new shape or included (in case of inflation)/excluded (in case of deflation) in the new shape
    let on_edge = |p: &Point| {
        simpl
            .segment_iter()
            .any(|segment| segment.distance_to(p) < simpl.diameter * 1e-6)
    };

    let orig_boundary_points = orig.boundary_points();
    let simpl_boundary_points = simpl.boundary_points();

    for p in orig_boundary_points
        .iter()
        .filter(|p| !simpl_boundary_points.contains(p))
    {
        let vertex_on_edge = on_edge(p);
        let vertex_in_simpl = simpl.collides_with(p);

        let error = match mode {
            ShapeModifyMode::Inflate => !vertex_in_simpl && !vertex_on_edge,
            ShapeModifyMode::Deflate => vertex_in_simpl && !vertex_on_edge,
        };

        if error {
            error!(
                "[PS] point {:?} from original shape is incorrect in simplified shape (original vertices: {:?}, simplified vertices: {:?})",
                p,
                orig.vertices.iter().map(|p| (p.0, p.1)).collect_vec(),
                simpl.vertices.iter().map(|p| (p.0, p.1)).collect_vec()
            );
            return false; //point is not in the new shape and does not collide with it
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn simplify_and_concavity_closing_skip_arc_bearing_shapes() {
        let shape = rounded_top_shape();

        let simplified = simplify_shape(&shape, ShapeModifyMode::Inflate, 0.5);
        let closed = close_narrow_concavities(&shape, ShapeModifyMode::Inflate, (0.5, 0.5));

        assert_eq!(simplified.vertices, shape.vertices);
        assert_eq!(simplified.bulges, shape.bulges);
        assert_eq!(closed.vertices, shape.vertices);
        assert_eq!(closed.bulges, shape.bulges);
    }

    #[cfg(feature = "curves")]
    #[test]
    fn curved_offset_uses_cavalier_and_preserves_arcs() {
        let shape = rounded_top_shape();

        let inflated = offset_shape(&shape, ShapeModifyMode::Inflate, 0.1).unwrap();
        let deflated = offset_shape(&shape, ShapeModifyMode::Deflate, 0.1).unwrap();

        assert!(inflated.has_arcs());
        assert!(deflated.has_arcs());
        assert!(inflated.area > shape.area);
        assert!(deflated.area < shape.area);
    }

    #[cfg(not(feature = "curves"))]
    #[test]
    fn curved_offset_requires_curves_feature() {
        let shape = rounded_top_shape();

        assert!(offset_shape(&shape, ShapeModifyMode::Inflate, 0.1).is_err());
    }
}
