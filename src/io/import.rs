use crate::collision_detection::CDEConfig;
use crate::entities::Item;
use crate::entities::{Container, InferiorQualityZone, N_QUALITIES};
use crate::geometry::OriginalShape;
use crate::geometry::geo_enums::RotationRange;
use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::Arc;
use crate::geometry::primitives::Point;
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::SPolygon;
use crate::geometry::shape_modification::{ShapeModifyConfig, ShapeModifyMode};
use crate::geometry::{DTransformation, Transformation};
use crate::io::ext_repr::{
    ExtArcSeg, ExtBulgedPolygon, ExtContainer, ExtItem, ExtSPolygon, ExtShape,
};
use anyhow::{Result, bail, ensure};
use float_cmp::approx_eq;
use itertools::Itertools;
use log::debug;
use std::f64::consts::{FRAC_PI_2, PI};

const BULGE_EPSILON: f64 = 1.0e-6;
const POINT_MATCH_EPSILON: f64 = 1.0e-4;
const MAX_CAVALIER_SWEEP: f64 = PI;
const MAX_EXPLICIT_ARC_SWEEP: f64 = FRAC_PI_2;

/// Controls whether imported curved external shapes stay native or are converted to straight
/// segments for legacy/regression benchmarking.
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum CurvedGeometryMode {
    Native,
    Tessellated { tolerance: f64 },
}

/// Converts external representations of items and containers into internal ones.
#[derive(Clone, Debug, Copy)]
pub struct Importer {
    /// Modification config for item outer boundaries.
    /// Driven by `min_item_separation` (half the requested item-item clearance goes to each item).
    pub shape_modify_config: ShapeModifyConfig,
    /// Modification config for bin/container boundaries and static bin-internal hazards.
    /// Driven by `min_bin_separation`, compensated for item inflation so item-bin clearance can
    /// differ from item-item clearance.
    pub bin_modify_config: ShapeModifyConfig,
    /// Modification config for item-internal holes (inner rings of items). Applied with
    /// [`ShapeModifyMode::Deflate`] so a positive offset *shrinks* the free pocket, forcing items
    /// nested inside a hole to keep their distance from the surrounding parent's hole boundary.
    /// Driven by `min_hole_separation`.
    pub hole_modify_config: ShapeModifyConfig,
    pub cde_config: CDEConfig,
    pub curved_geometry_mode: CurvedGeometryMode,
}

impl Importer {
    /// Creates a new instance with the given configuration.
    ///
    /// * `cde_config` - Configuration for the CDE (Collision Detection Engine).
    /// * `simplify_tolerance` - See [`ShapeModifyConfig`].
    /// * `min_item_separation` - Optional minimum separation distance between items. Half is
    ///   applied to item outlines (Inflate). When `min_bin_separation` is not provided, this also
    ///   preserves the previous item-bin behavior by applying the other half to bins.
    /// * `min_hole_separation` - Optional minimum separation distance between an item nested inside
    ///   another item's hole and the surrounding hole's boundary. Applied as a Deflate offset on
    ///   the hole's free pocket.
    /// * `narrow_concavity_cutoff` - Optional definition for closing narrow concavities.
    pub fn new(
        cde_config: CDEConfig,
        simplify_tolerance: Option<f64>,
        min_item_separation: Option<f64>,
        min_hole_separation: Option<f64>,
        narrow_concavity_cutoff: Option<(f64, f64)>,
    ) -> Importer {
        Importer::new_with_separations(
            cde_config,
            simplify_tolerance,
            min_item_separation,
            min_item_separation,
            min_hole_separation,
            narrow_concavity_cutoff,
        )
    }

    /// Creates a new importer with independent item-item, item-bin, and item-hole clearances.
    pub fn new_with_separations(
        cde_config: CDEConfig,
        simplify_tolerance: Option<f64>,
        min_item_separation: Option<f64>,
        min_bin_separation: Option<f64>,
        min_hole_separation: Option<f64>,
        narrow_concavity_cutoff: Option<(f64, f64)>,
    ) -> Importer {
        let item_offset = min_item_separation.map(|sep| sep / 2.0);
        let item_offset_value = item_offset.unwrap_or(0.0);
        let bin_offset = min_bin_separation
            .map(|sep| sep - item_offset_value)
            .or(item_offset);
        let hole_offset = min_hole_separation.map(|sep| sep - item_offset_value);

        Importer {
            shape_modify_config: ShapeModifyConfig {
                offset: item_offset,
                simplify_tolerance,
                narrow_concavity_cutoff,
            },
            bin_modify_config: ShapeModifyConfig {
                offset: bin_offset,
                simplify_tolerance,
                narrow_concavity_cutoff,
            },
            hole_modify_config: ShapeModifyConfig {
                offset: hole_offset,
                simplify_tolerance,
                narrow_concavity_cutoff,
            },
            cde_config,
            curved_geometry_mode: CurvedGeometryMode::Native,
        }
    }

    /// Sets how curved input geometry should be represented internally.
    pub fn with_curved_geometry_mode(mut self, curved_geometry_mode: CurvedGeometryMode) -> Self {
        self.curved_geometry_mode = curved_geometry_mode;
        self
    }

    fn import_shape_as_simple_polygon(&self, shape: &ExtShape) -> Result<SPolygon> {
        let shape = import_shape_as_simple_polygon(shape)?;
        apply_curved_geometry_mode(shape, self.curved_geometry_mode)
    }

    pub fn import_item(&self, ext_item: &ExtItem) -> Result<Item> {
        debug!("[IMPORT] starting item {:?}", ext_item.id);

        // Holes (inner rings) are only meaningful for `ExtShape::Polygon`. For other
        // shape types they are conceptually empty.
        let mut hole_ext: Vec<ExtSPolygon> = vec![];

        let original_shape = {
            let shape = match &ext_item.shape {
                ExtShape::Polygon(ep) => {
                    // Preserve the inner rings; they will become hole pockets on the item.
                    hole_ext = ep.inner.clone();
                    import_simple_polygon(&ep.outer)?
                }
                shape => self.import_shape_as_simple_polygon(shape)?,
            };
            OriginalShape {
                pre_transform: centering_transformation(&shape),
                shape,
                modify_mode: ShapeModifyMode::Inflate,
                modify_config: self.shape_modify_config,
            }
        };

        // Use the SAME pre_transform as the outer so holes translate identically.
        let pre_transform = original_shape.pre_transform;
        let original_holes: Vec<OriginalShape> = hole_ext
            .into_iter()
            .map(|esp| -> Result<OriginalShape> {
                let shape = import_simple_polygon(&esp)?;
                Ok(OriginalShape {
                    shape,
                    pre_transform,
                    // Holes are free pockets carved out of the item. We Deflate them so a positive
                    // offset shrinks the pocket — items nested inside must keep their distance.
                    modify_mode: ShapeModifyMode::Deflate,
                    modify_config: self.hole_modify_config,
                })
            })
            .collect::<Result<_>>()?;

        let base_quality = ext_item.min_quality;

        let allowed_orientations = match ext_item.allowed_orientations.as_ref() {
            Some(a_o) => {
                if a_o.is_empty() || (a_o.len() == 1 && a_o[0] == 0.0) {
                    RotationRange::None
                } else {
                    RotationRange::Discrete(a_o.iter().map(|angle| angle.to_radians()).collect())
                }
            }
            None => RotationRange::Continuous,
        };

        Item::new_with_holes(
            ext_item.id as usize,
            original_shape,
            original_holes,
            allowed_orientations,
            base_quality,
            self.cde_config.item_surrogate_config,
        )
    }

    pub fn import_container(&self, ext_cont: &ExtContainer) -> Result<Container> {
        assert!(
            ext_cont.zones.iter().all(|zone| zone.quality < N_QUALITIES),
            "All quality zones must have lower quality than N_QUALITIES, set N_QUALITIES to a higher value if required"
        );

        let original_outer = {
            let outer = self.import_shape_as_simple_polygon(&ext_cont.shape)?;
            OriginalShape {
                shape: outer,
                pre_transform: DTransformation::empty(),
                modify_mode: ShapeModifyMode::Deflate,
                modify_config: self.bin_modify_config,
            }
        };

        let holes = match &ext_cont.shape {
            ExtShape::SimplePolygon(_)
            | ExtShape::Rectangle { .. }
            | ExtShape::Circle { .. }
            | ExtShape::BulgedPolygon(_)
            | ExtShape::Arcs(_) => vec![],
            ExtShape::Polygon(jp) => {
                let json_holes = &jp.inner;
                json_holes
                    .iter()
                    .map(import_simple_polygon)
                    .collect::<Result<Vec<SPolygon>>>()?
            }
            ExtShape::MultiPolygon(_) => {
                unimplemented!("No support for multipolygon shapes yet")
            }
        };

        let mut shapes_inferior_qzones = (0..N_QUALITIES)
            .map(|q| {
                ext_cont
                    .zones
                    .iter()
                    .filter(|zone| zone.quality == q)
                    .map(|zone| match &zone.shape {
                        ExtShape::Rectangle {
                            x_min,
                            y_min,
                            width,
                            height,
                        } => Rect::try_new(*x_min, *y_min, x_min + width, y_min + height)
                            .map(|r| r.into()),
                        ExtShape::SimplePolygon(esp) => import_simple_polygon(esp),
                        ExtShape::Circle { .. }
                        | ExtShape::BulgedPolygon(_)
                        | ExtShape::Arcs(_) => self.import_shape_as_simple_polygon(&zone.shape),
                        ExtShape::Polygon(_) => {
                            unimplemented!("No support for polygon to simplepolygon conversion yet")
                        }
                        ExtShape::MultiPolygon(_) => {
                            unimplemented!("No support for multipolygon shapes yet")
                        }
                    })
                    .collect::<Result<Vec<SPolygon>>>()
            })
            .collect::<Result<Vec<Vec<SPolygon>>>>()?;

        //merge the container holes with quality == 0
        shapes_inferior_qzones[0].extend(holes);

        //convert the shapes to inferior quality zones
        let quality_zones = shapes_inferior_qzones
            .into_iter()
            .enumerate()
            .map(|(q, zone_shapes)| {
                let original_shapes = zone_shapes
                    .into_iter()
                    .map(|s| OriginalShape {
                        shape: s,
                        pre_transform: DTransformation::empty(),
                        modify_mode: ShapeModifyMode::Inflate,
                        modify_config: self.bin_modify_config,
                    })
                    .collect_vec();
                InferiorQualityZone::new(q, original_shapes)
            })
            .collect::<Result<Vec<InferiorQualityZone>>>()?;

        Container::new(
            ext_cont.id as usize,
            original_outer,
            quality_zones,
            self.cde_config,
        )
    }
}

fn import_shape_as_simple_polygon(shape: &ExtShape) -> Result<SPolygon> {
    match shape {
        ExtShape::Rectangle {
            x_min,
            y_min,
            width,
            height,
        } => Rect::try_new(*x_min, *y_min, x_min + width, y_min + height).map(SPolygon::from),
        ExtShape::Circle { cx, cy, r } => import_circle(*cx, *cy, *r),
        ExtShape::SimplePolygon(esp) => import_simple_polygon(esp),
        ExtShape::BulgedPolygon(ebp) => import_bulged_polygon(ebp),
        ExtShape::Arcs(arcs) => import_arc_segments(arcs),
        ExtShape::Polygon(ep) => import_simple_polygon(&ep.outer),
        ExtShape::MultiPolygon(_) => bail!("No support for multipolygon shapes yet"),
    }
}

pub fn import_simple_polygon(sp: &ExtSPolygon) -> Result<SPolygon> {
    let vertices = sp.0.iter().map(|(x, y)| (*x, *y, 0.0)).collect_vec();
    import_bulged_vertices(vertices)
}

pub fn import_bulged_polygon(sp: &ExtBulgedPolygon) -> Result<SPolygon> {
    import_bulged_vertices(sp.0.clone())
}

fn apply_curved_geometry_mode(
    shape: SPolygon,
    curved_geometry_mode: CurvedGeometryMode,
) -> Result<SPolygon> {
    match curved_geometry_mode {
        CurvedGeometryMode::Native => Ok(shape),
        CurvedGeometryMode::Tessellated { tolerance } => {
            ensure!(
                tolerance.is_finite() && tolerance > 0.0,
                "curve tessellation tolerance must be positive and finite"
            );
            if !shape.has_arcs() {
                return Ok(shape);
            }
            let vertices = shape
                .tessellated_edge_iter(tolerance)
                .map(|edge| edge.start)
                .collect_vec();
            let tessellated = SPolygon::new(vertices)?;
            ensure_no_self_intersections(&tessellated)?;
            Ok(tessellated)
        }
    }
}

fn import_circle(cx: f64, cy: f64, r: f64) -> Result<SPolygon> {
    ensure!(
        cx.is_finite() && cy.is_finite(),
        "invalid circle center: ({cx}, {cy})"
    );
    ensure!(r.is_finite() && r > 0.0, "invalid circle radius: {r}");

    let quarter_bulge = (FRAC_PI_2 / 4.0).tan();
    import_bulged_vertices(vec![
        (cx + r, cy, quarter_bulge),
        (cx, cy + r, quarter_bulge),
        (cx - r, cy, quarter_bulge),
        (cx, cy - r, quarter_bulge),
    ])
}

fn import_arc_segments(arcs: &[ExtArcSeg]) -> Result<SPolygon> {
    ensure!(
        !arcs.is_empty(),
        "arc shape must contain at least one segment"
    );

    let mut points = Vec::new();
    let mut bulges = Vec::new();
    let mut first_point = None;
    let mut previous_end = None;

    for arc in arcs {
        ensure!(
            arc.center.0.is_finite()
                && arc.center.1.is_finite()
                && arc.radius.is_finite()
                && arc.start_angle.is_finite()
                && arc.sweep.is_finite(),
            "invalid arc segment: center={:?}, radius={}, start_angle={}, sweep={}",
            arc.center,
            arc.radius,
            arc.start_angle,
            arc.sweep
        );
        ensure!(arc.radius > 0.0, "invalid arc radius: {}", arc.radius);
        ensure!(
            arc.sweep.abs() > BULGE_EPSILON,
            "arc segment sweep must be non-zero"
        );

        let center = Point(arc.center.0, arc.center.1);
        let split_count = split_count_for_sweep(arc.sweep, MAX_EXPLICIT_ARC_SWEEP);
        let sweep_step = arc.sweep / split_count as f64;
        let bulge = (sweep_step / 4.0).tan();

        let start = point_on_circle(center, arc.radius, arc.start_angle);
        if let Some(end) = previous_end {
            ensure!(
                points_match(start, end),
                "arc segments are not contiguous: previous end {:?}, next start {:?}",
                end,
                start
            );
        } else {
            first_point = Some(start);
        }

        for part in 0..split_count {
            let angle = arc.start_angle + sweep_step * part as f64;
            points.push(point_on_circle(center, arc.radius, angle));
            bulges.push(bulge);
        }
        previous_end = Some(point_on_circle(
            center,
            arc.radius,
            arc.start_angle + arc.sweep,
        ));
    }

    ensure!(
        points_match(
            previous_end.expect("arc segments checked non-empty"),
            first_point.unwrap()
        ),
        "arc segments do not form a closed boundary"
    );

    import_boundary_parts(points, bulges)
}

fn import_bulged_vertices(vertices: Vec<(f64, f64, f64)>) -> Result<SPolygon> {
    let mut points = vertices.iter().map(|(x, y, _)| Point(*x, *y)).collect_vec();
    let mut bulges = vertices.iter().map(|(_, _, bulge)| *bulge).collect_vec();

    //Strip the last vertex if it is the same as the first one
    if points.len() > 1 && points_match(points[0], points[points.len() - 1]) {
        points.pop();
        bulges.pop();
    }
    //Remove duplicates that are consecutive (e.g. [1, 2, 2, 3] -> [1, 2, 3])
    eliminate_degenerate_boundary_vertices(&mut points, &mut bulges);
    //Bail if there are any non-consecutive duplicates.
    if points.len() != points.iter().unique().count() {
        bail!("Simple polygon has non-consecutive duplicate vertices");
    }

    import_boundary_parts(points, bulges)
}

fn import_boundary_parts(points: Vec<Point>, bulges: Vec<f64>) -> Result<SPolygon> {
    let (points, bulges) = split_large_bulges(points, bulges, MAX_CAVALIER_SWEEP)?;
    let shape = SPolygon::new_with_bulges(points, bulges)?;
    ensure_no_self_intersections(&shape)?;
    Ok(shape)
}

fn split_large_bulges(
    points: Vec<Point>,
    bulges: Vec<f64>,
    max_sweep: f64,
) -> Result<(Vec<Point>, Vec<f64>)> {
    ensure!(
        points.len() == bulges.len(),
        "one bulge is required per boundary vertex"
    );

    let mut split_points = Vec::new();
    let mut split_bulges = Vec::new();
    for i in 0..points.len() {
        let start = points[i];
        let end = points[(i + 1) % points.len()];
        let bulge = bulges[i];

        ensure!(
            start.0.is_finite() && start.1.is_finite(),
            "invalid polygon vertex: {start:?}"
        );
        ensure!(bulge.is_finite(), "invalid polygon bulge: {bulge}");

        if bulge.abs() <= (max_sweep / 4.0).tan() + BULGE_EPSILON {
            split_points.push(start);
            split_bulges.push(bulge);
            continue;
        }

        let arc = Arc::try_from_bulge(start, end, bulge)?;
        let split_count = split_count_for_sweep(arc.sweep, max_sweep);
        let sweep_step = arc.sweep / split_count as f64;
        let split_bulge = (sweep_step / 4.0).tan();
        for part in 0..split_count {
            let point = if part == 0 {
                start
            } else {
                arc.point_at_angle(arc.start_angle() + sweep_step * part as f64)
            };
            split_points.push(point);
            split_bulges.push(split_bulge);
        }
    }
    Ok((split_points, split_bulges))
}

fn split_count_for_sweep(sweep: f64, max_sweep: f64) -> usize {
    usize::max((sweep.abs() / max_sweep).ceil() as usize, 1)
}

fn ensure_no_self_intersections(shape: &SPolygon) -> Result<()> {
    let n_segments = shape.n_segments();
    for i in 0..n_segments {
        let segment_i = shape.segment(i);
        for j in (i + 1)..n_segments {
            if segments_are_adjacent(i, j, n_segments) {
                continue;
            }
            let segment_j = shape.segment(j);
            ensure!(
                !segment_i.collides_with(&segment_j),
                "Simple polygon boundary self-intersects between segments {i} and {j}"
            );
        }
    }
    Ok(())
}

fn segments_are_adjacent(i: usize, j: usize, n_segments: usize) -> bool {
    (i + 1) % n_segments == j || (j + 1) % n_segments == i
}

fn point_on_circle(center: Point, radius: f64, angle: f64) -> Point {
    Point(
        center.0 + radius * angle.cos(),
        center.1 + radius * angle.sin(),
    )
}

fn points_match(a: Point, b: Point) -> bool {
    approx_eq!(f64, a.0, b.0, epsilon = POINT_MATCH_EPSILON)
        && approx_eq!(f64, a.1, b.1, epsilon = POINT_MATCH_EPSILON)
}

/// Returns a transformation that translates the shape's centroid to the origin.
pub fn centering_transformation(shape: &SPolygon) -> DTransformation {
    let Point(cx, cy) = shape.centroid();
    DTransformation::new(0.0, (-cx, -cy))
}

/// Converts an external transformation (applicable to the original shapes) to an internal transformation (used within `jagua-rs`).
///
/// * `ext_transf` - The external transformation.
/// * `pre_transf` - The transformation that was applied to the original shape to derive the internal representation.
pub fn ext_to_int_transformation(
    ext_transf: &DTransformation,
    pre_transf: &DTransformation,
) -> DTransformation {
    //1. undo pre-transform
    //2. do the absolute transformation

    Transformation::empty()
        .transform(&pre_transf.compose().inverse())
        .transform_from_decomposed(ext_transf)
        .decompose()
}

pub fn eliminate_degenerate_vertices(points: &mut Vec<Point>) {
    let mut bulges = vec![0.0; points.len()];
    eliminate_degenerate_boundary_vertices(points, &mut bulges);
}

fn eliminate_degenerate_boundary_vertices(points: &mut Vec<Point>, bulges: &mut Vec<f64>) {
    let mut indices_to_remove = vec![];
    let n_points = points.len();
    for i in 0..n_points {
        let j = (i + 1) % n_points;
        let p_i = points[i];
        let p_j = points[j];
        if points_match(p_i, p_j) {
            //points are equal, mark for removal
            indices_to_remove.push(i);
        }
    }
    //remove points in reverse order to avoid shifting indices
    indices_to_remove.sort_unstable_by(|a, b| b.cmp(a));
    for index in indices_to_remove {
        if index < points.len() {
            let j = (index + 1) % points.len();
            debug!(
                "[IMPORT] degenerate vertex eliminated (idx: {}, {:?}, {:?})",
                index, points[index], points[j]
            );
            points.remove(index);
            bulges.remove(index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::ext_repr::ExtArcSeg;
    use std::f64::consts::PI;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-3,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn legacy_simple_polygon_import_stays_straight() {
        let shape = import_simple_polygon(&ExtSPolygon(vec![
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
        ]))
        .unwrap();

        assert!(!shape.has_arcs());
        assert_eq!(shape.bulges, vec![0.0; 4]);
        assert_close(shape.area, 2.0);
    }

    #[test]
    fn bulged_polygon_import_preserves_arc_segments() {
        let shape = import_bulged_polygon(&ExtBulgedPolygon(vec![
            (-1.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 2.0, 1.0),
            (-1.0, 2.0, 0.0),
        ]))
        .unwrap();

        assert!(shape.has_arcs());
        assert_eq!(shape.n_segments(), 4);
        assert_close(shape.area, 4.0 + PI / 2.0);
    }

    #[test]
    fn circle_shape_imports_as_exact_quarter_arcs() {
        let shape = import_shape_as_simple_polygon(&ExtShape::Circle {
            cx: 0.0,
            cy: 0.0,
            r: 2.0,
        })
        .unwrap();

        assert!(shape.has_arcs());
        assert_eq!(shape.n_segments(), 4);
        assert_close(shape.area, 4.0 * PI);
        assert!(shape.bulges.iter().all(|bulge| bulge.abs() <= 1.0));
    }

    #[test]
    fn large_bulges_are_split_for_cavalier_compatibility() {
        let shape = import_bulged_polygon(&ExtBulgedPolygon(vec![
            (0.0, 0.0, 2.0),
            (2.0, 0.0, 0.0),
            (1.0, -1.0, 0.0),
        ]))
        .unwrap();

        assert!(shape.has_arcs());
        assert!(shape.n_segments() > 3);
        assert!(shape.bulges.iter().all(|bulge| bulge.abs() <= 1.0));
    }

    #[test]
    fn explicit_arcs_import_closed_boundary() {
        let shape = import_shape_as_simple_polygon(&ExtShape::Arcs(vec![
            ExtArcSeg {
                center: (0.0, 0.0),
                radius: 1.0,
                start_angle: 0.0,
                sweep: PI,
            },
            ExtArcSeg {
                center: (0.0, 0.0),
                radius: 1.0,
                start_angle: PI,
                sweep: PI,
            },
        ]))
        .unwrap();

        assert!(shape.has_arcs());
        assert_eq!(shape.n_segments(), 4);
        assert_close(shape.area, PI);
        assert!(shape.bulges.iter().all(|bulge| bulge.abs() <= 1.0));
    }

    #[test]
    fn bulged_polygon_json_schema_deserializes() {
        let shape: ExtShape = serde_json::from_str(
            r#"{
                "type": "bulged_polygon",
                "data": [[-1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,2.0,1.0],[-1.0,2.0,0.0]]
            }"#,
        )
        .unwrap();

        let imported = import_shape_as_simple_polygon(&shape).unwrap();
        assert!(imported.has_arcs());
    }

    #[test]
    fn circle_json_schema_deserializes() {
        let shape: ExtShape = serde_json::from_str(
            r#"{
                "type": "circle",
                "data": { "cx": 0.0, "cy": 0.0, "r": 1.0 }
            }"#,
        )
        .unwrap();

        let imported = import_shape_as_simple_polygon(&shape).unwrap();
        assert_close(imported.area, PI);
    }

    #[test]
    fn tessellated_curved_geometry_mode_converts_curves_to_lines() {
        let native = import_shape_as_simple_polygon(&ExtShape::Circle {
            cx: 0.0,
            cy: 0.0,
            r: 1.0,
        })
        .unwrap();
        let tessellated =
            apply_curved_geometry_mode(native, CurvedGeometryMode::Tessellated { tolerance: 0.01 })
                .unwrap();

        assert!(!tessellated.has_arcs());
        assert!(tessellated.n_vertices() > 4);
        assert!((tessellated.area - PI).abs() < 0.05);
    }

    #[test]
    fn tessellated_curved_geometry_mode_leaves_straight_shapes_straight() {
        let straight = import_simple_polygon(&ExtSPolygon(vec![
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
        ]))
        .unwrap();
        let tessellated = apply_curved_geometry_mode(
            straight.clone(),
            CurvedGeometryMode::Tessellated { tolerance: 0.01 },
        )
        .unwrap();

        assert!(!tessellated.has_arcs());
        assert_eq!(tessellated.vertices, straight.vertices);
        assert_eq!(tessellated.bulges, straight.bulges);
    }

    #[test]
    fn self_intersecting_import_is_rejected() {
        let err = import_simple_polygon(&ExtSPolygon(vec![
            (0.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (2.0, 0.0),
            (3.0, 1.0),
        ]))
        .unwrap_err();

        assert!(err.to_string().contains("self-intersects"));
    }
}
