use crate::collision_detection::quadtree::qt_traits::QTQueryable;
use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::{BoundarySegment, Rect, SPolygon};
use std::sync::Arc;

/// Defines a set of boundary segments from a hazard that is partially active in the [`QTNode`](crate::collision_detection::quadtree::QTNode).
#[derive(Clone, Debug)]
pub struct QTHazPartial {
    /// The boundary segments that are active in the quadtree-node.
    pub segments: Vec<BoundarySegment>,
    /// A bounding box that guarantees all segments are contained within it. (used for fail fast)
    pub ff_bbox: Rect,
}

impl QTHazPartial {
    pub fn from_entire_shape(shape: &SPolygon) -> Self {
        let segments = shape.segment_iter().collect();
        let ff_bbox = shape.bbox;
        Self { segments, ff_bbox }
    }

    /// Like [`Self::from_entire_shape`] but additionally folds the boundary segments of any inner rings
    /// (holes) into the partial hazard. The bounding box stays as the outer's bbox since holes
    /// always lie inside the outer.
    pub fn from_entire_shape_with_holes(shape: &SPolygon, holes: &[Arc<SPolygon>]) -> Self {
        let mut segments: Vec<BoundarySegment> = shape.segment_iter().collect();
        for h in holes {
            segments.extend(h.segment_iter());
        }
        let ff_bbox = shape.bbox;
        Self { segments, ff_bbox }
    }

    pub fn from_parent(parent: &QTHazPartial, restricted_segments: Vec<BoundarySegment>) -> Self {
        debug_assert!(!restricted_segments.is_empty());
        debug_assert!(restricted_segments
            .iter()
            .all(|segment| parent.segments.contains(segment)));
        let ff_bbox = {
            //calculate a bounding box around the boundary segments
            if parent.segments.len() == restricted_segments.len() {
                // If the segments cover the entire shape, use the shape's bounding box
                parent.ff_bbox
            } else {
                // Otherwise, calculate the bounding box from the segments
                let (mut x_min, mut y_min, mut x_max, mut y_max) = (
                    f32::INFINITY,
                    f32::INFINITY,
                    f32::NEG_INFINITY,
                    f32::NEG_INFINITY,
                );
                for segment in &restricted_segments {
                    let bbox = segment.bbox();
                    x_min = x_min.min(bbox.x_min);
                    y_min = y_min.min(bbox.y_min);
                    x_max = x_max.max(bbox.x_max);
                    y_max = y_max.max(bbox.y_max);
                }
                if x_min < x_max && y_min < y_max {
                    Rect {
                        x_min,
                        y_min,
                        x_max,
                        y_max,
                    }
                } else {
                    // If the segments are all aligned to an axis, use the parent bounding box
                    parent.ff_bbox
                }
            }
        };

        Self {
            segments: restricted_segments,
            ff_bbox,
        }
    }

    pub fn n_segments(&self) -> usize {
        self.segments.len()
    }
}

impl<T: QTQueryable> CollidesWith<T> for QTHazPartial {
    fn collides_with(&self, entity: &T) -> bool {
        // If the entity does not collide with the bounding box of the hazard, it cannot collide with the hazard
        entity.collides_with(&self.ff_bbox)
            && self
                .segments
                .iter()
                .any(|segment| entity.collides_with(segment))
    }
}
