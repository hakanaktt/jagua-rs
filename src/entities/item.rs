use std::sync::Arc;

use crate::geometry::OriginalShape;
use crate::geometry::fail_fast::SPSurrogateConfig;
use crate::geometry::geo_enums::RotationRange;
use crate::geometry::primitives::SPolygon;

use anyhow::Result;

/// Item to be produced.
#[derive(Clone, Debug)]
pub struct Item {
    pub id: usize,
    /// Original contour of the item as defined in the input
    pub shape_orig: Arc<OriginalShape>,
    /// Contour of the item to be used for collision detection
    pub shape_cd: Arc<SPolygon>,
    /// Original holes (inner rings) of the item, in the same coordinate system as `shape_orig`.
    /// Empty for items without holes (the common case).
    pub holes_orig: Vec<Arc<OriginalShape>>,
    /// Holes of the item in the internal (CD) coordinate system, ready to be transformed alongside `shape_cd`.
    /// `holes_cd[i]` corresponds to `holes_orig[i]`. Holes are converted with `ShapeModifyMode::Deflate`,
    /// so when `min_hole_separation > 0` the free pocket shrinks (items nested inside must stay clear of the boundary).
    pub holes_cd: Vec<Arc<SPolygon>>,
    /// Allowed rotations in which the item can be placed
    pub allowed_rotation: RotationRange,
    /// The minimum quality the item should be produced out of, if `None` the item requires full quality
    pub min_quality: Option<usize>,
    /// Configuration for the surrogate generation
    pub surrogate_config: SPSurrogateConfig,
}

impl Item {
    pub fn new(
        id: usize,
        original_shape: OriginalShape,
        allowed_rotation: RotationRange,
        min_quality: Option<usize>,
        surrogate_config: SPSurrogateConfig,
    ) -> Result<Item> {
        Self::new_with_holes(
            id,
            original_shape,
            vec![],
            allowed_rotation,
            min_quality,
            surrogate_config,
        )
    }

    /// Like [`Item::new`] but accepts inner rings (holes) that other items may eventually nest inside.
    /// At Phase A this only affects area accounting and SVG rendering — collision detection still treats
    /// the outer ring as a solid region.
    pub fn new_with_holes(
        id: usize,
        original_shape: OriginalShape,
        original_holes: Vec<OriginalShape>,
        allowed_rotation: RotationRange,
        min_quality: Option<usize>,
        surrogate_config: SPSurrogateConfig,
    ) -> Result<Item> {
        let shape_orig = Arc::new(original_shape);
        let shape_int = {
            let mut shape_int = shape_orig.convert_to_internal()?;
            shape_int.generate_surrogate(surrogate_config)?;
            Arc::new(shape_int)
        };
        let mut holes_orig = Vec::with_capacity(original_holes.len());
        let mut holes_cd = Vec::with_capacity(original_holes.len());
        for h in original_holes {
            let h_orig = Arc::new(h);
            let h_int = Arc::new(h_orig.convert_to_internal()?);
            holes_orig.push(h_orig);
            holes_cd.push(h_int);
        }
        Ok(Item {
            id,
            shape_orig,
            shape_cd: shape_int,
            holes_orig,
            holes_cd,
            allowed_rotation,
            min_quality,
            surrogate_config,
        })
    }

    /// Net area of the item (outer minus the area of all holes).
    pub fn area(&self) -> f32 {
        let outer = self.shape_orig.area();
        let holes: f32 = self.holes_orig.iter().map(|h| h.area()).sum();
        (outer - holes).max(0.0)
    }
}
