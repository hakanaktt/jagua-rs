use crate::entities::{Instance, Item};
use crate::geometry::DTransformation;
use crate::io::import::{Importer, ext_to_int_transformation};
use crate::probs::bpp::entities::{
    BPInstance, BPLayoutType, BPPlacement, BPProblem, BPSolution, Bin,
};
use crate::probs::bpp::io::ext_repr::{ExtBPInstance, ExtBPSolution};
use itertools::Itertools;
use rayon::prelude::*;

use anyhow::{Result, ensure};

/// Imports an instance into the library
pub fn import_instance(importer: &Importer, ext_instance: &ExtBPInstance) -> Result<BPInstance> {
    let items = {
        let mut items = ext_instance
            .items
            .par_iter()
            .map(|ext_item| {
                let item = importer.import_item(&ext_item.base)?;
                let demand = ext_item.demand as usize;
                Ok((item, demand))
            })
            .collect::<Result<Vec<(Item, usize)>>>()?;

        items.sort_by_key(|(item, _)| item.id);
        items.retain(|(_, demand)| *demand > 0);

        ensure!(
            items.iter().enumerate().all(|(i, (item, _))| item.id == i),
            "All items should have consecutive IDs starting from 0. IDs: {:?}",
            items.iter().map(|(item, _)| item.id).sorted().collect_vec()
        );
        ensure!(
            !items.is_empty(),
            "ExtBPInstance must have at least one item with positive demand"
        );

        items
    };

    let bins = {
        let mut bins: Vec<Bin> = ext_instance
            .bins
            .par_iter()
            .map(|ext_bin| {
                let container = importer.import_container(&ext_bin.base)?;
                Ok(Bin::new(container, ext_bin.stock, ext_bin.cost))
            })
            .collect::<Result<Vec<Bin>>>()?;

        bins.sort_by_key(|bin| bin.id);
        bins.retain(|bin| bin.stock > 0);
        ensure!(
            bins.iter().enumerate().all(|(i, bin)| bin.id == i),
            "All bins should have consecutive IDs starting from 0. IDs: {:?}",
            bins.iter().map(|bin| bin.id).sorted().collect_vec()
        );
        ensure!(
            !bins.is_empty(),
            "ExtBPInstance must have at least one bin with positive stock"
        );

        bins
    };

    Ok(BPInstance::new(items, bins))
}

/// Imports a solution into the library.
pub fn import_solution(instance: &BPInstance, ext_solution: &ExtBPSolution) -> BPSolution {
    let mut prob = BPProblem::new(instance.clone());

    for ext_layout in ext_solution.layouts.iter() {
        let bin_id = ext_layout.container_id as usize;
        let mut layout_id = BPLayoutType::Closed { bin_id };
        for ext_placement in ext_layout.placed_items.iter().cloned() {
            let item_id = ext_placement.item_id as usize;
            let d_transf = {
                let ext_transf = DTransformation::from(ext_placement.transformation);
                let item = instance.item(item_id);
                ext_to_int_transformation(&ext_transf, &item.shape_orig.pre_transform)
            };
            let (lkey, _) = prob.place_item(BPPlacement {
                layout_id,
                item_id,
                d_transf,
            });
            layout_id = BPLayoutType::Open(lkey);
        }
    }

    prob.save()
}
