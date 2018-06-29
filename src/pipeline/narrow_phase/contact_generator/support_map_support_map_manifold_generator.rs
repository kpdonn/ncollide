use bounding_volume::PolyhedralCone;
use js_func::log;
use math::{Isometry, Point, Vector};
use na::{Real, Unit};
use pipeline::narrow_phase::{ContactDispatcher, ContactManifoldGenerator};
use query::algorithms::gjk::GJKResult;
use query::algorithms::VoronoiSimplex;
use query::contacts_internal;
use query::Contact;
use query::{ContactKinematic, ContactManifold, ContactPrediction};
use shape::{Ball, FeatureId, Shape};
use utils::IdAllocator;

/// Collision detector between two balls.
pub struct SupportMapSupportMapManifoldGenerator<N: Real> {
    manifold: ContactManifold<N>,
    simplex: VoronoiSimplex<N>,
    last_gjk_dir: Option<Unit<Vector<N>>>,
}

impl<N: Real> Clone for SupportMapSupportMapManifoldGenerator<N> {
    fn clone(&self) -> SupportMapSupportMapManifoldGenerator<N> {
        SupportMapSupportMapManifoldGenerator {
            manifold: self.manifold.clone(),
            simplex: self.simplex.clone(),
            last_gjk_dir: self.last_gjk_dir.clone(),
        }
    }
}

impl<N: Real> SupportMapSupportMapManifoldGenerator<N> {
    /// Creates a new persistent collision detector between two balls.
    #[inline]
    pub fn new() -> SupportMapSupportMapManifoldGenerator<N> {
        SupportMapSupportMapManifoldGenerator {
            manifold: ContactManifold::new(),
            simplex: VoronoiSimplex::new(),
            last_gjk_dir: None,
        }
    }
}

impl<N: Real> ContactManifoldGenerator<N> for SupportMapSupportMapManifoldGenerator<N> {
    fn update(
        &mut self,
        _: &ContactDispatcher<N>,
        ida: usize,
        ma: &Isometry<N>,
        a: &Shape<N>,
        idb: usize,
        mb: &Isometry<N>,
        b: &Shape<N>,
        prediction: &ContactPrediction<N>,
        id_alloc: &mut IdAllocator,
    ) -> bool {
        if let (Some(a), Some(b)) = (a.as_support_map(), b.as_support_map()) {
            self.manifold.set_subshape_id1(ida);
            self.manifold.set_subshape_id2(idb);
            self.manifold.save_cache_and_clear(id_alloc);

            let contact_result = contacts_internal::support_map_against_support_map_with_params(
                ma,
                a,
                mb,
                b,
                prediction.linear,
                &mut self.simplex,
                self.last_gjk_dir,
            );

            match contact_result {
                GJKResult::ClosestPoints(world1, world2, normal) => {
                    self.last_gjk_dir = Some(normal);

                    let contact = Contact::new_wo_depth(world1, world2, normal);
                    let mut kinematic = ContactKinematic::new();
                    let _ = self.manifold.push(contact, kinematic, id_alloc);
                }
                GJKResult::NoIntersection(dir) => {
                    self.last_gjk_dir = Some(dir);
                }
                GJKResult::Intersection => unreachable!(),
                GJKResult::Proximity(_) => unreachable!(),
            }

            true
        } else {
            false
        }
    }

    #[inline]
    fn num_contacts(&self) -> usize {
        self.manifold.len()
    }

    #[inline]
    fn contacts<'a: 'b, 'b>(&'a self, out: &'b mut Vec<&'a ContactManifold<N>>) {
        if self.manifold.len() != 0 {
            out.push(&self.manifold)
        }
    }
}
