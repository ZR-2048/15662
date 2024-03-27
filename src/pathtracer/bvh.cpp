
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.

	//TODO
    root_idx = buildHelper(0, primitives.size()-1, max_leaf_size);
}

// A helper function to recursively build BVH
template<typename Primitive>
size_t BVH<Primitive>::buildHelper(size_t start, size_t end, size_t max_leaf_size) {
    // termination condition, leaf
    size_t size = end - start + 1;
//    printf("\n start and end is %d %d, size is %d", (int)start, (int)end, (int)size);
    if (size <= max_leaf_size){
//        printf("\n in leaf, start and end index is %d, %d", (int)start, (int)end);
        BBox leaf_box = primitives[start].bbox();
        for (size_t i=start+1; i <= end; i++){
            leaf_box.enclose(primitives[i].bbox());
        }
        // leaf has no children l and r, thus, l==r and is a leaf
        return new_node(leaf_box, start, size);
    }

    // bucket partition
    // initialize best cost
    float best_cost = std::numeric_limits<float>::max();
    size_t best_partition = 0;
    uint32_t best_axis = 0;   // 0 for x, 1 for y, 2 for z1
    BBox best_left_box=bbox(), best_right_box=bbox();

    for (uint32_t axis=0; axis<3; axis++){
        // rearrange primitives on every axis
        std::sort(primitives.begin()+start, primitives.begin()+end+1, [axis](const Primitive& a, const Primitive& b) -> bool {
            return (float)a.bbox().center()[axis] < (float)b.bbox().center()[axis];
        });

        // initialize buckets
        size_t num_bucket = 2;
        std::vector<SAHBucketData> buckets(num_bucket);
        for (auto& bucket : buckets) {
            bucket.bb = BBox();
            bucket.num_prims = 0;
        }

        // allocate primitives to buckets
        for (size_t i=start; i<=end; i++){  // include end
            i = (uint32_t)i;
            size_t bucket_idx = compute_bucket(primitives[i].bbox().center()[axis], primitives[start].bbox().center()[axis], primitives[end].bbox().center()[axis], num_bucket);
            buckets[bucket_idx].bb.enclose(primitives[i].bbox());
            buckets[bucket_idx].num_prims++;
        }

        size_t num_partition = num_bucket-1;

        // note!! partition is 0 here, but is not the actual index, neec convertion
        for (size_t partition=0; partition<num_partition; partition++){  // n buckets have n-1 partition
            BBox left_bbox, right_bbox;
            size_t left_count = 0, right_count = 0;
            // left part
            for (size_t k=0; k<=partition; k++){
                left_bbox.enclose(buckets[k].bb);
                left_count += buckets[k].num_prims;
            }
            // right part
            for (size_t k=partition+1; k<num_bucket; k++){
                right_bbox.enclose(buckets[k].bb);
                right_count += buckets[k].num_prims;
            }
            // SAH
            float current_cost = left_bbox.surface_area() * left_count + right_bbox.surface_area() * right_count;
            if (current_cost < best_cost){
                best_cost = current_cost;
                best_axis = axis;
                best_partition = start + left_count-1;  // conversion here, count->index, need-1
                best_left_box = left_bbox;
                best_right_box = right_bbox;
            }
        }
    }
    // todo: recurse here? yes! out of for loop!
    // todo: sort primitives using best axis here?
    // rearrange primitives based on best axis
    std::sort(primitives.begin()+start, primitives.begin()+end+1, [best_axis](const Primitive& a, const Primitive& b) -> bool {
        return (size_t)a.bbox().center()[best_axis] < (size_t)b.bbox().center()[best_axis];
    });

    // create current node
    // recurse build left and right here
//    printf("\n start and end is %d %d, size is %d", (int)start, (int)end, (int)size);
//    printf("\n best partition is %d", (int)best_partition);
    size_t left_child_index = buildHelper(start, best_partition, max_leaf_size);
    size_t right_child_index = buildHelper(best_partition+1, end, max_leaf_size);
    BBox current_box;
    current_box.enclose(best_left_box);
    current_box.enclose(best_right_box);
    return new_node(current_box, start, end-start+1, left_child_index, right_child_index);
}

template<typename Primitive>
size_t BVH<Primitive>::compute_bucket(float center, float bounds_min, float bounds_max, size_t num_buckets){
    float normalized_position = (center - bounds_min) / (bounds_max - bounds_min);
    size_t bucket_index = static_cast<size_t>(num_buckets * normalized_position);
    // constrain index in a reasonable range
    bucket_index = std::min(bucket_index, num_buckets - 1);
    return bucket_index;
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:

    Trace closest;
    closest.distance = std::numeric_limits<float>::infinity();
    find_closest_hit(ray, root_idx, closest);
    return closest;

//    Trace ret;
//    for(const Primitive& prim : primitives) {
//        Trace hit = prim.hit(ray);
//        ret = Trace::min(ret, hit);
//    }
//    return ret;
}

template<typename Primitive>
void BVH<Primitive>::find_closest_hit(const Ray& ray, size_t node_idx, Trace& closest) const {
    // nodes empty
    if(nodes.size()==0){
        return;
    }

    const Node& node = nodes[node_idx];
    Vec2 hit_times = Vec2(0.f, std::numeric_limits<float>::infinity());
    if (!node.bbox.hit(ray, hit_times)) {
        return;
    }

    if (hit_times.x > closest.distance) {
        return;
    }

    if (node.is_leaf()) {
        for (size_t i = node.start; i < node.start + node.size; ++i) {
            Trace hit = primitives[i].hit(ray);
            if (hit.hit && hit.distance < closest.distance) {
                closest = hit;
            }
        }
    } else {
        // method 1
//        find_closest_hit(ray, node.l, closest);
//        find_closest_hit(ray, node.r, closest);

        // front-to-back
        Vec2 hit_times_child1 = Vec2(0.f, std::numeric_limits<float>::infinity()), hit_times_child2 = Vec2(0.f, std::numeric_limits<float>::infinity());
        bool hit_child1 = nodes[node.l].bbox.hit(ray, hit_times_child1);
        bool hit_child2 = nodes[node.r].bbox.hit(ray, hit_times_child2);

        size_t first_child_idx, second_child_idx;
        if (hit_child1 && (!hit_child2 || hit_times_child1.x < hit_times_child2.x)) {
            first_child_idx = node.l;
            second_child_idx = node.r;
        } else if (hit_child2) {
            first_child_idx = node.r;
            second_child_idx = node.l;
        } else {
            return;
        }

        find_closest_hit(ray, first_child_idx, closest);

        if ((hit_child1 && hit_times_child1.y < closest.distance) ||
            (hit_child2 && hit_times_child2.y < closest.distance)) {
            find_closest_hit(ray, second_child_idx, closest);
        }
    }
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
