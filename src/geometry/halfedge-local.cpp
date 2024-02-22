
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge

    bisect_edge(e);

    HalfedgeRef h = e->halfedge->next;
    HalfedgeRef t = e->halfedge->twin;

    // edge case 1: line
    if (h->face->boundary && t->face->boundary){
        return std::nullopt;
    }

//    std::cout << h->face->boundary << " " << t->face->boundary;
    bool is_boundary = false;
    if (!h->face->boundary && !t->face->boundary){  // not boundary
        // do nothing
    }
    else {
        // change the reference of h and t if needed. Do not split on t side
        is_boundary = true;
//        std::cout << "\nthere is boundary";
        if (h->face->boundary){
            HalfedgeRef mid_h = h;
            h = t;
            t = mid_h;
        }
    }

    VertexRef vbi = t->vertex;
    VertexRef vh = h->next->next->vertex;  // vertex at h side
    VertexRef vt = t->next->next->vertex;  // vertex at t side
    assert(t->vertex == h->vertex);

    HalfedgeRef for_h = t->twin;
    HalfedgeRef for_t = h->twin;
    HalfedgeRef lat_h = h->next;
    HalfedgeRef lat_t = t->next;

    // Create new elements
    EdgeRef e_h_n = emplace_edge();  // edge_h_new
    EdgeRef e_t_n = emplace_edge();  // edge_twin_new
    FaceRef f_h = emplace_face();  // face_h
    FaceRef f_t = emplace_face();  // face_twin
    FaceRef f_old_h = h->face;
    FaceRef f_old_t = t->face;
    HalfedgeRef h_n = emplace_halfedge();  // h_new
//    interpolate_data({})
    HalfedgeRef h_n_t = emplace_halfedge();  //h_new_twin
    HalfedgeRef t_n = emplace_halfedge();  // twin_new
    HalfedgeRef t_n_t = emplace_halfedge();  // twin_new_twin
    h_n->twin = h_n_t;
    h_n_t->twin = h_n;
    t_n->twin = t_n_t;
    t_n_t->twin = t_n;
    e_h_n->halfedge = h_n;
    e_t_n->halfedge = t_n;
    f_h->halfedge = h_n;
    f_t->halfedge = t_n;

    // Reassign connectivity
    h_n_t->next = lat_h->next;
    h_n_t->vertex = h->vertex;
    h_n_t->edge = e_h_n;
    lat_h->next = h_n;
    for_h->next = h_n_t;

    if (!is_boundary){
        t_n_t->next = lat_t->next;
        t_n_t->vertex = t->vertex;
        t_n_t->edge = e_t_n;
        lat_t->next = t_n;
        for_t->next = t_n_t;
    }

    h_n->next = h;
    h_n->vertex = vh;
    h_n->edge = e_h_n;

    if (!is_boundary) {
        t_n->next = t;
        t_n->vertex = vt;
        t_n->edge = e_t_n;
    }

    // update new face
    h->face = f_h;
    lat_h->face = f_h;
    h_n->face = f_h;
    h_n_t->face = for_h->face;

    // update old face
    f_old_h->halfedge = for_h;
    f_old_t->halfedge = for_t;

    if (!is_boundary) {
        t->face = f_t;
        lat_t->face = f_t;
        t_n->face = f_t;
        t_n_t->face = for_t->face;
    }

//    std::cout << "\n h id " << h->id << "\n twin id " << t->id << "\n";
    if (is_boundary){
        // delete useless element on boundary
        erase_edge(e_t_n);
        erase_face(f_t);
        erase_halfedge(t_n);
        erase_halfedge(t_n_t);
    }

//    std::cout << describe();

    return vbi;

//	(void)e; //this line avoids 'unused parameter' warnings. You can delete it as you fill in the function.
//    std::cout << describe();
//    return std::nullopt;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

    // f->halfedge remains unchanged
    HalfedgeRef h = f->halfedge;

    // todo: square is the same as other shape(include triangle)
    // 1 . create the inner face
    do {
        // create a repeated line, which means there're 3 lines starts from a vertex for the face
        // that need extrude
        HalfedgeRef temp_h = emplace_halfedge();
        HalfedgeRef temp_h_t = emplace_halfedge();  // temp_h_twin, parallel to h
        EdgeRef temp_e = emplace_edge();
        // update the line
        temp_h->vertex = h->next->vertex;
        temp_h_t->vertex = h->vertex;
        temp_h->twin = temp_h_t;
        temp_h_t->twin = temp_h;
        temp_h->edge = temp_e;
        temp_h_t->edge = temp_e;
        temp_e->halfedge = temp_h;
        // todo: face hasn't been set
        temp_h->next = temp_h_t;
        temp_h_t->next = h->next;
        h->next = temp_h;

        // update h to next(original next)
        h = temp_h_t->next;
    } while (h != f->halfedge);

    // 2. created the line (and the face?) between inner face and original face
    // create 1-2 first
    do {
        // create the link line, for l4, is 0-2 / 1-3 ...
        HalfedgeRef temp_h = emplace_halfedge();  // line 1->2
        HalfedgeRef temp_h_t = emplace_halfedge();  // line 2->1
        VertexRef temp_v = emplace_vertex();  // vertex 2, initialize, the same as v1
        temp_v->position = h->next->vertex->position;
        // for newly added edge, only update line1-2  one time
        EdgeRef temp_e = emplace_edge();  // line 1-2
        // todo: update face and edge
        // update line and vertex
        // update h1
        /*
         * 1
         * | \
         * |  2
         * |   |
         * |   |
         * |  3(not emplaced yet)
         * |
         * 0
         */
        // todo: twin not updated yet
        // todo: how to judge whether the edge of a halfedge exists?
        // update line 1->2 (edge added here)
        temp_h->vertex = h->next->vertex;  // vertex 1
        temp_h->next = h->next;  // h->next: line 2->3 for now
        h->next = temp_h;  // h->next: line 1->2 now

        // add edge and twin
        temp_h->edge = temp_e;
        temp_e->halfedge = temp_h;
        temp_h->twin = temp_h_t;

        // add vertex 2
        temp_h->next->vertex = temp_v;  // vertex of line 2->3
        temp_v->halfedge = temp_h->next;

        // update line 2->1
        temp_h_t->vertex = temp_h->next->vertex;  // vertex 2
        temp_h_t->next = temp_h->next->next->next;  // 1->2 => 2->3 => 3-2.next()
        temp_h->next->next->next = temp_h_t;  // 3->2.next() should be 2->1

        // add edge and twin
        temp_h_t->edge = temp_e;
        temp_h_t->twin = temp_h;

        // update h
        h = temp_h_t->next;
    } while (h != f->halfedge);

    // update connectivity line 2<->3
    // h's 3->2 can only be updated after h is updated
    do {
        // update 2->3 (of this polygon)
        HalfedgeRef lat_h_t = h->next->twin;
//        std::cout << "\nopposite line twin is " << h->next->next->twin->id;
//        std::cout << "\nopposite line twin next is " << lat_h_t->next->next->next->twin->id;
        h->next->next->twin->next = lat_h_t->next->next->next->twin; //line 2->3
//        std::cout << "\nopposite line twin next is " << h->next->next->twin->id;

        // update h
        h = h->next->twin->next;

        // update 3->2 (of next polygon)
        h->next->next->next = lat_h_t;
        h->next->next->twin->vertex = lat_h_t->vertex;
    } while (h!=f->halfedge);

    // 3. update face
    do {
        // update face within each polygon
        HalfedgeRef h_within = h;
        FaceRef temp_f = emplace_face();
        h->next->next->twin->face = h->face;  // opposite side.twin

//        std::cout << "\nid of h_face is " << h->face->id;

        do {
            h_within->face = temp_f;
            // update h_within
            h_within = h_within->next;
        } while (h_within != h);
        temp_f->halfedge = h_within;

        // update h
        h = h->next->twin->next;
    } while (h != f->halfedge);

    // update f's reference
    f->halfedge = f->halfedge->next->next->twin;

//    std::cout << describe();

    return f;
//    (void)f;
//    return std::nullopt;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	//A2L1: Flip Edge

    // edge case 1: on boundary
    if (e->on_boundary()){
        return std::nullopt;
    }

    // not on boundary
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;
    VertexRef vh = h->vertex;
    VertexRef vt = t->vertex;
//    std::cout << "\nhere is h " << h->id;
//    std::cout << "\nhere is h.twin " << h->twin->id;
//    std::cout << "\nhere is h.next " << h->next->id;
//    std::cout << "\nhere is h.twin.next " << h->twin->next->vertex->id;

    // update 4 vertex
    vh->halfedge = t->next;
    vt->halfedge = h->next;
    t->next->next->vertex->halfedge = h;
    h->next->next->vertex->halfedge = t;

    // update face?

    // update <next> of the halfedge before
    assert(e->halfedge == h);
    HalfedgeRef for_h = e->halfedge;
    HalfedgeRef for_t = e->halfedge->twin;
    while (for_h->next != h){
        for_h = for_h->next;
    }
    for_h->next = t->next;
    while (for_t->next != t){
        for_t = for_t->next;
    }
    for_t->next = h->next;

    // next's face will change
    h->next->face = t->face;
    t->next->face = h->face;

    // update h, t: twin, face, edge unchanged
    h->vertex = t->next->next->vertex;
    t->vertex = h->next->next->vertex;

    h->next = h->next->next;
    t->next = t->next->next;

    // updaate <next> of the halfedge after (should be modified after h, t has changed)
    assert(e->halfedge == h);
    HalfedgeRef lat_h = e->halfedge;
    HalfedgeRef lat_t = e->halfedge->twin;
    while (lat_h->next != t->next){
        lat_h = lat_h->next;
    }
    lat_h->next = h;
    while (lat_t->next != h->next){
        lat_t = lat_t->next;
    }
    lat_t->next = t;

//    std::cout << describe();

    return e;
    // for null operation
//    return std::nullopt;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    // edge case 1: line
    if (e->halfedge->face == e->halfedge->twin->face){
        return std::nullopt;
    }

    // Collect elements
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;

    HalfedgeRef lat_h = h->next;
    HalfedgeRef lat_t = t->next;
//    std::cout<<"\nlat_h_next id " << lat_h->next->id;
    HalfedgeRef for_h = h;
    HalfedgeRef for_t = t;
    while (for_h->next != h){
        for_h = for_h->next;
    }
    while (for_t->next != t){
        for_t = for_t->next;
    }
    VertexRef vh = h->vertex;
    VertexRef vt = t->vertex;
    VertexRef vm = emplace_vertex();
    vm->position = (vh->position + vt->position) / 2;
    interpolate_data({vh, vt}, vm);

//    std::cout<<"\nlat_h id " << lat_h->id;
//    std::cout<<"\nfor_h id " << for_h->id;

    // edge case: triangle: collpased into a line
    // erase face if is a triangle
    /*
     * 0---1
     * \   |
     *  \  |
     *   \ |
     *    \2
     */
    bool isTriangle = false;
    // h: line 0-1; lat_h: line 0-2; dup_h: line 1-2
    HalfedgeRef dup_h = lat_h->next;  // assume 0-1 collapse, and keep 0-2 instead of 1-2
    HalfedgeRef dup_h_t = dup_h->twin;  // 1->2
    if (h->next->next->next == h){
        isTriangle = true;
//        std::cout << "\n dup_h is " << dup_h->id;
//        std::cout << "\n dup_h_t is " << dup_h_t->id;
        erase_face(dup_h->face);
        lat_h->next = dup_h_t->next;
        lat_h->face = dup_h_t->face;
        HalfedgeRef for_dup_h_t = dup_h_t;
        while (for_dup_h_t->next != dup_h_t){
            for_dup_h_t = for_dup_h_t->next;
        }
        for_dup_h_t->next = lat_h;

        // todo: what if dup_h_t_vertex/face refer to dup_h_t?
        if (dup_h->vertex->halfedge == dup_h){
            dup_h->vertex->halfedge = lat_h->next;
        }
        if (dup_h_t->face->halfedge == dup_h_t){
            dup_h_t->face->halfedge = lat_h;
        }

        erase_edge(dup_h->edge);
        erase_halfedge(dup_h);
        erase_halfedge(dup_h_t);
    }

    // No new elements created except new vertex
    // Reassign connectivity

    // change the halfedge of face to make sure it is not the erased one
    // e.g. the boundary face not refer to the erased halfedge
    if (t->face->halfedge==t){
        t->face->halfedge = lat_t;
    }

    for_h->next = lat_h;
    for_t->next = lat_t;

    // vertex concerning vm
    lat_h->vertex = vm;
    lat_t->vertex = vm;
//    lat_h->twin->vertex = vm;
//    lat_t->twin->vertex = vm;
    for_t->twin->vertex = vm;
    for_h->twin->vertex = vm;

    lat_t->twin->next->vertex = vm;  // if boundary, next is for_h_twin
                                     // if not boundary, next is a side of triangle with vertex vm
    if (!lat_t->edge->on_boundary()){  // not boundary
        for_h->twin->next->vertex = vm;
    }
    lat_h->twin->next->vertex = vm;
    if (lat_h->edge->on_boundary()){
        for_t->twin->next->vertex = vm;
    }

    vm->halfedge = lat_h;

    // delete extra elements
    erase_edge(e);
    erase_halfedge(h);
    erase_halfedge(t);
    erase_vertex(vh);
    erase_vertex(vt);
//    if (isTriangle){
//        erase_edge(lat_h_next->edge);
////        std::cout<<"\nlat_h_next_twin id " << lat_h_next_twin->id;
//        erase_halfedge(lat_h_next->twin);
//        erase_halfedge(lat_h_next);
//    }

//    std::cout<<describe();

    return vm;
//    return std::nullopt;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

    // edge case 1: shrink > 1 face is a point. < 0 can always work?
    // no other edge case
    if (shrink > 1){
        shrink = 1;
    }
    else{
        HalfedgeRef h = face->halfedge;
        VertexRef v = face->halfedge->vertex;
        VertexRef v_centroid = emplace_vertex();

//        std::cout << "\n move is " << move;
//        std::cout << "\n shrink is " << shrink;

        float v_size = 0.0f;
        do {
            v = h->vertex;
            v_centroid->position += v->position;

            v_size += 1.0f;

            // update h
            h = h->next;
        } while (h != face->halfedge);
        v_centroid->position /= v_size;
//        std::cout << "\n Centre v position " << v_centroid->position;

        assert(h == face->halfedge);
        do {
            VertexRef temp_v = h->vertex;
            VertexRef v0 = h->twin->next->vertex;

            temp_v->position = v0->position + shrink*(v_centroid->position-v0->position);
            temp_v->position += move;
            interpolate_data({v0, v_centroid}, temp_v);
//            std::cout << "\n changing v's position " << temp_v->id << " to " << temp_v->position;

            // update h
            h = h->next;
        } while (h != face->halfedge);

        erase_vertex(v_centroid);
    }


//    std::cout << "\n" << describe();
}

