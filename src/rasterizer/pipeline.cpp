// clang-format off
#include "pipeline.h"

#include <iostream>
#include <algorithm>

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include "framebuffer.h"
#include "sample_pattern.h"
template<PrimitiveType primitive_type, class Program, uint32_t flags>
void Pipeline<primitive_type, Program, flags>::run(std::vector<Vertex> const& vertices,
                                                   typename Program::Parameters const& parameters,
                                                   Framebuffer* framebuffer_) {
	// Framebuffer must be non-null:
	assert(framebuffer_);
	auto& framebuffer = *framebuffer_;

	// A1T7: sample loop
	// TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	//  	 This will probably involve inserting a loop of the form:
	// 		 	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//      	for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//   	 around some subset of the code.
	// 		 You will also need to transform the input and output of the rasterize_* functions to
	// 	     account for the fact they deal with pixels centered at (0.5,0.5).
    std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
    for (uint32_t s = 0; s < samples.size(); ++s) {
        Vec3 const &sample = samples[s];
        float offsetX = sample.x - 0.5f;
        float offsetY = sample.y - 0.5f;

        std::vector<ShadedVertex> shaded_vertices;
        shaded_vertices.reserve(vertices.size());

        //--------------------------
        // shade vertices:
        for (auto const &v: vertices) {
            ShadedVertex sv;
            Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
            shaded_vertices.emplace_back(sv);
        }

        //--------------------------
        // assemble + clip + homogeneous divide vertices:
        std::vector<ClippedVertex> clipped_vertices;

        // reserve some space to avoid reallocations later:
        if constexpr (primitive_type == PrimitiveType::Lines) {
            // clipping lines can never produce more than one vertex per input vertex:
            clipped_vertices.reserve(shaded_vertices.size());
        } else if constexpr (primitive_type == PrimitiveType::Triangles) {
            // clipping triangles can produce up to 8 vertices per input vertex:
            clipped_vertices.reserve(shaded_vertices.size() * 8);
        }
        // clang-format off

        //coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
        //x: [-1,1] -> [0,width]
        //y: [-1,1] -> [0,height]
        //z: [-1,1] -> [0,1] (OpenGL-style depth range)
        Vec3 const clip_to_fb_scale = Vec3{
                framebuffer.width / 2.0f,
                framebuffer.height / 2.0f,
                0.5f
        };
        Vec3 const clip_to_fb_offset = Vec3{
                0.5f * framebuffer.width,
                0.5f * framebuffer.height,
                0.5f
        };

        // helper used to put output of clipping functions into clipped_vertices:
        auto emit_vertex = [&](ShadedVertex const &sv) {
            ClippedVertex cv;
            float inv_w = 1.0f / sv.clip_position.w;
            cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;


            cv.inv_w = inv_w;
            cv.attributes = sv.attributes;
            clipped_vertices.emplace_back(cv);
        };

        // offset
//        std::vector<ClippedVertex> adjusted_vertices;
//        for (auto const v : clipped_vertices) {
//            ClippedVertex adjusted_v = v;
//            adjusted_v.fb_position.x += offsetX;
//            adjusted_v.fb_position.y += offsetY;
//            adjusted_vertices.push_back(adjusted_v);
//        }

        // actually do clipping:
        if constexpr (primitive_type == PrimitiveType::Lines) {
            for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
                clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
            }
        } else if constexpr (primitive_type == PrimitiveType::Triangles) {
            for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
                clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
            }
        } else {
            static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
        }

        //--------------------------
        // rasterize primitives:

        std::vector<Fragment> fragments;

        // helper used to put output of rasterization functions into fragments:
        auto emit_fragment = [&](Fragment const &f) { fragments.emplace_back(f); };

        // actually do rasterization:
        if constexpr (primitive_type == PrimitiveType::Lines) {
            for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
                ClippedVertex vertices1 = clipped_vertices[i];
                vertices1.fb_position.x += offsetX;
                vertices1.fb_position.y += offsetY;
                ClippedVertex vertices2 = clipped_vertices[i+1];
                vertices2.fb_position.x += offsetX;
                vertices2.fb_position.y += offsetY;
                rasterize_line(vertices1, vertices2, emit_fragment);
            }
        } else if constexpr (primitive_type == PrimitiveType::Triangles) {
            for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
//                std::printf("\nin sample %d\nadjusted vertices %f %f", s,adjusted_vertices[s].fb_position.x, adjusted_vertices[s].fb_position.y);
//                rasterize_triangle(clipped_vertices[i], clipped_vertices[i + 1], clipped_vertices[i + 2],
//                                   emit_fragment);
                ClippedVertex vertices1 = clipped_vertices[i];
                vertices1.fb_position.x += offsetX;
                vertices1.fb_position.y += offsetY;
                ClippedVertex vertices2 = clipped_vertices[i+1];
                vertices2.fb_position.x += offsetX;
                vertices2.fb_position.y += offsetY;
                ClippedVertex vertices3 = clipped_vertices[i+2];
                vertices3.fb_position.x += offsetX;
                vertices3.fb_position.y += offsetY;
                rasterize_triangle(vertices1, vertices2, vertices3,
                                   emit_fragment);
            }
        } else {
            static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
        }

        //--------------------------
        // depth test + shade + blend fragments:
        uint32_t out_of_range = 0; // check if rasterization produced fragments outside framebuffer
        // (indicates something is wrong with clipping)
        for (auto const &f: fragments) {

            // fragment location (in pixels):
            int32_t x = (int32_t) std::floor(f.fb_position.x);
            int32_t y = (int32_t) std::floor(f.fb_position.y);

            // if clipping is working properly, this condition shouldn't be needed;
            // however, it prevents crashes while you are working on your clipping functions,
            // so we suggest leaving it in place:
            if (x < 0 || (uint32_t) x >= framebuffer.width ||
                y < 0 || (uint32_t) y >= framebuffer.height) {
                ++out_of_range;
                continue;
            }

            // local names that refer to destination sample in framebuffer:
            float &fb_depth = framebuffer.depth_at(x, y, s);
            Spectrum &fb_color = framebuffer.color_at(x, y, s);

            // depth test:
            if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
                // "Always" means the depth test always passes.
            } else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
                // "Never" means the depth test never passes.
                continue; //discard this fragment
            } else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
                // "Less" means the depth test passes when the new fragment has depth less than the stored depth.
                // A1T4: Depth_Less
                // TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less".
                if (f.fb_position.z < fb_depth) {
                    fb_depth = f.fb_position.z;
                } else {
                    continue;
                }
            } else {
                static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
            }

            // if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
            if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
                fb_depth = f.fb_position.z;
            }

            // shade fragment:
            ShadedFragment sf;
            sf.fb_position = f.fb_position;
            Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

            // write color to framebuffer if color writes aren't disabled:
            if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {
                // blend fragment:
                if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
                    fb_color = sf.color;
                } else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
                    // A1T4: Blend_Add
                    // TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
                    fb_color += sf.color * sf.opacity; //<-- replace this line
                } else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
                    // A1T4: Blend_Over
                    // TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
                    // 		 You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
                    fb_color = sf.color + fb_color * (1 - sf.opacity); //<-- replace this line
                } else {
                    static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
                }
            }
        }
        if (out_of_range > 0) {
            if constexpr (primitive_type == PrimitiveType::Lines) {
                warn("Produced %d fragments outside framebuffer; this indicates something is likely "
                     "wrong with the clip_line function.",
                     out_of_range);
            } else if constexpr (primitive_type == PrimitiveType::Triangles) {
                warn("Produced %d fragments outside framebuffer; this indicates something is likely "
                     "wrong with the clip_triangle function.",
                     out_of_range);
            }
        }
    }
}

// -------------------------------------------------------------------------
// clipping functions

// helper to interpolate between vertices:
template<PrimitiveType p, class P, uint32_t F>
auto Pipeline<p, P, F>::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  	va, vb: endpoints of line
 *  	emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_line(ShadedVertex const& va, ShadedVertex const& vb,
                                      std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// Determine portion of line over which:
	// 		pt = (b-a) * t + a
	//  	-pt.w <= pt.x <= pt.w
	//  	-pt.w <= pt.y <= pt.w
	//  	-pt.w <= pt.z <= pt.w
	// ... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		// restrict range such that:
		// l + t * dl <= r + t * dr
		// re-arranging:
		//  l - r <= t * (dr - dl)
		if (dr == dl) {
			// want: l - r <= 0
			if (l - r > 0.0f) {
				// works for none of range, so make range empty:
				min_t = 1.0f;
				max_t = 0.0f;
			}
		} else if (dr > dl) {
			// since dr - dl is positive:
			// want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { // dr < dl
			// since dr - dl is negative:
			// want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	// local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
	}
}

/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  	va, vb, vc: vertices of triangle
 *  	emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// A1EC: clip_triangle
	// TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

// -------------------------------------------------------------------------
// rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result. 
 *    We will not explicitly check for 45 degree lines along the diamond edges (this will be extra credit),
 *    but you should be able to handle 45 degree lines in every other case (such as starting from pixel centers)
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function<void(Fragment const&)> const& emit_fragment) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	// A1T2: rasterize_line

	// TODO: Check out the block comment above this function for more information on how to fill in
	// this function!
	// The OpenGL specification section 3.5 may also come in handy.

	{
        // check which point is lower andw which axis is the main axis
        // Todo: do this later, now only consider a normal line √
        // Todo: hasn't considered the main axis

//        const ClippedVertex *pa = nullptr;  // should be nullptr
//        const ClippedVertex *pb = nullptr;

        int main_axis = (abs(va.fb_position.x-vb.fb_position.x) > abs(va.fb_position.y-vb.fb_position.y)) ? 1 : 0;

        float main_x1 = 0;
        float main_x2 = 0;
        float main_y1 = 0;
        float main_y2 = 0;
        float z1 = 0;
        float z2 = 0;
        // main axis: x
        if (main_axis==1 && va.fb_position.x<vb.fb_position.x){
            main_x1 = va.fb_position.x;
            main_x2 = vb.fb_position.x;
            main_y1 = va.fb_position.y;
            main_y2 = vb.fb_position.y;
            z1 = va.fb_position.z;
            z2 = vb.fb_position.z;
        }
        else if (main_axis==1 && va.fb_position.x>vb.fb_position.x) {
            main_x2 = va.fb_position.x;
            main_x1 = vb.fb_position.x;
            main_y2 = va.fb_position.y;
            main_y1 = vb.fb_position.y;
            z2 = va.fb_position.z;
            z1 = vb.fb_position.z;
        }
        else if (main_axis==0 && va.fb_position.y<vb.fb_position.y) {
            main_x1 = va.fb_position.y;
            main_x2 = vb.fb_position.y;
            main_y1 = va.fb_position.x;
            main_y2 = vb.fb_position.x;
            z1 = va.fb_position.z;
            z2 = vb.fb_position.z;
        }
        else if (main_axis==0 && va.fb_position.y>vb.fb_position.y) {
            main_x2 = va.fb_position.y;
            main_x1 = vb.fb_position.y;
            main_y2 = va.fb_position.x;
            main_y1 = vb.fb_position.x;
            z2 = va.fb_position.z;
            z1 = vb.fb_position.z;
        }

        float k = (main_y2 - main_y1) / (main_x2 - main_x1);
        float b = main_y1 - k*main_x1;

        // get the largest integer smaller than x
        uint32_t t1 = (uint32_t)floor(main_x1);
        uint32_t t2 = (uint32_t)floor(main_x2);

        for (uint32_t u = t1; u <= t2; u++) {
            float w = (float) (u + 0.5 - main_x1) / (main_x2 - main_x1);  // a percentage
            float v = (float) w * (main_y2 - main_y1) + main_y1;

            float z = (float) w * (z2 - z1) + z1;
            // raster a pixel (floor(u)+0.5, floor(v)+0.5)
            Fragment mid;  //middle of a pixel
            mid.fb_position = (main_axis==1) ? Vec3((float) floor(u) + 0.5f, (float) floor(v) + 0.5f, z) :  Vec3((float) floor(v) + 0.5f, (float) floor(u) + 0.5f, z);
            mid.attributes = va.attributes;
            mid.derivatives.fill(Vec2(0.0f, 0.0f));

            // diamoned exit rule
            float x, y;  // to get the bottom left point
            // deal with start and end point
            // special case: line in the same pixel
            std::cout << std::endl << "this is u,v " << u << " " << v;
            if (u==t1 && u==t2){
                main_x1 = (main_x1 - floor(main_x1));
                main_y1 = (main_y1 - floor(main_y1));
                main_x2 = (main_x2 - floor(main_x2));
                main_y2 = (main_y2 - floor(main_y2));
                int isEntering = main_x1+main_y1<1.5 && main_y1-main_x1>=-0.5 && !((main_x1==0.5) && (main_y1==1)) && !((main_x1==1) && (main_y1==0.5));
                int isExiting  = (main_x2+main_y2>1.5) || (main_y2-main_x2<-0.5) || ((main_x2==0.5) && (main_y2==1)) || ((main_x2==1) && (main_y2==0.5));

                if (isEntering && isExiting) {
                    emit_fragment(mid);
                }
                continue;
            }
            // start point
            else if (u == t1){
                std::cout << std::endl << "in start point";
                x = (float) (main_x1-floor(u));
                y = (float) (main_y1-floor(v));
//                b = b-floor(b); //translate to [0,1]
                if (b >= 0){
                    b = b-floor(b); //translate to [0,1]
                }
                else if (b >= -1 && b < 0){
                    // not sure what to do, just keep it
                }
                int isEntering = x+y<1.5 && y-x >= -0.5 && !((x==0.5) && (y==1)) && !((x==1) && (y==0.5));
                int isExiting = ((1.5-b)/(1+k) >= 0.5 && (1.5-b)/(1+k) <= 1) || ((b+0.5)/(1-k) > 0.5 && (b+0.5)/(1-k)<1);

                if (isEntering && isExiting) {
                    emit_fragment(mid);
                }
                // horizontal line cross the bottom of the diamond
                else if (x==0 && k == 0) {
                    emit_fragment(mid);
                }
                continue;
            }
            // end point
            else if (u == t2) {
                std::cout << std::endl << "in end point";
                x = (float) (main_x2-floor(u));
                y = (float) (main_y2-floor(v));

                int isExiting = x+y>1.5 || y-x<-0.5 || ((x==0.5) && (y==1)) || ((x==1) && (y==0.5));
                if (isExiting) {
                    emit_fragment(mid);
                }
                continue;
            }
//          shade all other points
            emit_fragment(mid);
        }
        // to make the program run
        // should be deleted
//        Fragment mid;  //middle of a pixel
//        mid.fb_position = va.fb_position;
//        mid.attributes = va.attributes;
//        mid.derivatives.fill(Vec2(0.0f, 0.0f));
//        emit_fragment(mid);
    }
}


/*
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  	(x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Smooth: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 * 	The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function<void(Fragment const&)> const& emit_fragment) {
	// NOTE: it is okay to restructure this function to allow these tasks to use the
	//  same code paths. Be aware, however, that all of them need to remain working!
	//  (e.g., if you break Flat while implementing Correct, you won't get points
	//   for Flat.)
    auto barycentric_coordinates = [&] (Vec2 a, Vec2 b, Vec2 c, Vec2 q) {
        // get barycentric coordinates
        float Sabc = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        float Sqbc = (b.x - q.x) * (c.y - q.y) - (c.x - q.x) * (b.y - q.y);
        float Sqac = (c.x - q.x) * (a.y - q.y) - (a.x - q.x) * (c.y - q.y);
        float Sqab = (a.x - q.x) * (b.y - q.y) - (b.x - q.x) * (a.y - q.y);

        float alpha = Sqbc / Sabc;
        float beta = Sqac / Sabc;
        float theta = Sqab / Sabc;

        return Vec3(alpha, beta, theta);
    };

    auto VectorFromPoints = [&](const Vec2& a, const Vec2& b) {
        return Vec2{b.x - a.x, b.y - a.y};
    };

    auto Cross = [&](Vec2 const& va, Vec2 const& vb){
        return va.x*vb.y - va.y*vb.x;
    };

    auto isHorizontal = [&](Vec2 & p1, Vec2 & p2, Vec2 & p3) {
        return p1.y == p2.y;
    };

    // Cross product
     auto onEdge = [&](const Vec2& p1, const Vec2& p2, const Vec2& q) {
         Vec2 pq = Vec2{q.x - p1.x, q.y - p1.y};
         Vec2 p2p1 = Vec2{p2.x - p1.x, p2.y - p1.y};

        float cross = pq.x * p2p1.y - pq.y * p2p1.x;

        return fabs(cross) < 1e-6;
    };

//     printf("\nva.fb.position.x %f va.fb.position.y %f", va.fb_position.x, va.fb_position.y);
    // derivatives
    Vec2 a_a = Vec2{va.fb_position.x, va.fb_position.y};
    Vec2 b_b = Vec2{vb.fb_position.x, vb.fb_position.y};
    Vec2 c_c = Vec2{vc.fb_position.x, vc.fb_position.y};
    Vec2 q_a = a_a;
    Vec2 q_a_up = Vec2{a_a.x, a_a.y+1};
    Vec2 q_a_right = Vec2{a_a.x+1, a_a.y};
    Vec3 q_a_barycentric = barycentric_coordinates(a_a,b_b,c_c,q_a);
    Vec3 q_a_up_barycentric = barycentric_coordinates(a_a,b_b,c_c,q_a_up);
    Vec3 q_a_right_barycentric = barycentric_coordinates(a_a,b_b,c_c,q_a_right);
    std::array<float, 5> q_a_attributes = {0,0,0,0,0};
    std::array<float, 5> q_a_up_attributes = {0,0,0,0,0};
    std::array<float, 5> q_a_right_attributes = {0,0,0,0,0};
    for (int i = 0; i < (int)(sizeof(va.attributes)/sizeof(va.attributes[0])); ++i) {
        q_a_attributes[i] = q_a_barycentric.x * va.attributes[i] + q_a_barycentric.y * vb.attributes[i] + q_a_barycentric.z * vc.attributes[i];
        q_a_up_attributes[i] = q_a_up_barycentric.x * va.attributes[i] + q_a_up_barycentric.y * vb.attributes[i] + q_a_up_barycentric.z * vc.attributes[i];
        q_a_right_attributes[i] = q_a_right_barycentric.x * va.attributes[i] + q_a_right_barycentric.y * vb.attributes[i] + q_a_right_barycentric.z * vc.attributes[i];
    }
    Vec2 tri_derivatives = Vec2{q_a_right_attributes[0]-q_a_attributes[0], q_a_up_attributes[0]-q_a_attributes[0]};
    // end derivatives

    if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		// A1T3: flat triangles
		// TODO: rasterize triangle (see block comment above this function).
        Vec2 a = Vec2{va.fb_position.x, va.fb_position.y};
        Vec2 b = Vec2{vb.fb_position.x, vb.fb_position.y};
        Vec2 c = Vec2{vc.fb_position.x, vc.fb_position.y};
        Vec2 ac = VectorFromPoints(a, c);
        Vec2 ab = VectorFromPoints(a, b);
        Vec2 bc = VectorFromPoints(b, c);
        Vec2 ba = VectorFromPoints(b, a);
        Vec2 cb = VectorFromPoints(c, b);
        Vec2 ca = VectorFromPoints(c, a);

        // get the leftmost/upmost edge
        Vec2 leftEdgeStart = a;
        Vec2 leftEdgeEnd = b;
        Vec2 topEdgeStart = a;
        Vec2 topEdgeEnd = b;
        bool hasTopEdge = false;
        // there must have a left edge, but may have a top edge
        if (isHorizontal(a, b, c)) {
            if (c.y < a.y){
                hasTopEdge = true;
                topEdgeStart = a; topEdgeEnd = b;
            }
            leftEdgeStart = a.x < b.x ? a : b;
            leftEdgeEnd = c;
        } else if (isHorizontal(a, c, b)) {
            if (b.y < a.y){
                hasTopEdge = true;
                topEdgeStart = a; topEdgeEnd = c;
            }
            leftEdgeStart = a.x < c.x ? a : c;
            leftEdgeEnd = b;
        } else if (isHorizontal(b, c, a)) {
            if (a.y < b.y){
                hasTopEdge = true;
                topEdgeStart = b; topEdgeEnd = c;
            }
            leftEdgeStart = b.x < c.x ? b : c;
            leftEdgeEnd = a;
        }
        // no horizontal line
        else {
            Vec2* leftMost = &a;
            if (b.x < leftMost->x || (b.x == leftMost->x && b.y < leftMost->y)) {
                leftMost = &b;
            }
            if (c.x < leftMost->x || (c.x == leftMost->x && c.y < leftMost->y)) {
                leftMost = &c;
            }

            if (leftMost == &a) {
                leftEdgeStart = a;
                leftEdgeEnd = (b.x == a.x && b.y > a.y) ? b : c;
            } else if (leftMost == &b) {
                leftEdgeStart = b;
                leftEdgeEnd = (c.x == b.x && c.y > b.y) ? c : a;
            } else { // leftMost == &c
                leftEdgeStart = c;
                leftEdgeEnd = (a.x == c.x && a.y > c.y) ? a : b;
            }
        }

        float minX = std::min({va.fb_position.x, vb.fb_position.x, vc.fb_position.x});
        float maxX = std::max({va.fb_position.x, vb.fb_position.x, vc.fb_position.x});
        float minY = std::min({va.fb_position.y, vb.fb_position.y, vc.fb_position.y});
        float maxY = std::max({va.fb_position.y, vb.fb_position.y, vc.fb_position.y});

        // judge whether a point is inside the triangle
        for (float x = floor(minX); x <= floor(maxX); ++x) {
            for (float y = floor(minY); y <= floor(maxY); ++y) {
                float centerX = x + 0.5f;
                float centerY = y + 0.5f;
                Vec2 q = Vec2{centerX, centerY};
                Vec2 cq = VectorFromPoints(c, q);
                Vec2 aq = VectorFromPoints(a, q);
                Vec2 bq = VectorFromPoints(b, q);
                float crossa = Cross(ab, aq);
                float crossb = Cross(bc, bq);
                float crossc = Cross(ca, cq);

                Fragment frag;
                frag.fb_position = Vec3(centerX, centerY, va.fb_position.z);
                frag.attributes = va.attributes;
                // top edge
                if (hasTopEdge && onEdge(topEdgeStart, topEdgeEnd, q)) {
                    emit_fragment(frag);
                }
                // left edge
                else if(onEdge(leftEdgeStart, leftEdgeEnd, q)){
                    emit_fragment(frag);
                }
                // in the triangle
                else if ((crossa > 0 && crossb > 0 && crossc > 0) ||
                    (crossa < 0 && crossb < 0 && crossc < 0)) {
                    emit_fragment(frag);
                }
                // else: do nothing
            }
        }
        // As a placeholder, here's code that draws some lines:
		//(remove this and replace it with a real solution)
//		Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(va, vb, emit_fragment);
//		Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(vb, vc, emit_fragment);
//		Pipeline<PrimitiveType::Lines, P, flags>::rasterize_line(vc, va, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Smooth) {
		// A1T5: screen-space smooth triangles
		// TODO: rasterize triangle (see block comment above this function).

        Vec2 a = Vec2{va.fb_position.x, va.fb_position.y};
        Vec2 b = Vec2{vb.fb_position.x, vb.fb_position.y};
        Vec2 c = Vec2{vc.fb_position.x, vc.fb_position.y};

        float minX = std::min({a.x, b.x, c.x});
        float maxX = std::max({a.x, b.x, c.x});
        float minY = std::min({a.y, b.y, c.y});
        float maxY = std::max({a.y, b.y, c.y});

        for (float x = floor(minX); x <= floor(maxX); ++x) {
            for (float y = floor(minY); y <= floor(maxY); ++y) {
                float centerX = x + 0.5f;
                float centerY = y + 0.5f;
                Vec2 q = Vec2{centerX, centerY};
                Vec2 q_up = Vec2{centerX, centerY+1};
                Vec2 q_right = Vec2{centerX+1, centerY};
                Vec3 q_barycentric = barycentric_coordinates(a,b,c,q);
                Vec3 q_up_barycentric = barycentric_coordinates(a,b,c,q_up);
                Vec3 q_right_barycentric = barycentric_coordinates(a,b,c,q_right);

                if (q_barycentric.x > 0 && q_barycentric.y > 0 && q_barycentric.z > 0){
                    std::array<float, 5> q_attributes = {0,0,0,0,0};
                    std::array<float, 5> q_up_attributes = {0,0,0,0,0};
                    std::array<float, 5> q_right_attributes = {0,0,0,0,0};

                    for (int i = 0; i < (int)sizeof(va.attributes)/(int)sizeof(va.attributes[0]); ++i) {
                        q_attributes[i] = q_barycentric.x * va.attributes[i] + q_barycentric.y * vb.attributes[i] + q_barycentric.z * vc.attributes[i];
                        q_up_attributes[i] = q_up_barycentric.x * va.attributes[i] + q_up_barycentric.y * vb.attributes[i] + q_up_barycentric.z * vc.attributes[i];
                        q_right_attributes[i] = q_right_barycentric.x * va.attributes[i] + q_right_barycentric.y * vb.attributes[i] + q_right_barycentric.z * vc.attributes[i];
                    }

                    Fragment frag;
                    // todo: z interpolation
                    frag.fb_position = Vec3(centerX, centerY, q_barycentric.x * va.fb_position.z + q_barycentric.y * vb.fb_position.z + q_barycentric.z * vc.fb_position.z);
                    frag.attributes = q_attributes;
                    frag.derivatives = {tri_derivatives, Vec2(0.0f)};

                    emit_fragment(frag);
                }
            }
        }

		// As a placeholder, here's code that calls the Flat interpolation version of the function:
		//(remove this and replace it with a real solution)
//		Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Flat>::rasterize_triangle(va, vb, vc, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		// A1T5: perspective correct triangles
		// TODO: rasterize triangle (block comment above this function).
        Vec2 a = Vec2{va.fb_position.x, va.fb_position.y};
        Vec2 b = Vec2{vb.fb_position.x, vb.fb_position.y};
        Vec2 c = Vec2{vc.fb_position.x, vc.fb_position.y};

        float minX = std::min({a.x, b.x, c.x});
        float maxX = std::max({a.x, b.x, c.x});
        float minY = std::min({a.y, b.y, c.y});
        float maxY = std::max({a.y, b.y, c.y});

        for (float x = floor(minX); x <= floor(maxX); ++x) {
            for (float y = floor(minY); y <= floor(maxY); ++y) {
                float centerX = x + 0.5f;
                float centerY = y + 0.5f;
                Vec2 q = Vec2{centerX, centerY};

                Vec3 q_barycentric = barycentric_coordinates(a, b, c, q);

                float interpolated_inv_w = q_barycentric.x * va.inv_w + q_barycentric.y * vb.inv_w + q_barycentric.z * vc.inv_w;

                float w = 1.0f / interpolated_inv_w;

                std::array<float, 5> q_attributes = {0,0,0,0,0};

                for (int i = 0; i < (int)sizeof(va.attributes)/(int)sizeof(va.attributes[0]); ++i) {
                    q_attributes[i] = (q_barycentric.x * va.attributes[i] * va.inv_w +
                                       q_barycentric.y * vb.attributes[i] * vb.inv_w +
                                       q_barycentric.z * vc.attributes[i] * vc.inv_w) * w;
                }

                if (q_barycentric.x > 0 && q_barycentric.y > 0 && q_barycentric.z > 0) {
                    Fragment frag;
                    frag.fb_position = Vec3(centerX, centerY,
                                            (q_barycentric.x * va.fb_position.z * va.inv_w +
                                             q_barycentric.y * vb.fb_position.z * vb.inv_w +
                                             q_barycentric.z * vc.fb_position.z * vc.inv_w) * w);
                    frag.attributes = q_attributes;
                    frag.derivatives = {tri_derivatives, Vec2(0.0f)};

                    emit_fragment(frag);
                }
            }
        }

		// As a placeholder, here's code that calls the Screen-space interpolation function:
		//(remove this and replace it with a real solution)
//		Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Smooth>::rasterize_triangle(va, vb, vc, emit_fragment);
	}
}

//-------------------------------------------------------------------------
// compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;